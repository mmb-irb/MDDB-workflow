from mddb_workflow.utils.auxiliar import load_yaml, load_json, is_glob
from mddb_workflow.utils.cache import Cache
from mddb_workflow.utils.file import File
from mddb_workflow.utils.constants import CACHE_FILENAME
from mddb_workflow.utils.type_hints import *
import pandas as pd
import subprocess
import jinja2
import time
import os
import glob
import sqlite3
import sys
import threading
import shutil
from pathlib import Path
from enum import Enum
from contextlib import contextmanager
import importlib.util


def _stream_output(stream, log_file_handle, output_stream):
    """Stream output from subprocess to both log file and terminal."""
    for line in stream:
        output_stream.write(line)
        output_stream.flush()
        log_file_handle.write(line)
        log_file_handle.flush()


class State(Enum):
    """Enumeration of possible workflow states."""
    NEW = 'new'
    RUNNING = 'running'
    ERROR = 'error'
    DONE = 'done'


class DatabaseLock:
    """Directory-based locking for SQLite database access on distributed filesystems.

    Uses os.mkdir() which is atomic across distributed filesystems like NFS and BeeGFS.
    This provides reliable locking where fcntl.flock() fails on network filesystems.

    Note: This implementation uses exclusive locking for both read and write operations
    to ensure consistency on distributed filesystems. While this is more conservative
    than true reader-writer locks, it guarantees correctness across all nodes.
    """

    def __init__(self, db_path: str | Path, timeout: float = 30.0, retry_interval: float = 0.1):
        """Initialize the database lock.

        Args:
            db_path: Path to the SQLite database file.
            timeout: Maximum time to wait for lock acquisition (seconds).
            retry_interval: Time between lock acquisition attempts (seconds).

        """
        self.db_path = Path(db_path).resolve()
        self.lock_dir = self._path_to_lock_dir(self.db_path)
        self.timeout = timeout
        self.retry_interval = retry_interval
        self._lock_count = 0  # For reentrant locking
        self._current_lock_type = None

    @staticmethod
    def _path_to_lock_dir(db_path: str | Path) -> Path:
        """Get the lock directory path for a given database path."""
        db_path = Path(db_path).resolve()
        return db_path.with_name('.lock_' + db_path.name)

    def acquire(self) -> bool:
        """Acquire the database lock using atomic directory creation."""
        # Handle reentrant locking
        if self._lock_count > 0:
            self._lock_count += 1
            return True

        start_time = time.time()
        while True:
            try:
                # os.mkdir is atomic across distributed filesystems
                os.mkdir(self.lock_dir)
                self._lock_count = 1
                break
            except FileExistsError:
                if time.time() - start_time >= self.timeout:
                    raise TimeoutError(f"Could not acquire lock on {self.lock_dir} within {self.timeout} seconds")
                time.sleep(self.retry_interval)
            except OSError as e:
                # Handle other OS errors (permission denied, etc.)
                raise RuntimeError(f"Failed to acquire lock: {e}") from e

    def release(self):
        """Release the database lock by removing the lock directory."""
        if self._lock_count > 0:
            self._lock_count -= 1
            if self._lock_count == 0:
                try:
                    os.rmdir(self.lock_dir)
                except FileNotFoundError:
                    pass  # Lock was already released (shouldn't happen normally)
                except OSError as e:
                    # Log warning but don't raise - lock file might be stale
                    import warnings
                    warnings.warn(f"Failed to release lock directory {self.lock_dir}: {e}")
                self._current_lock_type = None

    def force_release(self):
        """Force release the lock, useful for cleaning up stale locks.

        Use with caution - only call this if you're sure no other process holds the lock.
        """
        try:
            os.rmdir(self.lock_dir)
        except FileNotFoundError:
            pass
        except OSError:
            pass
        self._lock_count = 0
        self._current_lock_type = None

    def is_locked(self) -> bool:
        """Check if the lock is currently held by any process."""
        return self.lock_dir.exists()

    def close(self):
        """Release any held locks."""
        while self._lock_count > 0:
            self.release()

    def __del__(self):
        """Ensure lock is released on garbage collection."""
        self.close()

    @contextmanager
    def lock_file(self):
        """Context manager for acquiring an exclusive lock."""
        self.acquire()
        try:
            yield
        finally:
            self.release()


class Dataset:
    """Class to manage and process a dataset of MDDB projects and their MDs (replicas/subprojects) using a central SQLite database."""
    common_columns = ['state', 'message', 'process_id', 'slurm_job_id', 'last_modified']
    project_columns = ['uuid', 'rel_path', 'num_mds'] + common_columns
    md_columns = ['uuid', 'project_uuid', 'rel_path'] + common_columns
    date_format = "%H:%M:%S %d-%m-%Y"

    def __init__(self, dataset_path: str = 'mwf_ds.sql', lock_timeout: float = 30.0):
        """Initialize the Dataset object and connect to SQLite DB.

        Args:
            dataset_path (str): Path to the dataset storage file, normally an .db file.
            lock_timeout (float): Maximum time to wait for database lock (seconds).

        """
        self.dataset_path = Path(dataset_path).resolve()
        # if not self.dataset_path.exists():
        #     raise ValueError(f"Dataset path {self.dataset_path} does not exist.")
        self.root_path = self.dataset_path.parent
        self._lock = DatabaseLock(self.dataset_path, timeout=lock_timeout)
        self.conn = sqlite3.connect(dataset_path, check_same_thread=False)
        # Enable foreign key constraints (disabled by default in SQLite)
        self.conn.execute("PRAGMA foreign_keys = ON")
        self._ensure_tables()

    @property
    def locked_storage_file(self):
        """Context manager for shared (read/write) database access."""
        return self._lock.lock_file()

    def close(self):
        """Close the database connection and release locks."""
        self._lock.close()
        if self.conn:
            self.conn.close()

    def _abs_to_rel(self, abs_path: str | Path) -> str:
        """Convert any path to relative path from root_path."""
        abs_path = Path(abs_path).resolve()
        return os.path.relpath(str(abs_path), start=str(self.root_path))

    def _rel_to_abs(self, rel_path: str | Path) -> str:
        """Convert any path to absolute path from root_path."""
        rel_path = Path(rel_path)
        return str((self.root_path / rel_path).resolve())

    def _ensure_tables(self):
        """Create tables if they do not exist."""
        with self.locked_storage_file:
            cur = self.conn.cursor()
            # Create projects table
            cur.execute('''
                CREATE TABLE IF NOT EXISTS projects (
                    uuid TEXT PRIMARY KEY NOT NULL,
                    rel_path TEXT UNIQUE NOT NULL,
                    num_mds INTEGER DEFAULT 0,
                    state TEXT,
                    message TEXT,
                    process_id INTEGER,
                    slurm_job_id TEXT,
                    last_modified TEXT
                )
            ''')
            # Create mds table
            # MDs are linked to projects via project_uuid foreign key
            # ON DELETE CASCADE ensures MDs are deleted when project is deleted
            cur.execute('''
                CREATE TABLE IF NOT EXISTS mds (
                    uuid TEXT PRIMARY KEY NOT NULL,
                    project_uuid TEXT NOT NULL,
                    rel_path TEXT UNIQUE NOT NULL,
                    state TEXT,
                    message TEXT,
                    process_id INTEGER,
                    slurm_job_id TEXT,
                    last_modified TEXT,
                    FOREIGN KEY (project_uuid) REFERENCES projects(uuid) ON DELETE CASCADE
                )
            ''')

            # Backward compatibility: add runtime tracking columns to existing databases.
            self._ensure_column_exists(cur, 'projects', 'process_id', 'INTEGER')
            self._ensure_column_exists(cur, 'projects', 'slurm_job_id', 'TEXT')
            self._ensure_column_exists(cur, 'mds', 'process_id', 'INTEGER')
            self._ensure_column_exists(cur, 'mds', 'slurm_job_id', 'TEXT')

            # Create triggers to automatically maintain num_mds count
            # Trigger on INSERT: increment num_mds when a new MD is added
            cur.execute('''
                CREATE TRIGGER IF NOT EXISTS increment_num_mds
                AFTER INSERT ON mds
                BEGIN
                    UPDATE projects
                    SET num_mds = num_mds + 1
                    WHERE uuid = NEW.project_uuid;
                END
            ''')

            # Trigger on DELETE: decrement num_mds when an MD is removed
            cur.execute('''
                CREATE TRIGGER IF NOT EXISTS decrement_num_mds
                AFTER DELETE ON mds
                BEGIN
                    UPDATE projects
                    SET num_mds = num_mds - 1
                    WHERE uuid = OLD.project_uuid;
                END
            ''')

            self.conn.commit()

    @staticmethod
    def _ensure_column_exists(cur: sqlite3.Cursor, table: str, column: str, column_type: str):
        """Add a column to an existing table if it doesn't exist."""
        cur.execute(f"PRAGMA table_info({table})")
        columns = {row[1] for row in cur.fetchall()}
        if column not in columns:
            cur.execute(f"ALTER TABLE {table} ADD COLUMN {column} {column_type}")

    @staticmethod
    def _is_pid_running(process_id: int | str | None) -> bool | None:
        """Return whether a PID is alive in the current host namespace."""
        if process_id in (None, ''):
            return None
        try:
            pid = int(process_id)
        except (TypeError, ValueError):
            return None
        if pid <= 0:
            return None
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return False
        except PermissionError:
            return True
        except OSError:
            return False
        return True

    @staticmethod
    def _check_active_slurm_jobs(job_ids: list[str]) -> dict[str, bool]:
        """Check active SLURM jobs with a single squeue call.

        Returns an empty dict if squeue is unavailable or fails.
        """
        unique_job_ids = sorted({str(job_id).strip() for job_id in job_ids if str(job_id).strip()})
        if not unique_job_ids:
            return {}
        if shutil.which('squeue') is None:
            return {}

        cmd = ['squeue', '-h', '-o', '%A', '-j', ','.join(unique_job_ids)]
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=False, timeout=10)
        except (OSError, subprocess.SubprocessError):
            return {}

        if result.returncode != 0:
            return {}

        active_jobs = {line.strip() for line in result.stdout.splitlines() if line.strip()}
        return {job_id: job_id in active_jobs for job_id in unique_job_ids}

    def _is_running_entry_active(self, process_id: int | str | None, slurm_job_id: str | None, active_slurm_jobs: dict[str, bool]) -> bool | None:
        """Determine whether a RUNNING entry is still active."""
        slurm_value = str(slurm_job_id).strip() if slurm_job_id else ''
        if slurm_value and slurm_value in active_slurm_jobs:
            return active_slurm_jobs[slurm_value]

        pid_state = self._is_pid_running(process_id)
        if pid_state is not None:
            return pid_state

        return None

    @staticmethod
    def _stale_running_message(process_id: int | str | None, slurm_job_id: str | None) -> str:
        """Build a clear status message for externally terminated jobs."""
        details = []
        if process_id not in (None, ''):
            details.append(f"pid={process_id}")
        if slurm_job_id not in (None, ''):
            details.append(f"slurm_job_id={slurm_job_id}")
        details_text = f" ({', '.join(details)})" if details else ''
        return (
            f"Workflow process is no longer running{details_text}. "
            "Likely terminated externally (e.g. SLURM time limit, OOM, or manual cancellation)."
        )

    def _refresh_running_entries(self, scopes: list[str] | None = None):
        """Mark stale RUNNING entries as ERROR when backing process/job no longer exists."""
        scopes = scopes or ['projects', 'mds']
        running_entries: list[tuple[str, dict]] = []

        for table in scopes:
            cols = self.project_columns if table == 'projects' else self.md_columns
            for row in self.query_table(table, query_state=[State.RUNNING.value]):
                running_entries.append((table, dict(zip(cols, row))))

        if not running_entries:
            return

        slurm_job_ids = [entry['slurm_job_id'] for _, entry in running_entries if entry.get('slurm_job_id')]
        active_slurm_jobs = self._check_active_slurm_jobs(slurm_job_ids)

        for table, entry in running_entries:
            is_active = self._is_running_entry_active(
                entry.get('process_id'),
                entry.get('slurm_job_id'),
                active_slurm_jobs,
            )
            if is_active is False:
                self.update_status(
                    uuid=entry['uuid'],
                    project_uuid=entry.get('project_uuid') if table == 'mds' else None,
                    state=State.ERROR,
                    message=self._stale_running_message(entry.get('process_id'), entry.get('slurm_job_id')),
                )

    def _read_uuid_from_cache(self, directory: str, make_uuid: bool = False, project_uuid: str = None) -> tuple[str, str | None]:
        """Read UUID and project_uuid from cache file in the given directory."""
        directory = Path(directory).resolve()
        if not directory.exists():
            raise ValueError(f"Directory {directory} does not exist.")
        cache_file = directory / CACHE_FILENAME
        if cache_file.exists():
            cache = load_json(cache_file)
            uuid = cache.get('uuid')
            project_uuid = cache.get('project_uuid')
            if uuid:
                return uuid, project_uuid
        if make_uuid:
            cache = Cache(File(str(cache_file)), project_uuid=project_uuid)
            uuid = cache.retrieve('uuid')
            project_uuid = cache.retrieve('project_uuid')
            return uuid, project_uuid
        return None, None

    def add_project(self, directory: str, make_uuid: bool = False, verbose: bool = False):
        """Add a single project entry to the database."""
        uuid, _ = self._read_uuid_from_cache(directory, make_uuid=make_uuid)
        rel_path = self._abs_to_rel(directory)
        if not self.get_uuid_status(uuid):
            if verbose:
                print(f"Adding project: {rel_path} (UUID: {uuid})")
            self.update_status(uuid, state=State.NEW, message='No information recorded yet.', rel_path=rel_path)

    def add_md(self, directory: str, make_uuid: bool = False, verbose: bool = False):
        """Add a single MD entry to the database."""
        # Get the project UUID from the parent directory in case is not already cached
        project_uuid, _ = self._read_uuid_from_cache(Path(directory).parent.as_posix())
        if not project_uuid:
            raise ValueError(f"Directory '{Path(directory).name}' is not inside a valid project directory")
        # Ensure the parent project exists in the DB before adding MD
        if not self.get_uuid_status(project_uuid):
            raise ValueError(f"Project with UUID '{project_uuid}' does not exist in the database. Add the project first.")
        # Get MD UUID and add the project_uuid
        uuid, _ = self._read_uuid_from_cache(directory, make_uuid, project_uuid)
        rel_path = self._abs_to_rel(directory)
        if not project_uuid:
            raise ValueError(f"Directory '{directory}' does not appear to be an MD")
        if not self.get_uuid_status(uuid, project_uuid):
            if verbose:
                print(f"Adding MD: {rel_path} (UUID: {uuid}, Project UUID: {project_uuid})")
            self.update_status(uuid, state=State.NEW, message='No information recorded yet.', project_uuid=project_uuid, rel_path=rel_path)

    def add_entries(self,
        paths_or_globs: list[str] | str,
        ignore_dirs: list[str] | str = [],
        md_dirs: list[str] | str = [],
        verbose: bool = True,
    ):
        """Add multiple project and MD entries to the database from given paths or glob patterns.

        Args:
            paths_or_globs (list[str]): List of directory paths or glob patterns to search for projects.
            ignore_dirs (list[str]): List of directory paths or glob patterns to ignore.
            md_dirs (list[str]): List of directory paths or glob patterns to identify MDs within each project.
            verbose (bool): Whether to print verbose output.

        """
        # Normalize str inputs to lists
        paths_or_globs = _type_check_dir_list(paths_or_globs)
        ignore_dirs = _type_check_dir_list(ignore_dirs)
        md_dirs = _type_check_dir_list(md_dirs)

        project_directories = _resolve_directory_patterns(paths_or_globs)
        ignore_dirs = _resolve_directory_patterns(ignore_dirs)
        for project_dir in project_directories:
            rel_path = self._abs_to_rel(project_dir)
            if project_dir in ignore_dirs:
                if verbose: print(f"Ignoring project: {rel_path}")
                continue
            if not project_dir.exists():
                if verbose: print(f"Warning: Project directory {rel_path} does not exist. Skipping.")
                continue
            # Add project row (num_mds will be set automatically by triggers when MDs are added)
            self.add_project(project_dir, make_uuid=True, verbose=verbose)
            # Count MDs in this project
            project_md_dirs = [md for md in _resolve_directory_patterns(md_dirs, root_path=project_dir) if md.is_dir()]
            # Add MDs (replicas/subprojects) rows (triggers will auto-increment num_mds)
            for md_dir_path in project_md_dirs:
                self.add_md(md_dir_path, make_uuid=True, verbose=verbose)

    def get_uuid_status(self, uuid: str, project_uuid: str = None):
        """Retrieve a project or MD's status from the database by UUID."""
        status = None
        with self.locked_storage_file:
            cur = self.conn.cursor()
            if project_uuid:
                # This is an MD entry
                cur.execute("SELECT * FROM mds WHERE uuid=?", (uuid,))
                row = cur.fetchone()
                if row:
                    columns = [desc[0] for desc in cur.description]
                    status = dict(zip(columns, row))
                    status['scope'] = 'MD'
            else:
                # This is a project entry
                cur.execute("SELECT * FROM projects WHERE uuid=?", (uuid,))
                row = cur.fetchone()
                if row:
                    columns = [desc[0] for desc in cur.description]
                    status = dict(zip(columns, row))
                    status['scope'] = 'Project'

        if status and status.get('state') == State.RUNNING.value:
            active_slurm_jobs = self._check_active_slurm_jobs([status.get('slurm_job_id')])
            is_active = self._is_running_entry_active(
                status.get('process_id'),
                status.get('slurm_job_id'),
                active_slurm_jobs,
            )
            if is_active is False:
                self.update_status(
                    uuid=uuid,
                    project_uuid=project_uuid,
                    state=State.ERROR,
                    message=self._stale_running_message(status.get('process_id'), status.get('slurm_job_id')),
                )
                return self.get_uuid_status(uuid, project_uuid)
        return status

    def get_status(self, directory: str | Path):
        """Retrieve a project or MD's status from the database by directory path."""
        if isinstance(directory, Path):
            directory = directory.resolve().as_posix()
        uuid, project_uuid = self._read_uuid_from_cache(directory)
        return self.get_uuid_status(uuid, project_uuid)

    def update_status(self,
        uuid: str,
        state: State | str,
        message: str,
        project_uuid: str = None,
        rel_path: str = None
    ):
        """Update or insert (upsert) a project or MD's status in the database.

        Args:
            uuid: UUID of the project or MD
            state: State enum or string
            message: Status message
            project_uuid: If provided, this is an MD entry (requires rel_path)
            rel_path: Relative path to the directory (required for new entries)
            overwrite: Whether to replace existing entries with the same rel_path.

        """
        with self.locked_storage_file:
            cur = self.conn.cursor()
            if isinstance(state, State):
                state = state.value
            last_modified = time.strftime(self.date_format, time.localtime())

            if state == State.RUNNING.value:
                process_id = os.getpid()
                slurm_job_id = os.environ.get('SLURM_JOB_ID')
            else:
                process_id = None
                slurm_job_id = None

            if project_uuid:
                # This is an MD entry
                if rel_path is None:
                    # Try to get existing rel_path (lock already held, call internal)
                    cur.execute("SELECT rel_path FROM mds WHERE uuid=?", (uuid,))
                    row = cur.fetchone()
                    if row:
                        rel_path = row[0]
                    else:
                        raise ValueError("rel_path is required for new MD entries")

                try:
                    cur.execute('''
                        INSERT INTO mds (uuid, project_uuid, rel_path, state, message, process_id, slurm_job_id, last_modified)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                        ON CONFLICT(uuid) DO UPDATE SET
                            state=excluded.state,
                            message=excluded.message,
                            process_id=excluded.process_id,
                            slurm_job_id=excluded.slurm_job_id,
                            last_modified=excluded.last_modified,
                            rel_path=excluded.rel_path
                    ''', (uuid, project_uuid, rel_path, state, message, process_id, slurm_job_id, last_modified))
                except sqlite3.IntegrityError as e:
                    if 'UNIQUE constraint failed: mds.rel_path' in str(e):
                        # In case of rel_path collision, remove the existing entry and insert the new one
                        # This is the same as saying we maintain the rel_path and updating the uuid
                        cur.execute("SELECT uuid FROM mds WHERE rel_path=?", (rel_path,))
                        row = cur.fetchone()
                        if row:
                            existing_uuid = row[0]
                            cur.execute("DELETE FROM mds WHERE uuid=?", (existing_uuid,))
                            cur.execute('''
                                INSERT INTO mds (uuid, project_uuid, rel_path, state, message, process_id, slurm_job_id, last_modified)
                                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                            ''', (uuid, project_uuid, rel_path, state, message, process_id, slurm_job_id, last_modified))
                    else:
                        raise
            else:
                # This is a project entry
                if rel_path is None:
                    # Try to get existing rel_path (lock already held, call internal)
                    cur.execute("SELECT rel_path FROM projects WHERE uuid=?", (uuid,))
                    row = cur.fetchone()
                    if row:
                        rel_path = row[0]
                    else:
                        raise ValueError("rel_path is required for new project entries")
                try:
                    cur.execute('''
                        INSERT INTO projects (uuid, rel_path, state, message, process_id, slurm_job_id, last_modified, num_mds)
                        VALUES (?, ?, ?, ?, ?, ?, ?, 0)
                        ON CONFLICT(uuid) DO UPDATE SET
                            state=excluded.state,
                            message=excluded.message,
                            process_id=excluded.process_id,
                            slurm_job_id=excluded.slurm_job_id,
                            last_modified=excluded.last_modified,
                            rel_path=excluded.rel_path
                    ''', (uuid, rel_path, state, message, process_id, slurm_job_id, last_modified))
                except sqlite3.IntegrityError as e:
                    # Handling rel_path collision like for MDs
                    if 'UNIQUE constraint failed: projects.rel_path' in str(e):
                        # Get the UUID of the existing project with this rel_path
                        cur.execute("SELECT uuid FROM projects WHERE rel_path=?", (rel_path,))
                        row = cur.fetchone()
                        if row:
                            existing_uuid = row[0]
                            # Replacing project automatically removes linked MD rows via ON DELETE CASCADE.
                            cur.execute("DELETE FROM projects WHERE uuid=?", (existing_uuid,))
                            cur.execute('''
                                INSERT INTO projects (uuid, rel_path, state, message, process_id, slurm_job_id, last_modified, num_mds)
                                VALUES (?, ?, ?, ?, ?, ?, ?, 0)
                            ''', (uuid, rel_path, state, message, process_id, slurm_job_id, last_modified))

                    else:
                        raise

            self.conn.commit()

    def remove_entry(self, query_path: str | list[str], query_scope: str = "projects", verbose: bool = True):
        """Remove a single project or MD entry from the database by directory."""
        query_path = _type_check_dir_list(query_path)
        entries = self.query_table(query_scope, query_path=query_path)
        for entry in entries:
            uuid = entry[0]  # uuid is the first column in projects
            project_uuid = entry[1] if query_scope == 'mds' else None  # project_uuid is the second column in mds
            self.remove_entry_by_uuid(uuid, project_uuid=project_uuid, verbose=verbose)

    def remove_entry_by_uuid(self, uuid: str, project_uuid: str = None, verbose: bool = False):
        """Remove a single project or MD entry from the database by UUID."""
        with self.locked_storage_file:
            cur = self.conn.cursor()
            if project_uuid:
                # Remove from MDs
                cur.execute("DELETE FROM mds WHERE uuid=?", (uuid,))
            else:
                # Remove from projects
                cur.execute("DELETE FROM projects WHERE uuid=?", (uuid,))
            if verbose:
                if cur.rowcount > 0:
                    print(f"Deleted {'MD' if project_uuid else 'project'} with UUID '{uuid}'")
                else:
                    print(f"No entry found for UUID '{uuid}' to delete")
            self.conn.commit()

    def generate_inputs_yaml(self,
            inputs_template_path: str,
            inputs_generator: Optional[Callable | str] = None,
            overwrite: bool = False,
            inputs_filename: str = 'inputs.yaml',
            query_path: list[str] = ['*'],
            query_state: list[str] = [],
            query_message: list[str] = [],
        ):
        """Generate an inputs.yaml file for the project directory based on a Jinja2 template.

        Args:
            inputs_template_path (str):
                Path to the inputs Jinja2 template file to be used for generating the inputs files.
            inputs_generator (callable):
                A callable function intended for generating input values or a string path
                to a generator file with a inputs_generator(project_dir) function.
                Accepts the project directory name (DIR) as an argument and returns a
                dictionary of the key-value pairs to be used in the template.
            overwrite (bool):
                Whether to overwrite existing inputs.yaml files.
            inputs_filename (str):
                The name of the inputs YAML file to be generated.
            query_path (list[str]):
                List of directory glob patterns to filter the projects from the dataset.
            query_state (list[str]):
                List of states to filter the projects from the dataset.
            query_message (list[str]):
                List of glob patterns to filter the projects by message field (e.g., 'URLError*').


        """
        # Query projects directly from database with filters
        query_path = _type_check_dir_list(query_path)
        query_state = _type_check_dir_list(query_state)
        query_message = _type_check_dir_list(query_message)
        filtered_projects = self.query_table('projects', query_path, query_state, query_message)
        if len(filtered_projects) == 0:
            print("No projects found matching the specified query_path and query_state filters.")
            return

        # Load the template
        with open(inputs_template_path, 'r') as f:
            template_str = f.read()

        # Generate input values using the provided generator function if any
        if type(inputs_generator) is str:
            print(f"Loading inputs generator from file: {inputs_generator}")
            # Load the generator module
            spec = importlib.util.spec_from_file_location("generator", inputs_generator)
            generator_module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(generator_module)
            # Call the inputs_generator function
            if not hasattr(generator_module, 'inputs_generator'):
                raise ValueError(f"Generator file '{inputs_generator}' must define a inputs_generator(project_dir) function")

            inputs_generator = generator_module.inputs_generator

        # Generate inputs.yaml for each project directory
        for project_entry in filtered_projects:
            project_dict = dict(zip(self.project_columns, project_entry))
            project_dir = self._rel_to_abs(project_dict['rel_path'])

            template = jinja2.Template(template_str)
            inputs_yaml_path = os.path.join(project_dir, inputs_filename)
            if os.path.exists(inputs_yaml_path) and not overwrite:
                print(f"{inputs_yaml_path} already exists and overwrite is set to False. Skipping.")
                continue
            generated_input = {}
            if inputs_generator:
                # If generator returned None, use empty dict
                generated_input = inputs_generator(project_dir) or {}
            # Render the template with project defaults
            rendered_yaml = template.render(DIR=Path(project_dir).name,
                                            DATASET=os.path.relpath(self.dataset_path, start=project_dir),
                                            project_status=project_dict,
                                            **generated_input)
            # Write the rendered YAML to inputs.yaml
            print(f"Generating {inputs_yaml_path} for project {project_dict['rel_path']}")
            with open(inputs_yaml_path, 'w') as f:
                f.write(rendered_yaml)

    def scan(self,
        root_dir: str = None,
        ignore_dirs: list[str] = [],
        verbose: bool = False,
    ):
        """Scan directory tree for projects and MDs by finding CACHE_FILENAME files.

        This function walks the directory tree starting from root_dir (or self.root_path),
        looks for CACHE_FILENAME files, reads their UUIDs and project_uuids, and adds entries to the database.
        Projects are identified by the absence of project_uuid in their cache.
        MDs are identified by the presence of project_uuid in their cache.

        Args:
            root_dir: Root directory to scan (defaults to self.root_path)
            ignore_dirs: List of directory names to ignore
            verbose: Whether to print verbose output

        """
        root_path = Path(root_dir) if root_dir else self.root_path
        root_dir = root_path.resolve()
        ignore_dirs_set = set(ignore_dirs)

        # Separate projects and MDs based on presence of project_uuid in cache
        projects = {}  # {rel_path: uuid}
        mds = {}  # {rel_path: (uuid, project_uuid)}

        # Walk directory tree and find all cache files
        for dirpath, dirnames, filenames in os.walk(root_dir):
            # Skip ignored directories
            dirnames[:] = [d for d in dirnames if d not in ignore_dirs_set]

            if CACHE_FILENAME in filenames:
                dir_path = Path(dirpath)
                uuid, project_uuid = self._read_uuid_from_cache(dir_path)

                if not uuid:
                    if verbose:
                        print(f"Warning: No UUID found in {dir_path}. Skipping.")
                    continue

                rel_path = self._abs_to_rel(dir_path)
                if project_uuid:
                    # This is an MD (has project_uuid in cache)
                    mds[rel_path] = (uuid, project_uuid)
                else:
                    # This is a Project (no project_uuid in cache)
                    projects[rel_path] = uuid

        # Add all projects first
        for rel_path, uuid in projects.items():
            status = self.get_uuid_status(uuid)
            if not status:
                if verbose:
                    print(f"Adding project: {rel_path} (UUID: {uuid})")
                self.update_status(uuid, state=State.NEW, message='No information recorded yet.', rel_path=rel_path)
            elif status['rel_path'] != rel_path:
                if verbose:
                    print(f"Updating project path: {rel_path} (UUID: {uuid}) from {status['rel_path']} to {rel_path}")
                self.update_status(uuid, state=status['state'], message=status['message'], rel_path=rel_path)
            elif verbose:
                print(f"Project already exists: {rel_path} (UUID: {uuid})")

        # Then add all MDs
        for rel_path, (uuid, project_uuid) in mds.items():
            if not self.get_uuid_status(uuid, project_uuid):
                if verbose:
                    print(f"  Adding MD: {rel_path} (UUID: {uuid}, Project UUID: {project_uuid})")
                self.update_status(uuid, state=State.NEW, message='No information recorded yet.', project_uuid=project_uuid, rel_path=rel_path)
            elif verbose:
                print(f"  MD already exists: {rel_path} (UUID: {uuid})")

    def _build_where_clause(self, query_path: list[str] = None, query_state: list[str] = None, query_message: list[str] = None) -> tuple[str, list]:
        """Build a SQL WHERE clause from query parameters.

        Args:
            query_path: List of glob patterns to filter rel_path by (uses GLOB with OR)
            query_state: List of states to filter by (uses IN clause)
            query_message: List of glob patterns to filter message by (uses GLOB with OR)

        Returns:
            Tuple of (where_clause_string, parameters_list)

        """
        conditions = []
        params = []

        if query_state:
            placeholders = ','.join('?' * len(query_state))
            conditions.append(f"state IN ({placeholders})")
            params.extend(query_state)

        if query_path:
            # Filter out '*' patterns and convert glob patterns to SQL GLOB patterns
            # Python glob uses ** for recursive, SQL GLOB uses * for any chars
            glob_conditions = []
            for pattern in query_path:
                if pattern and pattern != '*':
                    sql_pattern = pattern.replace('**/', '*').replace('**', '*')
                    glob_conditions.append("rel_path GLOB ?")
                    params.append(sql_pattern)
            if glob_conditions:
                conditions.append(f"({' OR '.join(glob_conditions)})")

        if query_message:
            # Use GLOB pattern matching for message field
            # Supports wildcards: * (any chars), ? (single char)
            message_conditions = []
            for pattern in query_message:
                if pattern and pattern != '*':
                    message_conditions.append("message GLOB ?")
                    params.append(pattern)
            if message_conditions:
                conditions.append(f"({' OR '.join(message_conditions)})")

        where_clause = " AND ".join(conditions) if conditions else ""
        return where_clause, params

    def query_table(self, table: str, query_path: list[str] = None, query_state: list[str] = None, query_message: list[str] = None) -> list[tuple]:
        """Query an available table with optional filters.

        Args:
            table: Table name to query ('projects' or 'mds')
            query_path: List of glob patterns to filter rel_path by
            query_state: List of states to filter by
            query_message: List of glob patterns to filter message by (e.g., 'URLError*')

        Returns:
            List of tuples with project data

        """
        if table not in ['projects', 'mds']:
            raise ValueError("Table must be either 'projects' or 'mds'")
        with self.locked_storage_file:
            cur = self.conn.cursor()
            where_clause, params = self._build_where_clause(query_path, query_state, query_message)
            cols = self.project_columns if table == 'projects' else self.md_columns
            query = f"SELECT {', '.join(cols)} FROM {table}"
            if where_clause:
                query += f" WHERE {where_clause}"
            query += " ORDER BY rel_path"
            cur.execute(query, params)
            return cur.fetchall()

    @property
    def projects_table(self) -> list[tuple]:
        """Retrieve all project status from the SQLite database as a list of tuples."""
        return self.query_table('projects')

    @property
    def mds_table(self) -> list[tuple]:
        """Retrieve all MD status from the SQLite database as a list of tuples."""
        return self.query_table('mds')

    def get_dataframe(self,
            uuid_length=None,
            root_path=None,
            sort_by='last_modified',
            asc=False,
            include_logs: bool = False,
            query_path: list[str] | str = [],
            query_state: list[str] | str = [],
            query_message: list[str] | str = [],
            query_scope: str = None,
    ) -> 'pd.DataFrame':
        """Retrieve a joined view of projects and MDs as a single DataFrame
        with empty values for the not matching columns.
        Adds a 'scope' column to indicate if the row is from a project or an MD.

        Args:
            uuid_length (int, optional): If provided, truncates 'uuid' and 'project_uuid' columns to this length for display.
            root_path (str, optional): If provided, show the absolute paths relative to this root path.
            sort_by (str, optional): Column name to sort the final DataFrame by. Defaults to 'rel_path'.
            asc (bool, optional): Whether to sort in ascending order. Defaults to True.
            include_logs (bool, optional): If True, adds 'log_file' and 'err_file' columns with HTML links to the latest log files.
            query_path (list[str] | str): If provided, filters rows whose 'rel_path' matches these glob patterns.
            query_state (list[str] | str): If provided, filters rows whose 'state' matches this value/list of values.
            query_message (list[str] | str): If provided, filters rows whose 'message' matches these glob patterns (e.g., 'URLError*').
            query_scope (str, optional): If provided, filters rows whose 'scope' matches this value ('project'/'p' or 'md'/'m').

        """
        query_path = _type_check_dir_list(query_path)
        query_state = _type_check_dir_list(query_state)
        query_message = _type_check_dir_list(query_message)

        if query_scope:
            query_scope = query_scope.lower()
            if query_scope in ('p', 'project', 'projects'):
                scopes_to_query = ['projects']
            elif query_scope in ('m', 'md', 'mds'):
                scopes_to_query = ['mds']
            else:
                raise ValueError("query_scope must be either 'project'/'p' or 'md'/'m'")
        else:
            scopes_to_query = ['projects', 'mds']

        # Keep persisted status in sync with real process state before displaying.
        self._refresh_running_entries(scopes=scopes_to_query)

        # Query data directly from SQLite with filters applied
        dataframes = []
        for scope in scopes_to_query:
            cols = self.project_columns if scope == 'projects' else self.md_columns
            data = self.query_table(scope, query_path, query_state, query_message)
            df = pd.DataFrame(data, columns=cols)
            df['scope'] = scope
            dataframes.append(df)

        df_joined = pd.concat(dataframes, ignore_index=True) if len(dataframes) > 1 else dataframes[0]

        # Sort the DataFrame
        if sort_by == 'last_modified':
            df_joined['last_modified'] = pd.to_datetime(df_joined['last_modified'], format=self.date_format, errors='coerce')
        df_joined = df_joined.sort_values([sort_by], ascending=asc, na_position='first')
        if sort_by == 'last_modified':
            df_joined['last_modified'] = df_joined['last_modified'].dt.strftime(self.date_format)

        # Replace NaN with empty string for display
        df_joined = df_joined.fillna('')
        df_joined.set_index('uuid', inplace=True)
        # Convert num_mds from float to int if present
        if 'num_mds' in df_joined.columns:
            df_joined['num_mds'] = df_joined['num_mds'].apply(lambda x: int(x) if isinstance(x, float) and not pd.isna(x) else x)
        # Optionally make rel_path relative to root_path
        # if root_path is not None:
        #     root_path = os.path.abspath(root_path)
        #     df_joined['rel_path'] = df_joined['rel_path'].apply(
        #         lambda x: os.path.relpath(Path(x).resolve(), root_path)
        #     )

        # Optionally truncate uuid and project_uuid columns for display
        if uuid_length is not None and uuid_length > 0:
            # Truncate index (uuid)
            df_joined.index = df_joined.index.map(lambda x: x[:uuid_length] if isinstance(x, str) else x)
            # Truncate project_uuid column if present
            if 'project_uuid' in df_joined.columns:
                df_joined['project_uuid'] = df_joined['project_uuid'].apply(
                    lambda x: x[:uuid_length] if isinstance(x, str) and x else x
                )

        # Optionally add log file columns with HTML links
        if include_logs:
            def make_log_link(log_path: str | None, link_text: str) -> str:
                """Create an HTML link for a log file."""
                if log_path and Path(log_path).exists():
                    file_url = f"file://{log_path}"
                    return f'<a href="{file_url}" target="_blank">{link_text}</a>'
                return ''

            log_files = []
            err_files = []
            for idx, row in df_joined.iterrows():
                rel_path = row['rel_path']
                latest_out, latest_err = self._get_latest_log_files(rel_path)

                out_name = Path(latest_out).name if latest_out else ''
                err_name = Path(latest_err).name if latest_err else ''

                log_files.append(make_log_link(latest_out, out_name))
                err_files.append(make_log_link(latest_err, err_name))

            df_joined['log_file'] = log_files
            df_joined['err_file'] = err_files

        # Remove internal runtime tracking columns from display (kept in DB for reconciliation)
        cols_to_drop = ['process_id', 'slurm_job_id']
        df_joined = df_joined.drop(columns=[c for c in cols_to_drop if c in df_joined.columns])

        # Change the order of the columns to have 'scope' first
        cols = df_joined.columns.tolist()
        if not query_scope and 'project_uuid' in cols:
            cols.insert(0, cols.pop(cols.index('project_uuid')))
            cols.insert(1, cols.pop(cols.index('scope')))
        df_joined = df_joined[cols]
        return df_joined

    @property
    def dataframe(self) -> 'pd.DataFrame':
        """Retrieve the joined DataFrame view of projects and MDs with log file links."""
        return self.get_dataframe(uuid_length=8)

    def summary(self) -> pd.DataFrame:
        """Return a summary DataFrame with the count of projects in each state."""
        df = self.get_dataframe(
                    query_scope='project',
        )
        # Group by state and count
        summary = df.groupby("state").size().reset_index(name="count")
        # Sort by state order:
        state_order = {state: index for index, state in enumerate(State)}
        summary['state_order'] = summary['state'].apply(lambda x: state_order.get(State(x), -1))
        summary = summary.sort_values('state_order').drop(columns=['state_order']).reset_index(drop=True)
        return summary

    def error_summary(self) -> pd.DataFrame:
        """Return a summary of error messages for projects."""
        df = self.dataframe
        return (df.loc[(df['state'] == 'error') & (df['scope'] == 'projects')]
                .groupby('message').size()
                .reset_index(name='count')
                .sort_values('count', ascending=False)
                .reset_index(drop=True)
                )

    def _get_latest_log_files(self, rel_path: str) -> tuple[str | None, str | None]:
        """Find the most recent .out and .err log files in the rel_path/logs directory."""
        logs_dir = Path(self._rel_to_abs(rel_path)) / 'logs'
        if not logs_dir.exists():
            return None, None

        out_files = sorted(logs_dir.glob('*.out'), key=lambda f: f.stat().st_mtime, reverse=True)
        err_files = sorted(logs_dir.glob('*.err'), key=lambda f: f.stat().st_mtime, reverse=True)

        latest_out = str(out_files[0]) if out_files else None
        latest_err = str(err_files[0]) if err_files else None

        return latest_out, latest_err

    def display(self, **kwargs):
        """Display the dataframe with clickable log links in Jupyter."""
        from IPython.display import HTML, display
        # Ensure include_logs is True for display
        kwargs.setdefault('include_logs', True)
        kwargs.setdefault('uuid_length', 8)
        df = self.get_dataframe(**kwargs)
        html = df.to_html(escape=False)
        display(HTML(html))

    def _repr_html_(self) -> str:
        """Return HTML representation for automatic Jupyter notebook display."""
        return self.dataframe.to_html(escape=False)

    def launch(self,
        cmd: str = None,
        n_jobs: int = 0,
        pool_size: int = 1,
        slurm: bool = False,
        job_template: str = None,
        query_path: list[str] | str = [],
        query_state: list[str] | str = [],
        query_message: list[str] | str = [],
        debug: bool = False
    ):
        """Launch a command for each project directory in the dataset.

        Args:
            cmd (str):
                Custom command to run the workflow or any other command/script.
            n_jobs (int):
                Number of jobs to launch. If 0, all jobs are launched.
            pool_size (int):
                Number of parallel jobs to run using subprocess. Default is 1.
            slurm (bool):
                Whether to submit the workflow to SLURM.
            job_template (str):
                Path to the bash script or SLURM job template file. You can use Jinja2
                templating to customize the job script using the fields of
                the input YAML, the project directory, and whether SLURM is used.
            query_path (list[str] | str):
                If provided, filters rows whose 'rel_path' matches these glob patterns.
            query_state (list[str] | str):
                If provided, filters rows whose 'state' matches this value/list of values.
            query_message (list[str] | str):
                If provided, filters rows whose 'message' matches these glob patterns (e.g., 'URLError*').
            debug (bool):
                Only print the commands without executing them.

        """
        # Validation: ensure proper parameters are provided
        if slurm and not job_template:
            raise ValueError("job_template must be provided when slurm is True")
        cmd = 'mwf run' if cmd is None else cmd

        # Load job template if provided
        template_str = None
        if job_template:
            with open(job_template, 'r') as f:
                template_str = f.read()

        # Determine pool size for parallel execution
        pool_size = min(pool_size, os.cpu_count())
        if pool_size <= 0 and not slurm:
            pool_size = os.cpu_count()
        parallel = pool_size > 1 and not slurm
        jobs_to_run = []

        # Query projects directly from database with filters
        query_path = _type_check_dir_list(query_path)
        query_state = _type_check_dir_list(query_state)
        query_message = _type_check_dir_list(query_message)
        filtered_projects = self.query_table('projects', query_path, query_state, query_message)

        n = 0
        if debug:
            print(f"DEBUG: Found {len(filtered_projects)} projects matching the query filters.")
            if len(filtered_projects) > 0 and slurm:
                print("DEBUG: Job will not be submitted to SLURM in debug mode, but the job script will be generated.")
        for project_entry in filtered_projects:
            project_dict = dict(zip(self.project_columns, project_entry))
            rel_path = project_dict['rel_path']
            project_dir = self._rel_to_abs(rel_path)

            n += 1
            if n_jobs > 0 and n > n_jobs:
                break
            # Generate job script from template if provided
            if template_str:
                inputs_yaml_path = os.path.join(project_dir, 'inputs.yaml')
                inputs_config = load_yaml(inputs_yaml_path) if os.path.exists(inputs_yaml_path) else {}
                template = jinja2.Template(template_str)
                rendered_script = template.render(**inputs_config,
                                                project_status=project_dict,
                                                DIR=Path(project_dir).name,
                                                PATH=Path(project_dir),
                                                slurm=slurm)

                job_script_path = os.path.join(project_dir, 'mwf_slurm_job.sh')
                with open(job_script_path, 'w') as f:
                    f.write(rendered_script)
                os.chmod(job_script_path, 0o755)
            command_to_run = cmd if cmd else f"bash {os.path.basename(job_script_path)}"
            # Create logs directory if it doesn't exist
            log_dir = os.path.join(project_dir, 'logs')
            os.makedirs(log_dir, exist_ok=True)
            log_file = os.path.join(log_dir, f'mwf_{int(time.time())}.out')
            err_file = os.path.join(log_dir, f'mwf_{int(time.time())}.err')
            # Launch workflow
            if slurm:
                # SLURM execution
                if debug:
                    print(f'Submitting SLURM job script: {job_script_path}')
                    # print(f"In directory {project_dir}")
                    # print("sbatch --output=logs/mwf_%j.out --error=logs/mwf_%j.err mwf_slurm_job.sh ")
                else:
                    print(f"Submitting SLURM job for {project_dir}")
                    subprocess.run(['sbatch',
                                    '--output=logs/mwf_%j.out',
                                    '--error=logs/mwf_%j.err',
                                    job_script_path],
                                cwd=project_dir)
            elif parallel:
                # Parallel execution - collect jobs to run later
                if debug:
                    print(f"[parallel] cd {project_dir}")
                    print(f"[parallel] bash {job_script_path} > {log_file} 2> {err_file}")
                else:
                    jobs_to_run.append({
                        'project_dir': project_dir,
                        'rel_path': rel_path,
                        'command_to_run': command_to_run,
                        'log_file': log_file,
                        'err_file': err_file
                    })
            else:
                # Normal Python execution (sequential)
                if debug:
                    print(f"cd {project_dir}")
                    print(f"{command_to_run}")
                    continue
                print(f"Running job for dataset entry {rel_path}")
                self._run_sequential_job(command_to_run, project_dir, log_file, err_file)

        # Execute parallel jobs if any
        if parallel and jobs_to_run and not debug:
            self._run_parallel_jobs(jobs_to_run, max_concurrent=pool_size)

    def _run_sequential_job(self, command_to_run: str, project_dir: str, log_file: str, err_file: str):
        """Run a single job sequentially with real-time output to both terminal and log files."""
        proc = subprocess.Popen(
            command_to_run,
            cwd=project_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=True,
            text=True,
            bufsize=1
        )

        with open(log_file, 'w') as out_f, open(err_file, 'w') as err_f:
            # Create threads to handle stdout and stderr
            stdout_thread = threading.Thread(
                target=_stream_output,
                args=(proc.stdout, out_f, sys.stdout)
            )
            stderr_thread = threading.Thread(
                target=_stream_output,
                args=(proc.stderr, err_f, sys.stderr)
            )
            stdout_thread.start()
            stderr_thread.start()
            # Wait for threads to complete
            stdout_thread.join()
            stderr_thread.join()

        proc.wait()

    def _run_parallel_jobs(self, jobs: list[dict], max_concurrent: int):
        """Run jobs in parallel and report progress.

        Args:
            jobs (list[dict]): List of job dictionaries with keys:
                - project_dir: Full path to project directory
                - rel_path: Relative path for display
                - command_to_run: Command to run for the job
                - log_file: Path to the log file
                - err_file: Path to the error file
            max_concurrent (int): Maximum number of concurrent jobs.

        """
        print(f"Launching {len(jobs)} jobs in parallel (max {max_concurrent} concurrent)...")

        running_processes = {}
        pending_jobs = list(jobs)
        completed = 0
        total = len(jobs)
        open_files = {}  # Track open file handles (tuple of log_f, err_f)

        while pending_jobs or running_processes:
            # Start new jobs if we have capacity
            while pending_jobs and len(running_processes) < max_concurrent:
                job = pending_jobs.pop(0)
                print(f"Starting job for {job['rel_path']} ({len(running_processes)+1}/{total})")
                log_f = open(job['log_file'], 'w')
                err_f = open(job['err_file'], 'w')
                proc = subprocess.Popen(
                    job['command_to_run'],
                    cwd=job['project_dir'],
                    stdout=log_f,
                    stderr=err_f,
                    shell=True
                )
                running_processes[proc] = job
                open_files[proc] = (log_f, err_f)

            # Check for completed jobs
            for proc in list(running_processes.keys()):
                retcode = proc.poll()
                if retcode is not None:
                    job = running_processes.pop(proc)
                    log_f, err_f = open_files.pop(proc)
                    log_f.close()
                    err_f.close()
                    completed += 1
                    status = "DONE" if retcode == 0 else f"FAILED (exit code {retcode})"
                    print(f"[{completed}/{total}] {job['rel_path']}: {status}")

            if running_processes:
                time.sleep(1)  # Avoid busy waiting

        print(f"All {total} jobs completed.")

    def recount_all_mds(self, verbose: bool = False):
        """Count all MDs for each project and update the num_mds field in the projects table."""
        with self.locked_storage_file:
            cur = self.conn.cursor()
            # Get all project UUIDs
            cur.execute("SELECT uuid FROM projects")
            project_uuids = [row[0] for row in cur.fetchall()]
            for project_uuid in project_uuids:
                # Count MDs for this project
                cur.execute("SELECT COUNT(*) FROM mds WHERE project_uuid=?", (project_uuid,))
                count = cur.fetchone()[0]
                # Update num_mds
                cur.execute("UPDATE projects SET num_mds=? WHERE uuid=?", (count, project_uuid))
                if verbose:
                    print(f"Project {project_uuid}: num_mds set to {count}")
            self.conn.commit()


def _resolve_directory_patterns(dir_patterns: list[str], root_path: Path = Path('.')) -> list[Path]:
    """Resolve directory patterns (glob or absolute/relative paths)."""
    directories = []
    for dir_pattern in dir_patterns:
        if is_glob(dir_pattern):
            if Path(dir_pattern).is_absolute():
                matched_dirs = [Path(p) for p in glob.glob(dir_pattern)]
            else:
                matched_dirs = list(root_path.glob(dir_pattern))
            directories.extend([p.absolute() for p in matched_dirs if p.is_dir()])
        elif not Path(dir_pattern).is_dir():
            continue
        else:
            p = Path(dir_pattern)
            p = p.absolute() if not p.is_absolute() else p
            directories.append(p)
    return directories


def _type_check_dir_list(dir_list: list[str] | str) -> list[str]:
    """Ensure the input is a list of strings."""
    if isinstance(dir_list, str):
        dir_list = [dir_list]
    return dir_list
