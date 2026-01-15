from mddb_workflow.utils.auxiliar import load_yaml, is_glob
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
from pathlib import Path
from enum import Enum
import importlib.util


class State(Enum):
    """Enumeration of possible workflow states."""
    NEW = 'new'
    RUNNING = 'running'
    ERROR = 'error'
    DONE = 'done'


class Dataset:
    """Class to manage and process a dataset of MDDB projects and their MDs (replicas/subprojects) using a central SQLite database."""
    # TODO: log_file, err_file, group_id
    common_columns = ['state', 'message', 'last_modified']
    project_columns = ['uuid', 'rel_path', 'num_mds'] + common_columns
    md_columns = ['uuid', 'project_uuid', 'rel_path'] + common_columns
    date_format = "%H:%M:%S %d-%m-%Y"

    def __init__(self, dataset_path: str):
        """Initialize the Dataset object and connect to SQLite DB.

        Args:
            dataset_path (str): Path to the dataset storage file, normally an .db file.

        """
        self.dataset_path = Path(dataset_path).resolve()
        # if not self.dataset_path.exists():
        #     raise ValueError(f"Dataset path {self.dataset_path} does not exist.")
        self.root_path = self.dataset_path.parent
        self.conn = sqlite3.connect(dataset_path, check_same_thread=False)
        # Enable foreign key constraints (disabled by default in SQLite)
        self.conn.execute("PRAGMA foreign_keys = ON")
        self._ensure_tables()

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
        cur = self.conn.cursor()
        # Create projects table
        cur.execute('''
            CREATE TABLE IF NOT EXISTS projects (
                uuid TEXT PRIMARY KEY NOT NULL,
                rel_path TEXT UNIQUE NOT NULL,
                num_mds INTEGER DEFAULT 0,
                state TEXT,
                message TEXT,
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
                last_modified TEXT,
                FOREIGN KEY (project_uuid) REFERENCES projects(uuid) ON DELETE CASCADE
            )
        ''')

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

    def _read_uuid_from_cache(self, directory: str, make_uuid: bool = False, project_uuid: str = None) -> tuple[str, str | None]:
        """Read UUID and project_uuid from cache file in the given directory."""
        directory = Path(directory).resolve()
        if not directory.exists():
            raise ValueError(f"Directory {directory} does not exist.")
        cache_file = directory / CACHE_FILENAME
        if not cache_file.exists() and not make_uuid:
            return None, None
        cache = Cache(File(str(cache_file)), project_uuid=project_uuid)
        uuid = cache.retrieve('uuid')
        project_uuid = cache.retrieve('project_uuid')
        if not uuid:
            raise ValueError(f"No 'uuid' found in cache file '{cache_file}'")
        return uuid, project_uuid

    def add_project(self, directory: str, make_uuid: bool = False, verbose: bool = False):
        """Add a single project mentry to the database."""
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
        cur = self.conn.cursor()
        if project_uuid:
            # This is an MD entry
            cur.execute("SELECT * FROM mds WHERE uuid=?", (uuid,))
            row = cur.fetchone()
            if row:
                columns = [desc[0] for desc in cur.description]
                result = dict(zip(columns, row))
                result['scope'] = 'MD'
                return result
        else:
            # This is a project entry
            cur.execute("SELECT * FROM projects WHERE uuid=?", (uuid,))
            row = cur.fetchone()
            if row:
                columns = [desc[0] for desc in cur.description]
                result = dict(zip(columns, row))
                result['scope'] = 'Project'
                return result

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
        rel_path: str = None,
    ):
        """Update or insert a project or MD's status in the database.

        Args:
            uuid: UUID of the project or MD
            state: State enum or string
            message: Status message
            project_uuid: If provided, this is an MD entry (requires rel_path)
            rel_path: Relative path to the directory (required for new entries)

        """
        cur = self.conn.cursor()
        if isinstance(state, State):
            state = state.value
        last_modified = time.strftime(self.date_format, time.localtime())

        if project_uuid:
            # This is an MD entry
            if rel_path is None:
                # Try to get existing rel_path
                existing = self.get_uuid_status(uuid)
                if existing:
                    rel_path = existing['rel_path']
                else:
                    raise ValueError("rel_path is required for new MD entries")

            cur.execute('''
                INSERT INTO mds (uuid, project_uuid, rel_path, state, message, last_modified)
                VALUES (?, ?, ?, ?, ?, ?)
                ON CONFLICT(uuid) DO UPDATE SET
                    state=excluded.state,
                    message=excluded.message,
                    last_modified=excluded.last_modified,
                    rel_path=excluded.rel_path
            ''', (uuid, project_uuid, rel_path, state, message, last_modified))
        else:
            # This is a project entry
            if rel_path is None:
                # Try to get existing rel_path
                existing = self.get_uuid_status(uuid)
                if existing:
                    rel_path = existing['rel_path']
                else:
                    raise ValueError("rel_path is required for new project entries")

            cur.execute('''
                INSERT INTO projects (uuid, rel_path, state, message, last_modified, num_mds)
                VALUES (?, ?, ?, ?, ?, 0)
                ON CONFLICT(uuid) DO UPDATE SET
                    state=excluded.state,
                    message=excluded.message,
                    last_modified=excluded.last_modified,
                    rel_path=excluded.rel_path
            ''', (uuid, rel_path, state, message, last_modified))
        self.conn.commit()

    def remove_entry(self, directory: str | Path, verbose: bool = False):
        """Remove a single project or MD entry from the database by directory."""
        if isinstance(directory, Path):
            directory = directory.resolve().as_posix()
        uuid, project_uuid = self._read_uuid_from_cache(directory)
        self.remove_entry_by_uuid(uuid, project_uuid, verbose=verbose)

    def remove_entry_by_uuid(self, uuid: str, project_uuid: str = None, verbose: bool = False):
        """Remove a single project or MD entry from the database by UUID."""
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


        """
        # Query projects directly from database with filters
        query_path = _type_check_dir_list(query_path)
        query_state = _type_check_dir_list(query_state)
        filtered_projects = self.query_table('projects', query_path, query_state)

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

    def _build_where_clause(self, query_path: list[str] = None, query_state: list[str] = None) -> tuple[str, list]:
        """Build a SQL WHERE clause from query parameters.

        Args:
            query_path: List of glob patterns to filter rel_path by (uses GLOB with OR)
            query_state: List of states to filter by (uses IN clause)

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

        where_clause = " AND ".join(conditions) if conditions else ""
        return where_clause, params

    def query_table(self, table: str, query_path: list[str] = None, query_state: list[str] = None) -> list[tuple]:
        """Query an available table with optional filters.

        Args:
            table: Table name to query ('projects' or 'mds')
            query_path: List of glob patterns to filter rel_path by
            query_state: List of states to filter by

        Returns:
            List of tuples with project data

        """
        if table not in ['projects', 'mds']:
            raise ValueError("Table must be either 'projects' or 'mds'")
        cur = self.conn.cursor()
        where_clause, params = self._build_where_clause(query_path, query_state)
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
            query_scope (str, optional): If provided, filters rows whose 'scope' matches this value ('project'/'p' or 'md'/'m').

        """
        query_path = _type_check_dir_list(query_path)
        query_state = _type_check_dir_list(query_state)

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

        # Query data directly from SQLite with filters applied
        dataframes = []
        for scope in scopes_to_query:
            cols = self.project_columns if scope == 'projects' else self.md_columns
            data = self.query_table(scope, query_path, query_state)
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
        if root_path is not None:
            root_path = os.path.abspath(root_path)
            df_joined['rel_path'] = df_joined['rel_path'].apply(
                lambda x: os.path.relpath(Path(x).resolve(), root_path)
            )

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

    def launch_workflow(self,
        query_path: list[str] | str = [],
        query_state: list[str] | str = [],
        n_jobs: int = 0,
        pool_size: int = 1,
        slurm: bool = False,
        job_template: str = None,
        debug: bool = False
    ):
        """Launch the workflow for each project directory in the dataset.

        Args:
            query_path (list[str] | str):
                If provided, filters rows whose 'rel_path' matches these glob patterns.
            query_state (list[str] | str):
                If provided, filters rows whose 'state' matches this value/list of values.
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
            debug (bool):
                Only print the commands without executing them.

        """
        if not slurm:
            raise NotImplementedError("Only slurm=True is currently implemented")
        if slurm and not job_template:
            raise ValueError("job_template must be provided when slurm is True")
        else:
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
        filtered_projects = self.query_table('projects', query_path, query_state)

        n = 0
        for project_entry in filtered_projects:
            project_dict = dict(zip(self.project_columns, project_entry))
            rel_path = project_dict['rel_path']
            project_dir = self._rel_to_abs(rel_path)

            n += 1
            if n_jobs > 0 and n > n_jobs:
                break
            inputs_yaml_path = os.path.join(project_dir, 'inputs.yaml')
            inputs_config = load_yaml(inputs_yaml_path) if os.path.exists(inputs_yaml_path) else {}

            template = jinja2.Template(template_str)
            rendered_script = template.render(**inputs_config,
                                              project_status=project_dict,
                                              DIR=Path(project_dir).name,
                                              slurm=slurm)

            job_script_path = os.path.join(project_dir, 'mwf_slurm_job.sh')
            log_dir = os.path.join(project_dir, 'logs')
            os.makedirs(log_dir, exist_ok=True)
            with open(job_script_path, 'w') as f:
                f.write(rendered_script)
            os.chmod(job_script_path, 0o755)
            log_file = os.path.join(log_dir, f'mwf_{int(time.time())}.out')
            err_file = os.path.join(log_dir, f'mwf_{int(time.time())}.err')
            # Launch workflow
            if slurm:
                # SLURM execution
                if debug:
                    print(f"cd {project_dir}")
                    print(f'# Job script: {job_script_path}')
                    print("sbatch --output=logs/mwf_%j.out --error=logs/mwf_%j.err mwf_slurm_job.sh ")
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
                        'job_script_path': job_script_path,
                        'log_file': log_file,
                        'err_file': err_file
                    })
            else:
                # Normal Python execution (sequential)
                if debug:
                    print(f"cd {project_dir}")
                    print(f"bash {project_dir}/{os.path.basename(job_script_path)}")
                    continue
                print(f"Running job for {project_dir}")
                with open(log_file, 'w') as out_f, open(err_file, 'w') as err_f:
                    subprocess.run(job_script_path, cwd=project_dir, stdout=out_f, stderr=err_f)
        # Execute parallel jobs if any
        if parallel and jobs_to_run and not debug:
            self._run_parallel_jobs(jobs_to_run, max_concurrent=pool_size)

    def _run_parallel_jobs(self, jobs: list[dict], max_concurrent: int):
        """Run jobs in parallel and report progress.

        Args:
            jobs (list[dict]): List of job dictionaries with keys:
                - project_dir: Full path to project directory
                - rel_path: Relative path for display
                - job_script_path: Path to the job script
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
                log_f = open(job['log_file'], 'w')
                err_f = open(job['err_file'], 'w')
                proc = subprocess.Popen(
                    job['job_script_path'],
                    cwd=job['project_dir'],
                    stdout=log_f,
                    stderr=err_f,
                    shell=False
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


def _resolve_directory_patterns(dir_patterns: list[str], root_path: Path = Path('.')) -> list[Path]:
    """Resolve directory patterns (glob or absolute/relative paths)."""
    directories = []
    for dir_pattern in dir_patterns:
        if not Path(dir_pattern).is_dir():
            continue
        if is_glob(dir_pattern):
            if Path(dir_pattern).is_absolute():
                matched_dirs = [Path(p) for p in glob.glob(dir_pattern)]
            else:
                matched_dirs = list(root_path.glob(dir_pattern))
            directories.extend([p.absolute() for p in matched_dirs if p.is_dir()])
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
