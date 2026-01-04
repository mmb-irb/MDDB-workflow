from mddb_workflow.utils.auxiliar import load_yaml, is_glob, warn
from mddb_workflow.utils.type_hints import *
import pandas as pd
import subprocess
import jinja2
import time
import os
import sqlite3
from pathlib import Path
from enum import Enum


class State(Enum):
    """Enumeration of possible workflow states."""
    NOT_RUN = 'not_run'
    RUNNING = 'running'
    DONE = 'done'
    ERROR = 'error'


class Dataset:
    """Class to manage and process a dataset of MDDB projects and their MDs (replicas/subprojects) using a central SQLite database."""
    # TODO: log_file, err_file, group_id
    common_columns = ['state', 'message', 'last_modified']
    project_columns = ['rel_path', 'num_mds'] + common_columns
    md_columns = ['project_rel_path', 'md_dir'] + common_columns

    def __init__(self, dataset_path: Optional[str] = None):
        """Initialize the Dataset object and connect to SQLite DB.

        Args:
            dataset_path (str): Path to the root directory of the dataset.

        """
        if dataset_path is None:
            raise ValueError("dataset_path must be provided")
        self.root_path = Path(dataset_path).parent.resolve()
        self.conn = sqlite3.connect(dataset_path, check_same_thread=False)
        # Enable foreign key constraints (disabled by default in SQLite)
        self.conn.execute("PRAGMA foreign_keys = ON")
        self._ensure_tables()

    def _ensure_tables(self):
        """Create tables if they do not exist."""
        cur = self.conn.cursor()
        # Create projects table
        cur.execute('''
            CREATE TABLE IF NOT EXISTS projects (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                rel_path TEXT UNIQUE NOT NULL,
                num_mds INTEGER DEFAULT 0,
                state TEXT,
                message TEXT,
                last_modified TEXT
            )
        ''')
        # Create mds table
        # MDs are linked to projects via project_rel_path foreign key
        # ON DELETE CASCADE ensures MDs are deleted when project is deleted
        cur.execute('''
            CREATE TABLE IF NOT EXISTS mds (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                project_rel_path TEXT NOT NULL,
                md_dir TEXT NOT NULL,
                state TEXT,
                message TEXT,
                last_modified TEXT,
                UNIQUE(project_rel_path, md_dir),
                FOREIGN KEY (project_rel_path) REFERENCES projects(rel_path) ON DELETE CASCADE
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
                WHERE rel_path = NEW.project_rel_path;
            END
        ''')

        # Trigger on DELETE: decrement num_mds when an MD is removed
        cur.execute('''
            CREATE TRIGGER IF NOT EXISTS decrement_num_mds
            AFTER DELETE ON mds
            BEGIN
                UPDATE projects
                SET num_mds = num_mds - 1
                WHERE rel_path = OLD.project_rel_path;
            END
        ''')

        self.conn.commit()

    def _resolve_directory_patterns(self, dir_patterns: list[str], root_path: Path = None) -> list[Path]:
        """Resolve directory patterns (glob or absolute/relative paths)."""
        root_path = self.root_path if root_path is None else root_path
        directories = []
        for dir_pattern in dir_patterns:
            if is_glob(dir_pattern):
                matched_dirs = list(root_path.glob(dir_pattern))
                directories.extend([p.absolute() for p in matched_dirs if p.is_dir()])
            else:
                p = Path(dir_pattern)
                p = p.absolute() if not p.is_absolute() else p
                directories.append(p)
        return directories

    def _type_check_dir_list(self, dir_list: list[str] | str) -> list[str]:
        """Ensure the input is a list of strings."""
        if isinstance(dir_list, str):
            dir_list = [dir_list]
        return dir_list

    def add_entries(self,
        paths_or_globs: list[str] | str,
        ignore_dirs: list[str] | str = [],
        md_dirs: list[str] | str = ["replica*"],
        verbose: bool = False,
    ):
        """Scan all project directories and their MDs (replicas/subprojects) and register them in the database if not present.
        This should be called once after creating the Dataset to ensure all projects and MDs are tracked in the DB.
        """
        # Normalize str inputs to lists
        paths_or_globs = self._type_check_dir_list(paths_or_globs)
        ignore_dirs = self._type_check_dir_list(ignore_dirs)
        md_dirs = self._type_check_dir_list(md_dirs)

        project_directories = self._resolve_directory_patterns(paths_or_globs)
        ignore_dirs = [Path(d).resolve() for d in ignore_dirs]
        for project_dir in project_directories:
            if project_dir in ignore_dirs:
                if verbose: print(f"Ignoring project: {project_dir}")
                continue
            rel_path = project_dir.relative_to(self.root_path).as_posix()
            if not project_dir.exists():
                if verbose:
                    print(f"Warning: Project directory {project_dir} does not exist. Skipping.")
                continue
            # Count MDs in this project
            project_md_dirs = [md for md in self._resolve_directory_patterns(md_dirs, root_path=project_dir) if md.is_dir()]
            num_mds = len(project_md_dirs)
            # Add project row (num_mds will be set automatically by triggers when MDs are added)
            if not self.get_status(rel_path):
                if verbose:
                    print(f"Adding project: {rel_path} (with {num_mds} MDs)")
                self.update_status(rel_path, state=State.NOT_RUN, message='No information have been recorded yet.')
            # Add MDs (replicas/subprojects) rows (triggers will auto-increment num_mds)
            for md_dir_path in project_md_dirs:
                md_dir = md_dir_path.name
                if not self.get_status(rel_path, md_dir=md_dir):
                    if verbose:
                        print(f"  Adding MD: {rel_path}/{md_dir}")
                    self.update_status(rel_path, md_dir=md_dir, state=State.NOT_RUN, message='No information have been recorded yet.')

    def remove_entries(self,
        paths_or_globs: list[str] | str,
        ignore_dirs: list[str] | str = [],
        md_dirs: list[str] | str = None,
        verbose: bool = False,
    ):
        """Remove specified project directories or MDs from the database.
        If md_dirs is provided, only remove that MD (or those MDs) for each project.

        Args:
            paths_or_globs: List of directory paths or glob patterns
            ignore_dirs: List of directories to ignore
            md_dirs: directory, list or glob pattern of MD directories to remove.
            verbose: Whether to print verbose output

        """
        # Normalize str inputs to lists
        paths_or_globs = self._type_check_dir_list(paths_or_globs)
        ignore_dirs = self._type_check_dir_list(ignore_dirs)
        md_dirs = self._type_check_dir_list(md_dirs) if md_dirs else None

        project_directories = self._resolve_directory_patterns(paths_or_globs)
        for project_dir in project_directories:
            if project_dir in ignore_dirs:
                if verbose: print(f"Ignoring project: {project_dir}")
                continue
            rel_path = project_dir.relative_to(self.root_path).as_posix()
            cur = self.conn.cursor()
            if md_dirs:
                project_md_dirs = self._resolve_directory_patterns(md_dirs, root_path=project_dir)
                # Remove specific MDs (triggers will auto-decrement num_mds)
                for md_dir in project_md_dirs:
                    cur.execute("DELETE FROM mds WHERE project_rel_path=? AND md_dir=?", (rel_path, md_dir.name))
                    if verbose and cur.rowcount > 0:
                        print(f"Deleted MD '{md_dir.name}' from '{rel_path}'")
            else:
                # Remove entire project and all its MDs (CASCADE will handle MDs)
                cur.execute("DELETE FROM projects WHERE rel_path=?", (rel_path,))
                if verbose and cur.rowcount > 0:
                    print(f"Deleted project '{rel_path}' and all its MDs")
        self.conn.commit()

    def update_status(self,
        rel_path: str,
        state: State | str,
        message: str,
        md_dir: str = None,
    ):
        """Update or insert a project or MD's status in the database."""
        cur = self.conn.cursor()
        if isinstance(state, State):
            state = state.value
        last_modified = time.strftime("%H:%M:%S %d-%m-%Y", time.localtime())

        if md_dir is not None:
            cur.execute('''
                INSERT INTO mds (project_rel_path, md_dir, state, message, last_modified)
                VALUES (?, ?, ?, ?, ?)
                ON CONFLICT(project_rel_path, md_dir) DO UPDATE SET
                    state=excluded.state,
                    message=excluded.message,
                    last_modified=excluded.last_modified
            ''', (rel_path, md_dir, state, message, last_modified))
        else:
            cur.execute('''
                INSERT INTO projects (rel_path, state, message, last_modified, num_mds)
                VALUES (?, ?, ?, ?, 0)
                ON CONFLICT(rel_path) DO UPDATE SET
                    state=excluded.state,
                    message=excluded.message,
                    last_modified=excluded.last_modified
            ''', (rel_path, state, message, last_modified))
        self.conn.commit()

    def get_status(self, rel_path: str, md_dir: str = None):
        """Retrieve a project's status from the database."""
        cur = self.conn.cursor()
        if md_dir is not None:
            cur.execute("SELECT * FROM mds WHERE project_rel_path=? AND md_dir=?", (rel_path, md_dir))
        else:
            cur.execute("SELECT * FROM projects WHERE rel_path=?", (rel_path,))
        row = cur.fetchone()
        if row:
            columns = [desc[0] for desc in cur.description]
            return dict(zip(columns, row))
        return None

    @property
    def status(self) -> dict:
        """Retrieve all project and MD status from the SQLite database as a dictionary with 'projects' and 'mds' keys."""
        cur = self.conn.cursor()
        # Get all projects
        cur.execute(f"SELECT {', '.join(self.project_columns)} FROM projects ORDER BY rel_path")
        projects = cur.fetchall()
        # Get all MDs
        cur.execute(f"SELECT {', '.join(self.md_columns)} FROM mds ORDER BY project_rel_path, md_dir")
        mds = cur.fetchall()
        return {'projects': projects, 'mds': mds}

    @property
    def dataframes(self) -> dict:
        """Retrieve project and MD status from the SQLite database as pandas DataFrames.

        Returns:
            dict: Dictionary with 'projects' and 'mds' keys, each containing a DataFrame.

        """
        status = self.status
        # Create DataFrames
        df_projects = pd.DataFrame(status['projects'], columns=self.project_columns)
        df_projects.set_index('rel_path', inplace=True)
        # Create MDs DataFrame
        df_mds = pd.DataFrame(status['mds'], columns=self.md_columns)
        df_mds.set_index(['project_rel_path', 'md_dir'], inplace=True)
        return {
            'projects': df_projects,
            'mds': df_mds
        }

    @property
    def dataframe(self) -> 'pd.DataFrame':
        """Retrieve a joined view of projects and MDs as a single DataFrame
        with empty values for the not matching columns.
        Adds a 'scope' column to indicate if the row is from a project or an MD.
        """
        dfs = self.dataframes
        df_projects = dfs['projects'].reset_index()
        df_mds = dfs['mds'].reset_index()

        # Set missing columns for alignment
        df_projects['md_dir'] = ''
        df_projects['scope'] = 'Project'
        df_mds['num_mds'] = ''
        df_mds['scope'] = 'MD'
        df_mds = df_mds.rename(columns={'project_rel_path': 'rel_path'})

        # Concatenate, sort, and set multi-index
        df_joined = pd.concat([df_projects, df_mds], ignore_index=True)
        df_joined = df_joined.sort_values(['rel_path', 'md_dir'], na_position='first')
        df_joined.set_index(['rel_path', 'md_dir'], inplace=True)
        # Change the order of the columns to have 'scope' first
        cols = df_joined.columns.tolist()
        cols.insert(0, cols.pop(cols.index('scope')))
        df_joined = df_joined[cols]
        return df_joined


class OldDataset:
    """Placeholder for the old Dataset class function that have to be migrated or removed."""

    def show_groups(self, cmd=False):
        """Display the groups of projects based on their status messages."""
        if cmd:
            print("Project groups based on status messages:\n")
            status = self.status
            grouped = status.groupby('group')
            for group_id, group_df in grouped:
                print(f"Group {group_id}:")
                print(f"Message: {group_df['message'].iloc[0]}")
                print("Projects:")
                for rel_path in group_df.index:
                    print(f"  - {rel_path}")
                print()
        else:
            grouped = self.status.groupby('group').agg({
                'message': 'first',
                'state': 'count'
            }).rename(columns={'state': 'count'})
            return grouped

    def status_with_links(self) -> pd.DataFrame:
        """Return the status DataFrame with clickable log file links."""
        self._status = None  # Force reload
        df = self.status.copy()

        # Create clickable links for log files
        def make_out_link(row):
            if row['log_file'] and row['state'] != 'not_run':
                project_dir = os.path.join(self.root_path, row.name)
                log_path = os.path.join(project_dir, row['log_file'])
                # Create a file:// URL for local files
                file_url = f"file://{log_path}"
                return f'<a href="{file_url}" target="_blank">{row["log_file"].split("/")[-1]}</a>'
            return row['log_file']

        # Create clickable links for error log files
        def make_error_link(row):
            if row['err_file'] and row['state'] != 'not_run':
                project_dir = os.path.join(self.root_path, row.name)
                error_log_path = os.path.join(project_dir, row['err_file'])
                # Create a file:// URL for local files
                file_url = f"file://{error_log_path}"
                return f'<a href="{file_url}" target="_blank">{row["err_file"].split("/")[-1]}</a>'
            return row['err_file']

        df['log_file_link'] = df.apply(make_out_link, axis=1)
        df['err_file_link'] = df.apply(make_error_link, axis=1)
        return df

    def display_status_with_links(self):
        """Display the status DataFrame with clickable links in Jupyter."""
        from IPython.display import HTML, display

        df = self.status_with_links()
        # Drop the original log_file columns and rename the link columns
        df_display = df.drop(['log_file', 'err_file'], axis=1)
        df_display = df_display.rename(columns={
            'log_file_link': 'log_file',
            'err_file_link': 'err_file'
        })

        # Convert to HTML and display
        html = df_display.to_html(escape=False)
        display(HTML(html))

    def _repr_html_(self):
        """Return HTML representation for Jupyter notebooks."""
        df = self.status_with_links()
        # Drop the original log_file columns and rename the link columns
        df_display = df.drop(['log_file', 'err_file'], axis=1)
        df_display = df_display.rename(columns={
            'log_file_link': 'log_file',
            'err_file_link': 'err_file'
        })
        return df_display.to_html(escape=False)

    def generate_inputs_yaml(self,
        inputs_template_path: str,
        overwrite: bool = False,
        input_generator: Optional[Callable] = None,
    ):
        """Generate an inputs.yaml file in each project directory based on the dataset configuration.

        Args:
            inputs_template_path (str): The file path to the Jinja2 template file that will be
                used to generate the inputs YAML files.
            input_generator (callable): A callable function intended for generating input values.
                Currently, it is called with the project directory name (DIR) as its argument
            overwrite (bool): Whether to overwrite existing inputs.yaml files. Default is False.

        """
        # Load the template
        with open(inputs_template_path, 'r') as f:
            template_str = f.read()

        template = jinja2.Template(template_str)
        skipped = 0
        generated = 0
        for project_dir in self.project_directories:
            inputs_yaml_path = os.path.join(project_dir, 'inputs.yaml')
            if os.path.exists(inputs_yaml_path) and not overwrite:
                skipped += 1
                continue
            generated += 1
            # Get the directory name
            DIR = os.path.basename(os.path.normpath(project_dir))
            # Render the template with project defaults
            generated_input = input_generator(DIR) if input_generator else {}
            rendered_yaml = template.render(DIR=DIR, **generated_input)
            # Write the rendered YAML to inputs.yaml
            with open(inputs_yaml_path, 'w') as f:
                f.write(rendered_yaml)
        print(f"Generated {generated} inputs.yaml files. Skipped {skipped} existing files.")

    def launch_workflow(self,
        include_groups: list[int] = [],
        exclude_groups: list[int] = [],
        n_jobs: int = 0,
        pool_size: int = 1,
        slurm: bool = False,
        job_template: str = None,
        debug: bool = False):
        """Launch the workflow for each project directory in the dataset.

        Args:
            include_groups (list[int]):
                List of group IDs to be run.
            exclude_groups (list[int]):
                List of group IDs to be excluded.
            n_jobs (int):
                Number of jobs to launch. If 0, all jobs are launched.
            pool_size (int):
                Number of parallel jobs to run using subprocess. Default is 1.
            slurm (bool):
                Whether to submit the workflow to SLURM.
            job_template (str):
                Path to the bash script or SLURM job template file. You can use Jinja2
                templating to customize the job script using the fields of
                the input YAML and the columns of project status dataframe.
            debug (bool):
                Only print the commands without executing them.

        """
        # Include/exclude groups should not intersect
        if include_groups and exclude_groups:
            intersection = set(include_groups).intersection(set(exclude_groups))
            if intersection:
                raise ValueError(f"include_groups and exclude_groups intersect: {intersection}")

        if slurm and not job_template:
            raise ValueError("job_template must be provided when slurm is True")

        # Determine pool size for parallel execution
        pool_size = min(pool_size, os.cpu_count())
        if pool_size <= 0 and not slurm:
            pool_size = os.cpu_count()
        parallel = pool_size > 1 and not slurm
        jobs_to_run = []

        n = 0
        for project_dir in self.project_directories:
            rel_path = os.path.relpath(project_dir, self.root_path)
            project_status = self.status.loc[rel_path].to_dict()
            project_status['rel_path'] = rel_path
            # Check group inclusion/exclusion
            if project_status['group'] in exclude_groups or \
                (include_groups and project_status['group'] not in include_groups):
                continue
            n += 1
            if n_jobs > 0 and n > n_jobs:
                break
            inputs_yaml_path = os.path.join(project_dir, 'inputs.yaml')
            if not os.path.exists(inputs_yaml_path):
                warn(f"{inputs_yaml_path} not found. Skipping {project_dir}")
                continue

            inputs_config = load_yaml(inputs_yaml_path)

            with open(job_template, 'r') as f:
                template_str = f.read()

            template = jinja2.Template(template_str)
            rendered_script = template.render(**inputs_config, **project_status,
                                                DIR=os.path.basename(os.path.normpath(project_dir)),
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
                    print(f'Job script: {job_script_path}')
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
