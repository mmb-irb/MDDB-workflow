from mddb_workflow.utils.auxiliar import load_yaml, is_glob
from mddb_workflow.mwf import workflow
import pandas as pd
import subprocess
import jinja2
import time
import glob
import os


class Dataset:
    """
    Class to manage and process a dataset of MDDB projects.
    """
    def __init__(self, dataset_yaml_path: str):
        """
        Initializes the Dataset object.

        Args:
            dataset_yaml_path (str): Path to the dataset YAML file.
        """
        self.root_path = os.path.dirname(os.path.abspath(dataset_yaml_path))
        self.config = load_yaml(dataset_yaml_path)
        # Properties cache
        self.project_directories = self.get_project_directories()
        self.groups = self.get_groups()
        self._status = None

    def _resolve_directory_patterns(self, dir_patterns: list[str]) -> list[str]:
        """
        Helper method to resolve directory patterns (glob or absolute/relative paths).
        Validates that resolved paths are under the dataset root.

        Args:
            dir_patterns (list[str]): List of directory patterns to resolve.

        Returns:
            list[str]: List of resolved directory paths.
        
        Raises:
            ValueError: If a resolved path is outside the dataset root.
        """
        directories = []
        for dir_pattern in dir_patterns:
            if is_glob(dir_pattern):
                matched_dirs = glob.glob(os.path.join(self.root_path, dir_pattern))
                # keep only directories
                directories.extend([p for p in matched_dirs if os.path.isdir(p)])
            else:
                p = dir_pattern
                if not os.path.isabs(p):
                    p = os.path.abspath(os.path.join(self.root_path, p))
                # ensure the resolved path is under the dataset root
                if os.path.commonprefix([p, self.root_path]) != self.root_path:
                    raise ValueError(f"Project directory '{p}' is outside the dataset root '{self.root_path}'")
                directories.append(p)
        return directories

    def get_project_directories(self) -> list[str]:
        """
        Retrieves the list of project directories from the dataset configuration.

        Returns:
            list: List of project directory patterns.
        """
        return self._resolve_directory_patterns(self.config.get('project_directories', []))

    def get_groups(self) -> dict[str, list[str]]:
        groups = {}
        for group in self.config.get('groups', []):
            groups[group] = self._resolve_directory_patterns(self.config['groups'][group])
        return groups

    def generate_inputs_yaml(self, inputs_template_path: str, input_generator: callable, overwrite: bool = False):
        """
        Generates an inputs.yaml file in each project directory based on the dataset configuration.
        
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
        
        for project_dir in self.project_directories:
            inputs_yaml_path = os.path.join(project_dir, 'inputs.yaml')
            if os.path.exists(inputs_yaml_path) and not overwrite:
                continue
            
             # Get the directory name
            DIR = os.path.basename(os.path.normpath(project_dir))
            
            # Render the template with project defaults
            rendered_yaml = template.render(DIR=DIR, title=input_generator(DIR))

            # Write the rendered YAML to inputs.yaml
            with open(inputs_yaml_path, 'w') as f:
                f.write(rendered_yaml)
            break

    @property
    def status(self) -> pd.DataFrame:
        """
        Retrieves last line from logs from all project directories as a pandas DataFrame.

        Returns:
            pd.DataFrame: Index is project directory; columns: state, message, log_file, error_log_file.
        """
        if self._status is not None:
            return self._status

        rows = []
        for project_dir in self.project_directories:
            # RUBEN: por ahora usamos el mismo patron de log que en launch_workflow
            # en un futuro se recuperar el estado a partir de .register/.mwf_cache
            log_files = glob.glob(os.path.join(project_dir, 'logs', 'mwf*[0-9].out'))
            err_files = glob.glob(os.path.join(project_dir, 'logs', 'mwf*[0-9].err'))

            if log_files:
                log_files.sort()
                log_files = [log_files[-1]]  # Take the most recent one
                with open(log_files[0], 'r') as f:
                    lines = f.read().splitlines()
                    if lines:
                        last_line = lines[-1].strip()
                        if len(last_line) > 80:
                            last_line = last_line[:80]+'...'
                    else:
                        last_line = ''

                if last_line == 'Done!':
                    state, message, log_file = 'done', last_line, log_files[0]
                else:
                    state, message, log_file = 'error', last_line if last_line else 'Empty log file', log_files[0]
            else:
                state, message, log_file = 'not_run', 'No output log available', None

            # Handle error log files
            if err_files:
                err_files.sort()
                err_file = err_files[-1]  # Take the most recent one
            else:
                err_file = None

            # Check if files were modified recently (within last 5 minutes) to detect running state
            current_time = time.time()
            recently_modified_threshold = 300  # 5 minutes in seconds

            last_modified = ''
            if (log_file and os.path.exists(log_file) and (current_time - os.path.getmtime(log_file)) < recently_modified_threshold) or \
                (err_file and os.path.exists(err_file) and (current_time - os.path.getmtime(err_file)) < recently_modified_threshold):
                 if state != 'done':
                      state = 'running'
            else:
                # Save last modification time if not running
                if log_file and os.path.exists(log_file):
                    last_modified = time.strftime('%H:%M:%S %d/%m/%y', time.localtime(os.path.getmtime(log_file)))
            
            rows.append({
                'rel_path': os.path.relpath(project_dir, self.root_path),
                'state': state,
                'message': message,
                'log_file': os.path.relpath(log_file, project_dir) if log_file else '',
                'err_file': os.path.relpath(err_file, project_dir) if err_file else '',
                'last_modified': last_modified
            })

        df = pd.DataFrame(rows).set_index('rel_path').sort_index()
        # Assign an integer group id for identical messages, with messages sorted first
        unique_messages = sorted(df['message'].unique())
        if 'Done!' in unique_messages:
            # Ensure 'Done!' is always group 0, then assign other messages
            unique_messages.remove('Done!')
            mapping = {'Done!': 0}
            mapping.update({msg: idx + 1 for idx, msg in enumerate(unique_messages)})
        else:
            mapping = {msg: idx for idx, msg in enumerate(unique_messages)}
        df['group'] = df['message'].map(mapping).astype(int)
        self._status = df
        return self._status
    
    def show_groups(self, cmd=False):
        """
        Displays the groups of projects based on their status messages.
        """
        if cmd:
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
        """
        Returns the status DataFrame with clickable log file links.
        """
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
        """
        Display the status DataFrame with clickable links in Jupyter.
        """
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

    def launch_workflow(self,
        include_groups: list[int]=[],
        exclude_groups: list[int]=[],
        slurm: bool=False,
        job_template: str=None,
        debug: bool=False):
        """
        Launches the workflow for each project directory in the dataset.
        Args:
            include_groups (list[int]):
                List of group IDs to be run.
            exclude_groups (list[int]):
                List of group IDs to be excluded.
            slurm (bool):
                Whether to submit the workflow to SLURM.
            job_template (str):
                Path to the SLURM job template file. You can use Jinja2
                templating to customize the job script using the fields of
                the input YAML and the columns of project status dataframe.
        """
        # Include/exclude groups should not intersect
        if include_groups and exclude_groups:
            intersection = set(include_groups).intersection(set(exclude_groups))
            if intersection:
                raise ValueError(f"include_groups and exclude_groups intersect: {intersection}")
            
        if slurm and not job_template:
            raise ValueError("job_template must be provided when slurm is True")
        for project_dir in self.project_directories:
            rel_path = os.path.relpath(project_dir, self.root_path)
            project_status = self.status.loc[rel_path].to_dict()
            project_status['rel_path'] = rel_path
            # Check group inclusion/exclusion
            if project_status['group'] in exclude_groups or \
                (include_groups and project_status['group'] not in include_groups):
                continue
            # Launch workflow
            if slurm:
                # SLURM execution
                inputs_yaml_path = os.path.join(project_dir, 'inputs.yaml')
                if not os.path.exists(inputs_yaml_path):
                    print(f"Warning: {inputs_yaml_path} not found. Skipping {project_dir}")
                    continue
                
                inputs_config = load_yaml(inputs_yaml_path)

                with open(job_template, 'r') as f:
                    template_str = f.read()
                
                template = jinja2.Template(template_str)
                rendered_script = template.render(**inputs_config, **project_status)

                job_script_path = os.path.join(project_dir, 'mwf_slurm_job.sh')
                log_dir = os.path.join(project_dir, 'logs')
                os.makedirs(log_dir, exist_ok=True)
                with open(job_script_path, 'w') as f:
                    f.write(rendered_script)
                
                os.chmod(job_script_path, 0o755)

                if debug:
                    print(f"cd {project_dir}")
                    print(f"sbatch --output=logs/mwf_%j.out --error=logs/mwf_%j.err mwf_slurm_job.sh ")
                else:
                    print(f"Submitting SLURM job for {project_dir}")
                    subprocess.run(['sbatch', 
                                    '--output=logs/mwf_%j.out',
                                    '--error=logs/mwf_%j.err',
                                    job_script_path],
                                cwd=project_dir)

            else:
                # Normal Python execution
                raise NotImplementedError("Python execution is not implemented yet.")
                workflow(working_directory=project_dir)
