from model_workflow.utils.auxiliar import load_yaml, is_glob
from model_workflow.mwf import workflow
import glob
import subprocess
import os
import pandas as pd


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
        self._project_directories = None
        self._status = None

    @property
    def project_directories(self) -> list[str]:
        """
        Retrieves the list of project directories from the dataset configuration.

        Returns:
            list: List of project directory patterns.
        """
        if self._project_directories is not None:
            return self._project_directories
        
        config_project_directories = self.config.get('global', {}).get('project_directories', [])
        self._project_directories = []
        for i, dir_pattern in enumerate(config_project_directories):
            if is_glob(dir_pattern):
                matched_dirs = glob.glob(os.path.join(self.root_path, dir_pattern))
                # keep only directories
                self._project_directories.extend([p for p in matched_dirs if os.path.isdir(p)])
            else:
                p = dir_pattern
                if not os.path.isabs(p):
                    p = os.path.abspath(os.path.join(self.root_path, p))
                # ensure the resolved path is under the dataset root
                if os.path.commonprefix([p, self.root_path]) != self.root_path:
                    raise ValueError(f"Project directory '{p}' is outside the dataset root '{self.root_path}'")
                self._project_directories.append(p)

        return self._project_directories

    @property
    def status(self) -> pd.DataFrame:
        """
        Retrieves last line from logs from all project directories as a pandas DataFrame.

        Returns:
            pd.DataFrame: Index is project directory; columns: state, message, log_file.
        """
        if self._status is not None:
            return self._status

        rows = []
        for project_dir in self.project_directories:
            # RUBEN: por ahora usamos el mismo patron de log que en launch_workflow
            # en un futuro se recuperar el estado a partir de .register/.mwf_cache
            log_files = glob.glob(os.path.join(project_dir, '*[0-9].out'))
            if log_files:
                log_files = [log_files[-1]]  # Take the most recent one
                with open(log_files[0], 'r') as f:
                    last_line = f.read().splitlines()[-1].strip()

                if last_line == 'Done!':
                    state, message, log_file = 'done', last_line, log_files[0]
                else:
                    state, message, log_file = 'error', last_line, log_files[0]
            else:
                state, message, log_file = 'not_run', None, 'not_run'

            rows.append({
                'rel_path': os.path.relpath(project_dir, self.root_path),
                'state': state,
                'message': message,
                'log_file': os.path.relpath(log_file, project_dir) if log_file else None
            })

        df = pd.DataFrame(rows).set_index('rel_path').sort_index()
        # Assign an integer group id for identical messages, with messages sorted first
        unique_messages = sorted(df['message'].unique())
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

    def launch_workflow(self,
        include_groups: list[int]=[],
        exclude_groups: list[int]=[],
        slurm: bool=False,
        job_template: str=None):
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
        if slurm:
            # Lazy import jinja2 for now to not break environments
            import jinja2
        for project_dir in self.project_directories:
            project_status = self.status.loc[os.path.relpath(project_dir, self.root_path)].to_dict()
            # Check group inclusion/exclusion
            group_id = project_status['group']
            if group_id in exclude_groups or (include_groups and group_id not in include_groups):
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
                with open(job_script_path, 'w') as f:
                    f.write(rendered_script)
                
                os.chmod(job_script_path, 0o755)

                print(f"Submitting SLURM job for {project_dir}")
                subprocess.run(['sbatch', 'mwf_slurm_job.sh'], cwd=project_dir)

            else:
                # Normal Python execution
                raise NotImplementedError("Python execution is not implemented yet.")
                workflow(working_directory=project_dir)
