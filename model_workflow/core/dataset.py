from model_workflow.utils.auxiliar import load_yaml, is_glob
from model_workflow.mwf import workflow
import glob
import subprocess
import os
from pathlib import Path


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
        self.root_path = Path(dataset_yaml_path).parent
        self.config = load_yaml(dataset_yaml_path)
        self.project_directories = self._get_project_directories()

    def _get_project_directories(self):
        """
        Retrieves the list of project directories from the dataset configuration.

        Returns:
            list: List of project directory patterns.
        """

        config_project_directories = self.config.get('global', {}).get('project_directories', [])
        project_directories = []
        for i, dir_pattern in enumerate(config_project_directories):
            if is_glob(dir_pattern):
                matched_dirs = glob.glob(str(self.root_path / dir_pattern))
                # keep only directories
                project_directories.extend([p for p in matched_dirs if Path(p).is_dir()])
            else:
                project_directories.append(str(self.root_path / dir_pattern))

        return project_directories

    def launch_workflow(self, slurm=False, job_template=None):
        """
        Launches the workflow for each project directory in the dataset.
        Args:
            slurm (bool): Whether to submit the workflow to SLURM.
            job_template (str): Path to the SLURM job template file. You can use Jinja2 templating to customize the job script.
        """
        if slurm and not job_template:
            raise ValueError("job_template must be provided when slurm is True")
        if slurm:
            # Lazy import jinja2 for now to not break environments
            import jinja2
        for project_dir in self.project_directories:
            if slurm:
                inputs_yaml_path = Path(project_dir) / 'inputs.yaml'
                if not inputs_yaml_path.exists():
                    print(f"Warning: {inputs_yaml_path} not found. Skipping {project_dir}")
                    continue
                
                inputs_config = load_yaml(str(inputs_yaml_path))

                with open(job_template, 'r') as f:
                    template_str = f.read()
                
                template = jinja2.Template(template_str)
                rendered_script = template.render(**inputs_config)

                job_script_path = Path(project_dir) / 'mwf_slurm_job.sh'
                with open(job_script_path, 'w') as f:
                    f.write(rendered_script)
                
                os.chmod(job_script_path, 0o755)

                print(f"Submitting SLURM job for {project_dir}")
                subprocess.run(['sbatch', str(job_script_path)], cwd=project_dir)

            else:
                # Normal Python execution
                raise NotImplementedError("Python execution is not implemented yet.")
                workflow(working_directory=project_dir)