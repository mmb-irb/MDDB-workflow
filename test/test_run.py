import os
import pytest
import shutil
import numpy as np
from conftest import get_analysis_file

from model_workflow.utils.constants import *
from model_workflow.utils.type_hints import *
from model_workflow.utils.auxiliar import load_json
from model_workflow.mwf import project_requestables, md_requestables, workflow
from model_workflow.console import run_parser

@pytest.mark.release
@pytest.mark.parametrize("test_accession", ["A0001", "A01IP" ], scope="class")
class TestMWFRun:
    """Test full workflow for different accessions"""
    
    # Argument to reduce long test execution time
    task_arguments = {
        'clusters': {'frames_limit': 10, 'desired_n_clusters': 2},
        'pockets': {'maximum_pockets_number': 2},
        'dist': {'frames_limit': 2},
        # Default arguments for tasks
        'inter': {'frames_limit': 10},
        'energies': {'frames_limit': 2},
        'hbonds': {'time_splits': 100 },
        'rmsds': {'frames_limit': 10},
        'rgyr': {'frames_limit': 10},
        'pairwise': {'frames_limit': 10},
        'pca': {'frames_limit': 10},
        'perres': {'frames_limit': 10},
        'tmscore': {'frames_limit': 10},
        'sas': {'frames_limit': 10},
        'markov': { 'rmsd_selection': PROTEIN_AND_NUCLEIC },

    }
    def _run_and_log_task(self, task_name: str, task, target_obj, project_dir: str, capsys):
        """Helper method to run a task and log its output"""
        try:
            task.args = self.task_arguments.get(task_name, {})
            task(target_obj)
        except Exception as e:
            pytest.fail(f"Task '{task_name}' failed with error: {str(e)}")
        finally:
            # Capture stdout and stderr
            out, err = capsys.readouterr()
            # Write the output to a log file
            log_file = f"{project_dir}/logs/{task_name}_output.log"
            os.makedirs(os.path.dirname(log_file), exist_ok=True)
            with open(log_file, 'w') as f:
                f.write("STDOUT:\n")
                f.write(out)
                f.write("\nSTDERR:\n")
                f.write(err)

    @pytest.mark.parametrize("project_task", project_requestables.keys())
    def test_project_task(self, project: 'Project', project_task: str, capsys):
        """Test that each project task runs without errors"""
        project.overwritables = {project_task}
        self._run_and_log_task(project_task, project_requestables[project_task], project, project.directory, capsys)

    @pytest.mark.parametrize("md_task", md_requestables.keys())
    def test_md_task(self, project: 'Project', md_task: str, capsys):
        """Test that each analysis runs without errors"""
        if md_task == 'dihedrals' or \
            (md_task == 'pockets' and project.accession == 'A01IP') or \
            (md_task == 'dist' and project.accession == 'A01IP'):
            pytest.skip(f"Skipping analysis '{md_task}' for now.")
            
        md: MD = project.mds[0]
        md.overwritables = {md_task}
        self._run_and_log_task(md_task, md_requestables[md_task], md, project.directory, capsys)

    @pytest.mark.skip(reason="Test is skipped")
    def test_TMscores_analysis(self, project : 'Project'):
        """Test that RMSD analysis runs and produces expected output"""
        # Run the analysis
        analysis = 'tmscores'
        self.test_analysis_execution(project, analysis)

        # Check that the output file was created
        output_file = f"{project.directory}/replica_1/{OUTPUT_TMSCORES_FILENAME}"
        assert os.path.exists(output_file), f"Output file '{output_file}' was not created"
        
        # Get the reference analysis file
        analysis_file = get_analysis_file(project, analysis)
        # Load the results
        results = load_json(output_file)
        reference = load_json(analysis_file.absolute_path)
        # Check that the results match the expected output
        now = np.array(results['data'][0]['values']).mean()
        ref_mn = np.array(reference['data'][0]['values']).mean()

        assert np.isclose(now, ref_mn, atol=0.1), f"Output values do not match expected values. Now: {now}, Ref: {ref_mn}"

@pytest.mark.CI
@pytest.mark.release
class TestWorkflow:
    """Test the workflow execution"""

    def test_workflow_execution(self, test_data_dir: str):
        """Test that the workflow runs without errors, simulating console execution"""
        
        cwd = os.getcwd()
        working_directory = os.path.join(test_data_dir, 'output/test_workflow')
        # Ensure the working directory doesn't exist and it is empty
        if os.path.exists(working_directory):
            # Remove the directory if it exists
            shutil.rmtree(working_directory)
        # Copy the inputs from raw_project
        shutil.copytree(
            os.path.join(test_data_dir, 'input/raw_project'),
            working_directory,
            dirs_exist_ok=True
        )

        # Use argparse to get default values for all arguments
        # Parse an empty argument list to get the default values
        args = run_parser.parse_args([])
        args_dict = vars(args)
        
        # Modify only the specific arguments we want to change
        args_dict.update({
            #'input_structure_filepath': '5ggr-rs1-310k-0ns.pdb',
            #'input_trajectory_filepaths': 'dynamic-nopbc-2',
            'working_directory': working_directory,
            'setup': True,
        })

        # Create a copy without common arguments for splitting
        run_specific_args = args_dict.copy()
        # Remove common arguments that should not be passed to workflow or project
        del run_specific_args['no_symlinks']
        
        # Split arguments between project_parameters and workflow_args
        # as done in console.py's main function
        project_parameters = {}
        workflow_kwargs = {}
        
        # Arguments that are direct named parameters of the workflow function
        workflow_direct_arg_names = {
            'working_directory', 
            'download', 
            'setup', 
            'include', 
            'exclude', 
            'overwrite'
        }
        
        for k, v in run_specific_args.items():
            if k in workflow_direct_arg_names:
                workflow_kwargs[k] = v
            else:
                # All other arguments go to project_parameters
                project_parameters[k] = v
        
        # Call workflow in a way that mimics console.py's main function
        workflow(project_parameters=project_parameters, **workflow_kwargs)
        os.chdir(cwd)