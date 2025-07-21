import os
import pytest
import shutil
import numpy as np
from conftest import get_analysis_file

from model_workflow.utils.constants import *
from model_workflow.utils.type_hints import *
from model_workflow.utils.auxiliar import load_json
from model_workflow.mwf import project_requestables, md_requestables


@pytest.fixture(scope="class", params=["A0001", "A01IP"])
def test_accession(request):
    return request.param

@pytest.mark.release
class TestMWFRun:
    """Test full workflow for different accessions"""
    
    # Arguments to reduce long test execution time
    task_arguments = {
        'clusters': {'frames_limit': 10, 'desired_n_clusters': 2},
        'pockets': {'maximum_pockets_number': 2},
        'dist': {'frames_limit': 2},
        'energies': {'frames_limit': 2},
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
