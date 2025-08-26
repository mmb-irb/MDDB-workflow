import numpy as np
import os, sys, shutil, pytest

from model_workflow.utils.constants import *
from model_workflow.utils.type_hints import *
from model_workflow.utils.auxiliar import load_json, InputError
from model_workflow.mwf import project_requestables, md_requestables
from model_workflow.console import  main


@pytest.mark.release
class TestRunAll:
    """Test all the tasks workflow for reference accessions:
      - A0001: base case
      - A01IP: for membrane analyses
      - A01V7: for only ligand
    """
    @pytest.fixture(scope="class", params=["A0001", "A01IP", "A01V7"])
    def test_accession(self, request):
        return request.param
    
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
        if md_task == 'dihedrals':
            pytest.skip(f"Skipping analysis.")
        elif project.accession == 'A01IP' and md_task in ['pockets', 'pockets', 'dist', 'energies']:
            pytest.skip(f"Skipping analysis to save time.")

        md: MD = project.mds[0]
        md.overwritables = {md_task}
        self._run_and_log_task(md_task, md_requestables[md_task], md, project.directory, capsys)

@pytest.mark.release
class TestRunFlags:
    """Test the run subcommand with different flags"""

    def test_top_no(self, test_data_dir: str):
        """Test the flag -top no. Add , 'protmap' for coverage"""

        working_directory = os.path.join(test_data_dir, 'output/test_top_no')
        # Remove the directory if it already exists to ensure a clean state
        if os.path.exists(working_directory):
            shutil.rmtree(working_directory)
        # Copy the inputs from raw_project
        shutil.copytree(
            os.path.join(test_data_dir, 'input/raw_project'),
            working_directory,
            ignore=shutil.ignore_patterns('topology.tpr'),
            dirs_exist_ok=True)

        sys.argv = ['model_workflow', 'run',
                    '-dir', working_directory,
                    '-stru', 'raw_structure.pdb',
                    '-traj', 'raw_trajectory.xtc',
                    '-top', 'no',
                    '-i', 'setup', 'rmsds', 'protmap']
        main()
        # Change back to the original directory
        os.chdir(test_data_dir)

    def test_no_inputs(self, test_data_dir: str):
        """Test the workflow without no inputs yaml"""

        working_directory = os.path.join(test_data_dir, 'output/test_no_inputs')
        # Remove the directory if it already exists to ensure a clean state
        if os.path.exists(working_directory):
            shutil.rmtree(working_directory)
        # Copy the inputs from raw_project
        shutil.copytree(
            os.path.join(test_data_dir, 'input/raw_project'),
            working_directory,
            ignore=shutil.ignore_patterns('inputs.yaml'),
            dirs_exist_ok=True)

        sys.argv = ['model_workflow', 'run',
                    '-dir', working_directory,
                    '-top', 'topology.tpr',
                    '-md', 'replica_1',
                    'raw_structure.pdb',
                    'raw_trajectory.xtc',
                    '-i', 'setup', 'rmsds']

        with pytest.raises(InputError, match='Missing inputs file "inputs.yaml"'):
            main()

        os.chdir(test_data_dir)

@pytest.mark.release
class TestRunSpecial:
    """Test for special cases and specific accessions"""

    @pytest.mark.parametrize("test_accession", ["cg_test"], scope="class")
    def test_CG(self, project: 'Project'):
        """Test coarse-grained (CG) model."""
        # Only two tasks. In the future "all" should be supported 
        md = project.mds[0]
        md.get_processed_interactions(md)
        md.run_rmsds_analysis(md)

    @pytest.mark.parametrize("test_accession", ["test_020"], scope="class")
    def test_dihedrals(self, project: 'Project'):
        """Test dihedrals energy calculation."""
        project.get_dihedrals()


# ========== Analysis Output Tests =========
# Test to see if the output of the analysis is the same as the reference file
def get_analysis_file(project: 'Project', analysis_type: str):
    """Download and provide the standard structure file"""
    output_path = os.path.join(project.directory, f"mda.{analysis_type}_REF.json")
    file_obj = File(output_path)
    # Only download if file doesn't exist yet
    if not file_obj.exists:
        project.remote.download_analysis_data(analysis_type,  file_obj)
    return file_obj

class TestAnalysisOutput:
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

class TestDensityAnalysis:
    @pytest.fixture(scope="class")
    def test_accession(self):
        """Override the default accession for this test"""
        return "A01JP.1"
    
    @pytest.fixture(scope="class")
    def output_file(self, test_proj_dir):
        """Create an output file path for the density analysis results"""
        return os.path.join(test_proj_dir, "density_output.json")
    
    def test_density_analysis(self, project : 'Project'):
        """Test that density analysis runs and produces expected output"""
        # Download the reference file and run the analysis
        analysis_file = get_analysis_file(project, 'density')
        md : MD = project.mds[0]
        md.overwritables = {'density'}
        md.run_density_analysis()

        # Check that the output file was created
        output_file = f"{project.directory}/replica_1/{OUTPUT_DENSITY_FILENAME}"
        assert os.path.exists(output_file), f"Output file '{output_file}' was not created"
        
        # Load the results
        results = load_json(output_file)
        reference = load_json(analysis_file.path)

        # Compare structure and keys
        assert set(results.keys()) == set(reference.keys()), \
            f"Result keys don't match: {set(results.keys())} vs {set(reference.keys())}"
        
        assert set(results.get('data', {}).keys()) == set(reference.get('data', {}).keys()), \
            "Data subkeys don't match"
        
        # Check components exist
        assert 'comps' in results.get('data', {}), "Missing 'comps' key in results"
        assert 'z' in results.get('data', {}), "Missing 'z' key in results"

        # Check z array length matches
        assert len(results['data']['z']) == len(reference['data']['z']), \
            f"Z array length mismatch: {len(results['data']['z'])} vs {len(reference['data']['z'])}"
        
        # Check same number of components
        assert len(results['data']['comps']) == len(reference['data']['comps']), \
            f"Component count mismatch: {len(results['data']['comps'])} vs {len(reference['data']['comps'])}"