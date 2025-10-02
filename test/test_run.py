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

    @pytest.mark.parametrize("project_task", project_requestables.keys())
    def test_project_task(self, project: 'Project', project_task: str):
        """Test that each project task runs without errors"""
        project.overwritables = {project_task}
        project_requestables[project_task](project)

    @pytest.mark.parametrize("md_task", md_requestables.keys())
    def test_md_task(self, project: 'Project', md_task: str, capsys):
        """Test that each analysis runs without errors"""
        if md_task == 'dihedrals':
            pytest.skip(f"Skipping analysis.")
        elif project.accession == 'A01IP' and md_task in ['pockets', 'pockets', 'dist', 'energies']:
            pytest.skip(f"Skipping analysis to save time.")
        elif project.accession == 'A01V7' and md_task in ['pockets']:
            pytest.skip(f"Skipping analysis to save time.")

        md: MD = project.mds[0]
        md.overwritables = {md_task}
        md_requestables[md_task](md)

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

@pytest.mark.CI
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
