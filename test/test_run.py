import os
import sys
import shutil
import pytest
import requests
from mddb_workflow.utils.constants import *
from mddb_workflow.utils.type_hints import *
from mddb_workflow.utils.auxiliar import InputError, load_json
from mddb_workflow.mwf import project_requestables, md_requestables
from mddb_workflow.console import main


@pytest.mark.release
class TestRunAll:
    """Test all the tasks workflow for reference accessions.

    - A0001: base case
    - A01IP: for membrane analyses
    - A025N: for only ligand
    - A02F9: only lipids

    """

    @pytest.fixture(scope='class', params=['A0001', 'A01IP', 'A025N', 'A02F9'])
    def test_accession(self, request):
        """Fixture to provide different test accessions."""
        return request.param

    # Arguments to reduce long test execution time
    task_arguments = {
        'clusters': {'frames_limit': 10, 'desired_n_clusters': 2},
        'pockets': {'maximum_pockets_number': 2},
        'dist': {'frames_limit': 2},
        'energies': {'frames_limit': 2},
    }

    @pytest.mark.parametrize('project_task', project_requestables.keys())
    def test_project_task(self, project: 'Project', project_task: str):
        """Test that each project task runs without errors."""
        project.overwritables = {project_task}
        project_requestables[project_task](project)
        if project_task == 'pmeta':
            # Check that the pmeta directory was created
            pmeta_dir = os.path.join(project.directory, 'metadata.json')
            metadata = load_json(pmeta_dir)
            syskeys = {
                'A0001': ['protein', 'protein only'],
                'A01IP': ['protein', 'ligand', 'lipid', 'carbohydrate', 'membrane'],
                'A025N': ['ligand', 'ligand only'],
                'A02F9': ['lipid', 'solvent', 'lipid only', 'membrane']
            }
            assert metadata['SYSKEYS'] == syskeys[project.accession]

    @pytest.mark.parametrize('md_task', md_requestables.keys())
    def test_md_task(self, project: 'Project', md_task: str, capsys):
        """Test that each analysis runs without errors."""
        if md_task == 'dihedrals':
            pytest.skip('Skipping analysis.')
        elif project.accession == 'A01IP' and md_task in [
            'pockets',
            'pockets',
            'dist',
            'energies',
        ]:
            pytest.skip('Skipping analysis to save time.')
        elif project.accession == 'A01V7' and md_task in ['pockets']:
            pytest.skip('Skipping analysis to save time.')

        md: MD = project.mds[0]
        md.overwritables = {md_task}
        md_requestables[md_task](md)


@pytest.mark.release
class TestRunFlags:
    """Test the run subcommand with different flags."""

    @pytest.mark.CI
    def test_top_no(self, test_data_dir: str):
        """Test the flag -top no. Add , 'protmap' for coverage."""
        working_directory = os.path.join(test_data_dir, 'output/test_top_no')
        # Remove the directory if it already exists to ensure a clean state
        if os.path.exists(working_directory):
            shutil.rmtree(working_directory)
        # Copy the inputs from raw_project
        shutil.copytree(
            os.path.join(test_data_dir, 'input/raw_project'),
            working_directory,
            ignore=shutil.ignore_patterns('raw_topology.tpr'),
            dirs_exist_ok=True,
        )

        sys.argv = [
            'mddb_workflow', 'run',
            '-dir', working_directory,
            '-stru', 'raw_structure.pdb',
            '-traj', 'raw_trajectory.xtc',
            '-top', 'no',
            '-i', 'setup', 'rmsds', 'protmap',
        ]
        main()
        # Change back to the original directory
        os.chdir(test_data_dir)

    def test_no_inputs(self, test_data_dir: str):
        """Test the workflow without no inputs yaml."""
        working_directory = os.path.join(test_data_dir, 'output/test_no_inputs')
        # Remove the directory if it already exists to ensure a clean state
        shutil.copytree(
            os.path.join(test_data_dir, 'input/raw_project'),
            working_directory,
            ignore=shutil.ignore_patterns('inputs.yaml'),
            dirs_exist_ok=True
        )

        sys.argv = [
            'mddb_workflow',
            'run',
            '-dir', working_directory,
            '-top', 'raw_topology.tpr',
            '-md', 'replica_1', 'raw_structure.pdb', 'raw_trajectory.xtc',
            '-i', 'setup', 'rmsds',
        ]

        with pytest.raises(InputError, match='Missing inputs file "inputs.yaml"'):
            main()

        os.chdir(test_data_dir)

    def test_no_internet(self, test_data_dir: str, monkeypatch):
        """Test the workflow with no internet connection."""

        def raise_connection_error(*args, **kwargs):
            raise requests.exceptions.ConnectionError('No internet connection')

        monkeypatch.setattr(requests, 'get', raise_connection_error)

        working_directory = os.path.join(test_data_dir, 'output/test_no_internet')
        # Remove the directory if it already exists to ensure a clean state
        if os.path.exists(working_directory):
            shutil.rmtree(working_directory)
        # Copy the inputs from raw_project
        shutil.copytree(
            os.path.join(test_data_dir, 'input/raw_project'),
            working_directory,
            dirs_exist_ok=True,
        )

        sys.argv = [
            'mddb_workflow', 'run',
            '-dir', working_directory,
            '-top', 'raw_topology.tpr',
            '-md', 'replica_1', 'raw_structure.pdb', 'raw_trajectory.xtc',
            '-e', 'network', 'interdeps', 'membs', 'clusters', 'pockets',]

        main()
        os.chdir(test_data_dir)


@pytest.mark.CI
class TestRunSpecial:
    """Test for special cases and specific accessions."""

    @pytest.mark.parametrize('test_accession', ['cg_test', 'cg_test2'], scope='class')
    def test_CG(self, project: 'Project'):
        """Test coarse-grained (CG) model."""
        # Only two tasks. In the future "all" should be supported
        md = project.mds[0]
        md.get_processed_interactions(md)
        md.run_rmsds_analysis(md)

    @pytest.mark.parametrize('test_accession', ['test_020'], scope='class')
    def test_dihedrals(self, project: 'Project'):
        """Test dihedrals energy calculation."""
        project.get_dihedrals()
