import pytest
import pathlib
import sys
import os
from mddb_workflow.mwf import workflow
from mddb_workflow.console import main


# Set up paths
data_dir = pathlib.Path(__file__).parent.parent / 'data'
test_fld = data_dir / 'output/test_inputs'

@pytest.mark.CI
class TestInputs:
    """Test the inputs subcommand and its integration with the workflow."""

    def test_console(self):
        """Test that the inputs subcommand creates an inputs.yaml file."""
        # We change the directory because the inputs command expects to be run in the output directory
        cwd = os.getcwd()
        os.chdir(test_fld)
        sys.argv = ['mddb_workflow', 'inputs', '-ed', 'none']
        main()
        os.chdir(cwd)
        assert os.path.exists(f'{test_fld}/inputs.yaml')

    def test_empty_mds(self, setup_dummy_project):
        """Test that inputs mds are added when they are not defined in the inputs.yaml."""
        md_config = setup_dummy_project(test_fld)
        workflow(
            project_parameters={
                'input_topology_filepath': 'ala_ala.tpr',
                'input_md_config': md_config,
                'directory': str(test_fld),
            },
            include=['mock_task'],
        )
