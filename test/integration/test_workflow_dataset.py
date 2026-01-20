import shutil
import pathlib
from contextlib import contextmanager
from unittest.mock import patch
import pytest

import mddb_workflow.utils.auxiliar as aux
from mddb_workflow.mwf import workflow
from mddb_workflow.core.dataset import Dataset, State
import mddb_workflow.mwf as mwf


# Set up paths
data_dir = pathlib.Path(__file__).parent.parent / 'data'
dummy_dir = data_dir / 'input/dummy'
test_dir = data_dir / 'output/test_workflow_handler'
project_dir = test_dir / 'project'


def setup_test_project(num_replicas=2):
    """Set up a test project with specified number of replicas."""
    # Remove old test folder and create a new one
    shutil.rmtree(test_dir, ignore_errors=True)
    test_dir.mkdir(parents=True, exist_ok=True)
    project_dir.mkdir(parents=True, exist_ok=True)

    # Copy necessary files
    shutil.copy(dummy_dir / "inputs.yaml", project_dir)
    shutil.copy(dummy_dir / "gromacs/ala_ala.tpr", project_dir)

    # Create replica directories
    md_config = []
    for i in range(1, num_replicas + 1):
        replica_dir = project_dir / f"replica_{i}"
        replica_dir.mkdir(exist_ok=True)
        shutil.copy(dummy_dir / "gromacs/raw_trajectory.xtc", replica_dir)
        md_config.append([f'replica_{i}', 'raw_trajectory.xtc'])

    return md_config


@contextmanager
def run_workflow_with_mock_task(mock_task, dataset_path=None):
    """Context manager to run workflow with a mock task."""
    with patch.dict(mwf.requestables, {'mock_task': mock_task}):
        with patch.dict(mwf.md_requestables, {'mock_task': mock_task}):
            md_config = setup_test_project(num_replicas=2)
            print(md_config)
            workflow(
                project_parameters={
                    'input_topology_filepath': 'ala_ala.tpr',
                    'md_config': md_config,
                    'directory': str(project_dir),
                },
                include=['mock_task'],
                dataset_path=dataset_path,
            )
            yield


@pytest.mark.unit_int
class TestDatasetIntegration:
    """Test the Dataset integration with the workflow."""

    def test_entries_are_added_to_dataset(self):
        """Test that project and MD entries are added to the dataset."""
        # Create a temporary database
        db_path = pathlib.Path(test_dir) / 'test_dataset.db'

        try:
            def mock_task(md):
                pass  # Do nothing

            with run_workflow_with_mock_task(mock_task, dataset_path=db_path):
                pass

            # Check that entries were added
            ds = Dataset(db_path)

            # Should have 1 project entry
            assert len(ds.projects_table) == 1
            # Should have 2 MD entries
            assert len(ds.mds_table) == 2

            # Check project entry
            project_status = ds.get_status(project_dir)

            assert project_status is not None
            assert project_status['state'] == State.DONE.value
            assert project_status['message'] == 'Done!'
            assert project_status['num_mds'] == 2

            # Check MD entries
            md1_status = ds.get_status(project_dir/'replica_1')
            assert md1_status is not None
            assert 'replica_1' in md1_status['rel_path']
            assert md1_status['state'] == State.DONE.value
            assert md1_status['message'] == 'Done!'

            md2_status = ds.get_status(project_dir/'replica_2')
            assert md2_status is not None
            assert 'replica_2' in md2_status['rel_path']
            assert md2_status['state'] == State.DONE.value
            assert md2_status['message'] == 'Done!'

        finally:
            # Clean up
            db_path.unlink(missing_ok=True)

    def test_running_state_during_execution(self):
        """Test that state is set to RUNNING during task execution."""
        db_path = pathlib.Path(test_dir) / 'test_dataset.db'

        try:
            # Track states seen during execution
            states_seen = {'project': [], 'md1': [], 'md2': []}

            def mock_task(md: 'mwf.MD'):
                # Check database state during execution
                ds = Dataset(db_path)

                # Check MD state
                md_status = ds.get_status(md.directory)
                if md_status:
                    states_seen[f'md{1 if "replica_1" in md.directory else 2}'].append(md_status['state'])
                    # Should be RUNNING during execution
                    assert md_status['rel_path'] == ds._abs_to_rel(md.directory)
                    assert md_status['state'] == State.RUNNING.value
                    assert md_status['message'] == 'Running workflow...'

                # Also check project state on first MD
                if md.directory.endswith('replica_1'):
                    project_status = ds.get_status(project_dir)
                    if project_status:
                        # If we are running the MD tasks, the project tasks have finished
                        states_seen['project'].append(project_status['state'])
                        assert project_status['state'] == State.RUNNING.value

            with run_workflow_with_mock_task(mock_task, dataset_path=db_path):
                pass

            # Verify that we actually checked the RUNNING states
            assert State.RUNNING.value in states_seen['project']
            assert State.RUNNING.value in states_seen['md1']
            assert State.RUNNING.value in states_seen['md2']

            # After workflow completes, check final states are DONE
            ds = Dataset(db_path)

            project_status = ds.get_status(project_dir)
            assert project_status['state'] == State.DONE.value

            md1_status = ds.get_status(project_dir/'replica_1')
            assert md1_status['state'] == State.DONE.value

            md2_status = ds.get_status(project_dir/'replica_2')
            assert md2_status['state'] == State.DONE.value

        finally:
            db_path.unlink(missing_ok=True)

    def test_state_transitions_on_error(self):
        """Test that states transition to ERROR when workflow fails."""
        db_path = pathlib.Path(test_dir) / 'test_dataset.db'

        try:
            call_count = [0]
            def mock_task(md):
                call_count[0] += 1
                if call_count[0] == 1:
                    raise aux.TestFailure("Simulated error in replica 1")

            with pytest.raises(aux.TestFailure):
                with run_workflow_with_mock_task(mock_task, dataset_path=db_path):
                    pass

            # Check final states
            ds = Dataset(db_path)
            # Check project entry
            project_status = ds.get_status(project_dir)

            assert project_status is not None
            assert project_status['state'] == State.ERROR.value

            # First MD should be ERROR
            md1_status = ds.get_status(project_dir/'replica_1')
            assert md1_status['state'] == State.ERROR.value
            assert 'TestFailure' in md1_status['message']
            assert 'Simulated error in replica 1' in md1_status['message']

            # Second MD should be DONE (keep_going=False)
            md2_status = ds.get_status(project_dir/'replica_2')
            assert md2_status is None

        finally:
            db_path.unlink(missing_ok=True)
