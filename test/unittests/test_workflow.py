import shutil
import pathlib
from contextlib import contextmanager
from unittest.mock import patch
import pytest

import mddb_workflow.utils.auxiliar as aux
from mddb_workflow.mwf import workflow
import mddb_workflow.mwf as mwf


# Set up paths
data_dir = pathlib.Path(__file__).parent.parent / 'data'
dummy_dir = data_dir / 'input/dummy'
test_fld = data_dir / 'output/test_workflow'


def setup_two_replica_project(traj1="raw_trajectory", traj2="raw_trajectory"):
    """Set up a test project with two replicas."""
    # Remove old test folder and create a new one
    shutil.rmtree(test_fld, ignore_errors=True)
    test_fld.mkdir(parents=True, exist_ok=True)

    # Copy necessary files
    shutil.copy(dummy_dir / "inputs.yaml", test_fld)
    shutil.copy(dummy_dir / "gromacs/ala_ala.tpr", test_fld)

    # Create two replica directories
    (test_fld / "replica_1").mkdir(exist_ok=True)
    (test_fld / "replica_2").mkdir(exist_ok=True)
    shutil.copy(dummy_dir / f"gromacs/{traj1}.xtc", test_fld / "replica_1")
    shutil.copy(dummy_dir / f"gromacs/{traj2}.xtc", test_fld / "replica_2")


@contextmanager
def run_workflow_with_mock_task(mock_task):
    """Context manager to run workflow with a mock task in the test folder."""
    with patch.dict(mwf.requestables, {'mock_task': mock_task}):
        with patch.dict(mwf.md_requestables, {'mock_task': mock_task}):
            workflow(
                project_parameters={
                    'input_topology_filepath': 'ala_ala.tpr',
                    'md_config': [['replica_1', 'raw_trajectory.xtc'],
                                  ['replica_2', 'raw_trajectory.xtc']],
                },
                working_directory=str(test_fld),
                include=['mock_task'],
                keep_going=True,
            )


class TestReplicaErrorHandling:
    """Test the error handling in the workflow for MD replicas."""

    def test_workflow_continues_after_replica_error(self, capsys):
        """Test that workflow continues processing other MDs when one fails."""
        setup_two_replica_project()

        # Nonlocal to allow modification in nested function
        call_count = [0]
        def mock_task(md):
            call_count[0] += 1
            if call_count[0] == 1:
                aux.fail("Simulated error in replica 1")
                raise aux.TestFailure("Simulated exception in replica 1")
            elif call_count[0] == 2:
                print("Replica 2 processed successfully")

        run_workflow_with_mock_task(mock_task)
        captured = capsys.readouterr()
        print(captured.out)
        assert "Simulated error in replica 1" in captured.out
        assert "Finished with errors in 1 MD replicas" in captured.out
        assert call_count[0] == 2, "Both replicas should have been processed"

    def test_workflow_collects_multiple_replica_errors(self, capsys):
        """Test that workflow collects errors from multiple failing replicas."""
        setup_two_replica_project()

        call_count = [0]
        def mock_task(md):
            call_count[0] += 1
            aux.fail(f"Simulated error in replica {call_count[0]}")
            raise aux.TestFailure(f"Simulated exception in replica {call_count[0]}")

        run_workflow_with_mock_task(mock_task)
        captured = capsys.readouterr()
        print(captured.out)
        assert "Simulated error in replica 1" in captured.out
        assert "Simulated error in replica 2" in captured.out
        assert "Finished with errors in 2 MD replicas" in captured.out
        assert call_count[0] == 2, "Both replicas should have been attempted"

    def test_workflow_succeeds_without_errors(self, capsys):
        """Test that workflow prints 'Done!' when no errors occur."""
        setup_two_replica_project()

        call_count = [0]
        def mock_task(md):
            call_count[0] += 1
            print(f"Replica {call_count[0]} processed successfully")

        run_workflow_with_mock_task(mock_task)
        captured = capsys.readouterr()
        print(captured.out)
        assert "Done!" in captured.out
        assert "Finished with errors" not in captured.out
        assert call_count[0] == 2, "Both replicas should have been processed"

    def test_workflow_raises_non_test_failure_exceptions(self, capsys):
        """Test that non-TestFailure exceptions are re-raised and not caught."""
        setup_two_replica_project()

        call_count = [0]
        def mock_task(md):
            call_count[0] += 1
            # Raise a regular exception, not a TestFailure
            raise RuntimeError("Unexpected runtime error")

        # The RuntimeError should propagate up and not be caught
        with pytest.raises(RuntimeError, match="Unexpected runtime error") as exc_info:
            run_workflow_with_mock_task(mock_task)

        # Print the caught exception for debugging
        print(f"Caught exception: {exc_info.value}")
        assert call_count[0] == 1, "Only first replica should have been attempted before exception"
