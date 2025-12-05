import os
import shutil
import pathlib
from contextlib import contextmanager
from unittest.mock import patch
import pytest

from mddb_workflow.utils.auxiliar import fail, TestFailure
from mddb_workflow.mwf import workflow
import mddb_workflow.mwf as mwf


# Set up paths
data_dir = pathlib.Path(__file__).parent.parent / 'data'
dummy_dir = data_dir / 'input/dummy'
test_fld = data_dir / 'output/test_workflow'


def regenerate_test_fld():
    """Remove old test folder and create a new one."""
    shutil.rmtree(test_fld, ignore_errors=True)
    test_fld.mkdir(parents=True, exist_ok=True)


def setup_two_replica_project():
    """Set up a test project with two replicas."""
    regenerate_test_fld()

    # Copy necessary files
    shutil.copy(dummy_dir / "inputs.yaml", test_fld)
    shutil.copy(dummy_dir / "gromacs/ala_ala.tpr", test_fld)

    # Create two replica directories
    (test_fld / "replica_1").mkdir(exist_ok=True)
    (test_fld / "replica_2").mkdir(exist_ok=True)
    shutil.copy(dummy_dir / "gromacs/trajectory.xtc", test_fld / "replica_1")
    shutil.copy(dummy_dir / "gromacs/trajectory.xtc", test_fld / "replica_2")


@contextmanager
def run_workflow_with_mock_task(mock_task):
    """Context manager to run workflow with a mock task in the test folder."""
    cwd = pathlib.Path.cwd()
    os.chdir(test_fld)

    try:
        with patch.dict(mwf.requestables, {'mock_task': mock_task}):
            with patch.dict(mwf.md_requestables, {'mock_task': mock_task}):
                workflow(
                    project_parameters={
                        'input_topology_filepath': 'ala_ala.tpr',
                        'md_config': [['replica_1', 'trajectory.xtc'],
                                      ['replica_2', 'trajectory.xtc']],
                    },
                    include=['mock_task'],
                    keep_going=True,
                )
        yield
    finally:
        os.chdir(cwd)


class TestReplicaErrorHandling:
    """Test the error handling in the workflow for MD replicas."""

    def test_workflow_continues_after_replica_error(self, capsys):
        """Test that workflow continues processing other MDs when one fails."""
        setup_two_replica_project()

        call_count = [0]

        def mock_task(md):
            call_count[0] += 1
            if call_count[0] == 1:
                fail("Simulated error in replica 1")
                raise TestFailure("Simulated exception in replica 1")

        with run_workflow_with_mock_task(mock_task):
            pass

        assert call_count[0] == 2, "Both replicas should have been processed"

        captured = capsys.readouterr()
        print(captured.out)
        assert "Simulated error in replica 1" in captured.out
        assert "Finished with errors in 1 MD replicas" in captured.out

    def test_workflow_collects_multiple_replica_errors(self, capsys):
        """Test that workflow collects errors from multiple failing replicas."""
        setup_two_replica_project()

        call_count = [0]

        def mock_task(md):
            call_count[0] += 1
            fail(f"Simulated error in replica {call_count[0]}")
            raise TestFailure(f"Simulated exception in replica {call_count[0]}")

        with run_workflow_with_mock_task(mock_task):
            pass

        assert call_count[0] == 2, "Both replicas should have been attempted"

        captured = capsys.readouterr()
        print(captured.out)
        assert "Simulated error in replica 1" in captured.out
        assert "Simulated error in replica 2" in captured.out
        assert "Finished with errors in 2 MD replicas" in captured.out

    def test_workflow_succeeds_without_errors(self, capsys):
        """Test that workflow prints 'Done!' when no errors occur."""
        setup_two_replica_project()

        call_count = [0]

        def mock_task(md):
            call_count[0] += 1

        with run_workflow_with_mock_task(mock_task):
            pass

        assert call_count[0] == 2, "Both replicas should have been processed"

        captured = capsys.readouterr()
        print(captured.out)
        assert "Done!" in captured.out
        assert "Finished with errors" not in captured.out

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
            with run_workflow_with_mock_task(mock_task):
                pass

        # Print the caught exception for debugging
        print(f"Caught exception: {exc_info.value}")

        # Only the first replica should have been attempted before the exception propagated
        assert call_count[0] == 1, "Only first replica should have been attempted before exception"
