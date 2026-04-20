import pathlib
from contextlib import contextmanager
from unittest.mock import patch
import pytest

import mddb_workflow.utils.auxiliar as aux
from mddb_workflow.mwf import workflow
import mddb_workflow.mwf as mwf


# Set up paths
data_dir = pathlib.Path(__file__).parent.parent / 'data'
test_fld = data_dir / 'output/test_workflow'

@contextmanager
def run_workflow_with_mock_task(mock_task, keep_going=True):
    """Context manager to run workflow with a mock task in the test folder."""
    with patch.dict(mwf.requestables, {'mock_task': mock_task}):
        with patch.dict(mwf.md_requestables, {'mock_task': mock_task}):
            # mwf run -top ala_ala.tpr -md replica_1 raw_trajectory.xtc -md replica_2 raw_trajectory.xtc -dir test_fld -k
            workflow(
                project_parameters={
                    'input_topology_filepath': 'ala_ala.tpr',
                    'input_md_config': [['replica_1', 'raw_trajectory.xtc'],
                                  ['replica_2', 'raw_trajectory.xtc']],
                    'directory': str(test_fld),
                },
                include=['mock_task'],
                keep_going=keep_going,
            )


@pytest.mark.unit_int
class TestReplicaErrorHandling:
    """Test the error handling in the workflow for MD replicas."""

    def test_workflow_no_keep_going(self, capsys, setup_dummy_project):
        """Test that non-TestFailure exceptions are re-raised and not caught."""
        setup_dummy_project(test_fld)

        call_count = [0]
        def mock_task(md):
            call_count[0] += 1
            # Raise a regular exception, not a TestFailure
            raise RuntimeError("Unexpected runtime error")

        # The RuntimeError should propagate up and not be caught
        with pytest.raises(RuntimeError, match="Unexpected runtime error") as exc_info:
            run_workflow_with_mock_task(mock_task, keep_going=False)

        # Print the caught exception for debugging
        print(f"Caught exception: {exc_info.value}")
        assert call_count[0] == 1, "Only first replica should have been attempted before exception"

    def test_workflow_continues_after_replica_error(self, capsys, setup_dummy_project):
        """Test that workflow continues processing other MDs when one fails."""
        setup_dummy_project(test_fld)

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
        assert "Finished with errors in MD replicas" in captured.out
        assert call_count[0] == 2, "Both replicas should have been processed"

    def test_workflow_collects_multiple_replica_errors(self, capsys, setup_dummy_project):
        """Test that workflow collects errors from multiple failing replicas."""
        setup_dummy_project(test_fld)

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
        assert "Finished with errors in MD replicas" in captured.out
        assert call_count[0] == 2, "Both replicas should have been attempted"

    def test_workflow_succeeds_without_errors(self, capsys, setup_dummy_project):
        """Test that workflow prints 'Done!' when no errors occur."""
        setup_dummy_project(test_fld)

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
