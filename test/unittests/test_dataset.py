import os
import tempfile
from mddb_workflow.core.dataset import Dataset


def test_add_remove_entries():
    """Test that adding projects works correctly."""
    # Create a temporary directory to act as the dataset root
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create some fake project directories
        project_dirs = [os.path.join(tmpdir, f"proj{i}") for i in range(3)]
        # Add 2 replicas to the first project and 3 to the second
        replica_dirs = [os.path.join(project_dirs[0], f"replica{i}") for i in range(2)]
        replica_dirs += [os.path.join(project_dirs[1], f"replica{i}") for i in range(3)]
        project_dirs += [os.path.join(tmpdir, f"proj_glob{i}") for i in range(3)]
        project_dirs += [os.path.join(tmpdir, "wrong")]
        for d in project_dirs + replica_dirs:
            os.makedirs(d)
        db_path = os.path.join(tmpdir, "dataset.db")
        ds = Dataset(dataset_path=db_path)
        ds.add_entries(project_dirs,
                       ignore_dirs=[os.path.join(tmpdir, "proj_glob2")],
                       md_dirs="replica*",
                       verbose=True)
        for d in project_dirs:
            rel_path = os.path.relpath(d, tmpdir)
            status = ds.get_status(rel_path)
            if "glob2" in d:
                assert status is None  # Ignored directory should not be added
            else:
                assert status is not None
                assert status['state'] == 'not_run'
                assert status['message'] == 'No information have been recorded yet.'
            if "proj0" in d:
                assert status['num_mds'] == 2, f"Expected 2 MDs for {d}, got {status['num_mds']}"
            elif "proj1" in d:
                assert status['num_mds'] == 3
        # Now remove the 'wrong' project
        ds.remove_entries([os.path.join(tmpdir, "wrong")], verbose=True)
        assert ds.get_status("wrong") is None  # 'wrong' project should be removed
        # Remove specific replicas from proj0
        ds.remove_entries([os.path.join(tmpdir, "proj0")], md_dirs=['replica0'], verbose=True)
        assert ds.get_status("proj0")['num_mds'] == 1
        ds.remove_entries([os.path.join(tmpdir, "proj1")], md_dirs=['replica0', 'replica1'], verbose=True)
        assert ds.get_status("proj1")['num_mds'] == 1  # One replica should remain
