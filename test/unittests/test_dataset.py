import os
import tempfile
from mddb_workflow.core.dataset import Dataset


def test_adds_projects():
    """Test that adding projects works correctly."""
    # Create a temporary directory to act as the dataset root
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create some fake project directories
        project_dirs = [os.path.join(tmpdir, f"proj{i}") for i in range(3)]
        project_dirs += [os.path.join(tmpdir, f"proj_glob{i}") for i in range(3)]
        project_dirs += [os.path.join(tmpdir, "wrong")]
        for d in project_dirs:
            os.makedirs(d)
        db_path = os.path.join(tmpdir, "dataset.db")
        ds = Dataset(dataset_path=db_path)
        ds.add_projects(project_dirs, ignore_dirs=[os.path.join(tmpdir, "proj_glob2")], verbose=True)
        for d in project_dirs:
            rel_path = os.path.relpath(d, tmpdir)
            status = ds.get_status(rel_path)
            if "glob2" in d:
                assert status is None  # Ignored directory should not be added
            else:
                assert status is not None
                assert status['state'] == 'not_run'
                assert status['message'] == 'No information have been recorded yet.'
        # Now remove the 'wrong' project
        ds.remove_projects([os.path.join(tmpdir, "wrong")], verbose=True)
        assert ds.get_status("wrong") is None  # 'wrong' project should be removed
