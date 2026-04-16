import os, shutil, pytest
from pathlib import Path
from mddb_workflow.mwf import Project
from mddb_workflow.utils.constants import *


def pytest_addoption(parser):
    parser.addoption(
        "--all",
        action="store_true",
        default=False,
        help="Run all marked tests (CI, unit_int, release)"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "CI: tests related to continuous integration")
    config.addinivalue_line("markers", "release: tests related to release processes")
    config.addinivalue_line("markers", "unit_int: unit and integration tests")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--all"):
        # When --all-marks is used, select tests with any of these marks
        mark_names = {"CI", "unit_int", "release"}
        selected = []
        for item in items:
            item_marks = {mark.name for mark in item.iter_markers()}
            if item_marks & mark_names:  # If any mark matches
                selected.append(item)
        items[:] = selected


@pytest.fixture(scope="session")
def test_data_dir():
    test_data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    output_dir = os.path.join(test_data_dir, 'output')
    os.makedirs(output_dir, exist_ok=True)
    return test_data_dir


@pytest.fixture(scope="class")
def test_accession(request):
    """Get accession ID from parametrization or use default."""
    if hasattr(request, 'param'):
        return request.param
    return "A0001"


@pytest.fixture(scope="class")
def test_proj_dir(test_accession: str, test_data_dir: str):
    """Create a persistent directory for test data that remains across test runs."""
    # Reset the current working to avoid issues with workflow changes
    os.chdir(test_data_dir)
    # Create a project-specific test data directory
    test_proj_dir = os.path.join(test_data_dir, 'output', test_accession)
    # Remove the directory if it already exists to ensure a clean state
    if os.path.exists(test_proj_dir):
        shutil.rmtree(test_proj_dir)
    # Create the directory if it doesn't exist
    os.makedirs(test_proj_dir, exist_ok=True)
    return test_proj_dir


@pytest.fixture(scope="class")
def project(test_proj_dir: str, test_accession: str):
    # Nassa fails with less than 8 frames
    n_frames = 8 if test_accession == "seq001-1" else 5
    project = Project(directory=test_proj_dir, accession=test_accession, sample_trajectory=n_frames)
    return project


@pytest.fixture
def setup_dummy_project(test_data_dir):
    """Factory fixture to create a two-replica project layout for tests."""
    def _setup(test_fld: Path, traj="raw_trajectory", n_replicas=2, copy_inputs=True):
        dummy_dir = Path(test_data_dir) / 'input/dummy'
        # Remove old test folder and create a new one
        shutil.rmtree(test_fld, ignore_errors=True)
        test_fld.mkdir(parents=True, exist_ok=True)
        # Copy required input files
        if copy_inputs:
            shutil.copy(dummy_dir / "inputs.yaml", test_fld)
        shutil.copy(dummy_dir / "gromacs/ala_ala.tpr", test_fld)

        # Create two replica directories
        md_config = []
        for i in range(1, n_replicas + 1):
            (test_fld / f"replica_{i}").mkdir(exist_ok=True)
            shutil.copy(dummy_dir / f"gromacs/{traj}.xtc", test_fld / f"replica_{i}")
            md_config.append([f"replica_{i}", "raw_trajectory.xtc"])
        return md_config

    return _setup
