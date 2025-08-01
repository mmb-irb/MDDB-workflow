import os, shutil, pytest
from model_workflow.mwf import Project
from model_workflow.utils.constants import *


def pytest_configure(config):
    config.addinivalue_line( "markers", "CI: tests related to continuous integration")
    config.addinivalue_line( "markers", "release: tests related to release processes")

@pytest.fixture(scope="session")
def test_data_dir():
    test_data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    output_dir = os.path.join(test_data_dir, 'output')
    os.makedirs(output_dir, exist_ok=True)
    return test_data_dir

@pytest.fixture(scope="class")
def test_accession(request):
    """Get accession ID from parametrization or use default"""
    if hasattr(request, 'param'):
        return request.param
    return "A0001"

@pytest.fixture(scope="class")
def test_proj_dir(test_accession: str, test_data_dir: str):
    """Create a persistent directory for test data that remains across test runs"""
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
def project(test_proj_dir : str, test_accession: str):
    # Nassa fails with less than 8 frames
    n_frames = 8 if test_accession == "seq001-1" else 5
    project = Project(directory=test_proj_dir, accession=test_accession, sample_trajectory=n_frames)
    return project
