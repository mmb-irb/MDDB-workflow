import os
import pytest
from model_workflow.mwf import Project
from model_workflow.utils.constants import *
from model_workflow.utils.file import File
from unittest.mock import patch


def pytest_configure(config):
    config.addinivalue_line( "markers", "CI: tests related to continuous integration")
    config.addinivalue_line( "markers", "release: tests related to release processes")

@pytest.fixture(scope="session")
def test_data_dir():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

@pytest.fixture(scope="class")
def test_accession(request):
    """Get accession ID from parametrization or use default"""
    if hasattr(request, 'param'):
        return request.param
    return "A0001"

@pytest.fixture(scope="class")
def test_proj_dir(test_accession: str, test_data_dir: str):
    """Create a persistent directory for test data that remains across test runs"""
    # Create a project-specific test data directory
    test_proj_dir = os.path.join(test_data_dir, 'output', test_accession)
    # Create the directory if it doesn't exist
    os.makedirs(test_proj_dir, exist_ok=True)
    return test_proj_dir

@pytest.fixture(scope="class")
def project(test_proj_dir : str, test_accession: str):
    project = Project(directory=test_proj_dir, accession=test_accession, sample_trajectory=10)
    return project

def get_analysis_file(project: 'Project', analysis_type: str):
    """Download and provide the standard structure file"""
    output_path = os.path.join(project.directory, f"mda.{analysis_type}_REF.json")
    file_obj = File(output_path)
    # Only download if file doesn't exist yet
    if not file_obj.exists:
        project.remote.download_analysis_data(analysis_type,  file_obj)
    return file_obj