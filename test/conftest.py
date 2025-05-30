import os
import pytest
from model_workflow.mwf import Project
from model_workflow.utils.constants import *
from model_workflow.utils.file import File
from unittest.mock import patch
from io import StringIO

# Constants
DATABASE_URL = "https://irb-dev.mddbr.eu/api/"
TEST_DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_data")

def pytest_configure(config):
    config.addinivalue_line( "markers", "CI: tests related to continuous integration")
    config.addinivalue_line( "markers", "release: tests related to release processes")

@pytest.fixture
def capture_stdout():
    """Capture stdout for testing console output"""
    buffer = StringIO()
    with patch('sys.stdout', buffer):
        yield buffer

# Fixtures for analysis type and test accession
@pytest.fixture(scope="class")
def test_accession():
    """Default accession ID for tests"""
    return "A01M9.1"  # Default value

@pytest.fixture(scope="class")
def test_data_dir(test_accession: str):
    """Create a persistent directory for test data that remains across test runs"""
    # Create a project-specific test data directory
    test_data_dir = os.path.join(TEST_DATA_ROOT, test_accession)
    # Create the directory if it doesn't exist
    os.makedirs(test_data_dir, exist_ok=True)
    return test_data_dir

@pytest.fixture(scope="class")
def project(test_data_dir : str, test_accession: str):
    project = Project(directory=test_data_dir, accession=test_accession)
    return project

def get_analysis_file(project: 'Project', analysis_type: str):
    """Download and provide the standard structure file"""
    output_path = os.path.join(project.directory, f"mda.{analysis_type}_REF.json")
    file_obj = File(output_path)
    # Only download if file doesn't exist yet
    if not file_obj.exists:
        project.remote.download_analysis_data(analysis_type,  file_obj)
    return file_obj