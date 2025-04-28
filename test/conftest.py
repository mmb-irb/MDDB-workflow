import os
import pytest
from model_workflow.utils.constants import *
from model_workflow.utils.file import File
from model_workflow.utils.remote import Remote
from model_workflow.tools.topology_manager import setup_structure
from model_workflow.utils.auxiliar import load_json, load_yaml
from unittest.mock import patch
from io import StringIO

# Constants
DATABASE_URL = "https://irb-dev.mddbr.eu/api/"
TEST_DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_data")

@pytest.fixture(scope="class")
def test_accession():
    """Default accession ID for tests"""
    return "A01M9.1"  # Default value

@pytest.fixture(scope="class")
def test_data_dir(test_accession):
    """Create a persistent directory for test data that remains across test runs"""
    # Create a project-specific test data directory
    test_data_dir = os.path.join(TEST_DATA_ROOT, test_accession)
    # Create the directory if it doesn't exist
    os.makedirs(test_data_dir, exist_ok=True)
    return test_data_dir

@pytest.fixture(scope="class")
def remote_client(test_accession):
    """Create a Remote client for the test accession"""
    return Remote(database_url=DATABASE_URL, accession=test_accession)

@pytest.fixture(scope="class")
def structure_file(remote_client, test_data_dir, test_accession):
    """Download and provide the standard structure file"""
    output_path = os.path.join(test_data_dir, f"structure.pdb")
    file_obj = File(output_path)
    
    # Only download if file doesn't exist yet
    if not file_obj.exists:
        remote_client.download_standard_structure(file_obj)
    
    return file_obj


@pytest.fixture(scope="class")
def structure(structure_file):
    """Load the structure into a Structure object"""
    return setup_structure(structure_file.path)


@pytest.fixture(scope="class")
def topology_file(remote_client, test_data_dir):
    """Download and provide the standard topology file"""
    output_path = os.path.join(test_data_dir, STANDARD_TOPOLOGY_FILENAME)
    file_obj = File(output_path)
    
    # Only download if file doesn't exist yet
    if not file_obj.exists:
        remote_client.download_standard_topology(file_obj)
    
    return file_obj

@pytest.fixture(scope="class")
def trajectory_file(remote_client, test_data_dir):
    """Download and provide a trajectory file with a limited frame selection for testing"""
    # For testing, we only need a small subset of frames
    output_path = os.path.join(test_data_dir, TRAJECTORY_FILENAME)
    file_obj = File(output_path)
    
    # Only download if file doesn't exist yet
    if not file_obj.exists:
        # Download only 10 frames for faster testing
        remote_client.download_trajectory(
            file_obj,
            # frame_selection="1:10:1",  # First 10 frames
            format="xtc"
        )
    return file_obj

@pytest.fixture(scope="class")
def inputs_file(remote_client, test_data_dir):
    """Download and provide a trajectory file with a limited frame selection for testing"""
    # For testing, we only need a small subset of frames
    output_path = os.path.join(test_data_dir, DEFAULT_INPUTS_FILENAME)
    file_obj = File(output_path)
    
    # Only download if file doesn't exist yet
    if not file_obj.exists:
        # Download only 10 frames for faster testing
        remote_client.download_inputs_file(file_obj)
    return file_obj

@pytest.fixture(scope="class")
def analysis_file(remote_client, test_data_dir, analysis_type):
    """Download and provide the standard structure file"""
    output_path = os.path.join(test_data_dir, f"mda.{analysis_type}_REF.json")
    file_obj = File(output_path)
    
    # Only download if file doesn't exist yet
    if not file_obj.exists:
        remote_client.download_analysis_data(analysis_type,  file_obj)
    
    return file_obj

@pytest.fixture(scope="class")
def membrane_map(remote_client, test_data_dir):
    """Download and provide the standard structure file"""
    output_path = os.path.join(test_data_dir, MEMBRANE_MAPPING_FILENAME)
    file_obj = File(output_path)
    
    # Only download if file doesn't exist yet
    if not file_obj.exists:
        remote_client.download_analysis_data('mem-map',  file_obj)
    file_obj = load_json(file_obj.path)
    return file_obj


@pytest.fixture
def capture_stdout():
    """Capture stdout for testing console output"""
    buffer = StringIO()
    with patch('sys.stdout', buffer):
        yield buffer