import os
import pytest
import subprocess
from model_workflow.utils.constants import *
from model_workflow.utils.type_hints import *

@pytest.mark.release
class TestNassa:
    @pytest.fixture(scope="class")
    def test_accession(self):
        """Override the default accession for this test"""
        return "seq001-1"
    
    def test_helical(self, project: 'Project'):
        project.mds[0].run_helical_analysis(project)
    
    
    def test_nassa_helical(self, test_data_dir: str):
        os.chdir(os.path.join(test_data_dir, 'output'))
        command = "mwf nassa -hp -stru *.prmtop -traj *.xtc -top *.prmtop  -pdirs seq001-1 -mdir replica_1 -m"
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        assert result.returncode == 0
        # Add more assertions based on the expected output or side effects

    def test_nassa_create_config(self, test_data_dir: str):
        os.chdir(os.path.join(test_data_dir, 'output'))
        command = f"mwf nassa -w seq001-1/replica_1 -seq {test_data_dir}/input"
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        assert result.returncode == 0
        # Add more assertions based on the expected output or side effects

    def test_nassa_all(self, test_data_dir: str):
        os.chdir(os.path.join(test_data_dir, 'output'))
        command = "mwf nassa -c nassa.yml -all"
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        assert result.returncode == 0
        # Add more assertions based on the expected output or side effects