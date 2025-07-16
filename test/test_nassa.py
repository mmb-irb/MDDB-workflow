import os
import pytest
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
    
    # Test Workflow
    # mwf nassa -hp -pdirs . -stru replica_1/structure.pdb -traj replica_1/trajectory.xtc -top topology.prmtop -all