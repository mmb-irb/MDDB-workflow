import os
import pytest
import numpy as np
from conftest import get_analysis_file

from model_workflow.utils.constants import *
from model_workflow.utils.type_hints import *
from model_workflow.utils.auxiliar import load_json


@pytest.mark.release
class TestMWFRun:
    """Test full workflow for A0001: mwf run -proj A0001 -dir test/test_data/A0001.1"""

    @pytest.fixture(scope="class")
    def test_accession(self):
        """Override the default accession for this test"""
        return "A0001.1"
    
    def test_TMscores_analysis(self, project : 'Project'):
        """Test that RMSD analysis runs and produces expected output"""
        # Download the reference file and run the analysis
        analysis_file = get_analysis_file(project, 'tmscores')
        md : MD = project.mds[0]
        md.overwritables = {'tmscore'}
        md.run_tmscores_analysis()

        # Check that the output file was created
        output_file = f"{project.directory}/replica_1/{OUTPUT_TMSCORES_FILENAME}"
        assert os.path.exists(output_file), f"Output file '{output_file}' was not created"
        
        # Load the results
        results = load_json(output_file)
        reference = load_json(analysis_file.absolute_path)
        # Check that the results match the expected output
        now = np.array(results['data'][0]['values']).mean()
        ref_mn = np.array(reference['data'][0]['values']).mean()

        assert np.isclose(now, ref_mn, atol=0.1), f"Output values do not match expected values. Now: {now}, Ref: {ref_mn}"
