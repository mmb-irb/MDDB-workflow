import os
import pytest
import numpy as np
from model_workflow.mwf import Project
from model_workflow.utils.auxiliar import load_json
from model_workflow.utils.constants import OUTPUT_TMSCORES_FILENAME

class TestTMscores:
    @pytest.fixture(scope="class", autouse=True)
    def analysis_type(self):
        return "tmscores"
    
    @pytest.fixture(scope="class")
    def test_accession(self):
        """Override the default accession for this test"""
        return "A0001.1"
    
    def test_TMscores_analysis(self, analysis_file, test_data_dir):
        """Test that RMSD analysis runs and produces expected output"""

        os.chdir(test_data_dir) # Project brekas if we do not change
        proj = Project(directory=test_data_dir,
                       accession='A0001')
        proj.mds[0].overwritables = {'tmscore'}
        proj.mds[0].run_tmscores_analysis()

        # Check that the output file was created
        output_file = f"{test_data_dir}/replica_1/{OUTPUT_TMSCORES_FILENAME}"
        assert os.path.exists(output_file), f"Output file '{output_file}' was not created"
        
        # Load the results
        results = load_json(output_file)
        reference = load_json(analysis_file.absolute_path)
        # Check that the results match the expected output
        now = np.array(results['data'][0]['values']).mean()
        ref_mn = np.array(reference['data'][0]['values']).mean()

        assert np.isclose(now, ref_mn, atol=0.1), f"Output values do not match expected values. Now: {now}, Ref: {ref_mn}"