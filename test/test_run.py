import os
import pytest
import numpy as np
from model_workflow.mwf import Project, MD
from model_workflow.utils.file import File
from model_workflow.utils.auxiliar import load_json
from model_workflow.utils.constants import *

@pytest.mark.release
class TestMWFRun:
    @classmethod
    def setup_class(cls, test_data_dir):
        """Initialize shared resources for all tests in this class"""
        cls.test_data_dir = test_data_dir
        cls.project = Project(directory=test_data_dir, accession='A0001')
        cls.MD : MD = cls.project.mds[0]
        cls.MD.overwritables = {'tmscore'}

    @classmethod
    def get_analysis_file(cls, analysis_type: str):
        """Download and provide the standard structure file"""
        output_path = os.path.join(cls.project.directory, f"mda.{analysis_type}_REF.json")
        file_obj = File(output_path)
        # Only download if file doesn't exist yet
        if not file_obj.exists:
            cls.project.remote.download_analysis_data(analysis_type,  file_obj)
        return file_obj
    
    def test_TMscores_analysis(self):
        """Test that RMSD analysis runs and produces expected output"""
        # Download the reference file and run the analysis
        analysis_file = self.get_analysis_file('tmscores')
        self.MD.run_tmscores_analysis()

        # Check that the output file was created
        output_file = f"{self.test_data_dir}/replica_1/{OUTPUT_TMSCORES_FILENAME}"
        assert os.path.exists(output_file), f"Output file '{output_file}' was not created"
        
        # Load the results
        results = load_json(output_file)
        reference = load_json(analysis_file.absolute_path)
        # Check that the results match the expected output
        now = np.array(results['data'][0]['values']).mean()
        ref_mn = np.array(reference['data'][0]['values']).mean()

        assert np.isclose(now, ref_mn, atol=0.1), f"Output values do not match expected values. Now: {now}, Ref: {ref_mn}"