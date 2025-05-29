import os
import pytest
import numpy as np
from model_workflow.mwf import Project, MD
from model_workflow.utils.file import File
from model_workflow.utils.auxiliar import load_json
from model_workflow.utils.constants import *

@pytest.mark.release
class TestMWFRun:
    """Test full workflow for A0001: mwf run -proj A0001 -dir test/test_data/A0001.1"""
    @pytest.fixture(scope="class")
    def test_accession(self):
        """Override the default accession for this test"""
        return "A0001.1"
    
    def get_analysis_file(self, project: Project, analysis_type: str):
        """Download and provide the standard structure file"""
        output_path = os.path.join(project.directory, f"mda.{analysis_type}_REF.json")
        file_obj = File(output_path)
        # Only download if file doesn't exist yet
        if not file_obj.exists:
            project.remote.download_analysis_data(analysis_type,  file_obj)
        return file_obj
    
    def test_TMscores_analysis(self, project):
        """Test that RMSD analysis runs and produces expected output"""
        # Download the reference file and run the analysis
        analysis_file = self.get_analysis_file(project, 'tmscores')
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