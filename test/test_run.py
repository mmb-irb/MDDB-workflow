import os
import pytest
import numpy as np
from conftest import get_analysis_file

from model_workflow.utils.constants import *
from model_workflow.utils.type_hints import *
from model_workflow.utils.auxiliar import load_json
from model_workflow.mwf import default_analyses, requestables

@pytest.mark.release
class TestMWFRun:
    """Test full workflow for A0001: mwf run -proj A0001 -dir test/test_data/A0001.1"""

    @pytest.fixture(scope="class")
    def test_accession(self):
        """Override the default accession for this test"""
        return "A0001.1"
    
    
    @pytest.mark.parametrize("analysis_name", default_analyses.keys())
    def test_analysis_execution(self, project: 'Project', analysis_name: str, capsys):
        """Test that each analysis runs without errors"""
        if analysis_name in ['energies', 'pockets','clusters','tmscores']:
            pytest.skip(f"Skipping analysis '{analysis_name}' for now.")
        
        md: MD = project.mds[0]
        md.overwritables = {analysis_name}
        try:
            requestables[analysis_name](md)
        except Exception as e:
            pytest.fail(f"Analysis '{analysis_name}' failed with error: {str(e)}")
        finally:
            # Capture stdout and stderr
            out, err = capsys.readouterr()
            # Write the output to a log file
            log_file = f"{project.directory}/logs/{analysis_name}_output.log"
            os.makedirs(os.path.dirname(log_file), exist_ok=True)
            with open(log_file, 'w') as f:
                f.write("STDOUT:\n")
                f.write(out)
                f.write("\nSTDERR:\n")
                f.write(err)

    def test_TMscores_analysis(self, project : 'Project'):
        """Test that RMSD analysis runs and produces expected output"""
        # Run the analysis
        analysis = 'tmscores'
        md : MD = project.mds[0]
        md.overwritables = {analysis}
        requestables[analysis](md)

        # Check that the output file was created
        output_file = f"{project.directory}/replica_1/{OUTPUT_TMSCORES_FILENAME}"
        assert os.path.exists(output_file), f"Output file '{output_file}' was not created"
        
        # Get the reference analysis file
        analysis_file = get_analysis_file(project, 'analysis')
        # Load the results
        results = load_json(output_file)
        reference = load_json(analysis_file.absolute_path)
        # Check that the results match the expected output
        now = np.array(results['data'][0]['values']).mean()
        ref_mn = np.array(reference['data'][0]['values']).mean()

        assert np.isclose(now, ref_mn, atol=0.1), f"Output values do not match expected values. Now: {now}, Ref: {ref_mn}"
