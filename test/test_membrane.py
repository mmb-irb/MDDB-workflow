import os
import pytest
from conftest import get_analysis_file

from model_workflow.utils.constants import *
from model_workflow.utils.type_hints import *
from model_workflow.utils.auxiliar import load_json


class TestDensityAnalysis:

    @pytest.fixture(scope="class")
    def test_accession(self):
        """Override the default accession for this test"""
        return "A01JP.1"
    
    @pytest.fixture(scope="class")
    def output_file(self, test_proj_dir):
        """Create an output file path for the density analysis results"""
        return os.path.join(test_proj_dir, "density_output.json")
    
    def test_density_analysis(self, project : 'Project'):
        """Test that density analysis runs and produces expected output"""
        # Download the reference file and run the analysis
        analysis_file = get_analysis_file(project, 'density')
        md : MD = project.mds[0]
        md.overwritables = {'density'}
        md.run_density_analysis()

        # Check that the output file was created
        output_file = f"{project.directory}/replica_1/{OUTPUT_DENSITY_FILENAME}"
        assert os.path.exists(output_file), f"Output file '{output_file}' was not created"
        
        # Load the results
        results = load_json(output_file)
        reference = load_json(analysis_file.path)

        # Compare structure and keys
        assert set(results.keys()) == set(reference.keys()), \
            f"Result keys don't match: {set(results.keys())} vs {set(reference.keys())}"
        
        assert set(results.get('data', {}).keys()) == set(reference.get('data', {}).keys()), \
            "Data subkeys don't match"
        
        # Check components exist
        assert 'comps' in results.get('data', {}), "Missing 'comps' key in results"
        assert 'z' in results.get('data', {}), "Missing 'z' key in results"

        # Check z array length matches
        assert len(results['data']['z']) == len(reference['data']['z']), \
            f"Z array length mismatch: {len(results['data']['z'])} vs {len(reference['data']['z'])}"
        
        # Check same number of components
        assert len(results['data']['comps']) == len(reference['data']['comps']), \
            f"Component count mismatch: {len(results['data']['comps'])} vs {len(reference['data']['comps'])}"