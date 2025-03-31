import pytest
import os
from model_workflow.utils.auxiliar import load_json
from model_workflow.analyses.density import density


class TestDensityAnalysis:
    @pytest.fixture(scope="session")
    def analysis_type(self):
        return "density"
    
    @pytest.fixture(scope="function")
    def output_file(self, test_data_dir):
        """Create an output file path for the density analysis results"""
        return os.path.join(test_data_dir, "density_output.json")
    
    def test_density_analysis(self, structure_file, trajectory_file, structure, 
                             output_file, membrane_map, analysis_file):
        """Test that density analysis runs and produces expected output"""
        # Run the density analysis
        density(
            input_structure_filepath=structure_file.path,
            input_trajectory_filepath=trajectory_file.path,
            output_analysis_filepath=output_file,
            membrane_map=membrane_map,
            structure=structure,
            snapshots=10001,
        )
        
        # Check that the output file was created
        assert os.path.exists(output_file), "Output file was not created"
        
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