import pytest
import os
from model_workflow.utils.auxiliar import load_json
from model_workflow.analyses.density import density
import hashlib
import json


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
        # Compare results with reference data using checksums
        # Convert JSON to string in a deterministic way (sorted keys)
        result_str = json.dumps(results, sort_keys=True)
        reference_str = json.dumps(reference, sort_keys=True)

        # Compute checksums
        result_checksum = hashlib.md5(result_str.encode()).hexdigest()
        reference_checksum = hashlib.md5(reference_str.encode()).hexdigest()

        # Assert checksums match
        assert result_checksum == reference_checksum, (
            f"Results don't match reference. "
            f"Got checksum: {result_checksum}, expected: {reference_checksum}"
        )