import pytest
import inspect
import model_workflow.utils.vmd_spells as vmd_spells
from model_workflow.utils.auxiliar import load_json, load_yaml
from model_workflow.tools.process_interactions import process_interactions, load_interactions
from model_workflow.tools.get_reduced_trajectory import get_reduced_trajectory


class TestVMD:
    @pytest.fixture(scope="session")
    def test_accession(self):
        """Override the default accession with a test-specific one"""
        return "A01JP.1"  # Different accession for this test file

    @pytest.fixture(scope="session")
    def analysis_type(self):
        return "interactions"

    # Test for escape_tcl_selection
    def test_escape_tcl_selection(self):
        # Test with string containing TCL reserved characters
        selection = 'name "CA" and [resid 1 to 10]'
        escaped = vmd_spells.escape_tcl_selection(selection)
        assert escaped == 'name \\"CA\\" and \\[resid 1 to 10\\]'
        
        # Test with string not containing TCL reserved characters
        selection = 'name CA and resid 1 to 10'
        escaped = vmd_spells.escape_tcl_selection(selection)
        assert escaped == selection

    # Test for get_covalent_bonds
    def test_get_covalent_bonds(self, structure_file, topology_file):
        # Call the function
        result = vmd_spells.get_covalent_bonds(structure_file.path)
        result = [sorted(bonds) for bonds in result]

        # Load the standard topology file and compare the results
        standard_topology = load_json(topology_file.path)
        js_bonds = standard_topology.get('atom_bonds', None)
        js_bonds = [sorted(bonds) for bonds in js_bonds]
        assert result == js_bonds


    # Test for get_interface_atom_indices and get_covalent_bonds_between
    def test_get_interface(self, analysis_file, structure, structure_file, trajectory_file, inputs_file):
        ref_inter = load_interactions (analysis_file, structure)

        inputs_file = load_yaml(inputs_file.path)
        # Get the default value using inspect
        signature = inspect.signature(process_interactions)
        distance_cutoff = signature.parameters['distance_cutoff'].default
        
        reduced_trajectory_filepath, step, frames = get_reduced_trajectory(
        structure_file,
        trajectory_file,
        1214,
        1000,
        )
        
        out_inter = vmd_spells.get_interface_atom_indices(
            structure_file.path,
            reduced_trajectory_filepath, # we should be using reduced_trajectory_filepath, but snapshots < 50
            inputs_file['interactions'][0]['selection_1'],
            inputs_file['interactions'][0]['selection_2'],
            distance_cutoff)
        

        residue_indices_1 = sorted(list(set([ structure.atoms[atom_index].residue_index for atom_index in out_inter['selection_1_interface_atom_indices'] ])))
        residue_indices_2 = sorted(list(set([ structure.atoms[atom_index].residue_index for atom_index in out_inter['selection_2_interface_atom_indices'] ])))
        
        assert residue_indices_1 == ref_inter[0]['interface_indices_1']
        assert residue_indices_2 == ref_inter[0]['interface_indices_2']


        strong_bonds = vmd_spells.get_covalent_bonds_between(structure_file.path, 
                                                             inputs_file['interactions'][0]['selection_1'],
                                                             inputs_file['interactions'][0]['selection_2'])
        
        assert ref_inter[0]['strong_bonds'] == [list(bond) for bond in strong_bonds]
