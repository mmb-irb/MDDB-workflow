import os
import pytest
import model_workflow.utils.structures as structures
from model_workflow.utils.type_hints import *


@pytest.mark.CI
@pytest.mark.release
class TestStructure:
    """Test structure-related functionalities"""

    @pytest.mark.parametrize(
            "structure_file",
            ["A0001.pdb", "A0001.prmtop", "A0224.tpr", "3rbn.cif"])
    def test_load_structure(self, structure_file: str, test_data_dir: str):
        """Test loading a structure from a file"""
        structure_path = os.path.join(test_data_dir, 'input/raw_structures', structure_file)
        structure = structures.Structure.from_file(structure_path)
        assert structure is not None, "Failed to load the structure"
        assert hasattr(structure, 'atoms'), "Structure should have atoms attribute"
        assert len(structure.atoms) > 0, "Structure should contain atoms"

        # Check structure properties to see if something fails
        structure.check_incoherent_bonds(), "Structure should not have incoherent bonds"
        structure.check_repeated_atoms(fix_atoms=True, display_summary=True), "Structure should not have repeated atoms"
        structure.check_repeated_residues(fix_residues=True, display_summary=True), "Structure should not have repeated residues"
        structure.auto_chainer()
        
        # prmtop has no coords
        # MDAnalysis does not read coordinates from tpr:
        # https://userguide.mdanalysis.org/stable/formats/reference/tpr.html#tpr-gromacs-run-topology-files
        if not structure_file.endswith('.prmtop') and not structure_file.endswith('.tpr'):
            structure.filter('protein')
    
    @pytest.fixture(scope="function")
    def structure(self, test_data_dir: str):
        """Fixture to load a sample structure for testing."""
        structure_path = os.path.join(test_data_dir, 'input/raw_structures', "A0001.pdb")
        structure = structures.Structure.from_file(structure_path)
        return structure

    def test_structure_functions(self, structure: 'Structure'):
        """Test structure-related functions in the project"""
        print(structure)
        structure.display_summary()
        sel = structure.select('resname ARG')
        structure.name_selection(sel)
        structure.get_selection_classification(sel)
        structure.get_next_available_chain_name('B')
        structure.merge(structure)
        structure.raw_protein_chainer()
        structure.select_cartoon()
        structure.select_pbc_guess()

        # Make a repeated residue
        structure.residues.extend([structure.residues[-1],structure.residues[0]])
        structure.check_repeated_residues(fix_residues=True, display_summary=True)

        # Make a repeated chain
        structure.chains.append(structure.chains[0])
        structure.check_repeated_chains(fix_chains=True, display_summary=True)

    def test_chain_functions(self, structure: 'Structure'):
        """Test chain-related functions in the project"""
        chain: structures.Chain = structure.chains[0]
        print(chain)
        chain.set_residues(structure.residues[:5])
        chain.classification
        chain.set_classification('protein')
        chain.get_selection()
        chain.has_cg()

    def test_residue_functions(self, structure: 'Structure'):
        """Test residue-related functions in the project"""
        res: 'Residue' = structure.residues[0]
        print(res)
        res.split([0,1],list(set(res.atom_indices))[2:])
        residue_atom_names = list(set([ atom.name for atom in res.atoms ]))
        res.split_by_atom_names(residue_atom_names[:2],residue_atom_names[2:])
        res.get_formula()
        res.get_classification_by_name()
        res.is_bonded_with_residue(res)
        res.get_bonded_residues()
        res.set_atoms(structure.atoms[:5])

    def test_atom_functions(self, structure: 'Structure'):
        """Test atom-related functions in the project"""
        atoms: list['Atom'] = structure.atoms[:5]
        atom1, atom2, atom3, atom4 = atoms[:4]
        print(atom1)
        atom1.set_residue_index(42)
        atom1.set_residue(structure.residues[-1])
        atom1.get_selection()
        for atom in atoms:
            atom.is_carbohydrate_candidate()
        structures.calculate_angle(atom1, atom2, atom3)
        structures.calculate_torsion(atom1, atom2, atom3, atom4)
