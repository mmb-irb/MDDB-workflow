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
    
    def test_structure_functions(self, project: 'Project'):
        """Test structure-related functions in the project"""
        structure = project.structure
        res = structure.residues[0]
        res.split([0,1],list(set(res.atom_indices))[2:])

        residue_atom_names = list(set([ atom.name for atom in res.atoms ]))
        res.split_by_atom_names(residue_atom_names[:2],residue_atom_names[2:])