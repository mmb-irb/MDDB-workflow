import os
import pytest
import mddb_workflow.utils.structures as structures
from mddb_workflow.utils.type_hints import *


@pytest.mark.CI
class TestStructure:
    """Test structure-related functionalities."""

    @pytest.mark.parametrize(
            "structure_file",
            ["A0001.pdb", "A0001.prmtop", "A0224.tpr", "3rbn.cif"])
    def test_load_structure(self, structure_file: str, test_data_dir: str):
        """Test loading a structure from a file."""
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
        """Test structure-related functions in the project."""
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
        structure.residues.extend([structure.residues[-1], structure.residues[0]])
        structure.check_repeated_residues(fix_residues=True, display_summary=True)

        # Make a repeated chain
        structure.chains.append(structure.chains[0])
        structure.check_repeated_chains(fix_chains=True, display_summary=True)

    def test_chain_functions(self, structure: 'Structure'):
        """Test chain-related functions in the project."""
        chain: structures.Chain = structure.chains[0]
        print(chain)
        chain.set_residues(structure.residues)
        chain.classification
        chain.set_classification('protein')
        chain.get_selection()
        chain.has_cg()

    def test_residue_functions(self, structure: 'Structure'):
        """Test residue-related functions in the project."""
        res: 'Residue' = structure.residues[0]
        print(res)
        res.split([0, 1], list(set(res.atom_indices))[2:])
        residue_atom_names = list(set([atom.name for atom in res.atoms]))
        res.split_by_atom_names(residue_atom_names[:2], residue_atom_names[2:])
        res.get_formula()
        res.get_classification_by_name()
        res.is_bonded_with_residue(res)
        res.get_bonded_residues()
        res.set_atoms(structure.atoms[:5])

    def test_atom_functions(self, structure: 'Structure'):
        """Test atom-related functions in the project."""
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


def cg_structure(residues: list[tuple[str, str, str]]) -> 'Structure':
    """Build a coarse grained structure from (chain, residue name, bead names) tuples."""
    from mddb_workflow.utils.constants import CG_ATOM_ELEMENT
    lines, index = [], 1
    for number, (chain, resname, beads) in enumerate(residues, 1):
        for bead in beads.split():
            lines.append(f'ATOM  {index:5d} {bead:<4s} {resname:<4s}{chain}{number:4d}    {0:8.3f}{0:8.3f}{0:8.3f}  1.00  0.00')
            index += 1
    structure = structures.Structure.from_pdb('\n'.join(lines))
    for atom in structure.atoms:
        atom.element = CG_ATOM_ELEMENT
    return structure


class TestCoarseGrainClassification:
    """Coarse grained residues are classified by their residue names, and lipids also by their bead names."""

    @pytest.fixture(scope='class')
    def structure(self) -> 'Structure':
        return cg_structure([
            # Lipids whose residue names are not in the lipid residue names
            ('M', 'DIPE', 'NH3 PO4 GL1 GL2 C1A D2A D3A C4A C1B D2B D3B C4B'),
            ('M', 'ERGO', 'ROH R1 R2 R3 R4 R5 C1 C2'),
            ('M', 'BENZ', 'R1 R2 R3'),
            ('M', 'W', 'W'),
            # Nucleic residue names used for both DNA and RNA, in chains with RNA, DNA and no other names
            # The last one stays nucleic since the structure has both DNA and RNA names
            ('R', 'A', 'BB BB1 BB2 SC1 SC2 SC3 SC4'),
            ('R', 'U', 'BB BB1 BB2 SC1 SC2 SC3'),
            ('D', 'DT', 'BB1 BB2 BB3 SC1 SC2 SC3 SC4'),
            ('D', 'G', 'BB1 BB2 BB3 SC1 SC2 SC3 SC4'),
            ('N', 'C', 'BB BB1 BB2 SC1 SC2 SC3'),
        ])

    def test_classification(self, structure: 'Structure'):
        assert [residue.classification for residue in structure.residues] == \
            ['fatty', 'steroid', 'unknown', 'solvent', 'rna', 'rna', 'dna', 'dna', 'nucleic']

    @pytest.mark.parametrize('specific_name, classification', [('U', 'rna'), ('DT', 'dna')])
    def test_nucleic_classification_by_structure(self, specific_name: str, classification: str):
        """Nucleic residue names used for both DNA and RNA, in chains split from their strand with no other names,
        are classified by the names of the whole structure (e.g. coarse grained siRNA strands whose chains were assigned again).
        """
        structure = cg_structure([
            ('A', specific_name, 'BB BB1 BB2 SC1 SC2 SC3'),
            ('B', 'C', 'BB BB1 BB2 SC1 SC2 SC3'),
            ('C', 'G', 'BB BB1 BB2 SC1 SC2 SC3 SC4'),
            ('C', 'A', 'BB BB1 BB2 SC1 SC2 SC3 SC4'),
        ])
        assert [residue.classification for residue in structure.residues] == [classification] * 4

    def test_lipids(self, structure: 'Structure'):
        assert [structure.residues[index].name for index in structure.get_selection_residue_indices(structure.select_lipids())] == \
            ['DIPE', 'ERGO']

    def test_dummy_atoms(self):
        """Dummy atom names do not include beads of unsaturated tails, like D2A."""
        structure = cg_structure([('A', 'DUM', 'DM1 D2A C3A')])
        assert [atom.guess_element() for atom in structure.atoms] == ['Dm', 'Cg', 'Cg']
