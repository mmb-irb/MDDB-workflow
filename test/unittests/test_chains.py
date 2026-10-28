"""Tests for chain names beyond the letters in the alphabet.

Chains are named in PDB style (A, B, ..., Z, AA, AB, ...) with no limit in the number of chains.
The PDB format supports only one character per chain, so the mmCIF format is used instead when needed.
"""
import pytest

from mddb_workflow.utils.structures import Structure, get_chain_name, get_chain_name_index
from mddb_workflow.utils.auxiliar import InputError, ToolError


def build_waters_pdb(residue_count: int) -> str:
    """Build a PDB with water residues, all in chain A."""
    lines = []
    serial = 0
    for residue_number in range(1, residue_count + 1):
        for atom_name in ['O', 'H1', 'H2']:
            serial += 1
            x, y, z = [f'{coord:.3f}'.rjust(8) for coord in (float(serial), float(residue_number), 0.0)]
            lines.append((f'ATOM  {str(serial).rjust(5)}  {atom_name.ljust(3)} HOH A{str(residue_number).rjust(4)}'
                f'    {x}{y}{z}  1.00  0.00           {atom_name[0]}').ljust(80))
    return '\n'.join(lines) + '\n'


def set_chain_per_residue(structure: Structure):
    """Set a new chain for every residue."""
    for residue in structure.residues:
        structure.set_selection_chain_name(residue.get_selection(), structure.get_available_chain_name())


@pytest.mark.unit_int
@pytest.mark.parametrize('index, name', [
    (0, 'A'), (25, 'Z'), (26, 'AA'), (51, 'AZ'), (52, 'BA'), (701, 'ZZ'), (702, 'AAA'),
])
def test_chain_names(index: int, name: str):
    """Chain names follow the PDB style sequence."""
    assert get_chain_name(index) == name
    assert get_chain_name_index(name) == index


@pytest.mark.unit_int
@pytest.mark.parametrize('name', ['a', '1', 'A1', ''])
def test_chain_names_out_of_sequence(name: str):
    """Names with lower case letters or numbers are not part of the sequence."""
    assert get_chain_name_index(name) is None


@pytest.mark.unit_int
def test_available_chain_names():
    """There is no limit in the number of chains and lower case letters are never used."""
    structure = Structure.from_pdb(build_waters_pdb(60))
    set_chain_per_residue(structure)
    chain_names = [chain.name for chain in structure.chains]
    # Note that the first residue gets chain B since chain A is in use when it is set
    assert chain_names == [get_chain_name(index) for index in range(1, 61)]
    assert structure.get_next_available_chain_name('Z') == 'BJ'
    assert structure.get_next_available_chain_name('a') == 'A'


@pytest.mark.unit_int
def test_repeated_chains_beyond_alphabet():
    """Repeated chains are renamed even when there are more chains than letters in the alphabet."""
    structure = Structure.from_pdb(build_waters_pdb(60))
    set_chain_per_residue(structure)
    structure.chains[1].name = structure.chains[0].name
    assert structure.check_repeated_chains(fix_chains=True)
    chain_names = [chain.name for chain in structure.chains]
    assert len(set(chain_names)) == 60


@pytest.mark.unit_int
def test_pdb_renames_long_chains():
    """Long chain names are renamed in the PDB when there are enough letters."""
    structure = Structure.from_pdb(build_waters_pdb(30))
    set_chain_per_residue(structure)
    pdb_chains = [line[21] for line in structure.generate_pdb().split('\n') if line.startswith('ATOM')]
    assert len(set(pdb_chains)) == 30


@pytest.mark.unit_int
def test_pdb_chain_limit():
    """The PDB can not be written with more chains than letters unless limitations are accepted."""
    structure = Structure.from_pdb(build_waters_pdb(60))
    set_chain_per_residue(structure)
    with pytest.raises(InputError):
        structure.generate_pdb()
    pdb_content = structure.generate_pdb(accept_limitations=True)
    pdb_chains = set(line[21] for line in pdb_content.split('\n') if line.startswith('ATOM'))
    assert pdb_chains == {'X'}


@pytest.mark.unit_int
def test_vmd_chain_selection_with_pdb():
    """VMD gets a PDB when elements are not standard, so selections by chain are not allowed."""
    structure = Structure.from_pdb(build_waters_pdb(3))
    assert structure.get_vmd_format() == 'cif'
    assert len(structure.select('chain A')) == 9
    for atom in structure.atoms:
        atom.element = 'Cg'
    assert structure.get_vmd_format() == 'pdb'
    assert len(structure.select('resname HOH')) == 9
    with pytest.raises(ToolError):
        structure.select('resname HOH and chain A')


@pytest.mark.unit_int
def test_mmcif_round_trip():
    """Writing and reading back a mmCIF keeps atoms, residues and long chain names."""
    structure = Structure.from_pdb(build_waters_pdb(60))
    set_chain_per_residue(structure)
    structure.residues[0].icode = 'B'
    reparsed = Structure.from_mmcif(structure.generate_mmcif())
    def summary(structure: Structure) -> list:
        return [(atom.name, atom.element, atom.coords, atom.residue.name, atom.residue.number,
            atom.residue.icode, atom.chain.name) for atom in structure.atoms]
    assert summary(reparsed) == summary(structure)


@pytest.mark.unit_int
def test_mmcif_blank_chain():
    """Blank chains are kept when writing and reading back a mmCIF."""
    structure = Structure.from_pdb(build_waters_pdb(3))
    structure.chains[0].name = ' '
    reparsed = Structure.from_mmcif(structure.generate_mmcif())
    assert [chain.name for chain in reparsed.chains] == [' ']


@pytest.mark.unit_int
def test_suspicious_chains(capsys):
    """A warning is shown when many single residue chains share the residue name."""
    structure = Structure.from_pdb(build_waters_pdb(60))
    structure.check_suspicious_chains()
    assert 'Chains may be wrongly assigned' not in capsys.readouterr().out
    set_chain_per_residue(structure)
    structure.check_suspicious_chains()
    output = capsys.readouterr().out
    assert 'There are 60 chains made of a single HOH residue' in output
    assert 'There are 60 chains for 60 residues' in output


MMCIF_WITH_LABEL_NOTATION = '''data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.type_symbol
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N ALA A 1 0.0 0.0 0.0 N 10 ALA P N 1
ATOM 2 CA ALA A 1 1.0 0.0 0.0 C 10 ALA P CA 1
HETATM 3 O HOH B . 5.0 0.0 0.0 O 101 HOH P O 1
HETATM 4 O HOH B . 9.0 0.0 0.0 O 102 HOH P O 1
ATOM 5 N ALA A 1 0.0 0.0 1.0 N 10 ALA P N 2
'''


@pytest.mark.unit_int
def test_mmcif_author_notation():
    """The author notation is read by default, whatever the order of the fields."""
    structure = Structure.from_mmcif(MMCIF_WITH_LABEL_NOTATION)
    assert [chain.name for chain in structure.chains] == ['P']
    assert [(residue.name, residue.number) for residue in structure.residues] == [('ALA', 10), ('HOH', 101), ('HOH', 102)]
    assert [atom.coords for atom in structure.atoms] == [(0, 0, 0), (1, 0, 0), (5, 0, 0), (9, 0, 0)]
    # In the label notation, residues out of polymers have no number so waters are joined
    structure = Structure.from_mmcif(MMCIF_WITH_LABEL_NOTATION, author_notation=False)
    assert [chain.name for chain in structure.chains] == ['A', 'B']
    assert [(residue.name, residue.number) for residue in structure.residues] == [('ALA', 1), ('HOH', 1)]
