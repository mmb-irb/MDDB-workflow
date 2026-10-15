"""Tests for the residue splitting in Structure.from_pdb when PDB residue numbers overflow.

VMD writes residue numbers in decimal up to 9999, then in hexadecimal up to ffff,
and then a '****' placeholder for every further residue (e.g. MDBind systems with > 65535 residues).
Consecutive residues sharing the placeholder must still be parsed as separated residues.
"""
import pytest

from mddb_workflow.utils.structures import Structure, is_residue_number_placeholder


def pdb_line(serial: str, atom_name: str, residue_name: str, chain: str, residue_number: str,
             coords: tuple[float, float, float], element: str = '') -> str:
    """Build a PDB atom line respecting the PDB columns."""
    name = ' ' + atom_name.ljust(3) if len(atom_name) < 4 else atom_name
    x, y, z = [f'{coord:.3f}'.rjust(8) for coord in coords]
    return (f'ATOM  {serial.rjust(5)} {name} {residue_name.ljust(4)}{chain}{residue_number.rjust(4)}'
        f'    {x}{y}{z}  0.00  0.00          {element.rjust(2)}').ljust(80)


def build_pdb(residues: list[tuple[str, str, str, list[str]]]) -> str:
    """Build a PDB from a list of (chain, residue number, residue name, atom names)."""
    lines = []
    serial = 0
    for chain, residue_number, residue_name, atom_names in residues:
        for atom_name in atom_names:
            serial += 1
            # Like VMD, write atom serials above 99999 in hexadecimal
            serial_text = str(serial) if serial <= 99999 else hex(serial)[2:]
            coords = (float(serial), float(serial) * 0.5, -float(serial))
            lines.append(pdb_line(serial_text, atom_name, residue_name, chain, residue_number, coords))
    lines.append('END')
    return '\n'.join(lines) + '\n'


ALA = ['N', 'CA', 'C', 'O', 'CB']
WAT = ['O', 'H1', 'H2']

# Mimic the tail of a VMD converted.pdb: decimal and hexadecimal numbers followed by overflowed ones
OVERFLOW_RESIDUES = [
    ('X', '1', 'ALA', ALA),
    ('X', '2', 'Na+', ['Na+']),
    ('X', '3', 'WAT', WAT),
    ('X', 'fffe', 'WAT', WAT),
    ('X', 'ffff', 'WAT', WAT),
    ('X', '****', 'WAT', WAT),
    ('X', '****', 'WAT', WAT),
    ('X', '****', 'WAT', WAT),
    ('X', '****', 'Cl-', ['Cl-']),
    ('X', '****', 'Cl-', ['Cl-']),
    ('X', '****', 'WAT', WAT),
    ('X', '****', 'WAT', WAT),
]


@pytest.mark.unit_int
@pytest.mark.parametrize('residue_number, expected', [
    ('****', True), ('  **', True), ('9999', False), ('ffff', False), ('  -1', False), ('    ', False), ('A0B1', False),
])
def test_residue_number_placeholder(residue_number: str, expected: bool):
    """Only numbers made exclusively of non-alphanumerical characters are placeholders."""
    assert is_residue_number_placeholder(residue_number) == expected


@pytest.mark.unit_int
def test_overflowed_residue_numbers_are_split():
    """Every '****' residue becomes its own residue, numbered after the last one in the chain."""
    structure = Structure.from_pdb(build_pdb(OVERFLOW_RESIDUES))
    assert structure.atom_count == 5 + 1 + 8 * 3 + 2
    assert len(structure.residues) == len(OVERFLOW_RESIDUES)
    assert len(structure.chains) == 1
    # Residue names and atoms are conserved, one residue per input residue
    for residue, (_, _, residue_name, atom_names) in zip(structure.residues, OVERFLOW_RESIDUES):
        assert residue.name == residue_name
        assert [atom.name for atom in residue.atoms] == atom_names
    # Every overflowed water is a 3-atom residue
    overflowed_waters = [residue for residue in structure.residues[5:] if residue.name == 'WAT']
    assert len(overflowed_waters) == 5
    assert all(residue.atom_count == 3 for residue in overflowed_waters)
    # The placeholder makes the numeration be ignored, so residues are numbered sequentially in the chain
    # Overflowed residues continue from the last residue number and they are unique
    assert [residue.number for residue in structure.residues] == list(range(1, len(OVERFLOW_RESIDUES) + 1))
    # All waters are classified as solvent and thus selected
    water_atom_count = 8 * 3
    assert len(structure.select_water()) == water_atom_count
    assert all(residue.classification == 'solvent' for residue in structure.residues if residue.name == 'WAT')
    # Counter ions are also selected
    assert len(structure.select_water_and_counter_ions()) == water_atom_count + 1 + 2


@pytest.mark.unit_int
def test_overflowed_residue_numbers_round_trip():
    """Writing and reading back the parsed structure keeps the same residues."""
    structure = Structure.from_pdb(build_pdb(OVERFLOW_RESIDUES))
    reparsed = Structure.from_pdb(structure.generate_pdb())
    assert len(reparsed.residues) == len(structure.residues)
    assert [residue.name for residue in reparsed.residues] == [residue.name for residue in structure.residues]
    assert [residue.atom_count for residue in reparsed.residues] == [residue.atom_count for residue in structure.residues]
    assert len(reparsed.select_water()) == len(structure.select_water())


@pytest.mark.unit_int
def test_decimal_residue_numbers_unchanged():
    """A regular decimal PDB keeps its residue numbers and grouping."""
    residues = [
        ('A', '10', 'ALA', ALA),
        ('A', '11', 'ALA', ALA),
        ('A', '13', 'ALA', ALA),
        ('B', '1', 'WAT', WAT),
        ('B', '2', 'WAT', WAT),
        ('B', '2', 'WAT', ['O']),  # Same label as the previous residue: still merged as before
    ]
    structure = Structure.from_pdb(build_pdb(residues))
    assert structure.residue_numeration_base == 10
    assert [residue.number for residue in structure.residues] == [10, 11, 13, 1, 2]
    assert [residue.atom_count for residue in structure.residues] == [5, 5, 5, 3, 4]
    assert [chain.name for chain in structure.chains] == ['A', 'B']


@pytest.mark.unit_int
def test_hexadecimal_residue_numbers_unchanged():
    """A PDB with hexadecimal residue numbers and no placeholders is parsed in base 16."""
    residues = [
        ('X', 'fffe', 'WAT', WAT),
        ('X', 'ffff', 'WAT', WAT),
    ]
    structure = Structure.from_pdb(build_pdb(residues))
    assert structure.residue_numeration_base == 16
    assert [residue.number for residue in structure.residues] == [0xfffe, 0xffff]
