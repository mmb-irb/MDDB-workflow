"""Tests for the residue splitting in Structure.from_pdb when PDB residue numbers overflow.

VMD writes residue numbers in decimal up to 9999, then in hexadecimal up to ffff,
and then a '****' placeholder for every further residue (e.g. MDBind systems with > 65535 residues).
Consecutive residues sharing the placeholder must still be parsed as separated residues.
"""
import pytest

from mddb_workflow.utils.structures import Structure, is_residue_number_placeholder, parse_hybrid_residue_numbers


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
            # Keep coordinates within the 8 PDB columns (below 1000) for big structures
            coordinate = float(serial % 1000)
            coords = (coordinate, coordinate * 0.5, -coordinate)
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
def test_repeated_chain_blocks_are_not_merged():
    """Atoms in a later block of the same chain letter go to a new residue, not to the earlier one with the same label."""
    residues = [
        ('A', '0', 'WAT', WAT),
        ('B', '0', 'WAT', WAT),
        ('A', '0', 'WAT', WAT),
        ('B', '0', 'WAT', WAT),
    ]
    structure = Structure.from_pdb(build_pdb(residues))
    assert [residue.atom_indices for residue in structure.residues] == [[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]]
    assert [chain.name for chain in structure.chains] == ['A', 'B', 'A', 'B']
    assert [chain.residue_indices for chain in structure.chains] == [[0], [1], [2], [3]]


@pytest.mark.unit_int
def test_same_number_written_differently_is_not_merged():
    """Only the last residue is reused, even if an earlier one has the same number (e.g. '  10' and ' 010')."""
    residues = [
        ('A', '  10', 'WAT', WAT),
        ('A', ' 010', 'WAT', WAT),
        ('A', '  10', 'WAT', WAT),
    ]
    structure = Structure.from_pdb(build_pdb(residues))
    assert [residue.number for residue in structure.residues] == [10, 10, 10]
    assert [residue.atom_indices for residue in structure.residues] == [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
    assert [chain.residue_indices for chain in structure.chains] == [[0, 1, 2]]


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


def vmd_residue_number(number: int) -> str:
    """Write a residue number like VMD: decimal up to 9999, then hexadecimal."""
    return str(number) if number <= 9999 else hex(number)[2:]


@pytest.mark.unit_int
def test_hybrid_decimal_hexadecimal_residue_numbers():
    """Decimal numbers up to 9999 followed by hexadecimal ones (VMD) keep their values, e.g. '10' is 10, '2710' is 10000."""
    numbers = list(range(1, 10051))
    residues = [('X', vmd_residue_number(number), 'ALA', ['CA']) for number in numbers]
    structure = Structure.from_pdb(build_pdb(residues))
    assert [residue.number for residue in structure.residues] == numbers
    assert structure.residue_numeration_base == 16


@pytest.mark.unit_int
def test_hybrid_residue_numbers_without_letters():
    """Up to 10009 residues the hexadecimal numbers have no letters ('2710' to '2719'), but they are still hexadecimal."""
    numbers = list(range(9990, 10006))
    residues = [('X', vmd_residue_number(number), 'ALA', ['CA']) for number in numbers]
    structure = Structure.from_pdb(build_pdb(residues))
    assert [residue.number for residue in structure.residues] == numbers


@pytest.mark.unit_int
def test_hybrid_residue_numbers_restart_in_new_chain():
    """A new chain which starts again at 1 is decimal again after a hexadecimal run."""
    residues = [
        ('X', '9999', 'ALA', ['CA']),
        ('X', '2710', 'ALA', ['CA']),
        ('X', '271a', 'ALA', ['CA']),
        ('Y', '1', 'ALA', ['CA']),
        ('Y', '10', 'ALA', ['CA']),
        ('Y', '100', 'ALA', ['CA']),
    ]
    structure = Structure.from_pdb(build_pdb(residues))
    assert [residue.number for residue in structure.residues] == [9999, 10000, 10010, 1, 10, 100]


@pytest.mark.unit_int
def test_fully_hexadecimal_residue_numbers():
    """Numbers with hexadecimal letters below 10000 mean the whole numeration is hexadecimal."""
    residues = [('X', number, 'ALA', ['CA']) for number in ['1', 'a', 'f', '10', '1f']]
    structure = Structure.from_pdb(build_pdb(residues))
    assert structure.residue_numeration_base == 16
    assert [residue.number for residue in structure.residues] == [1, 10, 15, 16, 31]


@pytest.mark.unit_int
def test_decimal_residue_numbers_going_back_stay_decimal():
    """A decimal numeration which goes back is not mistaken for a hexadecimal one."""
    residues = [
        ('A', '9500', 'ALA', ['CA']),
        ('A', '3000', 'LIG', ['C1']),  # 0x3000 = 12288, too far after 9500
        ('B', '9999', 'ALA', ['CA']),
        ('C', '2710', 'ALA', ['CA']),  # Another chain: decimal
    ]
    structure = Structure.from_pdb(build_pdb(residues))
    assert structure.residue_numeration_base == 10
    assert [residue.number for residue in structure.residues] == [9500, 3000, 9999, 2710]


@pytest.mark.unit_int
def test_hybrid_residue_numbers_round_trip():
    """generate_pdb writes numbers above 9999 in hexadecimal and reading them back gives the same numbers."""
    numbers = list(range(1, 10051))
    residues = [('X', vmd_residue_number(number), 'ALA', ['CA']) for number in numbers]
    structure = Structure.from_pdb(build_pdb(residues))
    reparsed = Structure.from_pdb(structure.generate_pdb())
    assert [residue.number for residue in reparsed.residues] == numbers


def chain_residues(chain: str, numbers: list[str]) -> list[tuple[str, str, str, list[str]]]:
    """One single-atom residue per residue number, in a chain."""
    return [(chain, number, 'ALA', ['CA']) for number in numbers]


def parsed_numbers(residues: list[tuple[str, str, str, list[str]]]) -> list[int]:
    return [residue.number for residue in Structure.from_pdb(build_pdb(residues)).residues]


@pytest.mark.unit_int
def test_hybrid_residue_numbers_starting_above_9999():
    """A chain starting above 9999 is hexadecimal from its first number, before any letter ('2711' is 10001)."""
    numbers = list(range(10001, 10040))
    residues = [('X', vmd_residue_number(number), 'ALA', ['CA']) for number in numbers]
    structure = Structure.from_pdb(build_pdb(residues))
    assert [residue.number for residue in structure.residues] == numbers
    assert structure.residue_numeration_base == 16
    # The same once written by generate_pdb and read back
    reparsed = Structure.from_pdb(structure.generate_pdb())
    assert [residue.number for residue in reparsed.residues] == numbers


@pytest.mark.unit_int
@pytest.mark.parametrize('numbers', [['9500', '2750'], ['9998', '9999', '2800'], ['9999', '2711']])
def test_decimal_residue_numbers_going_back_to_2710_or_more_stay_decimal(numbers: list[str]):
    """With no hexadecimal letters, only '2710' right after 9999 is hexadecimal: other numbers going back are decimal."""
    structure = Structure.from_pdb(build_pdb(chain_residues('A', numbers)))
    assert [residue.number for residue in structure.residues] == [int(number) for number in numbers]
    assert structure.residue_numeration_base == 10


@pytest.mark.unit_int
@pytest.mark.parametrize('next_chain_numbers', [['9001', '9002'], ['3000', '3001'], ['2712', '2713']])
def test_hexadecimal_run_ends_at_chain_change(next_chain_numbers: list[str]):
    """A decimal chain after a hexadecimal run is decimal, even if its numbers could continue the run."""
    residues = chain_residues('X', ['9999', '2710', '2711']) + chain_residues('L', next_chain_numbers)
    assert parsed_numbers(residues) == [9999, 10000, 10001] + [int(number) for number in next_chain_numbers]


@pytest.mark.unit_int
def test_hexadecimal_run_at_the_start_of_a_new_chain():
    """A new chain whose numbers continue after 9999 is hexadecimal when it has hexadecimal letters further."""
    residues = chain_residues('A', ['9998', '9999']) + chain_residues('B', ['2710', '2711', '271a'])
    assert parsed_numbers(residues) == [9998, 9999, 10000, 10001, 10010]


@pytest.mark.unit_int
def test_hybrid_residue_numbers_with_gaps():
    """Numbers without letters just before hexadecimal ones are hexadecimal, also after a gap in the numeration."""
    # From 8999 to 10000: the run starts after the gap
    assert parsed_numbers(chain_residues('X', ['8999', '2710', '2711', '271a'])) == [8999, 10000, 10001, 10010]
    # From 10010 to 12288, more than the maximum gap: '3000' is decimal at first, but '300a' (12298) follows it
    residues = chain_residues('X', ['9999', '2710', '271a', '3000', '3001', '300a'])
    assert parsed_numbers(residues) == [9999, 10000, 10010, 12288, 12289, 12298]


@pytest.mark.unit_int
@pytest.mark.parametrize('texts, expected, any_hexadecimal', [
    ([], [], False),
    (['1', '2', '3'], [1, 2, 3], False),
    (['-2', '-1', '0', '1'], [-2, -1, 0, 1], False),
    (['9998', '9999', '0', '1'], [9998, 9999, 0, 1], False),  # Decimal numeration cycling after 9999
    (['9999', '2710', '2711'], [9999, 10000, 10001], True),
    (['2711', '271A', '271B'], [10001, 10010, 10011], True),  # Uppercase letters
    (['9999', '2710', '271a', '1000', '1001'], [9999, 10000, 10010, 1000, 1001], True),  # Going back to decimal
    ([' 2719', '271a '], [10009, 10010], True),  # Padding
    (['5', '-', '6'], [5, 6, 6], False),  # Placeholder: after the previous residue
    (['-', '--', '1'], [1, 2, 1], False),
    (['9999', '2710', '-', '271a'], [9999, 10000, 10001, 10010], True),
])
def test_parse_hybrid_residue_numbers(texts: list[str], expected: list[int], any_hexadecimal: bool):
    """Parse one chain of hybrid decimal and hexadecimal residue numbers."""
    assert parse_hybrid_residue_numbers(texts) == (expected, any_hexadecimal)


@pytest.mark.unit_int
def test_placeholder_residue_numbers_in_decimal_numeration():
    """A residue number made only of '-' is not a number: it follows the previous residue, as before the hybrid numeration."""
    structure = Structure.from_pdb(build_pdb(chain_residues('A', ['5', '-', '6'])))
    assert [residue.number for residue in structure.residues] == [5, 6, 6]
    assert structure.residue_numeration_base == 10
