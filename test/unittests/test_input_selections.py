"""Tests for the parsing of input atom selections (PBC, CG, UA and dummy atoms).

Every selection input may be empty (no atoms), a VMD selection, or the automatic flag to guess atoms.
The automatic flag may be followed by an additional VMD selection (e.g. auto + resname LIG).
"""
import pytest

from mddb_workflow.mwf import MD
from mddb_workflow.utils.auxiliar import InputError
from mddb_workflow.utils.selections import Selection
from mddb_workflow.utils.structures import Structure


def build_structure() -> Structure:
    """Build a structure with 2 waters (atoms 0-5) and a sodium ion (atom 6)."""
    atoms = [ ('O', 'HOH', 1), ('H1', 'HOH', 1), ('H2', 'HOH', 1), ('O', 'HOH', 2), ('H1', 'HOH', 2), ('H2', 'HOH', 2), ('NA', 'NA', 3) ]
    lines = []
    for serial, (atom_name, residue_name, residue_number) in enumerate(atoms, 1):
        x, y, z = [f'{coord:.3f}'.rjust(8) for coord in (float(serial), float(residue_number), 0.0)]
        lines.append((f'ATOM  {str(serial).rjust(5)}  {atom_name.ljust(3)} {residue_name.ljust(3)} A{str(residue_number).rjust(4)}'
            f'    {x}{y}{z}  1.00  0.00          {atom_name[0].rjust(2)}').ljust(80))
    return Structure.from_pdb('\n'.join(lines) + '\n')


def parse(selection_string, must_select: bool = False) -> list[int]:
    """Parse a selection string where the automatic guess is the first water."""
    selection = MD._parse_input_selection('test', 'TEST', selection_string, build_structure(),
        lambda: Selection([0, 1, 2]), must_select=must_select)
    return list(selection.atom_indices)


@pytest.mark.unit_int
@pytest.mark.parametrize('selection_string', [None, ''])
def test_empty_selection(selection_string):
    """An empty input selects no atoms."""
    assert parse(selection_string) == []


@pytest.mark.unit_int
def test_selection():
    """A VMD selection selects its atoms, and may select no atoms unless it is a must."""
    assert parse('resname NA') == [6]
    assert parse('resname LIG') == []
    with pytest.raises(InputError, match='TEST selection "resname LIG" selected no atoms'):
        parse('resname LIG', must_select=True)


@pytest.mark.unit_int
def test_automatic_selection():
    """The automatic flag guesses atoms, and an additional selection adds its atoms to the guess."""
    assert parse('auto') == [0, 1, 2]
    assert parse('auto + resname NA') == [0, 1, 2, 6]
    with pytest.raises(InputError, match='TEST selection "resname LIG" selected no atoms'):
        parse('auto + resname LIG')
