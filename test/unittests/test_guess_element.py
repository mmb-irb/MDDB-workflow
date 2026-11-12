"""Tests for the atom element guess from atom names and bonds."""
from mddb_workflow.utils.structures import Structure


def build_tail(ca2_coords: tuple[float, float, float]) -> Structure:
    """Build a united atom lipid tail end (C1-CA1-CA2) with the given CA2 coordinates."""
    atoms = [('C1', (0.0, 0.0, 0.0)), ('CA1', (1.53, 0.0, 0.0)), ('CA2', ca2_coords)]
    lines = []
    for serial, (name, coords) in enumerate(atoms, 1):
        x, y, z = [f'{coord:.3f}'.rjust(8) for coord in coords]
        lines.append((f'ATOM  {str(serial).rjust(5)}  {name.ljust(3)} POP A   1    {x}{y}{z}  1.00  0.00           C').ljust(80))
    return Structure.from_pdb('\n'.join(lines) + '\n')


def test_guess_element_after_new_coordinates():
    """A CA2 split by the periodic boundary looks like a calcium ion.
    Once coordinates are imaged the bonds must be guessed again and CA2 is a carbon.
    """
    structure = build_tail((1.53, 30.0, 0.0))
    ca2 = structure.atoms[2]
    assert ca2.guess_element() == 'Ca'
    structure.set_new_coordinates([(0.0, 0.0, 0.0), (1.53, 0.0, 0.0), (3.06, 0.0, 0.0)])
    assert ca2.guess_element() == 'C'
