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


def test_guess_element_bonded_cs_is_carbon():
    """A ligand methyl carbon named CS (e.g. MDBind 1ppl, 2qtg) must not be guessed as caesium.
    A lone CS atom is still a caesium ion.
    """
    lines = []
    for serial, (name, x) in enumerate([('C1', 0.0), ('CS', 1.53)], 1):
        lines.append((f'HETATM{str(serial).rjust(5)}  {name.ljust(3)} MTH A   1    {x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00').ljust(80))
    lines.append((f'HETATM    3  CS   CS A   2    {30.0:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00').ljust(80))
    structure = Structure.from_pdb('\n'.join(lines) + '\n')
    assert structure.atoms[1].guess_element() == 'C'
    assert structure.atoms[2].guess_element() == 'Cs'
