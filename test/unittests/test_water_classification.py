"""Water residues are classified as solvent whatever their name, also with virtual sites (e.g. TIP4P)."""
import pytest
from mddb_workflow.utils.structures import Structure


def build_residue(residue_name: str, atoms: list[tuple[str, str]], bonded_pairs: list[tuple[int, int]]) -> Structure:
    """Build a structure with one residue from atom names and elements, where an empty element is to be guessed."""
    lines = [f'ATOM  {index + 1:5d}  {name:<3} {residue_name:<4}A   1    {0:8.3f}{0:8.3f}{0:8.3f}  1.00  0.00          {element:>2}'
        for index, (name, element) in enumerate(atoms)]
    structure = Structure.from_pdb('\n'.join(lines) + '\n')
    bonds = [[] for _ in atoms]
    for first, second in bonded_pairs:
        bonds[first].append(second)
        bonds[second].append(first)
    structure.bonds = bonds
    return structure


@pytest.mark.parametrize('residue_name, atoms, bonded_pairs, classification', [
    # Water with a name which is not a standard solvent name
    ('wat3', [('OW', 'O'), ('HW1', 'H'), ('HW2', 'H')], [(0, 1), (0, 2)], 'solvent'),
    # TIP4P water whose virtual site (MW) has no element and no bonds
    ('wate', [('OW', 'O'), ('HW1', 'H'), ('HW2', 'H'), ('MW', '')], [(0, 1), (0, 2)], 'solvent'),
    # Ammonium is not water
    ('NH4', [('N', 'N'), ('H1', 'H'), ('H2', 'H'), ('H3', 'H'), ('H4', 'H')], [(0, 1), (0, 2), (0, 3), (0, 4)], 'other'),
])
def test_water_classification(residue_name, atoms, bonded_pairs, classification):
    structure = build_residue(residue_name, atoms, bonded_pairs)
    assert structure.residues[0].classification == classification
