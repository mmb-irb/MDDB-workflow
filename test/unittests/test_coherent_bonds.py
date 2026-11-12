"""Tests for the coherent bonds test with united atoms.

United atom residues (e.g. GROMOS, Berger lipids) have the hydrogens of their carbons absorbed.
Thus a terminal methyl is a single carbon with 1 bond, which is wrong for an atomistic carbon.
However these systems usually have explicit hydrogens elsewhere (e.g. water, polar hydrogens).
"""
import pytest

from mddb_workflow.utils.structures import Structure
from mddb_workflow.utils.structures import coherent_bonds_with_hydrogen, coherent_bonds_united_atom


# Every atom is a residue name, an atom name, an element and the names of the atoms it is bonded to in the residue
# A head with an explicit ammonium (NH3+) and a chain ending in a united atom methyl (C3)
UNITED_ATOM_LIPID = [
    ('N', 'N', ['H1', 'H2', 'H3', 'C1']),
    ('H1', 'H', ['N']), ('H2', 'H', ['N']), ('H3', 'H', ['N']),
    ('C1', 'C', ['N', 'C2']),
    ('C2', 'C', ['C1', 'C3']),
    ('C3', 'C', ['C2']),
]
WATER = [
    ('O', 'O', ['H1', 'H2']),
    ('H1', 'H', ['O']), ('H2', 'H', ['O']),
]
# An atomistic ethane which lost the hydrogens of one carbon (C2)
BROKEN_ETHANE = [
    ('C1', 'C', ['H1', 'H2', 'H3', 'C2']),
    ('H1', 'H', ['C1']), ('H2', 'H', ['C1']), ('H3', 'H', ['C1']),
    ('C2', 'C', ['C1']),
]


def build_structure(residues: list[tuple[str, list]]) -> Structure:
    """Build a structure with the given residues and set the bonds between their atoms."""
    lines = []
    atom_names = []
    for residue_number, (residue_name, atoms) in enumerate(residues, 1):
        for atom_name, element, _ in atoms:
            serial = len(atom_names) + 1
            x, y, z = [f'{coord:.3f}'.rjust(8) for coord in (float(serial), float(residue_number), 0.0)]
            lines.append((f'ATOM  {str(serial).rjust(5)}  {atom_name.ljust(3)} {residue_name} A{str(residue_number).rjust(4)}'
                f'    {x}{y}{z}  1.00  0.00          {element.rjust(2)}').ljust(80))
            atom_names.append((residue_number, atom_name))
    structure = Structure.from_pdb('\n'.join(lines) + '\n')
    bonds = []
    for residue_number, (_, atoms) in enumerate(residues, 1):
        for _, _, bonded_names in atoms:
            bonds.append([ atom_names.index((residue_number, name)) for name in bonded_names ])
    structure.bonds = bonds
    return structure


def get_incoherent_labels(structure: Structure, ua_selection: str = None) -> list[str]:
    """Get the residue and atom names of atoms with incoherent bonds."""
    selection = structure.select(ua_selection) if ua_selection else None
    return [ f'{structure.atoms[index].residue.name}.{structure.atoms[index].name}'
        for index in structure.check_incoherent_bonds(selection) ]


@pytest.mark.unit_int
def test_united_atom_limits():
    """United atoms may have less bonds than atomistic atoms, but never more."""
    assert coherent_bonds_united_atom['C'] == {'min': 1, 'max': 4}
    assert coherent_bonds_united_atom['O'] == {'min': 0, 'max': 2}
    assert coherent_bonds_united_atom['N'] == {'min': 1, 'max': 4}
    assert coherent_bonds_united_atom['H'] == coherent_bonds_with_hydrogen['H']
    for element, limits in coherent_bonds_united_atom.items():
        assert limits['min'] <= coherent_bonds_with_hydrogen[element]['min']
        assert limits['max'] == coherent_bonds_with_hydrogen[element]['max']


@pytest.mark.unit_int
def test_united_atom_methyls():
    """Terminal methyls in united atom residues pass only when they are in the united atom selection.
    Explicit hydrogens in the same residue and in water are still fine.
    """
    structure = build_structure([('UAL', UNITED_ATOM_LIPID), ('HOH', WATER)])
    assert get_incoherent_labels(structure) == ['UAL.C3']
    assert get_incoherent_labels(structure, 'resname UAL') == []


@pytest.mark.unit_int
def test_atomistic_residues_stay_strict():
    """Atoms out of the united atom selection keep the atomistic limits in the same structure."""
    structure = build_structure([('UAL', UNITED_ATOM_LIPID), ('ETH', BROKEN_ETHANE), ('HOH', WATER)])
    assert get_incoherent_labels(structure, 'resname UAL') == ['ETH.C2']


@pytest.mark.unit_int
def test_united_atoms_with_too_many_bonds():
    """United atoms with more bonds than the element allows still fail."""
    pentavalent_carbon = [
        ('C1', 'C', ['C2', 'C3', 'C4', 'C5', 'C6']),
        *[ (f'C{number}', 'C', ['C1']) for number in range(2, 7) ],
    ]
    structure = build_structure([('UAL', pentavalent_carbon), ('HOH', WATER)])
    assert get_incoherent_labels(structure, 'resname UAL') == ['UAL.C1']
