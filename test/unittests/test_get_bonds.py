"""Tests for the bonds mined from topology files."""
import pytest

from mddb_workflow.utils.file import File
from mddb_workflow.tools.get_bonds import mine_topology_bonds


def write_psf(path, atoms: list[tuple[str, str, float]], bonds: list[tuple[int, int]]):
    """Write a CHARMM PSF file with one residue, given its atom names, types and masses, and its bonds (1-based)."""
    lines = ['PSF EXT', '', '         1 !NTITLE', ' REMARKS test', '', f'{len(atoms):>10} !NATOM']
    for number, (name, atom_type, mass) in enumerate(atoms, 1):
        lines.append(f'{number:>10} SYS      1        RES      {name:<8} {atom_type:<6}   0.000000 {mass:>13.4f}           0')
    lines += ['', f'{len(bonds):>10} !NBOND: bonds', ''.join(f'{a:>10}{b:>10}' for a, b in bonds), '']
    for section in ['NTHETA: angles', 'NPHI: dihedrals', 'NIMPHI: impropers', 'NDON: donors', 'NACC: acceptors']:
        lines += [f'         0 !{section}', '', '']
    lines += ['         0 !NNB', '', '         0         0', '', '         1         0 !NGRP', '         0         0         0', '']
    path.write_text('\n'.join(lines))


@pytest.mark.parametrize('atoms, bonds, expected_bonds', [
    # A TIP3 water as in CHARMM PSF files, which also bond both hydrogens to keep the water rigid
    ([('OH2', 'OT', 15.9994), ('H1', 'HT', 1.008), ('H2', 'HT', 1.008)], [(1, 2), (1, 3), (2, 3)], [[1, 2], [0], [0]]),
    # A H2 molecule, whose hydrogens are really bonded
    ([('H1', 'HH', 1.008), ('H2', 'HH', 1.008)], [(1, 2)], [[1], [0]]),
], ids=['rigid water', 'H2'])
def test_bonds_between_hydrogens(tmp_path, atoms, bonds, expected_bonds):
    """Only the bond between hydrogens bonded to the same atom (rigid water) is not a covalent bond."""
    psf_path = tmp_path / 'topology.psf'
    write_psf(psf_path, atoms, bonds)
    atom_bonds = mine_topology_bonds(File(str(psf_path)))
    assert [sorted(atom_bonds) for atom_bonds in atom_bonds] == expected_bonds
