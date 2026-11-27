"""Tests for the bonds mined from topology files and the frames read to check them."""
import pathlib
import numpy as np
import pytest
import pytraj as pt

from mddb_workflow.utils.file import File
from mddb_workflow.utils.constants import GLOBALS
from mddb_workflow.utils.structures import Structure
from mddb_workflow.tools.get_bonds import mine_topology_bonds
from mddb_workflow.tools.get_pdb_frames import get_starting_pdb_frames

dummy_dir = pathlib.Path(__file__).parent.parent/'data/input/dummy'


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


@pytest.mark.unit_int
def test_starting_pdb_frames(tmp_path, monkeypatch):
    """The frames checked to find the reference frame are the first frames of the trajectory."""
    monkeypatch.setitem(GLOBALS, 'working_directory', str(tmp_path))
    structure = Structure.from_file(str(dummy_dir/'ala_ala.pdb'))
    # Set bonds as the workflow does, so they are not guessed when copying the structure
    structure.bonds = [[] for _ in structure.atoms]
    trajectory_path = str(dummy_dir/'gromacs/rmsd_jump.xtc')
    expected = pt.load(trajectory_path, top=str(dummy_dir/'ala_ala.pdb')).xyz[:3]
    frames = get_starting_pdb_frames(structure, trajectory_path, snapshots=201, frames_limit=3)
    coordinates = [[atom.coords for atom in Structure.from_pdb_file(frame).atoms] for frame in frames]
    assert np.allclose(coordinates, expected, atol=1e-3)
