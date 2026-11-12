"""Automatic translation of the imaging on synthetic systems, in nm."""
import mdtraj as mdt
import numpy as np
import pytest
from mddb_workflow.tools.image_and_fit import find_automatic_translation
from mddb_workflow.utils.auxiliar import InputError
from mddb_workflow.utils.file import File
from mddb_workflow.utils.selections import Selection

RECTANGULAR_BOX = np.diag([4, 4.5, 8]).astype(float)
# Gromacs triclinic box: the first vector along X and the second one in the XY plane
TRICLINIC_BOX = np.array([[4, 0, 0], [1.2, 4.2, 0], [-1, 1.5, 7.5]], dtype=float)
# Heights of the beads of every lipid from the bilayer midplane, so all beads are in contact with their neighbours
LIPID_HEIGHTS = np.array([0.2, 0.6, 1.0, 1.4, 1.8])


def wrap(positions: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Put every atom back in the box."""
    fractional = positions @ np.linalg.inv(box)
    return (fractional - np.floor(fractional)) @ box


def flat_bilayer(box: np.ndarray, height: float, spacing: float = 0.4) -> np.ndarray:
    """Get the beads of a flat bilayer in the plane of the first two box vectors, filling the box along them.
    Every lipid is a column of beads along Z, with leaflets in contact at the midplane, which is at the given height.
    """
    counts = np.ceil(np.linalg.norm(box[:2], axis=1) / spacing).astype(int)
    grid = np.array([(i / counts[0], j / counts[1], 0) for i in range(counts[0]) for j in range(counts[1])]) @ box
    heights = height + np.concatenate([-LIPID_HEIGHTS, LIPID_HEIGHTS])
    return wrap(np.concatenate([grid + [0, 0, z] for z in heights]), box)


def blob(center: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Get the beads of a small globular object (e.g. a ligand) around the given center, put back in the box."""
    offsets = np.array([(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]) * 0.3
    return wrap(center + offsets, box)


def sphere(radius: float, box: np.ndarray, spacing: float = 0.4) -> np.ndarray:
    """Get the beads of a filled sphere (e.g. a lipid nanoparticle) in the middle of the box."""
    axis = np.arange(-radius, radius + spacing, spacing)
    grid = np.array([(x, y, z) for x in axis for y in axis for z in axis])
    return grid[np.linalg.norm(grid, axis=1) <= radius] + np.ones(3) @ box / 2


def translate(tmp_path, positions: np.ndarray, box: np.ndarray, ignore_periodic: bool = False) -> np.ndarray:
    """Write the positions as a single frame trajectory and find its automatic translation."""
    filepath = str(tmp_path / 'trajectory.xtc')
    with mdt.formats.XTCTrajectoryFile(filepath, 'w') as writer:
        writer.write(positions[np.newaxis].astype(np.float32), box=box[np.newaxis].astype(np.float32))
    selection = Selection(list(range(len(positions))))
    return np.array(find_automatic_translation(File(filepath), selection, ignore_periodic=ignore_periodic))


def fractional_spans(positions: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Get the span of the atoms along every box vector, in box vector units, once they are in the box."""
    fractional = wrap(positions, box) @ np.linalg.inv(box)
    return fractional.max(axis=0) - fractional.min(axis=0)


@pytest.mark.parametrize('box', [RECTANGULAR_BOX, TRICLINIC_BOX], ids=['rectangular', 'triclinic'])
@pytest.mark.parametrize('midplane', [0, 0.1, 0.85])
def test_split_bilayer(tmp_path, box, midplane):
    """A bilayer of PBC lipids split by the boundary along its normal is moved to the middle of the box, but not along its plane."""
    bilayer = flat_bilayer(box, midplane * box[2, 2])
    translation = translate(tmp_path, bilayer, box, ignore_periodic=True)
    fractional_translation = translation @ np.linalg.inv(box)
    assert np.allclose(fractional_translation[:2], 0, atol=1e-3)
    # The bilayer is whole along its normal and the boundary is in the middle of the water
    thickness = 2 * LIPID_HEIGHTS.max() / box[2, 2]
    fractional = wrap(bilayer + translation, box) @ np.linalg.inv(box)
    assert fractional[:, 2].max() - fractional[:, 2].min() == pytest.approx(thickness, abs=1e-3)
    assert fractional[:, 2].mean() == pytest.approx(0.5, abs=1e-3)


@pytest.mark.parametrize('box', [RECTANGULAR_BOX, TRICLINIC_BOX], ids=['rectangular', 'triclinic'])
def test_centered_bilayer(tmp_path, box):
    """A bilayer of PBC lipids in the middle of the box is not translated."""
    translation = translate(tmp_path, flat_bilayer(box, box[2, 2] / 2), box, ignore_periodic=True)
    assert np.allclose(translation, 0, atol=1e-3)


@pytest.mark.parametrize('box', [RECTANGULAR_BOX, TRICLINIC_BOX], ids=['rectangular', 'triclinic'])
def test_split_object(tmp_path, box):
    """A globular object split by the boundaries along all box vectors is moved to the middle of the box."""
    ligand = blob(np.zeros(3), box)
    translation = translate(tmp_path, ligand, box)
    assert np.all(fractional_spans(ligand + translation, box) < 0.5)
    assert np.allclose(wrap(ligand + translation, box).mean(axis=0), np.ones(3) @ box / 2, atol=1e-3)


@pytest.mark.parametrize('box', [RECTANGULAR_BOX, TRICLINIC_BOX], ids=['rectangular', 'triclinic'])
def test_bilayer_and_ligand(tmp_path, box):
    """The bilayer is ignored along its plane, where the boundary is placed away from a ligand in the water.
    Along the normal the boundary is placed away from both.
    """
    bilayer = flat_bilayer(box, 0.5 * box[2, 2])
    # The ligand is split by the boundary of the first box vector
    ligand = blob(np.array([0, 0.5, 0.9]) @ box, box)
    translation = translate(tmp_path, np.concatenate([bilayer, ligand]), box, ignore_periodic=True)
    assert np.all(fractional_spans(ligand + translation, box) < 0.5)
    assert fractional_spans(bilayer + translation, box)[2] < 0.5


@pytest.mark.parametrize('box', [RECTANGULAR_BOX, TRICLINIC_BOX], ids=['rectangular', 'triclinic'])
def test_periodic_non_pbc_atoms(tmp_path, box):
    """Non-PBC atoms in contact with their own periodic image can not be kept whole, so it fails (e.g. a bilayer)."""
    with pytest.raises(InputError, match='in contact with their own periodic image'):
        translate(tmp_path, flat_bilayer(box, 0), box)


def test_nanoparticle_touching_its_image(tmp_path):
    """A lipid nanoparticle of non-PBC atoms which touches its own image along a single box vector fails."""
    box = np.diag([6, 9, 9]).astype(float)
    with pytest.raises(InputError, match='in contact with their own periodic image'):
        translate(tmp_path, wrap(sphere(3.1, box), box), box)


def test_filled_box(tmp_path):
    """PBC atoms in contact with their own periodic image along all box vectors are not translated."""
    counts = np.ceil(np.diag(RECTANGULAR_BOX) / 0.4).astype(int)
    lattice = np.array([(i, j, k) for i in range(counts[0]) for j in range(counts[1]) for k in range(counts[2])]) / counts
    translation = translate(tmp_path, lattice @ RECTANGULAR_BOX, RECTANGULAR_BOX, ignore_periodic=True)
    assert np.array_equal(translation, [0, 0, 0])


def test_object_longer_than_box(tmp_path, capsys):
    """An object longer than the box along the first box vector, but not periodic, fills it.
    Non-PBC atoms fail, while PBC atoms are not translated along this box vector, but they are along the others.
    """
    direction = np.array([1, 0.3, 0]) / np.linalg.norm([1, 0.3, 0])
    rod = np.arange(13)[:, np.newaxis] * 0.4 * direction + [0, 1, 2]
    box = np.diag([4, 8, 4]).astype(float)
    with pytest.raises(InputError):
        translate(tmp_path, wrap(rod, box), box)
    translation = translate(tmp_path, wrap(rod, box), box, ignore_periodic=True)
    assert translation[0] == 0
    assert 'fill the whole box along the box vector 1' in capsys.readouterr().out
    assert np.all(fractional_spans(rod + translation, box)[1:] < 0.5)
