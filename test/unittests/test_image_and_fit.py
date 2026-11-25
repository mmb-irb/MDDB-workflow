"""Automatic translation and fragments in the box of the imaging on synthetic systems, in nm."""
from types import SimpleNamespace
from typing import Optional
import mdtraj as mdt
import numpy as np
import pytest
from mddb_workflow.tools.image_and_fit import find_automatic_translation, find_layer_axis, put_fragments_in_box, XTC_PRECISION
from mddb_workflow.utils.auxiliar import InputError
from mddb_workflow.utils.file import File
from mddb_workflow.utils.selections import Selection

RECTANGULAR_BOX = np.diag([4, 4.5, 8]).astype(float)
# Gromacs triclinic box: the first vector along X and the second one in the XY plane
TRICLINIC_BOX = np.array([[4, 0, 0], [1.2, 4.2, 0], [-1, 1.5, 7.5]], dtype=float)
# Heights of the beads of every lipid from the bilayer midplane, so all beads are in contact with their neighbours
LIPID_HEIGHTS = np.array([0.2, 0.6, 1.0, 1.4, 1.8])
# Boxes where the leaflets of a bilayer touch their image across a thin water layer (0.3 nm)
THIN_RECTANGULAR_BOX = np.diag([4, 4.5, 3.9]).astype(float)
THIN_TRICLINIC_BOX = np.array([[4, 0, 0], [1.2, 4.2, 0], [-1, 1.5, 3.9]], dtype=float)


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


def write_frames(filepath: str, frames: list[np.ndarray], box: np.ndarray) -> File:
    """Write the positions of every frame as a trajectory with the same box in all frames."""
    with mdt.formats.XTCTrajectoryFile(filepath, 'w') as writer:
        writer.write(np.array(frames, dtype=np.float32), box=np.repeat(box[np.newaxis], len(frames), axis=0).astype(np.float32))
    return File(filepath)


def translate(tmp_path, positions: np.ndarray, box: np.ndarray, ignore_periodic: bool = False) -> np.ndarray:
    """Write the positions as a single frame trajectory and find its automatic translation."""
    trajectory_file = write_frames(str(tmp_path / 'trajectory.xtc'), [positions], box)
    selection = Selection(list(range(len(positions))))
    return np.array(find_automatic_translation(trajectory_file, selection, ignore_periodic=ignore_periodic))


def bilayer_and_water(box: np.ndarray, midplane: float) -> tuple[np.ndarray, np.ndarray, list[list[int]]]:
    """Get a flat bilayer, a water layer in the middle of the water between its images and the bonds of all atoms.
    Every lipid is a column of bonded beads and every water molecule is a single bead.
    """
    lipids = flat_bilayer(box, midplane)
    grid_count = len(lipids) // (2 * len(LIPID_HEIGHTS))
    water = wrap(lipids[:grid_count] + [0, 0, LIPID_HEIGHTS[0]] + box[2] / 2, box)
    atom_bonds = [[] for _ in range(len(lipids) + len(water))]
    # Beads of the same lipid are at the same grid point, in consecutive heights of the same leaflet
    for height in [*range(len(LIPID_HEIGHTS) - 1), *range(len(LIPID_HEIGHTS), 2 * len(LIPID_HEIGHTS) - 1)]:
        for point in range(grid_count):
            atom, bonded = height * grid_count + point, (height + 1) * grid_count + point
            atom_bonds[atom].append(bonded)
            atom_bonds[bonded].append(atom)
    return lipids, water, atom_bonds


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


def test_fragments_in_box_keep_bonds(tmp_path):
    """Molecules far from the box (e.g. water after the no-jump) are put back in the box with the same bonds.
    The box length is half a step of the xtc grid, so atoms of a molecule could be rounded in opposite directions.
    """
    box = np.diag([5.9095, 5.9095, 9])
    water = np.array([[0, 0, 0], [0.096, 0, 0], [-0.024, 0.093, 0]])
    centers = np.random.default_rng(0).uniform(0, 5.5, (1000, 3)) + [301, 0, 0]
    positions = np.round(np.concatenate([center + water for center in centers]), 3)
    input_file, output_file = File(str(tmp_path / 'input.xtc')), File(str(tmp_path / 'output.xtc'))
    with mdt.formats.XTCTrajectoryFile(input_file.path, 'w') as writer:
        writer.write(positions[np.newaxis].astype(np.float32), box=box[np.newaxis].astype(np.float32))
    atom_bonds = [ [atom + 1, atom + 2] if atom % 3 == 0 else [atom - atom % 3] for atom in range(len(positions)) ]
    put_fragments_in_box(input_file, output_file, atom_bonds, SimpleNamespace(atom_count=len(positions)), Selection([]))
    oxygens = np.arange(0, len(positions), 3)
    def bond_steps(trajectory_file: File) -> np.ndarray:
        with mdt.formats.XTCTrajectoryFile(trajectory_file.path) as reader:
            steps = np.rint(reader.read()[0][0].astype(np.float64) / XTC_PRECISION).astype(int)
        return np.concatenate([steps[oxygens + 1] - steps[oxygens], steps[oxygens + 2] - steps[oxygens]])
    assert np.array_equal(bond_steps(output_file), bond_steps(input_file))


LAYER_BOXES = [THIN_RECTANGULAR_BOX, THIN_TRICLINIC_BOX, RECTANGULAR_BOX, TRICLINIC_BOX]
LAYER_BOX_IDS = ['thin-rectangular', 'thin-triclinic', 'rectangular', 'triclinic']


@pytest.mark.parametrize('box', LAYER_BOXES, ids=LAYER_BOX_IDS)
def test_water_layer_axis(tmp_path, box):
    """Water between bilayers forms a layer along their normal, whether the bilayer touches its image or not.
    Water filling the box, in a drop or in a layer covering part of the plane only forms no layer.
    """
    _, water, _ = bilayer_and_water(box, 0)
    drop = sphere(1, box)
    filling = np.array([(i, j, k) for i in range(10) for j in range(10) for k in range(10)]) / 10 @ box
    partial = water[(water @ np.linalg.inv(box))[:, 0] < 0.7]
    for positions, layer_axis in [(water, 2), (drop, None), (filling, None), (partial, None)]:
        trajectory_file = write_frames(str(tmp_path / 'trajectory.xtc'), [positions], box)
        assert find_layer_axis(trajectory_file, Selection(list(range(len(positions))))) == layer_axis


@pytest.mark.parametrize('box', LAYER_BOXES, ids=LAYER_BOX_IDS)
def test_water_layer_at_boundary(tmp_path, box):
    """A bilayer drifting along its normal, from split by the boundary in the first frame, stays whole in the middle of
    the box in every frame when the water layer is kept at the boundary.
    Otherwise it is split whenever its midplane is close to the boundary.
    """
    frames = []
    for midplane in np.linspace(0, 1, 8, endpoint=False) * box[2, 2]:
        lipids, water, atom_bonds = bilayer_and_water(box, midplane)
        frames.append(np.concatenate([lipids, water]))
    input_file = write_frames(str(tmp_path / 'input.xtc'), frames, box)
    structure = SimpleNamespace(atom_count=len(frames[0]))
    water_selection = Selection(list(range(len(lipids), len(frames[0]))))
    # The lower leaflet is the first half of the lipid beads, and the distance between leaflets is fixed while whole
    lower = len(lipids) // 2
    separation = 2 * LIPID_HEIGHTS.mean() / box[2, 2]
    def leaflet_separations(boundary_selection: Optional[Selection], boundary_axis: Optional[int]) -> np.ndarray:
        output_file = File(str(tmp_path / 'output.xtc'))
        put_fragments_in_box(input_file, output_file, atom_bonds, structure, Selection([]), boundary_selection, boundary_axis)
        with mdt.formats.XTCTrajectoryFile(output_file.path) as reader:
            heights = (reader.read()[0][:, :len(lipids)].astype(np.float64) @ np.linalg.inv(box))[:, :, 2]
        if boundary_axis is not None: assert np.allclose(heights.mean(axis=1), 0.5, atol=1e-2)
        return heights[:, lower:].mean(axis=1) - heights[:, :lower].mean(axis=1)
    assert np.allclose(leaflet_separations(water_selection, 2), separation, atol=1e-2)
    assert np.any(np.abs(leaflet_separations(None, None) - separation) > 0.1)
