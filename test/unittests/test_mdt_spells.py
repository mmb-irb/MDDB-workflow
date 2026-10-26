"""Reading xtc trajectories frame by frame."""
import mdtraj as mdt
import pytest
import numpy as np
from mddb_workflow.utils.mdt_spells import iterate_xtc_frames


def write_xtc(filename: str, box: np.ndarray | None) -> np.ndarray:
    """Write a trajectory of 2 frames and 3 atoms with the given box vectors (nm), if any, and return its positions."""
    positions = np.arange(18, dtype=np.float32).reshape(2, 3, 3) / 10
    with mdt.formats.XTCTrajectoryFile(filename, 'w') as writer:
        writer.write(positions, box=None if box is None else np.repeat(box[np.newaxis], 2, axis=0))
    return positions


def test_box(tmp_path):
    box = np.diag([5, 6, 7]).astype(np.float32)
    positions = write_xtc(str(tmp_path / 'trajectory.xtc'), box)
    frames = list(iterate_xtc_frames(str(tmp_path / 'trajectory.xtc')))
    assert len(frames) == 2
    for (frame_positions, frame_box), expected_positions in zip(frames, positions):
        assert np.allclose(frame_positions, expected_positions, atol=1e-3)
        assert np.allclose(frame_box, box)


def test_no_box(tmp_path):
    """Trajectories with no box give a box full of zeros, which means there are no periodic boundaries."""
    write_xtc(str(tmp_path / 'trajectory.xtc'), None)
    frames = list(iterate_xtc_frames(str(tmp_path / 'trajectory.xtc')))
    assert len(frames) == 2
    assert all(frame_box.shape == (3, 3) and not frame_box.any() for _, frame_box in frames)


@pytest.mark.parametrize('box', [None, np.diag([5, 6, 7]).astype(np.float32)])
def test_jumping_atoms(tmp_path, box):
    """The atoms of the bond which changes the most in a jump are found, also in trajectories with no box."""
    from mddb_workflow.analyses.rmsd_check import find_jumping_atoms
    from mddb_workflow.utils.selections import Selection
    positions = np.array([[[0, 0, 0], [0.1, 0, 0], [0.2, 0, 0]], [[0, 0, 0], [0.1, 0, 0], [1.2, 0, 0]]], dtype=np.float32)
    with mdt.formats.XTCTrajectoryFile(str(tmp_path / 'trajectory.xtc'), 'w') as writer:
        writer.write(positions, box=None if box is None else np.repeat(box[np.newaxis], 2, axis=0))
    pairs = np.array([[0, 1], [1, 2]])
    jumping_atoms = find_jumping_atoms(str(tmp_path / 'trajectory.xtc'), Selection([0, 1, 2]), pairs, jump_frame=1)
    assert jumping_atoms.atom_indices == [1, 2]
