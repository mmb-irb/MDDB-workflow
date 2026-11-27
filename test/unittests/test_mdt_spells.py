"""Reading xtc trajectories frame by frame."""
from typing import Optional
import mdtraj as mdt
import pytest
import numpy as np
from scipy.spatial.transform import Rotation
from mddb_workflow.utils.mdt_spells import iterate_xtc_frames, guess_translation_fitted, guess_rotation_fitted, find_fitted_groups


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



def get_frames(fitted_residues: Optional[slice], weights: np.ndarray, drift: float = 0) -> tuple[list, list[list[int]]]:
    """Get 3 frames of a flexible chain of 40 residues (5 atoms each), fitted to some residues as a fitting tool does.
    Optionally, every frame is rotated and shifted afterwards (drift in degrees and tenths of nm) so it is not fitted anymore.
    Return the frames and the atom indices of every residue.
    """
    generator = np.random.default_rng(0)
    residue_indices = [ list(range(residue * 5, residue * 5 + 5)) for residue in range(40) ]
    chain = np.cumsum(generator.normal(size=(40, 3)), axis=0) * 0.3
    frames = []
    for frame in range(3):
        # Every frame moves a lot, and so does every residue
        coordinates = (chain[:, None, :] + generator.normal(size=(40, 5, 3)) * 0.3).reshape(-1, 3)
        coordinates = Rotation.random(random_state=frame).apply(coordinates) + generator.normal(size=3) * 3
        if fitted_residues is not None:
            group = np.array([ index for indices in residue_indices[fitted_residues] for index in indices ])
            group_weights = weights[group] / weights[group].sum()
            center = group_weights @ coordinates[group]
            if frames:
                reference = frames[0][0][group]
                reference_center = group_weights @ reference
                rotation, _ = Rotation.align_vectors(reference - reference_center, coordinates[group] - center, weights=group_weights)
                coordinates = rotation.apply(coordinates - center) + reference_center
        if drift and frames:
            coordinates = Rotation.from_rotvec([0, 0, np.radians(drift * frame)]).apply(coordinates) + drift * frame / 10
        frames.append((coordinates, None))
    # Coordinates are written in xtc files with this precision
    frames = [ (np.round(coordinates, 3), box) for coordinates, box in frames ]
    return frames, residue_indices


@pytest.mark.parametrize('fitted_residues', [slice(31, 34), slice(10, 11), slice(0, 40)])
@pytest.mark.parametrize('use_masses', [True, False])
def test_fitted_to_group(fitted_residues, use_masses):
    weights = np.random.default_rng(1).uniform(1, 16, 200) if use_masses else np.ones(200)
    frames, residue_indices = get_frames(fitted_residues, weights)
    groups = find_fitted_groups(frames, [(residue_indices, weights)])
    assert len(groups) == 1
    assert groups[0] == [ index for indices in residue_indices[fitted_residues] for index in indices ]
    atoms = [ index for indices in residue_indices for index in indices ]
    assert guess_translation_fitted(frames, atoms, groups)
    assert guess_rotation_fitted(frames, residue_indices, groups)


def test_fitted_with_other_weights():
    # The group is found only with the weights the fitting tool used
    masses = np.random.default_rng(1).uniform(1, 16, 200)
    frames, residue_indices = get_frames(slice(31, 34), masses)
    assert find_fitted_groups(frames, [(residue_indices, np.ones(200))]) == []
    assert len(find_fitted_groups(frames, [(residue_indices, np.ones(200)), (residue_indices, masses)])) == 1


@pytest.mark.parametrize('drift', [0.5, 2, 10])
def test_not_fitted(drift):
    # Fitted trajectories which are slightly rotated and shifted afterwards
    weights = np.ones(200)
    frames, residue_indices = get_frames(slice(31, 34), weights, drift)
    assert find_fitted_groups(frames, [(residue_indices, weights)]) == []
    # Neither is a trajectory which was never fitted
    frames, residue_indices = get_frames(None, weights)
    groups = find_fitted_groups(frames, [(residue_indices, weights)])
    assert groups == []
    atoms = [ index for indices in residue_indices for index in indices ]
    assert not guess_translation_fitted(frames, atoms, groups)
    assert not guess_rotation_fitted(frames, residue_indices, groups)
