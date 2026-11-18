"""Trajectory integrity test (sudden RMSD jumps)."""
from types import SimpleNamespace
import mdtraj as mdt
import numpy as np
import pytest
from mddb_workflow.analyses.rmsd_check import check_trajectory_integrity
from mddb_workflow.utils.auxiliar import TestFailure
from mddb_workflow.utils.selections import Selection
from mddb_workflow.utils.structures import Structure

WATERS = 1000
FRAMES = 400
# Odd frame in the middle of the trajectory
ODD_FRAME = 200


def write_rigid_waters(filename: str, odd_frame_change: str) -> Structure:
    """Write a trajectory of rigid waters which move and rotate at random, so their bonds change only by the xtc rounding.

    The odd frame has more rounding noise ('noise') or a hydrogen far from its oxygen ('broken').
    Return the structure with the water bonds.
    """
    rng = np.random.default_rng(0)
    angle = np.radians(104.52)
    water = np.array([[0, 0, 0], [0.09572, 0, 0], [0.09572 * np.cos(angle), 0.09572 * np.sin(angle), 0]])
    positions = np.empty((FRAMES, WATERS * 3, 3), dtype=np.float32)
    for frame in range(FRAMES):
        rotations = np.linalg.qr(rng.normal(size=(WATERS, 3, 3)))[0]
        centers = rng.uniform(1, 9, size=(WATERS, 1, 3))
        positions[frame] = (centers + water @ rotations.transpose(0, 2, 1)).reshape(-1, 3)
    if odd_frame_change == 'noise':
        positions[ODD_FRAME] += rng.uniform(-0.0005, 0.0005, size=positions[ODD_FRAME].shape)
    if odd_frame_change == 'broken':
        positions[ODD_FRAME, 1] += 0.5
    boxes = np.repeat(np.diag([10, 10, 10]).astype(np.float32)[np.newaxis], FRAMES, axis=0)
    with mdt.formats.XTCTrajectoryFile(filename, 'w') as writer:
        writer.write(positions, box=boxes)
    lines = []
    for atom_index in range(WATERS * 3):
        residue_number, atom_name = atom_index // 3 + 1, ['O', 'H1', 'H2'][atom_index % 3]
        lines.append(f'ATOM  {atom_index + 1:5d}  {atom_name:<3} HOH A{residue_number:4d}    '
            f'{0:8.3f}{0:8.3f}{0:8.3f}  1.00  0.00           {atom_name[0]}')
    structure = Structure.from_pdb('\n'.join(lines) + '\n')
    structure.bonds = [[index + 1, index + 2] if index % 3 == 0 else [index - index % 3] for index in range(WATERS * 3)]
    return structure


@pytest.mark.parametrize('odd_frame_change, fails', [('noise', False), ('broken', True)])
def test_rounding_noise_is_not_a_jump(tmp_path, monkeypatch, odd_frame_change, fails):
    """A frame whose bonds differ from the rest less than the xtc precision passes, though its z-score is high."""
    monkeypatch.chdir(tmp_path)
    structure = write_rigid_waters('trajectory.xtc', odd_frame_change)
    register = SimpleNamespace(tests={}, remove_warnings=lambda tag: None, update_test=lambda tag, value: None)
    arguments = dict(input_structure_filename='structure.pdb', input_trajectory_filename='trajectory.xtc',
        view_script_filename='jump.tcl', structure=structure, pbc_selection=Selection(list(range(WATERS * 3))),
        mercy=[], trust=[], register=register, check_selection='all', standard_deviations_cutoff=9,
        update_image=lambda: False, snapshots=FRAMES, md_index=0, skipped_frames=0)
    if fails:
        with pytest.raises(TestFailure, match='RMSD check has failed'):
            check_trajectory_integrity(**arguments)
    else:
        assert check_trajectory_integrity(**arguments)
