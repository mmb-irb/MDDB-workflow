"""Frames skipped from the beginning and from the end of trajectories, including an incomplete last frame."""
import mdtraj as mdt
import numpy as np
import pytest
from mddb_workflow.core.inputs import ProjectInputs
from mddb_workflow.mwf import Project
from mddb_workflow.utils.file import File
from mddb_workflow.utils.auxiliar import InputError
from mddb_workflow.utils.gmx_spells import copy_xtc_frames, count_complete_xtc_frames

ATOMS = 1000
FRAMES = 5


def write_xtc(path: str, frames: int = FRAMES) -> np.ndarray:
    """Write the first frames of an XTC file and return the coordinates of all of them."""
    rng = np.random.default_rng(0)
    coordinates = rng.uniform(0, 5, (FRAMES, ATOMS, 3)).astype(np.float32)
    with mdt.formats.XTCTrajectoryFile(path, 'w') as writer:
        writer.write(coordinates[:frames], box=np.tile(np.eye(3, dtype=np.float32) * 5, (frames, 1, 1)))
    return coordinates


@pytest.fixture
def truncated_xtc(tmp_path) -> tuple[str, np.ndarray]:
    """Write an XTC file whose last frame is cut in the middle (e.g. an interrupted file transfer)."""
    path = str(tmp_path / 'truncated.xtc')
    coordinates = write_xtc(path)
    with open(path, 'rb+') as xtc_file:
        size = xtc_file.seek(0, 2)
        xtc_file.truncate(size - 100)
    return path, coordinates


@pytest.fixture
def corrupted_xtc(tmp_path) -> str:
    """Write an XTC file followed by data which is not a frame."""
    path = str(tmp_path / 'corrupted.xtc')
    write_xtc(path)
    with open(path, 'ab') as xtc_file:
        xtc_file.write(b'not a frame' * 10)
    return path


def test_count_complete_frames(truncated_xtc, corrupted_xtc, tmp_path):
    path, _ = truncated_xtc
    assert count_complete_xtc_frames(path) == (FRAMES - 1, True)
    # A file cut in the middle of the magic number of its last frame
    # The first frames are the same in both files, so the shorter one ends where the last frame starts
    write_xtc(str(tmp_path / 'complete.xtc'))
    write_xtc(str(tmp_path / 'shorter.xtc'), frames=FRAMES - 1)
    with open(tmp_path / 'complete.xtc', 'rb') as xtc_file:
        data = xtc_file.read()
    last_frame_start = (tmp_path / 'shorter.xtc').stat().st_size
    cut_path = str(tmp_path / 'cut_magic.xtc')
    with open(cut_path, 'wb') as xtc_file:
        xtc_file.write(data[:last_frame_start + 2])
    assert count_complete_xtc_frames(cut_path) == (FRAMES - 1, True)
    # Data after the complete frames which is not a frame is an error
    with pytest.raises(RuntimeError, match='not valid data'):
        count_complete_xtc_frames(corrupted_xtc)


def test_incomplete_frame_is_an_error(truncated_xtc, tmp_path):
    """Copying frames fails on the incomplete one unless it is allowed and not requested."""
    path, _ = truncated_xtc
    # Every frame (e.g. to skip the first ones), allowed or not
    with pytest.raises(RuntimeError, match='skip_frames'):
        copy_xtc_frames(path, str(tmp_path / 'all.xtc'), lambda frame: frame >= 1)
    with pytest.raises(RuntimeError, match='skip_frames'):
        copy_xtc_frames(path, str(tmp_path / 'all.xtc'), lambda frame: frame >= 1, allow_incomplete_last_frame=True)
    # Every complete frame, but not allowed
    with pytest.raises(RuntimeError, match='corrupted'):
        copy_xtc_frames(path, str(tmp_path / 'listed.xtc'), list(range(FRAMES - 1)))
    # The incomplete frame, allowed
    with pytest.raises(RuntimeError, match='corrupted'):
        copy_xtc_frames(path, str(tmp_path / 'listed.xtc'), list(range(FRAMES)), allow_incomplete_last_frame=True)


def test_not_valid_data_is_an_error(corrupted_xtc, tmp_path):
    """Data after the complete frames which is not a frame is never taken as an incomplete frame."""
    with pytest.raises(RuntimeError, match='corrupted'):
        copy_xtc_frames(corrupted_xtc, str(tmp_path / 'listed.xtc'), list(range(FRAMES)), allow_incomplete_last_frame=True)


def test_complete_frames_are_copied(truncated_xtc, tmp_path):
    """Frames before the incomplete one are copied, as when the last frame is skipped."""
    path, coordinates = truncated_xtc
    output_path = str(tmp_path / 'skipped.xtc')
    assert copy_xtc_frames(path, output_path, range(1, FRAMES - 1), allow_incomplete_last_frame=True) == FRAMES - 2
    assert count_complete_xtc_frames(output_path) == (FRAMES - 2, False)
    with mdt.formats.XTCTrajectoryFile(output_path) as reader:
        copied = reader.read()[0]
    np.testing.assert_allclose(copied, coordinates[1:FRAMES - 1], atol=0.001)


def write_inputs(tmp_path, skip_frames: str) -> str:
    inputs_path = tmp_path / 'inputs.yaml'
    inputs_path.write_text(f'mds:\n  - name: replica 1\n  - name: replica 2\nskip_frames: {skip_frames}\n', encoding='utf-8')
    return str(inputs_path)


def test_skip_frames_from_both_ends(tmp_path):
    """Positive numbers skip frames from the beginning and negative numbers from the end."""
    project = Project(directory=str(tmp_path), inputs_filepath=write_inputs(tmp_path, '{0: [100, -1], 1: -2}'))
    assert project.skip_frames == {0: (100, 1), 1: (0, 2)}
    assert project.mds[0].input_skip_frames == (100, 1)


def test_skip_frames_from_console(tmp_path):
    """Values from the console come as pairs, replace the ones of the MDs they name and are saved in the inputs file."""
    inputs_filepath = write_inputs(tmp_path, '{0: [5, -2], 1: 7}')
    project = Project(directory=str(tmp_path), inputs_filepath=inputs_filepath, skip_frames=[[0, 10], [0, -1]])
    assert project.skip_frames == {0: (10, 1), 1: (7, 0)}
    assert ProjectInputs(File(inputs_filepath), directory_name=tmp_path.name).get('skip_frames') == {0: [10, -1], 1: 7}


def test_skip_frames_twice_from_the_same_end(tmp_path):
    with pytest.raises(InputError, match='at most one positive number'):
        Project(directory=str(tmp_path), inputs_filepath=write_inputs(tmp_path, '{0: [100, 50]}'))


def test_skip_frames_of_missing_md(tmp_path):
    with pytest.raises(InputError, match='MD index 2'):
        Project(directory=str(tmp_path), inputs_filepath=write_inputs(tmp_path, '{2: -1}'))
