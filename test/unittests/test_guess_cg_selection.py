"""Tests for the coordinates used to guess coarse grain and united atom residues."""
import pathlib
import numpy as np
import pytest
import pytraj as pt

from mddb_workflow.utils.file import File
from mddb_workflow.utils.structures import Structure
from mddb_workflow.tools.guess_cg_selection import get_frames_coordinates, FRAMES_LIMIT

dummy_dir = pathlib.Path(__file__).parent.parent/'data/input/dummy'


@pytest.mark.unit_int
def test_frames_coordinates_with_few_frames():
    """Every frame keeps its own coordinates when the trajectory has less frames than the limit."""
    structure = Structure.from_file(str(dummy_dir/'ala_ala.pdb'))
    trajectory_path = str(dummy_dir/'gromacs/raw_trajectory.xtc')
    # Load the whole trajectory in memory, so every frame has its own coordinates
    expected = pt.load(trajectory_path, top=str(dummy_dir/'ala_ala.pdb')).xyz
    assert 1 < len(expected) <= FRAMES_LIMIT
    coordinates = get_frames_coordinates(structure, File(trajectory_path))
    assert np.allclose(coordinates, expected, atol=0.01)
