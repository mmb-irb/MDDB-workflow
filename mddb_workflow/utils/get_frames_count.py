from mddb_workflow.utils.gmx_spells import get_xtc_frame_count
from mddb_workflow.utils.pyt_spells import get_frames_count as pyt_get_frames_count
from mddb_workflow.utils.type_hints import *

# Generic smart frame counter
# Uses the best frame-counting function depending on the trajectory format
def get_frames_count(
    structure_file : 'File',
    trajectory_file : 'File',
    verbose : bool = True,
) -> int:
    """Get the trajectory frames count."""
    # For XTC format we have a specific function which is way faster
    if trajectory_file.format == 'xtc':
        return get_xtc_frame_count(trajectory_file.path, verbose)
    # For the rest of formats use the generic function from pytraj
    return pyt_get_frames_count(structure_file, trajectory_file, verbose)