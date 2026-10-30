from mddb_workflow.utils.gmx_spells import get_xtc_simulation_box
from mddb_workflow.utils.type_hints import *

# This is a wrapper to make the gromacs function compatible with the tasks format
# Note that the trajectory file must go first, so the input files are processed before checking if the box is ignored
def mine_simulation_box_data (trajectory_file : 'File', box_ignored : bool):
    """Mine simulation box data: its vectors, its size, whether it is constant or dynamic, whether it is orthogonal and its shape"""
    # If the box is to be ignored then act as if there was no box
    if box_ignored:
        print(' The simulation box is ignored')
        return None, None, None, None, None
    return get_xtc_simulation_box(xtc_path=trajectory_file.path)
