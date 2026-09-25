from mddb_workflow.utils.gmx_spells import get_xtc_simulation_box
from mddb_workflow.utils.type_hints import *

# This is a wrapper to make the gromacs function compatible with the tasks format
def mine_simulation_box_data (trajectory_file : 'File'):
    """Mine simulation box data: its vectors, its size, whether it is constant or dynamic, whether it is orthogonal and its shape"""
    return get_xtc_simulation_box(xtc_path=trajectory_file.path)