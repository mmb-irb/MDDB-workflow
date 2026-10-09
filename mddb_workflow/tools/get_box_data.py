import numpy as np

from mddb_workflow.utils.auxiliar import InputError, MISSING_TOPOLOGY
from mddb_workflow.utils.gmx_spells import get_xtc_simulation_box, get_box_shape
from mddb_workflow.utils.topologies import Topology
from mddb_workflow.utils.type_hints import *

# Maximum relative difference allowed between the topology and trajectory box volumes
# The topology box may come from before the equilibration (e.g. tleap) so some difference is expected
BOX_VOLUME_TOLERANCE = 0.1

# This is a wrapper to make the gromacs function compatible with the tasks format
def mine_simulation_box_data (trajectory_file : 'File', topology_file : 'File', ignore_box : bool):
    """Mine simulation box data: its vectors, its size, whether it is constant or dynamic, whether it is orthogonal and its shape"""
    # If the box is to be ignored then act as if there was no box
    if ignore_box:
        print(' The simulation box is ignored')
        return None, None, None, None, None
    box_data = get_xtc_simulation_box(xtc_path=trajectory_file.path)
    check_topology_box(topology_file, box_data)
    return box_data

def check_topology_box (topology_file : 'File', box_data : tuple):
    """Make sure the trajectory box matches the topology box, if both have a box.
    Box shapes must match and, if the trajectory box is constant, volumes must be similar as well."""
    box, box_size, is_constant, is_orthogonal, box_shape = box_data
    # If there is no box in the trajectory or no topology then there is nothing to compare
    if box is None or topology_file == MISSING_TOPOLOGY: return
    topology_box = Topology(topology_file).get_box()
    if topology_box is None: return
    topology_shape = get_box_shape(topology_box)
    # A dynamic box may change its shape, so compare the shape in the first frame
    trajectory_shape = box_shape.split(' / ')[0]
    problems = []
    if topology_shape != trajectory_shape:
        problems.append(f'Box shape in topology is {topology_shape} but it is {trajectory_shape} in trajectory')
    if is_constant:
        topology_volume = abs(np.linalg.det(topology_box))
        trajectory_volume = abs(np.linalg.det(np.array(box)))
        difference = abs(trajectory_volume - topology_volume) / topology_volume
        if difference > BOX_VOLUME_TOLERANCE:
            problems.append(f'Box volume in topology is {topology_volume / 1000:.1f} nm³ '
                f'but it is {trajectory_volume / 1000:.1f} nm³ in trajectory ({difference:.0%} difference)')
    if not problems: return
    raise InputError('Simulation boxes in topology and trajectory do not match:\n ' + '\n '.join(problems) + '\n'
        ' The trajectory box is probably wrong, which would break any box-dependent test or analysis.\n'
        ' If you want to proceed anyway then use the "--ignore_box" argument to ignore the simulation box.')
