# Get frames count

from os.path import exists

import pytraj as pt

# LORE: This was tried also with mdtraj's iterload but pytraj was way faster

from model_workflow.utils.auxiliar import InputError

# Get the trajectory frames number using pytraj
def get_frames_count (
    input_topology_filename : str,
    input_trajectory_filename : str) -> int:
    
    print('Counting number of frames...')

    if not exists(input_trajectory_filename):
        raise InputError('Missing trajectroy file when counting frames: ' + input_trajectory_filename)
    
    if not exists(input_topology_filename):
        raise InputError('Missing topology file when counting frames: ' + input_topology_filename)

    # Load the trajectory from pytraj
    pt_trajectory = pt.iterload(
        input_trajectory_filename,
        input_topology_filename)

    # Return the frames number
    frames = pt_trajectory.n_frames
    print(' Frames: ' + str(frames))

    # If 0 frames were counted then there is something wrong with the file
    if frames == 0:
        raise InputError('Something went wrong when reading the trajectory')

    return frames