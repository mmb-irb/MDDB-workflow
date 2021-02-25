# Get frames count

import pytraj as pt

# Get the trajectory frames number using pytraj
def get_frames_count (
    input_topology_filename : str,
    input_trajectory_filename : str) -> int:
    
    # Load the trajectory from pytraj
    pt_trajectory = pt.iterload(
        input_trajectory_filename,
        input_topology_filename)

    # Return the frames number
    frames = pt_trajectory.n_frames
    return frames