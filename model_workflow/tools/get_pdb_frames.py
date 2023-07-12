from model_workflow.tools.get_pytraj_trajectory import get_pytraj_trajectory, get_reduced_pytraj_trajectory

import os
from subprocess import run, PIPE, Popen
import math
from typing import Optional

import pytraj as pt

CURSOR_UP_ONE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'
ERASE_PREVIOUS_LINE = CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE

# Build a generator which returns frames from the trajectory in pdb format
# The frames limit is the maximum number of frames to be iterated
# Note that the number of frames iterated may be less than the specified number
def get_pdb_frames (
    topology_filename : str,
    trajectory_filename : str,
    snapshots : int,
    frames_limit : Optional[int] = None,
    output_frames_prefix : str = 'frame',
):

    # Load the trajectory using pytraj
    trajectory = get_pytraj_trajectory(topology_filename, trajectory_filename)
    # WARNING: Do not read trajectory.n_frames to get the number of snapshots or you will read the whole trajectory
    # WARNING: This may be a lot of time for a huge trajectory. Use the snapshots input instead
    # In case we are missing a frames limit set the limit as the number os snapshots
    if frames_limit == None:
        frames_limit = snapshots

    # Set a maximum number of frames
    # If trajectory has more frames than the limit create a reduced trajectory
    reduced_trajectory, frames_step, frames_count = get_reduced_pytraj_trajectory(
        topology_filename,
        trajectory_filename,
        snapshots,
        frames_limit
    )
    frames_list = range(0, frames_count)

    def frames_generator():
        # Get the current directory at this point and use it to delete old files, in case we change the directory
        cwd = os.getcwd()
        # Print an empty line for the first 'ERASE_PREVIOUS_LINE' to not delete a previous log
        print()
        # Extract each frame in pdb format
        for f in frames_list:
            # Update the current frame log
            print(ERASE_PREVIOUS_LINE)
            print('Frame ' + str(f+1) + ' / ' + str(frames_count))
            current_frame = cwd + '/' + output_frames_prefix + str(f) + '.pdb'
            single_frame_trajectory = reduced_trajectory[f:f+1]
            pt.write_traj(current_frame, single_frame_trajectory, overwrite=True)
            yield current_frame
            # Delete current frame file before going for the next frame
            os.remove(current_frame)

    return frames_generator(), frames_step, frames_count

# Get a specific trajectory frame in pdb format
# Return the name of the generated file
def get_pdb_frame (
    topology_filename : str,
    trajectory_filename : str,
    frame : int
) -> str:

    # Load the trajectory using pytraj
    trajectory = get_pytraj_trajectory(topology_filename, trajectory_filename)
    trajectory_frame = trajectory[frame:frame+1]
    trajectory_frame_filename = 'frame' + str(frame) + '.pdb'
    pt.write_traj(trajectory_frame_filename, trajectory_frame, overwrite=True)
    return trajectory_frame_filename