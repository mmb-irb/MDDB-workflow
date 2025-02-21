from model_workflow.tools.get_pytraj_trajectory import get_pytraj_trajectory, get_reduced_pytraj_trajectory
from model_workflow.utils.auxiliar import reprint
from tqdm import tqdm
import os
from typing import Optional

import pytraj as pt
# Build a generator which returns frames from the trajectory in pdb format
# The frames limit is the maximum number of frames to be iterated
# Note that the number of frames iterated may be less than the specified number
def get_pdb_frames (
    topology_filename : str,
    trajectory_filename : str,
    snapshots : int,
    frames_limit : Optional[int] = None,
    output_frames_prefix : str = 'frame',
    pbar_bool : bool = False,
):
    # WARNING: Do not set a pytraj iterload trajectory and read its 'n_frames' to get the snapshots
    # WARNING: Trying to read the number o frames of a xtc trajectory will read the whole trajectory
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

    def frames_generator():
        # Get the current directory at this point and use it to delete old files, in case we change the directory
        cwd = os.getcwd()
        # Create a progress bar
        if pbar_bool: pbar = tqdm(initial=0, desc=' Frames', total=frames_count, unit='frame')
        # Or print an empty line for the reprint to not delete a previous log
        else: print()
        # Extract each frame in pdb format
        for f in range(frames_count):
            # Get the actual frame number
            # We display latter the frame with a +1 to make it 1-based instead of 0-based
            frame_number = f * frames_step
            # Update the current frame log
            if pbar_bool: pbar.update(1); pbar.refresh()
            else: reprint(f'Frame {frame_number+1} ({f+1} / {frames_count})')
            current_frame = f'{cwd}/{output_frames_prefix}{frame_number+1}.pdb'
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
    trajectory_frame_filename = f'frame{frame}.pdb'
    pt.write_traj(trajectory_frame_filename, trajectory_frame, overwrite=True)
    return trajectory_frame_filename