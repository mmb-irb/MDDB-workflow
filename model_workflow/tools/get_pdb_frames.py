from model_workflow.tools.get_frames_count import get_frames_count
from model_workflow.tools.get_reduced_trajectory import get_reduced_trajectory
from model_workflow.tools.get_pytraj_trajectory import get_reduced_pytraj_trajectory

from subprocess import run, PIPE, Popen
import math

import pytraj as pt

# Build a generator which returns frames from the trajectory in pdb format
# The frames limit is the maximum number of frames to be iterated
# Note that the number of frames iterated may be less than the specified number
def get_pdb_frames (
    topology_filename : str,
    trajectory_filename : str,
    frames_limit : int
):

    # Load the trajectory using pytraj
    trajectory = pt.iterload(trajectory_filename, topology_filename)
    # Get the total number of snapshots in the trajectory
    snapshots = trajectory.n_frames

    # Set a maximum number of frames
    # If trajectory has more frames than the limit create a reduced trajectory
    reduced_trajectory = get_reduced_pytraj_trajectory(
        topology_filename,
        trajectory_filename,
        frames_limit
    )
    frames_step = reduced_trajectory.step
    frames_count = reduced_trajectory.n_frames
    frames_list = range(0, frames_count)

    def frames_generator():

        # Extract each frame in pdb format
        for f in frames_list:
            print('Frame ' + str(f+1) + ' / ' + str(frames_count))
            current_frame = 'frame' + str(f) + '.pdb'
            single_frame_trajectory = reduced_trajectory[f:f+1]
            pt.write_traj(current_frame, single_frame_trajectory, overwrite=True)
            yield current_frame
            # Delete current frame file before going for the next frame
            run([
                "rm",
                current_frame,
            ], stdout=PIPE).stdout.decode()

    return frames_generator(), frames_step, frames_count

# Build a generator which returns frames from the trajectory in pdb format
# The frames limit is the maximum number of frames to be iterated
# Note that the number of frames iterated may be less
# WARNING: Chains will be removed by Gromacs
# DEPRECATED
def get_pdb_frames_gromacs (
    topology_filename : str,
    trajectory_filename : str,
    frames_limit : int
):

    # get the total number of snapshots in the trajectory
    snapshots = get_frames_count(topology_filename, trajectory_filename)

    # Set the frames we must extract
    # WARNING: The gromacs '-fr' option counts frames starting at 1, not at 0
    frames_list = range(1, snapshots +1)

    # Set a maximum number of frames
    # If trajectory has more frames than the limit create a reduced trajectory
    reduced_trajectory_filename = trajectory_filename
    if snapshots > frames_limit:
        reduced_trajectory_filename, frames_step, frames_count = get_reduced_trajectory(
            topology_filename,
            trajectory_filename,
            frames_limit,
        )
        # WARNING: The gromacs '-fr' option counts frames starting at 1, not at 0
        frames_list = range(1, frames_count +1)  # if frames > 1 else [1]
    else:
        frames_count = snapshots

    def frames_generator():

        # Extract each frame in pdb format
        frames_ndx = 'frames.ndx'
        for f in frames_list:
            print('Frame ' + str(f) + ' / ' + str(frames_count))
            # Extract the current frame
            current_frame = 'frame' + str(f) + '.pdb'
            # The frame selection input in gromacs works with a 'ndx' file
            with open(frames_ndx, 'w') as file:
                file.write('[frames]\n' + str(f))
            p = Popen([
                "echo",
                "System",
            ], stdout=PIPE)
            logs = run([
                "gmx",
                "trjconv",
                "-s",
                topology_filename,
                "-f",
                reduced_trajectory_filename,
                '-o',
                current_frame,
                "-fr",
                frames_ndx,
                '-quiet'
            ], stdin=p.stdout, stdout=PIPE).stdout.decode()
            p.stdout.close()

            yield current_frame

            # Delete current frame files before going for the next frame
            run([
                "rm",
                frames_ndx,
                current_frame,
            ], stdout=PIPE).stdout.decode()

    return frames_generator(), frames_step, frames_count