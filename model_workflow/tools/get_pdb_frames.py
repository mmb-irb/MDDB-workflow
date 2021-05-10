from model_workflow.tools.get_frames_count import get_frames_count
from model_workflow.tools.get_reduced_trajectory import get_reduced_trajectory

from subprocess import run, PIPE, Popen

# The frames limit is the maximum number of frames to be iterated
# Note that the number of frames iterated may be less
def get_pdb_frames (
    topology_filename : str,
    trajectory_filename : str,
    frames_limit : int
):

    # get the total number of snapshots in the trajectory
    snapshots = get_frames_count(topology_filename, trajectory_filename)

    # Set the frames where we calculate the sasa
    # WARNING: The gromacs '-fr' option counts frames starting at 1, not at 0
    frames_list = range(1, snapshots +1)

    # Set a maximum number of frames
    # If trajectory has more frames than the limit create a reduced trajectory
    reduced_trajectory_filename = trajectory_filename
    frames = None
    if snapshots > frames_limit:
        reduced_trajectory_filename, frames_step, frames_count = get_reduced_trajectory(
            topology_filename,
            trajectory_filename,
            snapshots,
            frames_limit,
        )
        # WARNING: The gromacs '-fr' option counts frames starting at 1, not at 0
        frames_list = range(1, frames_count +1)  # if frames > 1 else [1]
    else:
        frames = snapshots

    def frames_generator():

        # Calculate the sasa for each frame
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