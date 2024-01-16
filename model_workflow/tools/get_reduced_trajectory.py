from subprocess import run, PIPE, Popen

from os import rename
from os.path import exists
from math import ceil

from model_workflow.utils.constants import INCOMPLETE_PREFIX, GREY_HEADER, COLOR_END

# Get a reduced version of the provided trajectory
# Frames are taken along the whole trajectory
# Several analyses use this function since they use a reduced number of frames to work
# The output trajectory filename is set here and returned
# This is because reduced trajectory names must be standard, since they are reused
# If the reduced trajectory already exists do not build it again but return its filename
# In addition returns always the step and the expected final number of frames
def get_reduced_trajectory (
    input_topology_file : 'File',
    input_trajectory_file : 'File',
    snapshots : int,
    reduced_trajectory_frames_limit : int,
    ) -> list:

    # If the trajectory already has the reduced number of frames or less return here
    if reduced_trajectory_frames_limit >= snapshots:
        output_trajectory_filepath = input_trajectory_file.path
        step = 1
        frames = snapshots
        return output_trajectory_filepath, step, frames

    # Set the reduced trajectory filename
    output_trajectory_filename = 'f' + str(reduced_trajectory_frames_limit) + '.trajectory.xtc'
    # Set the output path in the same directory than the input trajectory
    output_trajectory_filepath = input_trajectory_file.basepath + '/' + output_trajectory_filename

    # Set the incomplete reduced trajectory filename and path as well
    # This prevents the workflow from using an incomplete reduced trajectroy in case the workflow was suddenly interrupted
    incomplete_trajectory_filename = INCOMPLETE_PREFIX + output_trajectory_filename
    incomplete_trajectory_filepath = input_trajectory_file.basepath + '/' + incomplete_trajectory_filename

    # Calculate the step between frames in the reduced trajectory to match the final number of frames
    # WARNING: Since the step must be an integer the thorical step must be rounded
    # This means the specified final number of frames may not be accomplished, but it is okey
    # WARNING: Since the step is rounded with the math.ceil function it will always be rounded up
    # This means the final number of frames will be the specified or less
    # CRITICAL WARNING:
    # This formula is exactly the same that the client uses to request stepped frames to the API
    # This means that the client and the workflow are coordinated and these formulas must not change
    # If you decide to change this formula (in both workflow and client)...
    # You will have to run again all the database analyses with reduced trajectories
    step = ceil(snapshots / reduced_trajectory_frames_limit)

    # Create the reduced trajectory if it does not exist yet
    if not exists(output_trajectory_filepath):
        print('Reducing trajectory from ' + str(snapshots) + ' to less than ' + str(reduced_trajectory_frames_limit) + ' frames')
        print(GREY_HEADER)
        # Run Gromacs
        p = Popen([
            "echo",
            "System",
        ], stdout=PIPE)
        logs = run([
            "gmx",
            "trjconv",
            "-s",
            input_topology_file.path,
            "-f",
            input_trajectory_file.path,
            '-o',
            incomplete_trajectory_filepath,
            '-skip',
            str(step),
            '-quiet'
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()
        print(COLOR_END)
        # Check the output file exists at the end
        if not exists(incomplete_trajectory_filepath):
            print(logs)
            print('Something went wrong with GROMACS while reducing the trajectory')
        # Once the trajectory is complete we rename it as complete
        rename(incomplete_trajectory_filepath, output_trajectory_filepath)

    # Calculate also the final number of frames given the current step and return this value
    # WARNING: It may seem that the final number of frames is math.floor(snapshots / step)
    # HOWEVER, the frame 0 also counts so it would be math.floor() + 1
    # HOWEVER, when snapshots / step % 0, the last frame is not returned
    # For this reason, the final number of frames is equal to the ceiling of the division
    frames = ceil(snapshots / step)

    # Return gromacs logs
    return output_trajectory_filepath, step, frames