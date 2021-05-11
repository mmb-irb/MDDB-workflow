from subprocess import run, PIPE, Popen

import os
import math

from model_workflow.tools.get_frames_count import get_frames_count

# Get a reduced version of the provided trajectory
# Frames are taken along the whole trajectory
# Several analyses use this function since they use a reduced number of frames to work
# The output trajectory filename is set here and returned
# This is because reduced trajectory names must be standard, since they are reused
# If the reduced trajectory already exists do not build it again but return its filename
# In addition returns always the step and the expected final number of frames
def get_reduced_trajectory (
    input_topology_filename : str,
    input_trajectory_filename : str,
    reduced_trajectory_frames_limit : int,
    ) -> list:

    trajectory_frames = get_frames_count(input_topology_filename,input_trajectory_filename)

    # if the trajectory already has the reduced number of frames or less return here
    if reduced_trajectory_frames_limit >= trajectory_frames:
        print('The trajectory to be reduced has already the final number of frames or less')
        output_trajectory_filename = input_trajectory_filename
        step = 1
        frames = trajectory_frames
        return output_trajectory_filename, step, frames

    # Set the reduced trajectory filename
    output_trajectory_filename = 'f' + str(reduced_trajectory_frames_limit) + '.trajectory.xtc'

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
    step = math.ceil(trajectory_frames / reduced_trajectory_frames_limit)

    # Create the reduced trajectory if it does not exist yet
    if not os.path.exists(output_trajectory_filename):
        # Run Gromacs
        p = Popen([
            "echo",
            "System",
        ], stdout=PIPE)
        logs = run([
            "gmx",
            "trjconv",
            "-s",
            input_topology_filename,
            "-f",
            input_trajectory_filename,
            '-o',
            output_trajectory_filename,
            '-skip',
            str(step),
            '-quiet'
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()

    # Calculate also the final number of frames given the current step and return this value
    # WARNING: It may seem that the final number of frames is math.floor(trajectory_frames / step)
    # HOWEVER, the frame 0 also counts so it would be math.floor() + 1
    # HOWEVER, when trajectory_frames / step % 0, the last frame is not returned
    # For this reason, the final number of frames is equal to the ceiling of the division
    frames = math.ceil(trajectory_frames / step)

    # Return gromacs logs
    return output_trajectory_filename, step, frames