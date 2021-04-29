from subprocess import run, PIPE, Popen

import math

# Get a reduced version of the provided trajectory
# Frames are taken along the whole trajectory
# Several analyses use this function since they use a reduced number of frames to work
def get_reduced_trajectory (
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_trajectory_filename : str,
    original_frames : int,
    final_frames : int,
    ) -> list:

    if final_frames >= original_frames:
        raise SystemExit('ERROR: The trajectory to be reduced has already the final number of frames or less')

    # Calculate the step between frames in the reduced trajectory to match the final number of frames
    # DANI: Esto no se -> Since the frame 0 also counts, we must substract 1 to the final frames number

    # WARNING: Since the step must be an integer the thorical step must be rounded
    # This means the specified final number of frames may not be accomplished, but it is okey

    # WARNING: Since the step is rounded with the math.ceil function it will always be rounded up
    # This means the final number of frames will be the specified or less

    # CRITICAL WARNING:
    # This formula is exactly the same that the client uses to request stepped frames to the API
    # This means that the client and the workflow are coordinated and these formulas must not change
    # If you decide to change this formula (in both workflow and client)...
    # You will have to run again all the database analyses with reduced trajectories
    #step = int(math.ceil(original_frames / final_frames))
    step = int(math.floor(original_frames / (final_frames - 1)))

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

    # Return gromacs logs
    return step