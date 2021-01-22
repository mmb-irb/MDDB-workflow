from subprocess import run, PIPE, Popen

import math

# Get a reduced version of the provided trajectory
# This trajectory will have as much frames as specified
# Frames are taken along the whole trajectory
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
    # Since the frame 0 also counts, we must substract 1 to the final frames number
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