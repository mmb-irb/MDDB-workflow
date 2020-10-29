# This script is used to get a reduced version of a given trajectory
# The summarized trajectory frames are frames selected along the whole trajectory
# This process is carried by Gromacs

from subprocess import run, PIPE

# input_topology_filename - The name string of the input topology file (path)
# Tested supported formats are .pdb and .tpr
# input_trajectory_filename - The name string of the input trajectory file (path)
# Tested supported formats are .trr and .xtc
# summarized_trajectory_filename - The name string of the output summarized trajectory file (path)
# Tested supported format is .pdb 
def get_summarized_trajectory (
    input_topology_filename : str,
    input_trajectory_filename : str,
    summarized_trajectory_filename : str
    ) -> str:

    # Run Gromacs
    logs = run([
        "echo",
        "System",
        "|",
        "gmx",
        "trjconv",
        "-s",
        input_topology_filename,
        "-f",
        input_trajectory_filenames,
        '-o',
        summarized_trajectory_filename,
        '-dt',
        '10',
        '-quiet'
    ], stdout=PIPE).stdout.decode()

    # Return gromacs logs
    return logs