# This script is used to get the first trajectory frame
# This process is carried by Gromacs

from subprocess import run, PIPE, Popen

# input_topology_filename - The name string of the input topology file (path)
# Tested supported formats are .pdb and .tpr
# input_trajectory_filename - The name string of the input trajectory file (path)
# Tested supported formats are .trr and .xtc
# first_frame_filename - The name string of the output first frame file (path)
# Tested supported format is .pdb
def get_first_frame (
    input_topology_filename : str,
    input_trajectory_filename : str,
    first_frame_filename : str
    ) -> str:

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
        first_frame_filename,
        '-dump',
        '0',
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

    # Return gromacs logs
    return logs