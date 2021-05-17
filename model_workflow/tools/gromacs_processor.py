import os
from subprocess import run, PIPE, Popen

# Convert a tpr file to a pdb file using gromacs
def tpr2pdb (
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_topology_filename : str):

    # Although only fitting is required, we must make sure atoms won't jump during the fitting
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
        output_topology_filename,
        '-dump',
        '0',
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()
