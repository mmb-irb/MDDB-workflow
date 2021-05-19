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

# Merge multiple xtc files into 1 single xtc file
def merge_xtc_files(input_filenames : list, output_filename : str):

    # Run Gromacs
    logs = run([
        "gmx",
        "trjcat",
        "-f",
        *input_trajectory_filenames,
        "-o",
        output_trajectory_filename,
    ], stdout=PIPE).stdout.decode()

# Get the first frame from a trajectory
# Tested input file formats: xtc
# Tested output file formats: xtc
def get_first_frame(input_filename : str, output_filename : str):

    # Run Gromacs
    logs = run([
        "gmx",
        "trjconv",
        "-f",
        input_filename,
        "-o",
        output_filename,
        "-dump",
        "0",
    ], stdout=PIPE).stdout.decode()