# Generic analyses
# Easy and fast trajectory analyses carried by Gromacs

from subprocess import run, PIPE, Popen
from model_workflow.tools.get_reduced_trajectory import get_reduced_trajectory

# RMSD
# 
# Perform the RMSd analysis 
# Use the first trajectory frame in .pdb format as a reference
def rmsd (
    input_reference_filename : str,
    input_trajectory_filename : str,
    group_name : str,
    output_analysis_filename : str,
    frames_limit : int) -> str:

    # Use a reduced trajectory in case the original trajectory has many frames
    reduced_trajectory_filename, step, frames = get_reduced_trajectory(
        input_reference_filename,
        input_trajectory_filename,
        frames_limit,
    )
    
    # Run Gromacs
    p = Popen([
        "echo",
        group_name, # Select group for least squares fit
        group_name, # Select group for RMSD calculation
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "rms",
        "-s",
        input_reference_filename,
        "-f",
        reduced_trajectory_filename,
        '-o',
        output_analysis_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

    # Return gromacs logs
    return logs


# Fluctuation
# 
# Perform the fluctuation analysis
# This analysis also produces a 'noelem' file which is never used
def rmsf (
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_analysis_filename : str,
    output_noelem_filename : str = 'noelem.pdb') -> str:
    
    # Run Gromacs
    p = Popen([
        "echo",
        "Protein",
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "rmsf",
        "-s",
        input_topology_filename,
        "-f",
        input_trajectory_filename,
        '-o',
        output_analysis_filename,
        '-oq',
        output_noelem_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

    # Return gromacs logs
    return logs

# Radius of gyration (Rgyr)
# 
# Perform the RMSd analysis 
# Use the first trajectory frame in .pdb format as a reference
def rgyr (
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_analysis_filename : str,
    frames_limit : int) -> str:

    # Use a reduced trajectory in case the original trajectory has many frames
    reduced_trajectory_filename, step, frames = get_reduced_trajectory(
        input_topology_filename,
        input_trajectory_filename,
        frames_limit,
    )
    
    # Run Gromacs
    p = Popen([
        "echo",
        "Protein",
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "gyrate",
        "-s",
        input_topology_filename,
        "-f",
        reduced_trajectory_filename,
        '-o',
        output_analysis_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

    # Return gromacs logs
    return logs