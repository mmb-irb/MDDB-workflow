# Generic analyses
# Easy and fast trajectory analyses carried by Gromacs

from subprocess import run, PIPE, Popen
from model_workflow.tools.get_reduced_trajectory import get_reduced_trajectory
from model_workflow.tools.xvg_parse import xvg_parse

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

# Look for sudden raises of RMSd values from one frame to another
def rmsd_check (
    input_topology_filename : str,
    input_trajectory_filename : str
    ):

    print('Checking trajectory intergity')

    # Select the whole protein to check the RMSd
    test_group = 'Protein'

    # Set the name for the output of the test rmsd
    test_filename = 'test.rmsd.xvg'

    # Run Gromacs
    p = Popen([
        "echo",
        test_group, # Select group for least squares fit
        test_group, # Select group for RMSD calculation
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "rms",
        "-s",
        input_topology_filename,
        "-f",
        input_trajectory_filename,
        '-o',
        test_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

    # Read the output and do the check
    test = xvg_parse(test_filename, ['Times', 'Values'])
    values = test['Values']
    previous = values[0]
    for i, value in enumerate(values):
        if abs(value - previous) > 1:
            raise ValueError('There is something wrong with RMSd values. Check frame ' + str(i))
        previous = value

    # Remove the test xvg file since it is not required anymore
    run([
        "rm",
        test_filename,
    ], stdout=PIPE).stdout.decode()

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