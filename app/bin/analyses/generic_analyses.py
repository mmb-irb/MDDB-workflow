# Generic analyses
# Easy and fast trajectory analyses carried by Gromacs

from subprocess import run, PIPE, Popen

# RMSD
# 
# Perform the RMSd analysis 
# Use the first trajectory frame in .pdb format as a reference
def rmsd (
    input_first_frame_filename : str,
    input_trajectory_filename : str,
    output_analysis_filename : str) -> str:
    
    # Run Gromacs
    p = Popen([
        "echo",
        "Protein",
        "Protein",
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "rms",
        "-s",
        input_first_frame_filename,
        "-f",
        input_trajectory_filename,
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
    output_analysis_filename : str) -> str:
    
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
        input_trajectory_filename,
        '-o',
        output_analysis_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

    # Return gromacs logs
    return logs