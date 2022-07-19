# Generic analyses
# Easy and fast trajectory analyses carried by Gromacs

from subprocess import run, PIPE, Popen
import os

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
        "System",
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

    # If the output does not exist at this point it means something went wrong with gromacs
    if not os.path.exists(output_analysis_filename):
        print(logs)
        raise SystemExit('Something went wrong with GROMACS')

    # Return gromacs logs
    return logs