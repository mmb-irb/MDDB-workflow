# This script is used to merge trajectory files and convert both topology and trajectory files formats
# Everything is carried by VMD software
# Humphrey, W., Dalke, A. and Schulten, K., "VMD - Visual Molecular Dynamics", J. Molec. Graphics, 1996, vol. 14, pp. 33-38. 
# http://www.ks.uiuc.edu/Research/vmd/

import os
import glob

from subprocess import run, PIPE

# Set he path to a script with all commands needed for vmd to parse the topology file
commands_filename = 'parser.vmd'

# input_topology_filename - The name string of the input topology file (path)
# Tested supported formats are .mae (WARNING: May not be fully compatible with the workflow), .psf and .pdb
# input_trajectory_filenames - The name string of the input trajectory filename(s) (path)
# Multiple files may be selected with bash syntax (e.g. *.dcd)
# Tested supported formats are .dcd and .xtc
# output_topology_filename - The name string of the output topology file (path)
# Tested supported format is .pdb
# output_trajectory_filename - The name string of the output trajectory file (path)
# Tested supported format is .trr
def processor (
    input_topology_filename : str,
    input_trajectory_filenames : str,
    output_topology_filename : str,
    output_trajectory_filename : str
    ) -> str:

    # Write all lines but the last line: 'END'
    with open(commands_filename, "w") as file:
        file.write('set select [atomselect top "all"]\n')
        file.write('animate write trr ' + output_trajectory_filename + ' waitfor all sel $select\n')
        file.write('$select frame first\n')
        file.write('$select writepdb ' + output_topology_filename + '\n')
        file.write('exit\n')

    # The subprocess.run command does not deal with the '*' bash syntax to grab multiple files
    # Instead, we have to do it manually
    trajectory_filenames = [ input_trajectory_filenames ]
    if '*' in input_trajectory_filenames:
        trajectory_filenames_list = glob.glob(input_trajectory_filenames)
        trajectory_filenames = sorted(trajectory_filenames_list)

    # Run VMD
    logs = run([
        "vmd",
        input_topology_filename,
        *trajectory_filenames,
        "-e",
        commands_filename,
        "-dispdev",
        "none"
    ], stdout=PIPE).stdout.decode()

    # Check the output files have been created
    # Otherwise print logs and throw an error
    if not ( os.path.exists(output_topology_filename) and os.path.exists(output_trajectory_filename) ) :
        print(logs)
        raise SystemExit("ERROR: Something was wrong with VMD")

    # Return VMD logs
    return logs

