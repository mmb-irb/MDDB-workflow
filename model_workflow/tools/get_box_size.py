# This script is used to find out the box size (x, y and z) in a pdb topology
# This process is carried by Gromacs

import os
from subprocess import run, PIPE, Popen

from model_workflow.utils.constants import GROMACS_EXECUTABLE

# Set the box analysis filename
# This analysis is used here to mine the box size data
# It is never used further
box_analysis = 'box.xvg'

# input_topology_filename - The name string of the input topology file (path)
# Tested supported formats are .pdb and .tpr
# input_trajectory_filename - The name string of the input trajectory file (path)
# Tested supported formats are .trr and .xtc


def get_box_size(
    input_topology_filename: str,
    input_trajectory_filename: str,
) -> tuple:

    # Generate the box analysis
    # WARNING: Do not use the first_frame here instead of the trajectory
    # In modified topologies the first frame pdb may have lost box size data
    p = Popen([
        "echo",
        "System",
    ], stdout=PIPE)
    logs = run([
        GROMACS_EXECUTABLE,
        "traj",
        "-s",
        input_topology_filename,
        "-f",
        input_trajectory_filename,
        '-ob',
        box_analysis,
        '-b',
        '0',
        '-e',
        '0',
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE, stderr=PIPE).stdout.decode()
    p.stdout.close()

    if not os.path.exists(box_analysis):
        print(logs)
        raise SystemExit('ERROR: Something went wrong with Gromacs')

    # Read the box analysis and get the desired data
    boxsizex, boxsizey, boxsizez = "", "", ""
    with open(box_analysis, 'r') as file:
        for line in file:
            if line.startswith(("#", "@")) == False:

                # Simulation box 'x' size
                boxsizex = float(line.split()[1])

                # Simulation box 'y' size
                boxsizey = float(line.split()[2])

                # Simulation box 'z' size
                boxsizez = float(line.split()[3])

    # Display it
    print('Box size: (' + str(boxsizex) + ',' +
          str(boxsizey) + ',' + str(boxsizez) + ')')

    # Remove the box analysis
    os.remove(box_analysis)

    # Return gromacs logs
    return (boxsizex, boxsizey, boxsizez)
