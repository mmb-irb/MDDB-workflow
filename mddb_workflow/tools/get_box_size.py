# This script is used to find out the box size (x, y and z) in a pdb topology
# This process is carried by Gromacs

from os import remove

from mddb_workflow.utils.gmx_spells import run_gromacs

# Set the box analysis filename
# This analysis is used here to mine the box size data
# It is never used further
BOX_ANALYSIS = 'box.xvg'

# input_topology_filename - The name string of the input topology file (path)
# Tested supported formats are .pdb and .tpr
# input_trajectory_filename - The name string of the input trajectory file (path)
# Tested supported formats are .trr and .xtc
# DANI: This is very old code and probably will not work most times
def get_box_size(
    input_topology_filename: str,
    input_trajectory_filename: str,
) -> tuple:
    # Generate the box analysis
    # WARNING: Do not use the first_frame here instead of the trajectory
    # In modified topologies the first frame pdb may have lost box size data
    run_gromacs(f'traj -s {input_topology_filename} -f {input_trajectory_filename} \
                -ob {BOX_ANALYSIS} -b 0 -e 0', user_input = 'System',
                expected_output_filepath = BOX_ANALYSIS)
    # Read the box analysis and get the desired data
    boxsizex, boxsizey, boxsizez = "", "", ""
    with open(BOX_ANALYSIS, 'r') as file:
        for line in file:
            if line.startswith(("#", "@")) == False:
                # Simulation box 'x' size
                boxsizex = float(line.split()[1])
                # Simulation box 'y' size
                boxsizey = float(line.split()[2])
                # Simulation box 'z' size
                boxsizez = float(line.split()[3])
    # Remove the box analysis
    remove(BOX_ANALYSIS)
    # Display it
    print(f'Box size: ({boxsizex},{boxsizey},{boxsizez})')
    # Return gromacs logs
    return (boxsizex, boxsizey, boxsizez)
