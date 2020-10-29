# This script is used to find out the box size (x, y and z) in a pdb topology
# This process is carried by Gromacs

from subprocess import run, PIPE

# Set the box analysis filename
# This analysis is used here to mine the box size data
# It is never used further
box_analysis = 'box.xvg'

# input_topology_filename - The name string of the input topology file (path)
# Tested supported formats are .pdb and .tpr
# input_trajectory_filename - The name string of the input trajectory file (path)
# Tested supported formats are .trr and .xtc
def get_box_size (
    input_topology_filename : str,
    input_trajectory_filename : str,
    ) -> tuple:

    # Check if the box analysis is required
    if required(box_analysis):

        # Generate the box analysis
        # WARNING: Do not use the first_frame here instead of the trajectory
        # In modified topologies the first frame pdb may have lost box size data
        logs = run([
            "echo",
            "System",
            "|",
            "gmx",
            "traj",
            "-s",
            input_topology_filename,
            "-f",
            input_trajectory_filenames,
            '-o',
            box_analysis,
            '-b',
            '0',
            '-e',
            '0',
            '-quiet'
        ], stdout=PIPE).stdout.decode()

    # Read the box analysis and get the desired data
    with open(box_analysis,'r') as file:
        for line in file:
            if line.startswith(("#","@")) == False:

                # Simulation box 'x' size
                boxsizex = float(line.split()[1])

                # Simulation box 'y' size
                boxsizey = float(line.split()[2])

                # Simulation box 'z' size
                boxsizez = float(line.split()[3])

    # Display it
    print ('Box size: (' + str(boxsizex) + ',' + str(boxsizey) + ',' + str(boxsizez) + ')')

    # Remove the box analysis
    run([
        "rm",
        box_analysis,
    ], stdout=PIPE).stdout.decode()

    # Return gromacs logs
    return (boxsizex, boxsizey, boxsizez)

# Set a function to check if a process must be run (True) or skipped (False)
# i.e. check if the output file already exists and reapeated analyses must be skipped
def required (analysis_filename : str):
    if os.path.exists(analysis_filename) and skip_repeats:
        return False
    return True