# This script is used to fit a trajectory (i.e. eliminate translation and rotation)
# This process is carried by Gromacs

import os
from subprocess import run, PIPE

# input_topology_filename - The name string of the input topology file (path)
# Tested supported formats are .pdb and .tpr
# input_trajectory_filename - The name string of the input trajectory file (path)
# Tested supported formats are .trr and .xtc
# output_trajectory_filename - The name string of the output trajectory file (path)
# Tested supported format is .xtc
# preprocess_protocol - An int number indicating the image and fit protocols to be done
def image_and_fit (
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_trajectory_filename : str,
    preprocess_protocol : int,
    ) -> str:

    # Create indexes file to select only specific topology regions

    indexes = 'indexes.ndx'
    if required(indexes):
        logs = run([
            "echo",
            "'!\"Water_and_ions\"\nq'",
            "|",
            "gmx",
            "make_ndx",
            "-f",
            input_topology_filename,
            '-o',
            indexes,
            '-quiet'
        ], stdout=PIPE).stdout.decode()

    # Imaging --------------------------------------------------------------------------------------

    # Single protein
    if preprocess_protocol == 2:

        # FIRST STEP
        # Process the trajectory to center the protein in the simulation
        # Output only non water and non ion atoms
        # '-pbc nojump' sets atoms to dont jump across the box
        # Use this for a protein inside a membrane
        logs = run([
            "echo",
            "'Protein'",
            "'!Water_and_ions'",
            "|",
            "gmx",
            "trjconv",
            "-s",
            input_topology_filename,
            "-f",
            input_trajectory_filename,
            '-o',
            output_trajectory_filename,
            "-n",
            indexes,
            '-quiet'
        ], stdout=PIPE).stdout.decode()

        # SECOND STEP
        # Adapt the topology to the new reduced trayectory
        # Create the new reduced topology as the current 'topology'
        logs = run([
            "echo",
            "'!Water_and_ions'",
            "|",
            "gmx",
            "convert-tpr",
            "-s",
            input_topology_filename,
            '-o',
            input_topology_filename,
            "-n",
            indexes,
            '-quiet'
        ], stdout=PIPE).stdout.decode()
        
    # Protein inside membrane
    elif preprocess_protocol == 3:

        # FIRST STEP
        # Process the trajectory to center the protein in the simulation
        # Output only non water and non ion atoms
        # '-pbc nojump' sets atoms to dont jump across the box
        # Use this for a protein inside a membrane
        logs = run([
            "echo",
            "'Protein'",
            "'!Water_and_ions'",
            "|",
            "gmx",
            "trjconv",
            "-s",
            input_topology_filename,
            "-f",
            input_trajectory_filename,
            '-o',
            output_trajectory_filename,
            "-n",
            indexes,
            '-pbc',
            'nojump',
            '-center',
            '-quiet'
        ], stdout=PIPE).stdout.decode()

        # SECOND STEP
        # Adapt the topology to the new reduced trayectory
        # Create the new reduced topology as the current 'topology'
        logs = run([
            "echo",
            "'!Water_and_ions'",
            "|",
            "gmx",
            "convert-tpr",
            "-s",
            input_topology_filename,
            '-o',
            input_topology_filename,
            "-n",
            indexes,
            '-quiet'
        ], stdout=PIPE).stdout.decode()

        # FURTHER STEPS
        # '-pbc mol' sets the type of periodic boundary condition treatment and puts the center of mass of molecules in the box
        # '-ur compact' sets the unit cell representation and puts all atoms at the closest distance from the center of the box
        logs = run([
            "echo",
            "'Protein'",
            "'System'",
            "|",
            "gmx",
            "trjconv",
            "-s",
            input_topology_filename,
            "-f",
            output_trajectory_filename,
            '-o',
            output_trajectory_filename,
            "-n",
            indexes,
            '-pbc',
            'mol',
            '-center',
            '-ur',
            'compact',
            '-quiet'
        ], stdout=PIPE).stdout.decode()

    # Two interacting proteins
    elif preprocess_protocol == 4:

        # FIRST STEP
        # Translate all atoms inside the box so the contact zone between both proteins is in the box center
        # WARNING: The vector in the '-trans' option may be different among different trajectories
        logs = run([
            "echo",
            "'!Water_and_ions'",
            "|",
            "gmx",
            "trjconv",
            "-s",
            input_topology_filename,
            "-f",
            input_trajectory_filename,
            '-o',
            output_trajectory_filename,
            "-n",
            indexes,
            '-trans',
            '0',
            '4',
            '0',
            '-pbc',
            'atom', # 'residue' may work also
            '-ur',
            'compact',
            '-quiet'
        ], stdout=PIPE).stdout.decode()

        # SECOND STEP
        # Adapt the topology to the new reduced trayectory
        logs = run([
            "echo",
            "'!Water_and_ions'",
            "|",
            "gmx",
            "convert-tpr",
            "-s",
            input_topology_filename,
            '-o',
            input_topology_filename,
            "-n",
            indexes,
            '-quiet'
        ], stdout=PIPE).stdout.decode()

        # FURTHER STEPS
        # New topology must me created to conserve the atoms translation before uing the '-pbc nojump' option
        # Otherwise, the initial position of the topology would be used and the translation would be ignored
        
        # Set the new centered topolgy file name
        # The format of the new topology must be .gro, since .tpr is not allowed and .pdb won't work
        # This "temporal" topology is used only in the next steps and it is not required further
        pivot_topology = 'topology.gro'

        # Select the first frame of the recently centered trayectory as the new topology
        logs = run([
            "echo",
            "'System'",
            "|",
            "gmx",
            "trjconv",
            "-s",
            input_topology_filename,
            "-f",
            output_trajectory_filename,
            '-o',
            pivot_topology,
            '-dump',
            '0',
            '-quiet'
        ], stdout=PIPE).stdout.decode()

        # Finally, perform the '-pbc nojump' to avoid non-sense jumping of any protein or atom
        logs = run([
            "echo",
            "'System'",
            "|",
            "gmx",
            "trjconv",
            "-s",
            pivot_topology,
            "-f",
            output_trajectory_filename,
            '-o',
            output_trajectory_filename,
            '-pbc',
            'nojump',
            '-quiet'
        ], stdout=PIPE).stdout.decode()

    # Fitting --------------------------------------------------------------------------------------

    # The trajectory to fit is the already imaged trajectory
    # However, if there was no imaging, the trajectory to fit is the input trajectory
    trajectroy_to_fit = output_trajectory_filename
    if preprocess_protocol == 1:
        trajectroy_to_fit = input_trajectory_filename

    # Run Gromacs
    logs = run([
        "echo",
        "Protein",
        "System",
        "|",
        "gmx",
        "trjconv",
        "-s",
        input_topology_filename,
        "-f",
        trajectroy_to_fit,
        '-o',
        output_trajectory_filename,
        '-fit',
        'rot+trans',
        '-quiet'
    ], stdout=PIPE).stdout.decode()

    # ----------------------------------------------------------------------------------------------

    # Return the imaged and fitted trajectory filename
    return output_trajectory_filename

# Set a function to check if a process must be run (True) or skipped (False)
# i.e. check if the output file already exists and reapeated analyses must be skipped
def required (analysis_filename : str):
    if os.path.exists(analysis_filename) and skip_repeats:
        return False
    return True