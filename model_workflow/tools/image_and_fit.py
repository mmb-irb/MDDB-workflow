# This script is used to fit a trajectory (i.e. eliminate translation and rotation)
# This process is carried by Gromacs

import os
from subprocess import run, PIPE, Popen

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
    p = Popen([
        "echo",
        "!\"Water_and_ions\"\nq",
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "make_ndx",
        "-f",
        input_topology_filename,
        '-o',
        indexes,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

    # Imaging --------------------------------------------------------------------------------------

    # Only fitting
    if preprocess_protocol == 1:
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
            output_trajectory_filename,
            '-pbc',
            'nojump',
            '-quiet'
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()

    # Single protein
    if preprocess_protocol == 2:

        # FIRST STEP
        # Process the trajectory to center the protein in the simulation
        # Output only non water and non ion atoms
        # '-pbc nojump' sets atoms to dont jump across the box
        # Use this for a protein inside a membrane
        p = Popen([
            "echo",
            "Protein",
            "!Water_and_ions",
        ], stdout=PIPE)
        logs = run([
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
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()

        # SECOND STEP
        # Adapt the topology to the new reduced trayectory
        # Create the new reduced topology as the current 'topology'
        p = Popen([
            "echo",
            "!Water_and_ions",
        ], stdout=PIPE)
        logs = run([
            "gmx",
            "convert-tpr",
            "-s",
            input_topology_filename,
            '-o',
            input_topology_filename,
            "-n",
            indexes,
            '-quiet'
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()
        
    # Protein inside membrane
    elif preprocess_protocol == 3:

        # FIRST STEP
        # Process the trajectory to center the protein in the simulation
        # Output only non water and non ion atoms
        # '-pbc nojump' sets atoms to dont jump across the box
        # Use this for a protein inside a membrane
        p = Popen([
            "echo",
            "Protein",
            "!Water_and_ions",
        ], stdout=PIPE)
        logs = run([
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
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()

        # SECOND STEP
        # Adapt the topology to the new reduced trayectory
        # Create the new reduced topology as the current 'topology'
        p = Popen([
            "echo",
            "!Water_and_ions",
        ], stdout=PIPE)
        logs = run([
            "gmx",
            "convert-tpr",
            "-s",
            input_topology_filename,
            '-o',
            input_topology_filename,
            "-n",
            indexes,
            '-quiet'
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()

        # FURTHER STEPS
        # '-pbc mol' sets the type of periodic boundary condition treatment and puts the center of mass of molecules in the box
        # '-ur compact' sets the unit cell representation and puts all atoms at the closest distance from the center of the box
        p = Popen([
            "echo",
            "Protein",
            "System",
        ], stdout=PIPE)
        logs = run([
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
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()

    # Two interacting proteins
    elif preprocess_protocol == 4:

        # FIRST STEP
        # Translate all atoms inside the box so the contact zone between both proteins is in the box center
        # WARNING: The vector in the '-trans' option may be different among different trajectories
        p = Popen([
            "echo",
            "!Water_and_ions",
        ], stdout=PIPE)
        logs = run([
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
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()

        # SECOND STEP
        # Adapt the topology to the new reduced trayectory
        p = Popen([
            "echo",
            "!Water_and_ions",
        ], stdout=PIPE)
        logs = run([
            "gmx",
            "convert-tpr",
            "-s",
            input_topology_filename,
            '-o',
            input_topology_filename,
            "-n",
            indexes,
            '-quiet'
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()

        # FURTHER STEPS
        # New topology must me created to conserve the atoms translation before uing the '-pbc nojump' option
        # Otherwise, the initial position of the topology would be used and the translation would be ignored
        
        # Set the new centered topolgy file name
        # The format of the new topology must be .gro, since .tpr is not allowed and .pdb won't work
        # This "temporal" topology is used only in the next steps and it is not required further
        pivot_topology = 'topology.gro'

        # Select the first frame of the recently centered trayectory as the new topology
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
            output_trajectory_filename,
            '-o',
            pivot_topology,
            '-dump',
            '0',
            '-quiet'
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()

        # Finally, perform the '-pbc nojump' to avoid non-sense jumping of any protein or atom
        p = Popen([
            "echo",
            "System",
        ], stdout=PIPE)
        logs = run([
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
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()

    # Fitting --------------------------------------------------------------------------------------

    # The trajectory to fit is the already imaged trajectory
    # However, if there was no imaging, the trajectory to fit is the input trajectory
    trajectroy_to_fit = output_trajectory_filename
    if preprocess_protocol > 0:
        trajectroy_to_fit = input_trajectory_filename

    # Run Gromacs
    p = Popen([
        "echo",
        "Protein",
        "System",
    ], stdout=PIPE)
    logs = run([
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
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

    # ----------------------------------------------------------------------------------------------

    # Return the imaged and fitted trajectory filename
    return output_trajectory_filename