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
    output_topology_filename : str,
    output_trajectory_filename : str,
    preprocess_protocol : int,
    ) -> str:

    print(' Preprocess protocol: ' + str(preprocess_protocol))

    if preprocess_protocol == 0:
        return

    # Imaging --------------------------------------------------------------------------------------

    # Only fitting (protocol 1) has nothing to do here

    # Single protein
    if preprocess_protocol == 2:

        # Basic imaging
        # '-pbc mol' sets atoms in the same molecule to stay together
        p = Popen([
            "echo",
            "System",
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
            'atom',
            '-center',
            '-quiet'
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()

        # Select the first frame of the recently translated trayectory as the new topology
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
            output_topology_filename,
            '-dump',
            '0',
            '-quiet'
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()

    # Two interacting proteins, which may need to be translated
    elif preprocess_protocol == 3:

        # FIRST STEP
        # Translate all atoms inside the box so the contact zone between both proteins is in the box center
        # WARNING: The vector in the '-trans' option may be different among different trajectories
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

        # Select the first frame of the recently translated trayectory as the new topology
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
            output_topology_filename,
            '-dump',
            '0',
            '-quiet'
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()

    else:
        raise ValueError('There is no process protocol ' + str(preprocess_protocol))

    # -----------------------------------------------------------------------------------------------

    # Perform the '-pbc nojump' to avoid non-sense jumping of any protein or atom
    # This step is mandatory for all protocols
    p = Popen([
        "echo",
        "System",
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "trjconv",
        "-s",
        output_topology_filename,
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
    if preprocess_protocol == 1:
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