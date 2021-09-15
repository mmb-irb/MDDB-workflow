# This script is used to fit a trajectory (i.e. eliminate translation and rotation)
# This process is carried by Gromacs

import os
from subprocess import run, PIPE, Popen

from model_workflow.tools.topology_manager import get_chains, set_chains
from model_workflow.tools.formats import is_tpr

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
    input_tpr_filename : str, # This is optional for mosts protocols
    output_topology_filename : str,
    output_trajectory_filename : str,
    preprocess_protocol : int,
    translation : list
    ) -> str:

    print(' Preprocess protocol: ' + str(preprocess_protocol))

    if preprocess_protocol == 0:
        return

    # First of all save chains
    # Gromacs will delete chains so we need to recover them after
    chains_backup = get_chains(input_topology_filename)

    # In order to run the imaging protocol 4 we need a .tpr file, not just the .pdb file
    # This is because there is a '-pbc whole' step which only works with a .tpr file
    if preprocess_protocol == 4 and not is_tpr(input_tpr_filename):
        raise ValueError('In order to run protocol 4 it is mandatory to provide a .tpr file')

    # Imaging --------------------------------------------------------------------------------------

    # Only fitting (protocol 1) has nothing to do here
    if preprocess_protocol == 1:
        pass

    # Single protein
    elif preprocess_protocol == 2:

        # Basic imaging
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

        # Select the first frame of the recently imaged trayectory as the new topology
        reset_structure (input_topology_filename, output_trajectory_filename, output_topology_filename)


    # Two interacting proteins, which may need to be translated
    elif preprocess_protocol == 3 or preprocess_protocol == 4:

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
            str(translation[0]),
            str(translation[1]),
            str(translation[2]),
            '-pbc',
            'atom', # 'residue' may work also
            '-ur',
            'compact',
            '-quiet'
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()

        # Select the first frame of the recently translated and imaged trayectory as the new topology
        reset_structure (input_topology_filename, output_trajectory_filename, output_topology_filename)

    else:
        raise ValueError('There is no process protocol ' + str(preprocess_protocol))

    # -----------------------------------------------------------------------------------------------

    # Perform the '-pbc nojump' to avoid non-sense jumping of any protein or atom
    # This step is mandatory for all protocols but the protocol 4 (membranes)
    # Memebranes look better when their lipids jump across boundaries
    # So we make a '-pbc whole' to make sure it is the whole lipid who jumps and not individual atoms
    if preprocess_protocol == 4:

        p = Popen([
            "echo",
            "System",
        ], stdout=PIPE)
        logs = run([
            "gmx",
            "trjconv",
            "-s",
            input_tpr_filename,
            "-f",
            output_trajectory_filename,
            '-o',
            output_trajectory_filename,
            '-pbc',
            'whole',
            '-quiet'
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()
    
        # Select the first frame of the recently made whole trayectory as the new topology
        reset_structure (input_topology_filename, output_trajectory_filename, output_topology_filename)

    else:
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
        'Protein',
        "System",
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "trjconv",
        "-s",
        output_topology_filename,
        "-f",
        trajectroy_to_fit,
        '-o',
        output_trajectory_filename,
        '-fit',
        'rot+trans',
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

    # Recover chains
    set_chains(output_topology_filename, chains_backup)


# Get the first frame of a trajectory
def reset_structure (
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_topology_filename : str,
):
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
        output_topology_filename,
        '-dump',
        '0',
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()