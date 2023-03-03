# This script is used to fit a trajectory (i.e. eliminate translation and rotation)
# This process is carried by Gromacs

from os import remove
from os.path import exists
from subprocess import run, PIPE, Popen

from model_workflow.tools.topology_manager import get_chains, set_chains
from model_workflow.tools.formats import is_tpr

from mdtoolbelt.structures import Structure

# Set the deault centering selection (vmd syntax): protein and nucleic acids
center_selection = 'protein or nucleic'
center_selection_name = 'protein_and_nucleic_acids'
center_selection_filename = '.index.ndx'

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

    # Set a custom index file (.ndx) to select protein and nucleic acids
    # This is useful for some imaging steps
    structure = Structure.from_pdb_file(input_topology_filename)
    selection = structure.select(center_selection, syntax='vmd')
    # Convert the selection to a ndx file which gromacs can read
    selection_ndx = selection.to_ndx(center_selection_name)
    with open(center_selection_filename, 'w') as file:
        file.write(selection_ndx)

    # In order to run the imaging protocol 4 we need a .tpr file, not just the .pdb file
    # This is because there is a '-pbc res' step which only works with a .tpr file
    if preprocess_protocol == 4 and not ( input_tpr_filename and is_tpr(input_tpr_filename) ):
        raise SystemExit('In order to run protocol 4 it is mandatory to provide a .tpr file')

    # Imaging --------------------------------------------------------------------------------------

    # Only fitting (protocol 1) has nothing to do here
    if preprocess_protocol == 1:
        pass

    # Basic imaging
    elif preprocess_protocol == 2:

        # Run Gromacs
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


    # Translated imaging
    # The vector in the '-trans' option is passed as a console argument when launching the workflow
    elif preprocess_protocol == 3:

        # Run Gromacs
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
            'atom',
            #'-ur', # WARNING: This makes everything stay the same
            #'compact', # WARNING: This makes everything stay the same
            '-quiet'
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()

        # Select the first frame of the recently translated and imaged trayectory as the new topology
        reset_structure (input_topology_filename, output_trajectory_filename, output_topology_filename)

    # Residues centered imaging (then the 'nojump' and the fitting steps are skipped)
    elif preprocess_protocol == 4:

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
            input_tpr_filename,
            "-f",
            input_trajectory_filename,
            '-o',
            output_trajectory_filename,
            '-pbc',
            'res', # Note that the 'res' option requires a tpr to be passed
            '-center',
            '-ur',
            'compact',
            '-n',
            center_selection_filename,
            '-quiet'
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()

        # Select the first frame of the recently translated and imaged trayectory as the new topology
        reset_structure (input_topology_filename, output_trajectory_filename, output_topology_filename)

    else:
        raise ValueError('There is no process protocol ' + str(preprocess_protocol))

    # -----------------------------------------------------------------------------------------------

    # Perform the '-pbc nojump' to avoid non-sense jumping of any molecule
    # This step is mandatory for all protocols but the protocol 4 (membranes)
    # Memebranes look better when their lipids jump across boundaries
    if preprocess_protocol != 4:
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
            # Expanding the box may help in some situations
            # However there are secondary effects for the trajectory
            #'-box',
            #'999',
            #'999',
            #'999',
            '-pbc',
            'nojump',
            '-quiet'
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()

    # Fitting --------------------------------------------------------------------------------------

    # Remove translation and rotation in the trajectory using the Protein as reference
    # Simulations with membranes have the mebrane fixed so they do not require this step
    # Actually, fitting the protein usually leads to sudden membrane movements which do not look good
    if preprocess_protocol != 4:

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
            '-n',
            center_selection_filename,
            '-quiet'
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()

    # Recover chains
    set_chains(output_topology_filename, chains_backup)

    # Clean up the index files
    if exists(center_selection_filename):
        remove(center_selection_filename)


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