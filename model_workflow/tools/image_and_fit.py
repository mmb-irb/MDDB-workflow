# This script is used to fit a trajectory (i.e. eliminate translation and rotation)
# This process is carried by Gromacs

from os import remove
from os.path import exists
from subprocess import run, PIPE, Popen

from model_workflow.tools.topology_manager import get_chains, set_chains
from model_workflow.tools.formats import is_tpr

from mdtoolbelt.structures import Structure

# Set the default centering/fitting selection (vmd syntax): protein and nucleic acids
center_selection = 'protein or nucleic'
center_selection_name = 'protein_and_nucleic_acids'
center_selection_filename = '.index.ndx'

# input_structure_filename - The name string of the input topology file (path)
# Tested supported formats are .pdb and .tpr
# input_trajectory_filename - The name string of the input trajectory file (path)
# Tested supported formats are .trr and .xtc
# output_trajectory_filename - The name string of the output trajectory file (path)
# Tested supported format is .xtc
# image - Set if it must me imaged
# fit - Set if it must be fitted
# translation - Set how it must be translated during imaging
# pbc_selection - Set a selection to exclude PBC residues from both the centering and the fitting focuses
def image_and_fit (
    input_structure_filename : str,
    input_trajectory_filename : str,
    input_topology_filename : str, # This is optional if there are no PBC residues
    output_structure_filename : str,
    output_trajectory_filename : str,
    image : bool, fit : bool,
    translation : list,
    pbc_selection : str,
    ) -> str:

    print(' Image: ' + str(image) + ' | Fit: ' + str(fit))

    if not image and not fit:
        return

    # First of all save chains
    # Gromacs will delete chains so we need to recover them after
    chains_backup = get_chains(input_structure_filename)

    # Set a custom index file (.ndx) to select protein and nucleic acids
    # Also exclude PBC residues from this custom selection
    # This is useful for some imaging steps (centering and fitting)
    structure = Structure.from_pdb_file(input_structure_filename)
    system_selection = structure.select('all', syntax='vmd')
    custom_selection = structure.select(center_selection, syntax='vmd')
    if pbc_selection:
        parsed_pbc_selection = structure.select(pbc_selection, syntax='vmd')
        custom_selection -= parsed_pbc_selection
    # Convert both selections to a single ndx file which gromacs can read
    system_selection_ndx = system_selection.to_ndx('System')
    selection_ndx = custom_selection.to_ndx(center_selection_name)
    with open(center_selection_filename, 'w') as file:
        file.write(system_selection_ndx)
        file.write(selection_ndx)

    # In order to run the imaging with PBC residues we need a .tpr file, not just the .pdb file
    # This is because there is a '-pbc mol' step which only works with a .tpr file
    is_tpr_available = input_topology_filename and is_tpr(input_topology_filename)
    has_pbc_residues = bool(pbc_selection)
    if has_pbc_residues and not is_tpr_available:
        raise SystemExit('In order to image a simulation with PBC residues it is mandatory to provide a .tpr file')

    # Imaging --------------------------------------------------------------------------------------

    if image:

        # Check if coordinates are to be translated 
        must_translate = translation != [0, 0, 0]

        # If so run the imaging process without the '-center' flag
        if must_translate:

            # Run Gromacs
            p = Popen([
                "echo",
                "System",
            ], stdout=PIPE)
            logs = run([
                "gmx",
                "trjconv",
                "-s",
                input_structure_filename,
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

        # Otherwise, center the custom selection
        else:

            # Run Gromacs
            p = Popen([
                "echo",
                center_selection_name,
                "System",
            ], stdout=PIPE)
            logs = run([
                "gmx",
                "trjconv",
                "-s",
                input_structure_filename,
                "-f",
                input_trajectory_filename,
                '-o',
                output_trajectory_filename,
                '-pbc',
                'atom',
                '-center',
                '-n',
                center_selection_filename,
                '-quiet'
            ], stdin=p.stdout, stdout=PIPE).stdout.decode()
            p.stdout.close()

        # Select the first frame of the recently imaged trayectory as the new topology
        reset_structure (input_structure_filename, output_trajectory_filename, output_structure_filename)

        # If there are PBC residues then run a '-pbc mol' to make all residues stay inside the box anytime
        if has_pbc_residues:

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
                output_trajectory_filename,
                '-o',
                output_trajectory_filename,
                '-pbc',
                'mol', # Note that the 'mol' option requires a tpr to be passed
                '-n',
                center_selection_filename,
                '-quiet'
            ], stdin=p.stdout, stdout=PIPE).stdout.decode()
            p.stdout.close()

            # Select the first frame of the recently translated and imaged trayectory as the new topology
            reset_structure (input_structure_filename, output_trajectory_filename, output_structure_filename)


        # -----------------------------------------------------------------------------------------------

        # If there are no PBC residues then run a '-pbc nojump' to avoid non-sense jumping of any molecule
        else:
            p = Popen([
                "echo",
                "System",
            ], stdout=PIPE)
            logs = run([
                "gmx",
                "trjconv",
                "-s",
                output_structure_filename,
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

    # Remove translation and rotation in the trajectory using the custom selection as reference
    # WARNING: Sometimes, simulations with membranes do not require this step and they may get worse when fitted
    if fit:

        # The trajectory to fit is the already imaged trajectory
        # However, if there was no imaging, the trajectory to fit is the input trajectory
        trajectroy_to_fit = output_trajectory_filename if image else input_trajectory_filename

        # Run Gromacs
        p = Popen([
            "echo",
            center_selection_name,
            "System",
        ], stdout=PIPE)
        logs = run([
            "gmx",
            "trjconv",
            "-s",
            output_structure_filename,
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
    set_chains(output_structure_filename, chains_backup)

    # Clean up the index file
    # if exists(center_selection_filename):
    #     remove(center_selection_filename)


# Get the first frame of a trajectory
def reset_structure (
    input_structure_filename : str,
    input_trajectory_filename : str,
    output_structure_filename : str,
):
    p = Popen([
        "echo",
        "System",
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "trjconv",
        "-s",
        input_structure_filename,
        "-f",
        input_trajectory_filename,
        '-o',
        output_structure_filename,
        '-dump',
        '0',
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()