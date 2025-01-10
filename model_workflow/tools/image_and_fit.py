# This script is used to fit a trajectory (i.e. eliminate translation and rotation)
# This process is carried by Gromacs

from os import remove
from os.path import exists
from subprocess import run, PIPE, Popen

from model_workflow.tools.topology_manager import get_chains, set_chains
from model_workflow.utils.constants import GROMACS_EXECUTABLE, GREY_HEADER, COLOR_END
from model_workflow.utils.structures import Structure
from model_workflow.utils.auxiliar import InputError

# Set the default centering/fitting selection (vmd syntax): protein and nucleic acids
center_selection = 'protein or nucleic'
center_selection_name = 'protein_and_nucleic_acids'
center_selection_filename = '.index.ndx'

# Image and fit a trajectory
# image - Set if it must me imaged
# fit - Set if it must be fitted
# translation - Set how it must be translated during imaging
# input_pbc_selection - Set a selection to exclude PBC residues from both the centering and the fitting focuses
def image_and_fit (
    # Tested supported formats are .pdb and .tpr
    input_structure_file : 'File',
    # Tested supported formats are .trr and .xtc
    input_trajectory_file : 'File',
    # This must be in .tpr format
    # This is optional if there are no PBC residues
    input_topology_file : 'File', 
    output_structure_file : 'File',
    output_trajectory_file : 'File',
    image : bool, fit : bool,
    translation : list,
    input_pbc_selection : str,
    ) -> str:

    print(f' Image: {image} | Fit: {fit}')

    if not image and not fit:
        return

    # First of all save chains
    # Gromacs will delete chains so we need to recover them after
    chains_backup = get_chains(input_structure_file.path)

    # Set a custom index file (.ndx) to select protein and nucleic acids
    # Also exclude PBC residues from this custom selection
    # This is useful for some imaging steps (centering and fitting)
    structure = Structure.from_pdb_file(input_structure_file.path)
    system_selection = structure.select('all', syntax='vmd')
    custom_selection = structure.select(center_selection, syntax='vmd')
    if not custom_selection:
        raise SystemExit(f'The default selection to center ({center_selection}) is empty. Please image your simulation manually.')
    if input_pbc_selection:
        parsed_input_pbc_selection = structure.select(input_pbc_selection, syntax='vmd')
        custom_selection -= parsed_input_pbc_selection
    elif input_pbc_selection == 'guess':
        parsed_input_pbc_selection = structure.select_pbc_guess()
        custom_selection -= parsed_input_pbc_selection
    # Convert both selections to a single ndx file which gromacs can read
    system_selection_ndx = system_selection.to_ndx('System')
    selection_ndx = custom_selection.to_ndx(center_selection_name)
    with open(center_selection_filename, 'w') as file:
        file.write(system_selection_ndx)
        file.write(selection_ndx)

    # In order to run the imaging with PBC residues we need a .tpr file, not just the .pdb file
    # This is because there is a '-pbc mol' step which only works with a .tpr file
    is_tpr_available = input_topology_file and input_topology_file.format == 'tpr'
    has_pbc_atoms = bool(input_pbc_selection)
    if image and has_pbc_atoms and not is_tpr_available:
        raise InputError('In order to image a simulation with PBC residues it is mandatory to provide a .tpr file')

    # Imaging --------------------------------------------------------------------------------------

    if image:

        # Check if coordinates are to be translated 
        must_translate = translation != [0, 0, 0]

        # Set logs grey
        print(GREY_HEADER)

        # If so run the imaging process without the '-center' flag
        if must_translate:
            p = Popen([
                "echo",
                "System",
            ], stdout=PIPE)
            process = run([
                GROMACS_EXECUTABLE,
                "trjconv",
                "-s",
                input_structure_file.path,
                "-f",
                input_trajectory_file.path,
                '-o',
                output_trajectory_file.path,
                '-trans',
                str(translation[0]),
                str(translation[1]),
                str(translation[2]),
                '-pbc',
                'atom',
                #'-ur', # WARNING: This makes everything stay the same
                #'compact', # WARNING: This makes everything stay the same
                '-quiet'
            ], stdin=p.stdout, stdout=PIPE)
        # If no translation is required and we have a tpr then use it to make molecules whole
        elif is_tpr_available:
            p = Popen([
                "echo",
                center_selection_name,
                "System",
            ], stdout=PIPE)
            process = run([
                GROMACS_EXECUTABLE,
                "trjconv",
                "-s",
                input_topology_file.path,
                "-f",
                input_trajectory_file.path,
                '-o',
                output_trajectory_file.path,
                '-pbc',
                'whole',
                '-quiet'
            ], stdin=p.stdout, stdout=PIPE)
        # Otherwise, center the custom selection
        else:
            p = Popen([
                "echo",
                center_selection_name,
                "System",
            ], stdout=PIPE)
            process = run([
                GROMACS_EXECUTABLE,
                "trjconv",
                "-s",
                input_structure_file.path,
                "-f",
                input_trajectory_file.path,
                '-o',
                output_trajectory_file.path,
                '-pbc',
                'atom',
                '-center',
                '-n',
                center_selection_filename,
                '-quiet'
            ], stdin=p.stdout, stdout=PIPE)
        # Run the defined process
        logs = process.stdout.decode()
        p.stdout.close()
        print(COLOR_END)

        # Check the output exists
        if not output_trajectory_file.exists:
            print(logs)
            raise Exception('Something went wrong with Gromacs')

        # Select the first frame of the recently imaged trayectory as the new topology
        reset_structure (input_structure_file.path, output_trajectory_file.path, output_structure_file.path)

        # If there are PBC residues then run a '-pbc mol' to make all residues stay inside the box anytime
        if has_pbc_atoms:
            print(GREY_HEADER)
            # Run Gromacs
            p = Popen([
                "echo",
                "System",
            ], stdout=PIPE)
            logs = run([
                GROMACS_EXECUTABLE,
                "trjconv",
                "-s",
                input_topology_file.path,
                "-f",
                output_trajectory_file.path,
                '-o',
                output_trajectory_file.path,
                '-pbc',
                'mol', # Note that the 'mol' option requires a tpr to be passed
                '-n',
                center_selection_filename,
                '-quiet'
            ], stdin=p.stdout, stdout=PIPE).stdout.decode()
            p.stdout.close()
            print(COLOR_END)

            # Select the first frame of the recently translated and imaged trayectory as the new topology
            reset_structure (input_structure_file.path, output_trajectory_file.path, output_structure_file.path)


        # -----------------------------------------------------------------------------------------------

        # If there are no PBC residues then run a '-pbc nojump' to avoid non-sense jumping of any molecule
        else:
            print(GREY_HEADER)
            # Run Gromacs
            p = Popen([
                "echo",
                "System",
            ], stdout=PIPE)
            logs = run([
                GROMACS_EXECUTABLE,
                "trjconv",
                "-s",
                output_structure_file.path,
                "-f",
                output_trajectory_file.path,
                '-o',
                output_trajectory_file.path,
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
            print(COLOR_END)

    # Fitting --------------------------------------------------------------------------------------

    # Remove translation and rotation in the trajectory using the custom selection as reference
    # WARNING: Sometimes, simulations with membranes do not require this step and they may get worse when fitted
    if fit:

        # The trajectory to fit is the already imaged trajectory
        # However, if there was no imaging, the trajectory to fit is the input trajectory
        structure_to_fit = output_structure_file if image else input_structure_file
        trajectroy_to_fit = output_trajectory_file if image else input_trajectory_file

        print(GREY_HEADER)
        # Run Gromacs
        p = Popen([
            "echo",
            center_selection_name,
            "System",
        ], stdout=PIPE)
        logs = run([
            GROMACS_EXECUTABLE,
            "trjconv",
            "-s",
            structure_to_fit.path,
            "-f",
            trajectroy_to_fit.path,
            '-o',
            output_trajectory_file.path,
            '-fit',
            'rot+trans',
            '-n',
            center_selection_filename,
            '-quiet'
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()
        print(COLOR_END)
        
        # If there is no output structure at this time (i.e. there was no imaging) then set symlink here
        if not output_structure_file.exists:
            output_structure_file.set_symlink_to(input_structure_file)

    # Recover chains
    set_chains(output_structure_file.path, chains_backup)

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
        GROMACS_EXECUTABLE,
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