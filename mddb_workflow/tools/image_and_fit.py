# This script is used to fit a trajectory (i.e. eliminate translation and rotation)
# This process is carried by Gromacs

from mddb_workflow.utils.structures import Structure
from mddb_workflow.utils.auxiliar import MISSING_TOPOLOGY, warn
from mddb_workflow.utils.gmx_spells import run_gromacs
from mddb_workflow.utils.type_hints import *

# Set a pair of independent auxiliar functions to save and recover chains from a pdb file
# WARNING: These functions must be used only when the pdb has not changed in number of atoms


def get_chains(pdb_filename: str) -> list:
    """Get a list with each atom chain from a PDB."""
    structure = Structure.from_pdb_file(pdb_filename)
    chains = structure.chains
    return chains


def set_chains(pdb_filename: str, chains: list):
    """Set the chains for a PDB file."""
    structure = Structure.from_pdb_file(pdb_filename)
    structure.chains = chains
    structure.generate_pdb_file(pdb_filename)

# Set some constant values to handle NDX atom index selections
SYSTEM_SELECTION_NAME = 'system'
CENTER_SELECTION_NAME = 'center'
INDEX_FILEPATH = '.index.ndx'


def image_and_fit(
    input_structure_file: 'File',
    input_trajectory_file: 'File',
    input_topology_file: Union['File', Exception],
    output_structure_file: 'File',
    output_trajectory_file: 'File',
    image: bool,
    fit: bool,
    translation: list[float],
    structure: 'Structure',
    pbc_selection: 'Selection',
    ) -> str:
    """Image and fit a trajectory.

    Args:
        input_structure_file:
            The input structure file. Tested supported formats are .pdb and .tpr.
        input_trajectory_file:
            The input trajectory file. Tested supported formats are .xtc and .trr.
        input_topology_file:
            The input topology file. Tested supported format is .tpr.
            This is optional if there are no PBC residues.
        output_structure_file:
            The output structure file.
        output_trajectory_file:
            The output trajectory file.
        image:
            Set if it must me imaged
        fit:
            Set if it must be fitted
        translation:
            Set how it must be translated during imaging
        structure:
            The structure of the system. Note that this is an early provisional structure.
        pbc_selection:
            Set a selection to exclude PBC residues from both the centering and the fitting focuses.
            Selection of atoms which are to stay in PBC conditions after imaging.
            Note that this is an early provisional atom selection.

    """
    print(f' Image: {image} | Fit: {fit}')

    if not image and not fit:
        return

    # Using a TPR topology may provide some advantage in the imaging process
    # This includes connectivity and atom mass data
    # This is important to keep atoms together when using -pbc res, mol or whole
    is_tpr_available = input_topology_file != MISSING_TOPOLOGY and input_topology_file.format == 'tpr'
    if not is_tpr_available:
        warn('Since we are missing a TPR file the automatic imaging protocol is more prone to fail.')
    # Check if we have PBC atoms
    # LORE: We were running the '-pbc nojump' step only when no pbc atoms were present in the system
    # LORE: However now we always run a no-jump followed by a '-pbc res' step thus recovering PBC atoms
    has_pbc_atoms = bool(pbc_selection)

    # First of all save chains
    # Gromacs will delete chains so we need to recover them after
    chains_backup = get_chains(input_structure_file.path)

    # Set a custom index file (.ndx) to select reference atoms for centering and fitting steps
    # We will use as center all atoms which are not under PBC
    system_selection = structure.select('all', syntax='vmd')
    center_selection = structure.invert_selection(pbc_selection)
    if not center_selection:
        warn(f'The default selection to center and fit (non PBC) is empty. These steps will be skipped.')
    # Convert both selections to a single ndx file which gromacs can read
    with open(INDEX_FILEPATH, 'w') as file:
        system_selection_ndx = system_selection.to_ndx(SYSTEM_SELECTION_NAME)
        file.write(system_selection_ndx)
        if center_selection:
            center_selection_ndx = center_selection.to_ndx(CENTER_SELECTION_NAME)
            file.write(center_selection_ndx)

    # Imaging --------------------------------------------------------------------------------------

    if image:
        print(' Running imaging steps')

        # Set the latest structure and trajectory files
        # They are the ones to be used as inputs for the following steps
        latest_structure = input_structure_file
        latest_trajectory = input_trajectory_file

        # Check if coordinates are to be translated
        must_translate = translation != [0, 0, 0]

        # Translate the system if requested
        if must_translate:
            print(' Running translation')
            # WARNING: Adding '-ur compact' makes everything stay the same
            run_gromacs(f'trjconv -s {latest_structure.path} \
                -f {latest_trajectory.path} -o {output_trajectory_file.path} \
                -trans {translation[0]} {translation[1]} {translation[2]}',
                user_input=SYSTEM_SELECTION_NAME, show_error_logs=True)
            latest_trajectory = output_trajectory_file

            # Select the first frame of the recently imaged trayectory as the new structure
            reset_structure(
                latest_structure.path,
                latest_trajectory.path,
                output_structure_file.path)
            latest_structure = output_structure_file

        # Run an intial -pbc res step just to keep residue atoms together before the no-jump step
        # Note that this can be done only if we have the TPR file
        if is_tpr_available:
            print(' Running initial -pbc res')
            run_gromacs(f'trjconv -s {input_topology_file.path} \
                -f {latest_trajectory.path} -o {output_trajectory_file.path} \
                -pbc res -n {INDEX_FILEPATH}',
                user_input=SYSTEM_SELECTION_NAME, show_error_logs=True)
            latest_trajectory = output_trajectory_file
            # Select the first frame of the recently imaged trayectory as the new structure
            reset_structure(
                latest_structure.path,
                latest_trajectory.path,
                output_structure_file.path)
            latest_structure = output_structure_file

        # !!! We assume that at this point (after a possible translation) the first frame is "correct"
        # This means that the different molecules are placed as desired
        # e.g. two protein monomers are together, and not interacting across boundaries
        # To keep things like this we will run a no-jump step
        # DANI: The nojump step may cause artifacts, specially in already pre-imaged simulations
        print(' Running no-jump')
        run_gromacs(f'trjconv -s {latest_structure.path} \
            -f {latest_trajectory.path} -o {output_trajectory_file.path} \
            -pbc nojump', user_input=SYSTEM_SELECTION_NAME, show_error_logs=True)
        latest_trajectory = output_trajectory_file

        # Place all residues in the box again in case we have PBC residues
        # This is critical to recover all PBC residues which were diluted because of the no-jump step
        # Also center the non-PBC region if there is a center selection
        # We try to do both processes in a single step just to be more efficient
        if has_pbc_atoms and center_selection:
            print(' Running center -pbc res')
            run_gromacs(f'trjconv -s {latest_structure.path} \
                -f {latest_trajectory.path} -o {output_trajectory_file.path} \
                -pbc res -center -n {INDEX_FILEPATH}',
                user_input=f'{CENTER_SELECTION_NAME} {SYSTEM_SELECTION_NAME}',
                show_error_logs=True)
            latest_trajectory = output_trajectory_file
        # If there is no center selection then just place all residues in the box
        elif has_pbc_atoms:
            print(' Running -pbc res')
            run_gromacs(f'trjconv -s {latest_structure.path} \
                -f {latest_trajectory.path} -o {output_trajectory_file.path} \
                -pbc res -n {INDEX_FILEPATH}',
                user_input=SYSTEM_SELECTION_NAME, show_error_logs=True)
            latest_trajectory = output_trajectory_file
        # If there are no PBC residues then just center the system
        elif center_selection:
            print(' Running center')
            run_gromacs(f'trjconv -s {latest_structure.path} \
                -f {latest_trajectory.path} -o {output_trajectory_file.path} \
                -center -n {INDEX_FILEPATH}',
                user_input=f'{CENTER_SELECTION_NAME} {SYSTEM_SELECTION_NAME}',
                show_error_logs=True)
            latest_trajectory = output_trajectory_file
        # If there is no center and no PBC then do nothing, but this should never happen

        # Select the first frame of the recently imaged trayectory as the new structure
        reset_structure(latest_structure.path, latest_trajectory.path, output_structure_file.path)
        latest_structure = output_structure_file

    # Fitting --------------------------------------------------------------------------------------

    # Remove translation and rotation in the trajectory using the custom selection as reference
    # WARNING: Sometimes, simulations with membranes do not require this step and they may get worse when fitted
    if fit:
        print(' Running fitting step')

        # The trajectory to fit is the already imaged trajectory
        # However, if there was no imaging, the trajectory to fit is the input trajectory
        structure_to_fit = output_structure_file if image else input_structure_file
        trajectroy_to_fit = output_trajectory_file if image else input_trajectory_file

        # Here we use the structure, not the topology, even if it is a TPR
        # Note that using the topology would make useless the last 'nojump' process

        # Run Gromacs
        run_gromacs(f'trjconv -s {structure_to_fit.path} \
            -f {trajectroy_to_fit.path} -o {output_trajectory_file.path} \
            -fit rot+trans -n {INDEX_FILEPATH}',
            user_input=f'{CENTER_SELECTION_NAME} {SYSTEM_SELECTION_NAME}',
            show_error_logs=True)

        # If there is no output structure at this time (i.e. there was no imaging) then create it now
        # If there was imaging, we reset again to it its not far from the origin
        # Note that the input structure is not necessarily the output structure in this scenario
        # The fit may have been done using the topology and there may be an offset, so better dump it
        reset_structure(
            input_structure_file.path,
            output_trajectory_file.path,
            output_structure_file.path)

    # Recover chains
    set_chains(output_structure_file.path, chains_backup)

    # Clean up the index file
    # if exists(INDEX_FILEPATH):
    #     remove(INDEX_FILEPATH)


def reset_structure(
    input_structure_filepath: str,
    input_trajectory_filepath: str,
    output_structure_filepath: str,
):
    """Get the first frame of a trajectory and use it to reset the structure."""
    print(' Reseting structure after imaging step')
    run_gromacs(f'trjconv -s {input_structure_filepath} \
        -f {input_trajectory_filepath} -o {output_structure_filepath} \
        -dump 0', user_input=SYSTEM_SELECTION_NAME)
