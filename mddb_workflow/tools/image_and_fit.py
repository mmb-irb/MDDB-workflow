# This script is used to fit a trajectory (i.e. eliminate translation and rotation)
# This process is carried by Gromacs

from mddb_workflow.utils.structures import Structure
from mddb_workflow.utils.auxiliar import InputError, MISSING_TOPOLOGY
from mddb_workflow.utils.gmx_spells import run_gromacs
from mddb_workflow.utils.type_hints import *

# Set a pair of independent auxiliar functions to save and recover chains from a pdb file
# WARNING: These functions must be used only when the pdb has not changed in number of atoms

# Get a list with each atom chain from a pdb
def get_chains (pdb_filename : str) -> list:
    structure = Structure.from_pdb_file(pdb_filename)
    chains = structure.chains
    return chains

# Set each atom chain from a pdb from a list
def set_chains (pdb_filename : str, chains : list):
    structure = Structure.from_pdb_file(pdb_filename)
    structure.chains = chains
    structure.generate_pdb_file(pdb_filename)

# Set the default centering/fitting selection (vmd syntax): protein and nucleic acids
CENTER_SELECTION_NAME = 'protein_and_nucleic_acids'
CENTER_INDEX_FILEPATH = '.index.ndx'

# Image and fit a trajectory
# image - Set if it must me imaged
# fit - Set if it must be fitted
# translation - Set how it must be translated during imaging
# pbc_selection - Set a selection to exclude PBC residues from both the centering and the fitting focuses
def image_and_fit (
    # Tested supported formats are .pdb and .tpr
    input_structure_file : 'File',
    # Tested supported formats are .trr and .xtc
    input_trajectory_file : 'File',
    # This must be in .tpr format
    # This is optional if there are no PBC residues
    input_topology_file : Union['File', Exception],
    output_structure_file : 'File',
    output_trajectory_file : 'File',
    image : bool,
    fit : bool,
    translation : list,
    # Note that this is an early provisional structure
    structure : 'Structure',
    # Selection of atoms which are to stay in PBC conditions after imaging
    # Note that this is an early provisional atom selection
    pbc_selection : 'Selection',
    ) -> str:

    print(f' Image: {image} | Fit: {fit}')

    if not image and not fit:
        return

    # In order to run the imaging with PBC residues we need a .tpr file, not just the .pdb file
    # This is because there is a '-pbc mol' step which only works with a .tpr file
    is_tpr_available = input_topology_file != MISSING_TOPOLOGY and input_topology_file.format == 'tpr'
    has_pbc_atoms = bool(pbc_selection)
    if image and not is_tpr_available and has_pbc_atoms:
        raise InputError('In order to image a simulation with PBC residues using Gromacs it is mandatory to provide a .tpr file')

    # If we have coarse grain then we can only fit if we have a tpr file
    # It will not work otherwise since we need atom masses
    # For some image protocols we need to work with pdbs at some points so we directly surrender
    if structure.select_cg():
        if image: raise InputError('We cannot image a coarse grain simulation using Gromacs')
        if not is_tpr_available:
            raise InputError('We cannot fit a coarse grain simulation using Gromacs without a TPR file')

    # First of all save chains
    # Gromacs will delete chains so we need to recover them after
    chains_backup = get_chains(input_structure_file.path)

    # Set a custom index file (.ndx) to select protein and nucleic acids
    # This is useful for some imaging steps (centering and fitting)
    system_selection = structure.select('all', syntax='vmd')
    protein_selection = structure.select_protein()
    nucleic_selection = structure.select_nucleic()
    custom_selection = protein_selection + nucleic_selection
    if not custom_selection:
        raise InputError(f'The default selection to center (protein or nucleic) is empty. Please image your simulation manually.')
    # Exclude PBC residues from this custom selection
    custom_selection -= pbc_selection
    # Convert both selections to a single ndx file which gromacs can read
    system_selection_ndx = system_selection.to_ndx('System')
    selection_ndx = custom_selection.to_ndx(CENTER_SELECTION_NAME)
    with open(CENTER_INDEX_FILEPATH, 'w') as file:
        file.write(system_selection_ndx)
        file.write(selection_ndx)

    # Imaging --------------------------------------------------------------------------------------

    if image:
        print(' Running first imaging step')

        # Check if coordinates are to be translated
        must_translate = translation != [0, 0, 0]

        # If so run the imaging process without the '-center' flag
        if must_translate:
            # WARNING: Adding '-ur compact' makes everything stay the same
            run_gromacs(f'trjconv -s {input_structure_file.path} \
                -f {input_trajectory_file.path} -o {output_trajectory_file.path} \
                -trans {translation[0]} {translation[1]} {translation[2]} \
                -pbc atom', user_input = 'System', show_error_logs = True)
        # If no translation is required and we have a tpr then use it to make molecules whole
        elif is_tpr_available:
            run_gromacs(f'trjconv -s {input_topology_file.path} \
                -f {input_trajectory_file.path} -o {output_trajectory_file.path} \
                -pbc whole', user_input = 'System', show_error_logs = True)
        # Otherwise, center the custom selection
        else:
            run_gromacs(f'trjconv -s {input_structure_file.path} \
                -f {input_trajectory_file.path} -o {output_trajectory_file.path} \
                -pbc atom -center -n {CENTER_INDEX_FILEPATH}',
                user_input = f'{CENTER_SELECTION_NAME} System', show_error_logs = True)

        # Select the first frame of the recently imaged trayectory as the new topology
        reset_structure (input_structure_file.path, output_trajectory_file.path, output_structure_file.path)

        # If there are PBC residues then run a '-pbc mol' to make all residues stay inside the box anytime
        if has_pbc_atoms:
            print(' Fencing PBC molecules back to the box')
            # Run Gromacs
            run_gromacs(f'trjconv -s {input_topology_file.path} \
                -f {input_trajectory_file.path} -o {output_trajectory_file.path} \
                -pbc mol -n {CENTER_INDEX_FILEPATH}',
                user_input = 'System', show_error_logs = True)

            # Select the first frame of the recently translated and imaged trayectory as the new topology
            reset_structure (input_structure_file.path, output_trajectory_file.path, output_structure_file.path)


        # -----------------------------------------------------------------------------------------------

        # If there are no PBC residues then run a '-pbc nojump' to avoid non-sense jumping of any molecule
        else:
            print(' Running no-jump')
            # Run Gromacs
            # Expanding the box may help in some situations but there are secondary effects
            # i.e. -box 999 999 999
            run_gromacs(f'trjconv -s {output_structure_file.path} \
                -f {output_trajectory_file.path} -o {output_trajectory_file.path} \
                -pbc nojump',
                user_input = 'System', show_error_logs = True)

    # Fitting --------------------------------------------------------------------------------------

    # Remove translation and rotation in the trajectory using the custom selection as reference
    # WARNING: Sometimes, simulations with membranes do not require this step and they may get worse when fitted
    if fit:
        print(' Running fitting step')

        # The trajectory to fit is the already imaged trajectory
        # However, if there was no imaging, the trajectory to fit is the input trajectory
        structure_to_fit = output_structure_file if image else input_structure_file
        trajectroy_to_fit = output_trajectory_file if image else input_trajectory_file

        # If we have a TPR then use it here
        # DANI: No estoy seguro de si funcionaría bien después de un imageado
        # DANI: Esto siempre lo he hecho usando la estructura de la primera frame ya imageada
        if is_tpr_available: structure_to_fit = input_topology_file

        # Run Gromacs
        run_gromacs(f'trjconv -s {structure_to_fit.path} \
            -f {trajectroy_to_fit.path} -o {output_trajectory_file.path} \
            -fit rot+trans -n {CENTER_INDEX_FILEPATH}',
            user_input = f'{CENTER_SELECTION_NAME} System',
            show_error_logs = True)

        # If there is no output structure at this time (i.e. there was no imaging) then create it now
        # Note that the input structure is not necessarily the output structure in this scenario
        # The fit may have been done using the topology and there may be an offset, so better dump it
        if not output_structure_file.exists:
            reset_structure (input_structure_file.path, output_trajectory_file.path, output_structure_file.path)

    # Recover chains
    set_chains(output_structure_file.path, chains_backup)

    # Clean up the index file
    # if exists(CENTER_INDEX_FILEPATH):
    #     remove(CENTER_INDEX_FILEPATH)


# Get the first frame of a trajectory
def reset_structure (
    input_structure_filepath : str,
    input_trajectory_filepath : str,
    output_structure_filepath : str,
):
    print(' Reseting structure after imaging step')
    run_gromacs(f'trjconv -s {input_structure_filepath} \
        -f {input_trajectory_filepath} -o {output_structure_filepath} \
        -dump 0', user_input = 'System')
