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
    translation : list[float],
    # Note that this is an early provisional structure
    structure : 'Structure',
    # Selection of atoms which are to stay in PBC conditions after imaging
    # Note that this is an early provisional atom selection
    pbc_selection : 'Selection',
    ) -> str:

    print(f' Image: {image} | Fit: {fit}')

    if not image and not fit:
        return

    # Using a TPR topology may provide some advantage in the imaging process
    # This includes connectivity and atom mass data
    # LORE: This was important back in the day to use '-pbc mol' and '-pbc whole'
    # LORE: However we do not rely in these features anymore given they were arbitrary
    is_tpr_available = input_topology_file != MISSING_TOPOLOGY and input_topology_file.format == 'tpr'
    # Check if we have PBC atoms
    # LORE: We were running the '-pbc nojump' step only when no pbc atoms were present in the system
    # LORE: However now we always run a no-jump followed by a '-pbc res' step thus recovering PBC atoms
    has_pbc_atoms = bool(pbc_selection)
    
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
                user_input = 'System', show_error_logs = True)
            latest_trajectory = output_trajectory_file
            
            # Select the first frame of the recently imaged trayectory as the new topology
            reset_structure (
                latest_structure.path,
                latest_trajectory.path,
                output_structure_file.path)
            latest_structure = output_structure_file

        # We assume that at this point (after a possible trnaslation) the first frame is "correct"
        # This means that the different molecules are placed as desired
        # e.g. two protein monomers are together, and not interacting across boundaries
        # To keep things like this we will run a no-jump step
        # DANI: The nojump step may cause artifacts, specially in already pre-imaged simulations
        # DANI: I will disable this part by now
        # print(' Running no-jump')
        # run_gromacs(f'trjconv -s {latest_structure.path} \
        #     -f {latest_trajectory.path} -o {output_trajectory_file.path} \
        #     -pbc nojump', user_input = 'System', show_error_logs = True)
        # latest_trajectory = output_trajectory_file
            
        # Center the non-PBC region while we put all residues in the box
        # This is critical to recover all PBC regios which were diluted because of the no-jump step
        print(' Running ceneterd -pbc res')
        run_gromacs(f'trjconv -s {latest_structure.path} \
            -f {latest_trajectory.path} -o {output_trajectory_file.path} \
            -pbc res -center -n {CENTER_INDEX_FILEPATH}',
            user_input = f'{CENTER_SELECTION_NAME} System', show_error_logs = True)
        latest_trajectory = output_trajectory_file
        
        # Select the first frame of the recently imaged trayectory as the new topology
        reset_structure (latest_structure.path, latest_trajectory.path, output_structure_file.path)
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
            -fit rot+trans -n {CENTER_INDEX_FILEPATH}',
            user_input = f'{CENTER_SELECTION_NAME} System',
            show_error_logs = True)

        # If there is no output structure at this time (i.e. there was no imaging) then create it now
        # Note that the input structure is not necessarily the output structure in this scenario
        # The fit may have been done using the topology and there may be an offset, so better dump it
        if not output_structure_file.exists:
            reset_structure (
                input_structure_file.path,
                output_trajectory_file.path,
                output_structure_file.path)

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
