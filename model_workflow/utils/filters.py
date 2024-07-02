from model_workflow.utils.auxiliar import InputError
from model_workflow.utils.gmx_spells import filter_structure, filter_trajectory
from model_workflow.utils.structures import Structure

# Note tha this function only accepts pdb/gro structures and xtc/trr trajectories
accepted_structure_formats = ['pdb']
accepted_trajectory_formats = ['xtc', 'trr']

# Filter both structure and trajectory
# Note that this function is not format-smart
def filter_atoms (
    input_structure_file : 'File',
    input_trajectory_file : 'File',
    output_structure_file : 'File',
    output_trajectory_file : 'File',
    selection_string : str,
    selection_syntax : str
):  

    # Check formats are as expected
    # Check also input files exist
    if not input_structure_file:
        raise InputError('Missing input structure filename')
    if input_structure_file.format not in accepted_structure_formats:
        raise InputError(f'Not valid input structure format ({input_structure_file.format}). Accepted structure formats: ' + ','.join(accepted_structure_formats))
    if not input_structure_file.exists:
        raise InputError('Missing input file ' + input_structure_file.path)
    if input_trajectory_file:
        if input_trajectory_file.format not in accepted_trajectory_formats:
            raise InputError(f'Not valid input trajectory format ({input_trajectory_file.format}). Accepted trajectory formats: ' + ','.join(accepted_trajectory_formats))
        if not input_trajectory_file.exists:
            raise InputError('Missing input file ' + input_trajectory_file)
    if output_structure_file:
        if output_structure_file.format not in accepted_structure_formats:
            raise InputError(f'Not valid output structure format ({output_structure_file.format}). Accepted structure formats: ' + ','.join(accepted_structure_formats))
    if output_trajectory_file:
        if output_trajectory_file.format not in accepted_trajectory_formats:
            raise InputError(f'Not valid output trajectory format ({output_trajectory_file.format}). Accepted trajectory formats: ' + ','.join(accepted_trajectory_formats))
    if not output_structure_file and not output_trajectory_file:
        raise InputError('Missing output')
    # Check we are not missing additional inputs
    if not selection_string:
        print('Missing input selection string -> Water and counter ions will be filtered')
    elif not selection_syntax:
        raise InputError('Missing input selection syntax')

    # Parse the selection
    structure = Structure.from_pdb_file(input_structure_file.path)
    if selection_string:
        selection = structure.select(selection_string, syntax=selection_syntax)
        # If the selection is empty then war the user
        if not selection:
            raise InputError(f'Selection {selection_string} is empty')
    else:
        # If the selection is missing then filter out water and ions by default
        water_selection = structure.select_water()
        counter_ions_selection = structure.select_counter_ions()
        filter_selection = water_selection + counter_ions_selection
        selection = structure.invert_selection(filter_selection)
        # If the selection is empty then war the user
        if not selection:
            raise InputError('There are no water or counter ions')

    # Run the filters
    if output_structure_file:
        filter_structure(input_structure_file, output_structure_file, selection)
    if input_trajectory_file and output_trajectory_file:
        filter_trajectory(input_structure_file, input_trajectory_file, output_trajectory_file, selection)