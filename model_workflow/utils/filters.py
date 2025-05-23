from model_workflow.utils.auxiliar import InputError
from model_workflow.utils.formats import get_format_set_suitable_function
from model_workflow.utils.gmx_spells import filter_structure, filter_trajectory, filter_tpr
from model_workflow.utils.pyt_spells import filter_topology
from model_workflow.utils.structures import Structure
from model_workflow.utils.type_hints import *
from model_workflow.utils.conversions import convert
from model_workflow.utils.file import File
from inspect import getfullargspec

# Set functions to performe structure conversions
# These functions must have 'input_structure_file' and 'output_structure_file' keywords
# These functions must have the 'format_sets' property
# These functions may have the 'input_trajectory_file' keyword
structure_filtering_functions = [ filter_structure, filter_tpr, filter_topology ]

# Set functions to performe trajectory conversions
# These functions must have 'input_trajectory_file' and 'output_trajectory_file' keywords
# These functions must have the 'format_sets' property
trajectory_filtering_functions = [ filter_trajectory ]

# Auxiliar PDB file used to instatiate the reference structure sometimes
AUXILIAR_PDB_FILE = File('.auxiliar.pdb')

# Handle filterings of different structure and trajectory formats
def filter_atoms (
    input_structure_file :  Optional['File'] = None,
    output_structure_file : Optional['File'] = None,
    input_trajectory_file :  Optional['File'] = None,
    output_trajectory_file : Optional['File'] = None,
    selection_string : Optional[str] = None,
    selection_syntax : str = 'vmd'
):

    # If we have no output at all then there is nothing to do
    if not output_structure_file and not output_trajectory_file:
        raise InputError('No output structure neither output trajectory was specified')
    # If we have output but not input we must complain
    if output_structure_file and not input_structure_file:
        raise InputError('Missing input structure when output structure is requested')
    if output_trajectory_file and not input_trajectory_file:
        raise InputError('Missing input trajectory when output trajectory is requested')

    # Check input files to exist
    for input_file in [ input_structure_file, input_trajectory_file ]:
        if input_file and not input_file.exists:
            raise InputError('Could not find input file ' + input_file.path)

    # Check we are not missing additional inputs
    if not selection_string:
        print('Missing input selection string -> Water and counter ions will be filtered')
    elif not selection_syntax:
        raise InputError('Missing input selection syntax')
    
    # We need a PDB to instantiate the reference structure
    # The structure file may also be instantiated from some topology file formats but this is not much reliable
    # If the input structure is not a PDB then we make a conversion to generate it
    structure_filepath = input_structure_file.path
    if input_structure_file.format != 'pdb':
        # To do so we need coordinates so make sure they were passed
        if not input_trajectory_file:
            raise InputError('In order to filter we need coordinates since some atom selections may rely in close atoms.\n' +
                '  Your input structure is a topology thus not having coordinates.\n' +
                '  Please provide a trajectory file as well')
        # Generate a PDB file using both input topology and trajectory
        structure_filepath = AUXILIAR_PDB_FILE.path
        convert(input_structure_filepath=input_structure_file.path,
                output_structure_filepath=structure_filepath,
                input_trajectory_filepaths=input_trajectory_file.path)

    # Parse the selection
    selection = None
    # Hope the input structure is in a supported format
    structure = Structure.from_file(structure_filepath)
    if selection_string:
        # Parse the input selection
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
        
    # Remove the axuiliar PDB file, if exists
    if AUXILIAR_PDB_FILE.exists: AUXILIAR_PDB_FILE.remove()

    # Check if any input file has an non-standard extension of a supported format
    # If so then we create a symlink with the standard extension
    # Save created symlinks to remove them at then of the process
    symlink_files = []
    if input_structure_file and input_structure_file.extension != input_structure_file.format:
        input_structure_file = input_structure_file.get_standard_file()
        symlink_files.append(input_structure_file)
    if input_trajectory_file and input_trajectory_file.extension != input_trajectory_file.format:
        input_trajectory_file = input_trajectory_file.get_standard_file()
        symlink_files.append(input_trajectory_file)

    # Filter the structure if an output structure file was provided
    if output_structure_file:
        print(f'Filtering structure {input_structure_file.path} to {output_structure_file.path}')
        # Find a suitable function according to the formats
        # Note than given the nature of this logic we may encounter a function which converts the file format
        # This function is not intended for that but this is not a problem either
        # If the user specifies different input/output formats and the function can do it then go ahead
        request_format_set = {
            'inputs': {
                'input_structure_file': { input_structure_file.format },
                'input_trajectory_file': { input_trajectory_file.format }
            },
            'outputs': {
                'output_structure_file': { output_structure_file.format }
            }
        }
        suitable = next(get_format_set_suitable_function(
            available_functions=structure_filtering_functions,
            available_request_format_sets=[request_format_set],
        ), None)
        # If there is no function to handle this specific conversion we stop here
        if not suitable:
            if input_structure_file.format == output_structure_file.format:
                raise InputError(f'Filtering structure files in {input_structure_file.format} format is not supported')
            else:
                raise InputError(f'Filtering structure from {input_structure_file.format} to {output_structure_file.format} is not supported')
        filtering_function, formats = suitable
        # Find the function keywords
        # This is important since some functions may need a trajectory input in addition
        filtering_function_keywords = getfullargspec(filtering_function)[0]
        required_trajectory = 'input_trajectory_filename' in filtering_function_keywords
        # If we need a trajectory then pass it as argument as well
        if required_trajectory:
            # Make sure an input trajectory was passed
            if not input_trajectory_file:
                raise InputError('The structure input format ' + input_structure_file.format +
                ' is missing coordinates and the output format ' + output_structure_file.format +
                ' needs them. An input trajectory file is required.')
            filtering_function(
                input_structure_file=input_structure_file,
                input_trajectory_file=input_trajectory_file,
                output_structure_file=output_structure_file,
                input_selection=selection
            )
        # Otherwise use just the structure as input
        else:
            filtering_function(
                input_structure_file=input_structure_file,
                output_structure_file=output_structure_file,
                input_selection=selection
            )

    # Filter the trajectory if an output trajectory file was provided
    if output_trajectory_file:
        print(f'Filtering trajectory {output_trajectory_file.path} to {output_trajectory_file.path}')
        # Otherwise, we must convert
        # Choose the right conversion function according to input and output formats
        request_format_set = {
            'inputs': {
                'input_structure_file': { input_structure_file.format },
                'input_trajectory_file': { input_trajectory_file.format }
            },
            'outputs': {
                'output_trajectory_file': { output_trajectory_file.format }
            }
        }
        suitable = next(get_format_set_suitable_function(
            available_functions=trajectory_filtering_functions,
            available_request_format_sets=[request_format_set],
        ), None)
        # If there is no function to handle the filtering then we stop here
        if not suitable:
            if input_trajectory_file.format == output_trajectory_file.format:
                raise InputError(f'Filtering trajectory files in {input_trajectory_file.format} format is not supported')
            else:
                raise InputError(f'Filtering trajectory from {input_trajectory_file.format} to {output_trajectory_file.format} is not supported')
        filtering_function, formats = suitable
        # Get the input structure expected format
        expected_input_structure_formats = formats['inputs'].get('input_structure_file', False)
        # If the function expects an input structure in any format then pass the structure
        if expected_input_structure_formats:
            filtering_function(
                input_structure_file=input_structure_file,
                input_trajectory_file=input_trajectory_file,
                output_trajectory_file=output_trajectory_file,
                input_selection=selection
            )
        # If the function expects an input structure as None then pass None
        elif expected_input_structure_formats == None:
            filtering_function(
                input_structure_file=None,
                input_trajectory_file=input_trajectory_file,
                output_trajectory_file=output_trajectory_file,
                input_selection=selection
            )
        # If the function does not expect an input structure then do not pass it
        else:
            filtering_function(
                input_trajectory_file=input_trajectory_file,
                output_trajectory_file=output_trajectory_file,
                input_selection=selection
            )

    # Remove generated symlinks, if any
    for symlink_file in symlink_files:
        symlink_file.remove()