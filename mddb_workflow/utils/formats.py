import os
from typing import Optional, Callable, Generator
from inspect import getfullargspec

def get_format (filename : str) -> str:
    """ Get a filename format. """
    if not filename:
        return None
    return filename.split('.')[-1]

def get_format_set_suitable_function (
    available_functions : list[Callable],
    available_request_format_sets : list[dict],
) -> Generator[ Optional[ tuple[ Callable, dict ]], None, None ]:
    """
    Find a function which is suitable for any of the available request "format sets".
    All functions are checked for each request format set before jumping to another and they are evaluated in order.
    A function and new generated format set with formats in common are returned.
    None is returned when there is no suitable function.

    WARNING: Available functions must have the 'format_sets' property.
    Format sets are dictionaries which specify input and output formats.
    Consider function format sets as 'required input' and 'available output'.
    Consider request format sets as 'available input' and 'required output'.
    Both inputs and outputs are dictionaries where keys are function arguments and values are sets of supported formats.
    Alternatively, an argument may have None as value to represent unnecessary requirements or missing availabilities.
    
    An example is shown below:
    {
        'inputs': {
            'input_structure_filename': {'tpr'},
            'input_trajectory_filenames': {'xtc', 'trr'},
        },
        'outputs': {
            'output_trajectory_filename': {'pdb', 'gro'}
        },
    }

    Functions may have multiple format sets since different input formats may lead to different output formats.
    """
    # -------------------------------------------------------------
    # Use this to see what is going on
    # print('REQUEST')
    # print(available_request_format_sets)
    # print('AVAILABLE FUNCTIONS')
    # print([ func.format_sets for func in  available_functions ])
    # -------------------------------------------------------------
    # Try with each request format set
    for request_format_set in available_request_format_sets:
        # Search functions to match formats for every required argument
        for function in available_functions:
            # Test every function format set independently
            for function_format_set in function.format_sets:
                # Check format keys are compatible
                if not check_format_sets_compability(request_format_set, function_format_set):
                    raise RuntimeError(f'Format keys are not compatible with function {function.__name__}')
                # Check the function inputs to be fulfilled by the request inputs
                required_inputs = function_format_set.get('inputs', None)
                available_inputs = request_format_set.get('inputs', None)
                common_inputs = get_common_argument_formats(required_inputs, available_inputs)
                # If any of the common format sets was empty it means formats do not match
                if not common_inputs:
                    continue
                # Check the request outputs to be fulfilled by the function outputs
                required_outputs = request_format_set.get('outputs', None)
                available_outputs = function_format_set.get('outputs', None)
                common_outputs = get_common_argument_formats(required_outputs, available_outputs)
                # If any of the common format sets was empty it means formats do not match
                if not common_outputs:
                    continue
                # Generate a new format set with the common formats for every argument
                common_format_set = { 'inputs': common_inputs, 'outputs': common_outputs }
                # Otherwise we have the function
                yield function, common_format_set

def get_common_argument_formats (required_arguments : dict, available_arguments : dict):
    """
    Get compatible formats between two groups of arguments.
    All required argument formats must be fulfilled by the available argument formats.
    Arguments are defined as dictionaries with argument names as keys and sets of available formats as values.
    e.g.

    {
        'input_structure_filename': {'tpr'},
        'input_trajectory_filenames': {'xtc', 'trr'},
    }
    """
    # If there are not required arguments we return an empty dictionary
    if not required_arguments:
        return {'not_required'}
    # If there are not available arguments we return None
    if not available_arguments:
        return None
    # Set a dictionary with the same keys that the input arguments including only the common formats
    common_argument_formats = {}
    # Iterate over each required argument
    for required_argument, required_formats in required_arguments.items():
        # If there is not format required for this argument we set None for this argument common formats
        if required_formats == None:
            common_argument_formats[required_argument] = None
            continue
        # Get the available formats for the required input
        available_formats = available_arguments.get(required_argument, None) # DANI cuidao aquí, que igual el get no hacía falta
        # If the available arguments are missing this required argument then the function is not compatible
        if available_formats == None:
            return None
        # Find the formats in common between both the required and the available arguments
        common_formats = available_formats.intersection(required_formats)
        # If the common formats set is empty for this required argument then the function is not compatible
        if not common_formats:
            return None
        common_argument_formats[required_argument] = common_formats
    return common_argument_formats


def check_format_sets_compability (request_format_set : dict, function_format_set : dict) -> bool:
    """
    Check two format sets to be compatible.
    Both function and request format sets must match in their requirements.
    i.e. all function format set input arguments must be included in request format set input arguments.
    i.e. all request format set output arguments must be included in function format set output arguments.
    """
    # Check the function inputs keyowrds to exist in the request input arguments
    required_inputs = function_format_set.get('inputs', None)
    if required_inputs:
        required_input_arguments = required_inputs.keys()
        available_inputs = request_format_set.get('inputs', None)
        if not available_inputs:
            print('ERROR: Missing inputs')
            return False
        available_input_arguments = available_inputs.keys()
        for argument in required_input_arguments:
            if argument not in available_input_arguments:
                print('ERROR: Missing ' + argument + ' argument')
                return False
    # Check the request output keyowrds to exist in the function output arguments
    available_outputs = request_format_set.get('outputs', None)
    if available_outputs:
        required_output_arguments = available_outputs.keys()
        available_outputs = function_format_set.get('outputs', None)
        if not available_outputs:
            print('ERROR: Missing outputs')
            return False
        available_output_arguments = function_format_set['outputs'].keys()
        for argument in required_output_arguments:
            if argument not in available_output_arguments:
                print('ERROR: Missing ' + argument + ' argument')
                return False
    return True

# WARNING: This function makes only sets for those functions whose output can be reused as input of others
# WARNING: i.e. functions which return structure/trajectory files
# WARNING: This function should be called only when "get_format_set_suitable_function" has failed
def get_format_set_suitable_combination (
    available_functions : list[Callable],
    available_request_format_sets : list[dict],
) -> Generator[Optional[tuple[Callable, dict]], None, None]:

    # Try with each request format set
    for request_format_set in available_request_format_sets:

        # Get the required outputs
        required_outputs = request_format_set.get('outputs', None)

        # For each function + format set possibility which is compatible with the required inputs, return the available outputs
        def get_combinations (
            current_functions : list[Callable],
            current_function_common_inputs : list[dict],
            available_inputs : dict,
        ):
            # Search functions to match formats for every required argument
            for function in available_functions:
                # Test every function format set independently
                for function_format_set in function.format_sets:
                    # Check format keys are compatible
                    if not check_format_sets_compability(request_format_set, function_format_set):
                        raise ValueError('Format keys are not compatible with function ' + str(function.__name__))
                    # We add the current function to the list of functions to combine
                    new_functions = [ *current_functions, function ]
                    # Check the function inputs to be fulfilled by the request inputs
                    required_inputs = function_format_set.get('inputs', None)
                    common_inputs = get_common_argument_formats(required_inputs, available_inputs)
                    new_function_common_inputs = [ *current_function_common_inputs, common_inputs ]
                    # If any of the common format sets was empty it means formats do not match
                    if not common_inputs:
                        continue
                    # Check the request outputs to be fulfilled by the function outputs
                    available_outputs = function_format_set.get('outputs', None)
                    common_outputs = get_common_argument_formats(required_outputs, available_outputs)
                    # If any of the common format sets was not empty then we have found a successful combination
                    if common_outputs:
                        # Use this print to display which functions are selected
                        #print([ cf.__name__ + ' from ' + cf.__module__ for cf in new_functions ])
                        yield new_functions, new_function_common_inputs, common_outputs
                    # If we have no common outputs yet we must make another jump
                    # Merge current available inputs with the current function available outputs
                    # The using all these formats as available inputs, try to find another function
                    current_structure_inputs = available_inputs.get('input_structure_filename', set())
                    current_trajectory_inputs = available_inputs.get('input_trajectory_filenames', set())
                    new_structure_inputs = available_outputs.get('output_structure_filename', set())
                    new_trajectory_inputs = available_outputs.get('output_trajectory_filename', set())
                    next_structure_inputs = current_structure_inputs.union(new_structure_inputs)
                    next_trajectory_inputs = current_trajectory_inputs.union(new_trajectory_inputs)
                    # If current function did not add new formats to the already available formats then current step was useless
                    # In this case, stop here
                    if len(next_structure_inputs) == len(current_structure_inputs) and len(next_trajectory_inputs) == len(current_trajectory_inputs):
                        continue
                    # Build the new avaiable inputs dictionary
                    new_available_inputs = {
                        'input_structure_filename': next_structure_inputs,
                        'input_trajectory_filenames': next_trajectory_inputs
                    }
                    # In case we have new available input formats we find a a new function to get the final desired output format
                    for results in get_combinations(new_functions, new_function_common_inputs, new_available_inputs):
                        if results:
                            yield results

        # Get every possible combination after combining the corresponding functions
        first_available_inputs = request_format_set.get('inputs', None)
        for functions, function_common_inputs, last_common_outputs in get_combinations([], [], first_available_inputs):
            # Combine all functions into one single function
            def combined_function (
                input_structure_filename : Optional[str] = None,
                input_trajectory_filenames : Optional[list[str]] = None,
                output_structure_filename : Optional[str] = None,
                output_trajectory_filename : Optional[str] = None
            ):
                auxiliar_filenames = []
                available_structure_filenames = [ input_structure_filename ]
                available_trajectory_filenames = [ input_trajectory_filenames ] # This is a list of lists
                current_input_structure_filename = input_structure_filename
                current_input_trajectory_filenames = input_trajectory_filenames
                current_output_structure_filename = None
                current_output_trajectory_filename = None
                functions_count = len(functions)
                for i, function in enumerate(functions):
                    # Get the next function common inputs in order to know what format we must output
                    next_function_index = i + 1
                    already_existing_structure = None
                    already_existing_trajectories = None
                    if next_function_index < functions_count:
                        next_function_common_inputs = function_common_inputs[next_function_index]
                        # Find the formats for the outputs. Use the first common format to do so
                        # First select the structure format
                        next_function_common_structure_formats = next_function_common_inputs.get('input_structure_filename', None)
                        if next_function_common_structure_formats:
                            output_structure_format = list(next_function_common_structure_formats)[0]
                            # Set the output structure filename
                            # Set the ouput as None if there is a structure with the desired format already
                            # Otherwise, create it using an auxiliar filename
                            already_existing_structure = next(
                                ( structure for structure in available_structure_filenames if get_format(structure) == output_structure_format ),
                                None
                            )
                            if already_existing_structure:
                                current_output_structure_filename = None
                            else:
                                auxiliar_structure_filename = '.structure.' + output_structure_format
                                current_output_structure_filename = auxiliar_structure_filename
                                auxiliar_filenames.append(auxiliar_structure_filename)
                        else:
                            current_output_structure_filename = None
                        # Then select the trajectory format
                        next_function_common_trajectory_formats = next_function_common_inputs.get('input_trajectory_filenames', None)
                        if next_function_common_trajectory_formats:
                            output_trajectory_format = list(next_function_common_trajectory_formats)[0]
                            # Set the output trajectory filenames
                            # Set the ouput as None if there are trajectories with the desired format already
                            # Otherwise, create it using an auxiliar filename
                            already_existing_trajectories = next(
                                ( trajectories for trajectories in available_trajectory_filenames if get_format(trajectories[0]) == output_trajectory_format ),
                                None
                            )
                            if already_existing_trajectories:
                                current_output_trajectory_filename = None
                            else:
                                auxiliar_trajectory_filename = '.trajectory.' + output_trajectory_format
                                current_output_trajectory_filename = auxiliar_trajectory_filename
                                auxiliar_filenames.append(auxiliar_trajectory_filename)
                        else:
                            current_output_trajectory_filename = None
                    # In case this is the last function use the final output filenames
                    else:
                        current_output_structure_filename = output_structure_filename
                        current_output_trajectory_filename = output_trajectory_filename
                    # Set the arguments to be passed to the function
                    # This has to be anticipated since we cannot pass an argument the function does not expect
                    converting_function_arguments = getfullargspec(function)[0]
                    passing_arguments = {}
                    if 'input_structure_filename' in converting_function_arguments:
                        passing_arguments['input_structure_filename'] = current_input_structure_filename
                    if 'input_trajectory_filenames' in converting_function_arguments:
                        passing_arguments['input_trajectory_filenames'] = current_input_trajectory_filenames
                    if 'output_structure_filename' in converting_function_arguments:
                        passing_arguments['output_structure_filename'] = current_output_structure_filename
                    if 'output_trajectory_filename' in converting_function_arguments:
                        passing_arguments['output_trajectory_filename'] = current_output_trajectory_filename
                    # Excute the current function
                    function(**passing_arguments)
                    # Now set the inputs for the next function
                    # Also update the available structure/trajectory files in case a further function wants to reuse them
                    if already_existing_structure:
                        current_input_structure_filename = already_existing_structure
                    else:
                        current_input_structure_filename = current_output_structure_filename
                        available_structure_filenames.append(current_output_structure_filename)
                    if already_existing_trajectories:
                        current_input_trajectory_filenames = already_existing_trajectories
                    else:
                        current_input_trajectory_filenames = [ current_output_trajectory_filename ]
                        available_trajectory_filenames.append([ current_output_trajectory_filename ])
                # Remove auxililar files
                for auxiliar_filename in auxiliar_filenames:
                    os.remove(auxiliar_filename)

            # Set the combined function format set
            combined_format_set = {
                'inputs': first_available_inputs,
                'outputs': last_common_outputs
            }
            combined_function.format_sets = [combined_format_set]

            yield combined_function, combined_format_set


# Structure file formats
def is_pdb (filename : str) -> bool:
    return filename[-4:] == '.pdb'

def is_psf (filename : str) -> bool:
    return filename[-4:] == '.psf'

def is_tpr (filename : str) -> bool:
    return filename[-4:] == '.tpr'

def is_gro (filename : str) -> bool:
    return filename[-4:] == '.gro'

def is_prmtop (filename : str) -> bool:
    return filename[-7:] == '.prmtop'

def is_top (filename : str) -> bool:
    return filename[-4:] == '.top'

# Trajectory file formats

def is_xtc (filename : str) -> bool:
    return filename[-4:] == '.xtc'

def is_dcd (filename : str) -> bool:
    return filename[-4:] == '.dcd'

def is_netcdf (filename : str) -> bool:
    return filename[-3:] == '.nc'

def are_xtc (filenames : list) -> bool:
    return all([ is_xtc(filename) for filename in filenames ])

def are_dcd (filenames : list) -> bool:
    return all([ is_dcd(filename) for filename in filenames ])

def are_netcdf (filenames : list) -> bool:
    return all([ is_netcdf(filename) for filename in filenames ])

# Extra formats logic

# Check if a file may be read by pytraj according to its format
def is_pytraj_supported (filename : str) -> bool:
    return is_prmtop(filename) or is_top(filename) or is_psf(filename)

# From GitHub:
# ParmFormatDict = {
#     "AMBERPARM": AMBERPARM,
#     "PDBFILE": PDBFILEPARM,
#     "MOL2FILE": MOL2FILEPARM,
#     "CHARMMPSF": CHARMMPSF,
#     "CIFFILE": CIFFILE,
#     "GMXTOP": GMXTOP,
#     "SDFFILE": SDFFILE,
#     "TINKER": TINKERPARM,
#     "UNKNOWN_PARM": UNKNOWN_PARM,
# }

def get_pytraj_parm_format (filename : str) -> str:
    """ Get the pytraj format key for the write_parm function for a specific file according to its format. """
    if is_prmtop(filename):
        return 'AMBERPARM'
    if is_psf(filename):
        return 'CHARMMPSF'
    if is_top(filename):
        return 'GMXTOP'
    if is_pdb(filename):
        return 'PDBFILE'
    raise ValueError('The file ' + filename + ' format is not supported')
