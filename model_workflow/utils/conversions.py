import os
from shutil import copyfile
from typing import List, Optional
from inspect import getfullargspec

from model_workflow.utils.formats import get_format_set_suitable_function, get_format_set_suitable_combination
from model_workflow.utils.file import File
from model_workflow.utils.vmd_spells import vmd_to_pdb
from model_workflow.utils.gmx_spells import get_structure, get_structure_alone
from model_workflow.utils.gmx_spells import merge_and_convert_trajectories as gmx_merge_and_convert_trajectories
from model_workflow.utils.mdt_spells import merge_and_convert_trajectories as mdt_merge_and_convert_trajectories
from model_workflow.utils.mdt_spells import merge_and_convert_trajectories_alternative as mdt_merge_and_convert_trajectories_alternative
from model_workflow.utils.mdt_spells import merge_and_convert_trajectories_unefficient as mdt_merge_and_convert_trajectories_unefficient
from model_workflow.utils.vmd_spells import merge_and_convert_trajectories as vmd_merge_and_convert_trajectories
from model_workflow.utils.auxiliar import InputError

# Set functions to performe structure conversions
# These functions must have 'input_structure_filename' and 'output_structure_filename' keywords
# These functions must have the 'format_sets' property
# These functions may have the 'input_trajectory_filename' keyword
structure_converting_functions = [ get_structure, get_structure_alone, vmd_to_pdb ]

# Set functions to performe trajectory conversions
# These functions must have 'input_trajectory_filename' and 'output_trajectory_filepath' keywords
# These functions must have the 'format_sets' property
trajectory_converting_functions = [
    mdt_merge_and_convert_trajectories,
    gmx_merge_and_convert_trajectories,
    mdt_merge_and_convert_trajectories_alternative, # This should only be used in mdcrd to xtc/trr
    vmd_merge_and_convert_trajectories,
    mdt_merge_and_convert_trajectories_unefficient
]

# Handle conversions of different structure and trajectory formats
# Merge multiple input trajectories into one single output trajectory
# Inputs are the original strucutre and/or trajectory files and the list of possible output filenames
# Only one of each group of output filenames will be generated (if possible)
# Return the names of the generated output files
def convert (
    input_structure_filepath :  Optional[str] = '',
    output_structure_filepath : Optional[str] = '',
    input_trajectory_filepaths :  Optional[List[str]] = [],
    output_trajectory_filepath : Optional[str] = ''
):

    # If we have output but not input we must complain
    if output_structure_filepath and not input_structure_filepath:
        raise InputError('Missing input structure')
    if output_trajectory_filepath and not input_trajectory_filepaths or len(input_trajectory_filepaths) == 0:
        raise InputError('Missing input trajectory')

    # If the input trajectory filename is not a list but a single string (which should not happen) then fix it
    if type(input_trajectory_filepaths) == str:
        input_trajectory_filepaths = [input_trajectory_filepaths]

    # Parse input filepaths to actual files
    # Note that this step automatically raise input errors if any extension is not recognized
    input_structure_file = File(input_structure_filepath)
    output_structure_file = File(output_structure_filepath)
    input_trajectory_files = [ File(path) for path in input_trajectory_filepaths ]
    output_trajectory_file = File(output_trajectory_filepath)

    # Check input files to exist
    input_files = [ input_structure_file ] + input_trajectory_files
    for input_file in input_files:
        if input_file and not input_file.exists:
            raise InputError('Missing input file ' + input_file.path)

    # Check all input trajectory formats are the same
    input_trajectory_formats = set([ trajectory_file.format for trajectory_file in input_trajectory_files ])
    if len(input_trajectory_formats) > 1:
        raise InputError('Input trajectories must have the same format')
        
    # Get the first trajectory as a sample for those processes which do not require the whole trajectory
    trajectory_sample = input_trajectory_files[0] if len(input_trajectory_files) > 0 else File(None)

    # Check if any input file has an non-standard extension of a supported format
    # If so then we create a symlink with the standard extension
    # Save created symlinks to remove them at then of the process
    symlink_files = []
    if input_structure_file and input_structure_file.extension != input_structure_file.format:
        input_structure_file = input_structure_file.get_standard_file()
        symlink_files.append(input_structure_file)
    if trajectory_sample and trajectory_sample.extension != trajectory_sample.format:
        input_trajectory_files = [ trajectory_file.get_standard_file() for trajectory_file in input_trajectory_files ]
        symlink_files += input_trajectory_files
        trajectory_sample = input_trajectory_files[0]

    # Get file formats
    input_structure_format = input_structure_file.format
    output_structure_format = output_structure_file.format
    input_trajectory_format = trajectory_sample.format
    output_trajectory_format = output_trajectory_file.format

    # Convert the structure
    # Do it inside a function just to return as soon as we are done
    def convert_structure ():
        # If there is no output filename it means we have nothing to do here
        if not output_structure_file:
            return
        # If the input and output names match then we are done
        if input_structure_file.path == output_structure_file.path:
            return
        # If input and output formats are the same then just copy the file with the new name
        if input_structure_format == output_structure_format:
            copyfile(input_structure_file.path, output_structure_file.path)
            return
        print('Getting structure in ' + output_structure_format + ' format from ' + input_structure_format + ' file')
        # Otherwise, we must convert
        # Choose the right conversion function according to input and output formats
        request_format_set = {
            'inputs': {
                'input_structure_filename': { input_structure_format },
                'input_trajectory_filename': { input_trajectory_format }
            },
            'outputs': {
                'output_structure_filename': { output_structure_format }
            }
        }
        suitable = next(get_format_set_suitable_function(
            available_functions=structure_converting_functions,
            available_request_format_sets=[request_format_set],
        ), None)
        # If there is no function to handle this specific conversion we stop here
        if not suitable:
            raise InputError('Conversion from ' + input_structure_format +
                ' to ' + output_structure_format + ' is not supported')
        converting_function, formats = suitable
        # Find the function keywords
        # This is important since some functions may need a trajectory input in addition
        converting_function_keywords = getfullargspec(converting_function)[0]
        required_trajectory = 'input_trajectory_filename' in converting_function_keywords
        if required_trajectory:
            if len(input_trajectory_files) == 0:
                raise InputError('The structure input format ' + input_structure_format +
                ' is missing coordinates and the output format ' + output_structure_format +
                ' needs them. An input trajectory file is required.')
            converting_function(
                input_structure_filename=input_structure_file.path,
                input_trajectory_filename=trajectory_sample.path,
                output_structure_filename=output_structure_file.path
            )
        else:
            converting_function(
                input_structure_filename=input_structure_file.path,
                output_structure_filename=output_structure_file.path
            )
    convert_structure()

    def convert_trajectory ():
        # If there is no output filename it means we have nothing to do here
        if not output_trajectory_file:
            return
        # If the input and output names match then we are done
        trajectory_files_count = len(input_trajectory_files)
        if trajectory_files_count == 1 and trajectory_sample == output_trajectory_file:
            return
        # If there is only 1 input trajectory and it has the same format that the output then just copy the file with the new name
        if trajectory_files_count == 1 and input_trajectory_format == output_trajectory_format:
            copyfile(trajectory_sample.path, output_trajectory_file.path)
            return
        print('Converting trajectory format from ' + input_trajectory_format + ' to ' + output_trajectory_format)
        # Otherwise, we must convert
        # Choose the right conversion function according to input and output formats
        request_format_set = {
            'inputs': {
                'input_structure_filename': { input_structure_format },
                'input_trajectory_filenames': { input_trajectory_format }
            },
            'outputs': {
                'output_trajectory_filename': { output_trajectory_format }
            }
        }
        suitable = next(get_format_set_suitable_function(
            available_functions=trajectory_converting_functions,
            available_request_format_sets=[request_format_set],
        ), None)
        # If there is no function to handle this specific conversion we try to combine several functions in order to do it
        if not suitable:
            print('WARNING: There is no function to do the conversion directly. Trying to combine multiple functions...')
            suitable = next(get_format_set_suitable_combination(
                available_functions=trajectory_converting_functions,
                available_request_format_sets=[request_format_set],
            ), None)
        # If there is no function to handle this specific conversion we stop here
        if not suitable:
            raise InputError('Conversion from ' + input_trajectory_format +
                ' to ' + output_trajectory_format + ' is not supported')
        converting_function, formats = suitable
        # Get the input structure expected format
        expected_input_structure_formats = formats['inputs'].get('input_structure_filename', False)
        # Get the absolute paths of input trajectory files
        trajectory_filepaths = [ trajectory_file.path for trajectory_file in input_trajectory_files ]
        # If the function expects any fromat then pass the structure
        if expected_input_structure_formats:
            converting_function(
                input_structure_filename=input_structure_file.path,
                input_trajectory_filenames=trajectory_filepaths,
                output_trajectory_filename=output_trajectory_file.path
            )
        # If the function expects None then pass None
        elif expected_input_structure_formats == None:
            converting_function(
                input_structure_filename=None,
                input_trajectory_filenames=trajectory_filepaths,
                output_trajectory_filename=output_trajectory_file.path
            )
        # If the function has not the input structure argument then do not pass it
        else:
            converting_function(
                input_trajectory_filenames=trajectory_filepaths,
                output_trajectory_filename=output_trajectory_file.path
            )
    convert_trajectory()

    # Remove generated symlinks, if any
    for symlink_file in symlink_files:
        symlink_file.remove()