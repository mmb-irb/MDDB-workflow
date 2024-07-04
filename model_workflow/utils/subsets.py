from inspect import getfullargspec
from typing import List

from model_workflow.utils.formats import get_format_set_suitable_function
from model_workflow.utils.pyt_spells import get_frames_count
from model_workflow.utils.gmx_spells import get_trajectory_subset as gmx_get_trajectory_subset
from model_workflow.utils.mdt_spells import get_trajectory_subset as mdt_get_trajectory_subset

# Set the functions to perform single frame gettering
subset_functions = [ gmx_get_trajectory_subset, mdt_get_trajectory_subset ]

# Get specific frames from a trajectory
def get_trajectory_subset (
    input_structure_file : 'File',
    input_trajectory_file : 'File',
    output_trajectory_file : 'File',
    start : int = 0,
    end : int = None,
    step : int = 1,
    frames : List[int] = [],
    skip : List[int] = [],
):
    # If there is no end set then the end as the last frame of the simulation
    if end == None:
        end = get_frames_count(input_structure_file, input_trajectory_file)

    # How step and skip are combined is not intuitive and it could be missleading
    # For this reason both arguments are not allowed together
    if step and step != 1 and skip and len(skip) > 0:
        raise SystemExit("Arguments 'step' and 'skip' are not allowed together. Please do it in 2 separated calls.")

    # End must be grater than start
    if end != None and end < start:
        raise SystemExit('End frame must be posterior to start frame')

    # We need an output trajectory filename
    if not output_trajectory_file:
        raise SystemExit('Missing output trajectory filename')

    # In case the frames argument is passed
    if frames and len(frames) > 0:
        print('Specific frames were passed. Other arguments will be ignored (start, end, step and skip)')
        if frames != sorted(frames):
            print('WARNING: Note that the reduced trajectories will keep frames in their original order, not the input order')

    # Get the input formats
    input_structure_format = input_structure_file.format
    input_trajectory_format = input_trajectory_file.format
    output_trajectory_format = output_trajectory_file.format
    format_set = {
        'inputs': {
            'input_structure_filename': None if input_structure_format == None else { input_structure_format },
            'input_trajectory_filename': { input_trajectory_format }
        },
        'outputs': {
            'output_trajectory_filename': { output_trajectory_format }
        }
    }
    # Get a suitable function to do the job
    suitables = next(get_format_set_suitable_function(
        available_functions=subset_functions,
        available_request_format_sets=[format_set],
    ), None)
    if not suitables:
        raise SystemExit('There is no subset function which supports the requested formats')
    suitable_function, formats = suitables
    # Call the subset function
    # Check if the subset function requires the input structure argument before
    suitable_function_keywords = getfullargspec(suitable_function)[0]
    required_structure = 'input_structure_filename' in suitable_function_keywords
    if required_structure:
        suitable_function(input_structure_file.path, input_trajectory_file.path, output_trajectory_file.path, start, end, step, frames, skip)
    else:
        suitable_function(input_trajectory_file.path, output_trajectory_file.path, start, end, step, frames, skip)