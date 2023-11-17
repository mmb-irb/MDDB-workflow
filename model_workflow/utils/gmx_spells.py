from os import remove, rename
from os.path import exists
from shutil import copyfile
from subprocess import run, PIPE, Popen

from typing import List

from model_workflow.utils.file import File

# Get the first frame from a trajectory
def get_first_frame (input_structure_filename : str, input_trajectory_filename : str, output_frame_filename : str):
    # Run Gromacs
    if input_structure_filename:
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
            "-o",
            output_frame_filename,
            "-dump",
            "0",
            "-quiet"
        ], stdin=p.stdout, stderr=PIPE).stderr.decode()
    else:
        logs = run([
            "gmx",
            "trjconv",
            "-f",
            input_trajectory_filename,
            "-o",
            output_first_frame_filename,
            "-dump",
            "0",
            "-quiet"
        ], stderr=PIPE).stderr.decode()
    # If output has not been generated then warn the user
    if not exists(output_frame_filename):
        print(logs)
        raise SystemExit('Something went wrong with Gromacs')

# Set function supported formats
get_first_frame.format_sets = [
    {
        'inputs': {
            'input_structure_filename': {'tpr', 'pdb', 'gro'},
            'input_trajectory_filename': {'xtc', 'trr'}
        },
        'outputs': {
            'output_frame_filename': {'pdb', 'gro'}
        }
    },
    {
        'inputs': {
            'input_structure_filename': None,
            'input_trajectory_filename': {'xtc', 'trr'}
        },
        'outputs': {
            'output_frame_filename': {'xtc', 'trr'}
        }
    },
    {
        'inputs': {
            'input_structure_filename': None,
            'input_trajectory_filename': {'pdb'}
        },
        'outputs': {
            'output_frame_filename': {'pdb', 'xtc', 'trr'}
        }
    },
    {
        'inputs': {
            'input_structure_filename': None,
            'input_trajectory_filename': {'gro'}
        },
        'outputs': {
            'output_frame_filename': {'gro', 'xtc', 'trr'}
        }
    }
]

# Get the structure using the first frame getter function
def get_structure (input_structure_filename : str, input_trajectory_filename : str, output_structure_filename : str):
    get_first_frame(input_structure_filename, input_trajectory_filename, output_structure_filename)
get_structure.format_sets = [
    {
        'inputs': {
            'input_structure_filename': {'tpr', 'pdb', 'gro'},
            'input_trajectory_filename': {'xtc', 'trr'}
        },
        'outputs': {
            'output_structure_filename': {'pdb', 'gro'}
        }
    }
]

# Convert the structure using the first frame getter function (no trajectory is required)
def get_structure_alone (input_structure_filename : str, output_structure_filename : str):
    get_first_frame(input_structure_filename, input_structure_filename, output_structure_filename)
get_structure_alone.format_sets = [
    {
        'inputs': {
            'input_structure_filename': {'pdb', 'gro'},
        },
        'outputs': {
            'output_structure_filename': {'pdb', 'gro'}
        }
    }
]

# Get gromacs supported trajectories merged and converted to a different format
def merge_and_convert_trajectories (input_trajectory_filenames : List[str], output_trajectory_filename : str):
    # Get trajectory formats
    sample_trajectory = input_trajectory_filenames[0]
    sample_trajectory_file = File(sample_trajectory)
    input_trajectories_format = sample_trajectory_file.format
    output_trajectory_file = File(output_trajectory_filename)
    output_trajectory_format = output_trajectory_file.format
    auxiliar_single_trajectory_filename = '.single_trajectory.' + input_trajectories_format
    # If we have multiple trajectories then join them
    if len(input_trajectory_filenames) > 1:
        single_trajectory_filename = auxiliar_single_trajectory_filename
        logs = run([
            "gmx",
            "trjcat",
            "-f",
            *input_trajectory_filenames,
            "-o",
            single_trajectory_filename,
            "-quiet"
        ], stderr=PIPE).stderr.decode()
        # If output has not been generated then warn the user
        if not exists(single_trajectory_filename):
            print(logs)
            raise SystemExit('Something went wrong with Gromacs')
    else:
        single_trajectory_filename = sample_trajectory
    # In case input and output formats are different we must convert the trajectory
    if input_trajectories_format != output_trajectory_format:
        logs = run([
            "gmx",
            "trjconv",
            "-f",
            single_trajectory_filename,
            "-o",
            output_trajectory_filename,
            "-quiet"
        ], stderr=PIPE).stderr.decode()
        # If output has not been generated then warn the user
        if not exists(output_trajectory_filename):
            print(logs)
            raise SystemExit('Something went wrong with Gromacs')
    else:
        copyfile(single_trajectory_filename, output_trajectory_filename)
    # Remove residual files
    if exists(auxiliar_single_trajectory_filename):
        remove(auxiliar_single_trajectory_filename)

merge_and_convert_trajectories.format_sets = [
    {
        'inputs': {
            'input_trajectory_filenames': {'xtc', 'trr'}
        },
        'outputs': {
            'output_trajectory_filename': {'xtc', 'trr'}
        }
    },
    {
        'inputs': {
            'input_trajectory_filenames': {'pdb'}
        },
        'outputs': {
            'output_trajectory_filename': {'pdb', 'xtc', 'trr'}
        }
    },
    {
        'inputs': {
            'input_trajectory_filenames': {'gro'}
        },
        'outputs': {
            'output_trajectory_filename': {'gro', 'xtc', 'trr'}
        }
    }
]

# Get specific frames from a trajectory
def get_trajectory_subset (
    input_trajectory_filename : str,
    output_trajectory_filename : str,
    start : int = 0,
    end : int = None,
    step : int = 1,
    frames : List[int] = [],
    skip : List[int] = []
):
    # Set a list with frame indices from
    output_frames = frames if frames and len(frames) > 0 else [ frame for frame in range(start, end, step) if frame not in skip ]

    # Generate the ndx file to target the desired frames
    auxiliar_ndx_filename = '.frames.ndx'
    generate_frames_ndx(output_frames, auxiliar_ndx_filename)

    # Now run gromacs trjconv command in order to extract the desired frames
    p = Popen([
        "echo",
        "System",
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "trjconv",
        "-f",
        input_trajectory_filename,
        "-o",
        output_trajectory_filename,
        "-fr",
        auxiliar_ndx_filename,
        "-quiet"
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

    # If output has not been generated then warn the user
    if not exists(output_trajectory_filename):
        print(logs)
        raise SystemExit('Something went wrong with Gromacs (main conversion)')

    # Cleanup the auxiliar ndx file
    remove(auxiliar_ndx_filename)


get_trajectory_subset.format_sets = [
    {
        'inputs': {
            'input_trajectory_filename': {'xtc', 'trr'}
        },
        'outputs': {
            'output_trajectory_filename': {'xtc', 'trr'}
        }
    }
]

# Filter trajectory atoms
def filter_structure (
    input_structure_filename : str,
    output_structure_filename : str,
    input_selection : 'Selection'
):
    # Generate a ndx file with the desired selection
    filter_selection_name = 'filter'
    filter_index_content = input_selection.to_ndx(selection_name=filter_selection_name)
    filter_index_filename = '.filter.ndx'
    with open(filter_index_filename, 'w') as file:
        file.write(filter_index_content)

    # Filter the structure
    p = Popen([
        "echo",
        filter_selection_name,
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "editconf",
        "-f",
        input_structure_filename,
        '-o',
        output_structure_filename,
        '-n',
        filter_index_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

    # Check the output file exists at this point
    # If not then it means something went wrong with gromacs
    if not exists(output_structure_filename):
        print(logs)
        raise SystemExit('Something went wrong with Gromacs')

    # Cleanup the index file
    remove(filter_index_filename)

filter_structure.format_sets = [
    {
        'inputs': {
            'input_structure_filename': {'pdb', 'gro'},
        },
        'outputs': {
            'output_structure_filename': {'pdb', 'gro'}
        }
    }
]

# Filter trajectory atoms
def filter_trajectory (
    input_structure_filename : str,
    input_trajectory_filename : str,
    output_trajectory_filename : str,
    input_selection : 'Selection'
):
    # Generate a ndx file with the desired selection
    filter_selection_name = 'filter'
    filter_index_content = input_selection.to_ndx(selection_name=filter_selection_name)
    filter_index_filename = '.filter.ndx'
    with open(filter_index_filename, 'w') as file:
        file.write(filter_index_content)

    # Filter the trajectory
    p = Popen([
        "echo",
        filter_selection_name,
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
        '-n',
        filter_index_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

    # Check the output file exists at this point
    # If not then it means something went wrong with gromacs
    if not exists(output_trajectory_filename):
        print(logs)
        raise SystemExit('Something went wrong with Gromacs')

    # Cleanup the index file
    remove(filter_index_filename)

filter_trajectory.format_sets = [
    {
        'inputs': {
            'input_structure_filename': {'tpr', 'pdb', 'gro'},
            'input_trajectory_filename': {'xtc', 'trr'}
        },
        'outputs': {
            'output_trajectory_filename': {'xtc', 'trr'}
        }
    }
]

# Filter trajectory atoms
def filter_tpr (
    input_tpr_filename : str,
    output_tpr_filename : str,
    input_selection : 'Selection'
):
    # Generate a ndx file with the desired selection
    filter_selection_name = 'filter'
    filter_index_content = input_selection.to_ndx(selection_name=filter_selection_name)
    filter_index_filename = '.filter.ndx'
    with open(filter_index_filename, 'w') as file:
        file.write(filter_index_content)

    # Filter the tpr
    p = Popen([
        "echo",
        filter_selection_name,
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "convert-tpr",
        "-f",
        input_tpr_filename,
        '-o',
        output_tpr_filename,
        '-n',
        filter_index_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

    # Check the output file exists at this point
    # If not then it means something went wrong with gromacs
    if not exists(output_tpr_filename):
        print(logs)
        raise SystemExit('Something went wrong with Gromacs')

    # Cleanup the index file
    remove(filter_index_filename)

filter_tpr.format_sets = [
    {
        'inputs': {
            'input_tpr_filename': {'tpr'},
        },
        'outputs': {
            'output_tpr_filename': {'tpr'}
        }
    }
]

# Join xtc files
# This is a minimal implementation of 'gmx trjcat' used in loops
def merge_xtc_files (current_file : str, new_file : str):
    # If the current file does nt exist then set the new file as the current file
    if not exists(current_file):
        rename(new_file, current_file)
        return
    # Run trjcat
    logs = run([
        "gmx",
        "trjcat",
        "-f",
        new_file,
        current_file,
        '-o',
        current_file,
        '-quiet'
    ],
    stdout=PIPE,
    stderr=PIPE
    ).stdout.decode()

# Generate a ndx file with a selection of frames
def generate_frames_ndx (frames : List[int], filename : str):
    # Add a header 
    content = '[ frames ]\n'
    count = 0
    for frame in frames:
        # Add a breakline each 15 indices
        count += 1
        if count == 15:
            content += '\n'
            count = 0
        # Add a space between indices
        # Atom indices go from 0 to n-1
        # Add +1 to the index since gromacs counts from 1 to n
        content += str(frame + 1) + ' '
    content += '\n'
    # Write the file
    with open(filename, 'w') as file:
        file.write(content)