from os import remove, rename
from os.path import exists
from shutil import copyfile
from subprocess import run, PIPE, Popen
from re import search

from model_workflow.utils.constants import GROMACS_EXECUTABLE, GREY_HEADER, COLOR_END
from model_workflow.utils.file import File
from model_workflow.utils.type_hints import *

# Get the first frame from a trajectory
def get_first_frame (input_structure_filename : str, input_trajectory_filename : str, output_frame_filename : str):
    # Run Gromacs
    if input_structure_filename:
        p = Popen([
            "echo",
            "System",
        ], stdout=PIPE)
        process = run([
            GROMACS_EXECUTABLE,
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
        ], stdin=p.stdout, stdout=PIPE, stderr=PIPE)
    else:
        process = run([
            GROMACS_EXECUTABLE,
            "trjconv",
            "-f",
            input_trajectory_filename,
            "-o",
            output_frame_filename,
            "-dump",
            "0",
            "-quiet"
        ], stdout=PIPE, stderr=PIPE)
    # Make the process run as logs are saved and decoded
    logs = process.stdout.decode()
    # If output has not been generated then warn the user
    if not exists(output_frame_filename):
        print(logs)
        error_logs = process.stderr.decode()
        print(error_logs)
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
            GROMACS_EXECUTABLE,
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
            GROMACS_EXECUTABLE,
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
        GROMACS_EXECUTABLE,
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
    input_structure_file : 'File',
    output_structure_file : 'File',
    input_selection : 'Selection'
):
    # Generate a ndx file with the desired selection
    filter_selection_name = 'filter'
    filter_index_content = input_selection.to_ndx(selection_name=filter_selection_name)
    filter_index_filename = '.filter.ndx'
    with open(filter_index_filename, 'w') as file:
        file.write(filter_index_content)

    # Filter the structure
    print(GREY_HEADER, end='\r')
    p = Popen([
        "echo",
        filter_selection_name,
    ], stdout=PIPE)
    logs = run([
        GROMACS_EXECUTABLE,
        "editconf",
        "-f",
        input_structure_file.path,
        '-o',
        output_structure_file.path,
        '-n',
        filter_index_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()
    print(COLOR_END, end='\r')

    # Check the output file exists at this point
    # If not then it means something went wrong with gromacs
    if not output_structure_file.exists:
        print(logs)
        raise SystemExit('Something went wrong with Gromacs')

    # Cleanup the index file
    remove(filter_index_filename)

filter_structure.format_sets = [
    {
        'inputs': {
            'input_structure_file': {'pdb', 'gro'},
        },
        'outputs': {
            'output_structure_file': {'pdb', 'gro'}
        }
    }
]

# Filter trajectory atoms
def filter_trajectory (
    input_structure_file : 'File',
    input_trajectory_file : 'File',
    output_trajectory_file : 'File',
    input_selection : 'Selection'
):
    # Generate a ndx file with the desired selection
    filter_selection_name = 'filter'
    filter_index_content = input_selection.to_ndx(selection_name=filter_selection_name)
    filter_index_filename = '.filter.ndx'
    with open(filter_index_filename, 'w') as file:
        file.write(filter_index_content)

    # Filter the trajectory
    print(GREY_HEADER, end='\r')
    p = Popen([
        "echo",
        filter_selection_name,
    ], stdout=PIPE)
    logs = run([
        GROMACS_EXECUTABLE,
        "trjconv",
        "-s",
        input_structure_file.path,
        "-f",
        input_trajectory_file.path,
        '-o',
        output_trajectory_file.path,
        '-n',
        filter_index_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()
    print(COLOR_END, end='\r')

    # Check the output file exists at this point
    # If not then it means something went wrong with gromacs
    if not output_trajectory_file.exists:
        print(logs)
        raise SystemExit('Something went wrong with Gromacs')

    # Cleanup the index file
    remove(filter_index_filename)

filter_trajectory.format_sets = [
    {
        'inputs': {
            'input_structure_file': {'tpr', 'pdb', 'gro'},
            'input_trajectory_file': {'xtc', 'trr'}
        },
        'outputs': {
            'output_trajectory_file': {'xtc', 'trr'}
        }
    }
]

# Set a regular expression to further mine data from gromacs logs
GROMACS_SYSTEM_ATOMS_REGEX = r'System\) has[ ]+([0-9]*) elements'

# Mine system atoms count from gromacs logs
def mine_system_atoms_count (logs : str) -> int:
    system_atoms_match = search(GROMACS_SYSTEM_ATOMS_REGEX, logs)
    if not system_atoms_match:
        raise ValueError('Failed to mine Gromacs error logs')
    return int(system_atoms_match[1])

# Count TPR atoms
def get_tpr_atom_count (tpr_filepath : str) -> int:
    # Run Gromacs only to see the number of atoms in the TPR
    p = Popen([ "echo", "whatever" ], stdout=PIPE)
    process = run([
        GROMACS_EXECUTABLE,
        "convert-tpr",
        "-s",
        tpr_filepath,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE, stderr=PIPE)
    error_logs = process.stderr.decode()
    p.stdout.close()
    # Mine the number of atoms in the system from the logs
    atom_count = mine_system_atoms_count(error_logs)
    return atom_count

# Filter topology atoms
# DANI: Note that a TPR file is not a structure but a topology
# DANI: However it is important that the argument is called 'structure' for the format finder
def filter_tpr (
    input_structure_file : 'File',
    output_structure_file : 'File',
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
        GROMACS_EXECUTABLE,
        "convert-tpr",
        "-s",
        input_structure_file.path,
        '-o',
        output_structure_file.path,
        '-n',
        filter_index_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

    # Check the output file exists at this point
    # If not then it means something went wrong with gromacs
    if not output_structure_file.exists:
        print(logs)
        raise SystemExit('Something went wrong with Gromacs')

    # Cleanup the index file
    remove(filter_index_filename)

filter_tpr.format_sets = [
    {
        'inputs': {
            'input_structure_file': {'tpr'},
        },
        'outputs': {
            'output_structure_file': {'tpr'}
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
        GROMACS_EXECUTABLE,
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

# DANI: No se usa, pero me costó un rato ponerla a punto así que la conservo
# Set a function to read and parse xpm files with a single matrix
# Inspired in https://gromacswrapper.readthedocs.io/en/latest/_modules/gromacs/fileformats/xpm.html#XPM
def parse_xpm (filename : str) -> List[ List[float] ]:
    with open(filename) as file:
        # First lines include metadata such as the title, description and legend
        # Read lines until we find the start of the array
        metadata = [file.readline()]
        while not metadata[-1].startswith("static char *gromacs_xpm[]"):
            metadata.append(file.readline())
        # The next line will contain the dimensions of the matrix
        # e.g. "7 7   14 1",
        dimensions = file.readline().replace('"','').replace(',','').split()
        x_dimension, y_dimension, entity_count, x_stride = [ int(i) for i in dimensions ]
        # Next lines contain every entity definition
        # Every entity has the following:
        # - A letter which is used to refer this entity later in the matrix
        # - A color in #XXXXXX format
        # - The actual value of the entity
        # e.g. "A  c #FFFFFF " /* "0" */,
        # e.g. "B  c #EBEBFF " /* "1" */,
        entities = {}
        for i in range(entity_count):
            line = file.readline()
            entity_id = line[1]
            entity_color = line[6:13]
            entity_value = line[18:].split('"')[1]
            entities[entity_id] = { 'value': entity_value, 'color': entity_color }
        # Next lines are the matrix axis values
        x_axis = []
        y_axis = []
        # Every line has a maximum of 80 labels
        x_lines = math.ceil(x_dimension / 80)
        for l in range(x_lines):
            line = file.readline()[12:-3].split()
            x_axis += [ int(v) for v in line ]
        y_lines = math.ceil(y_dimension / 80)
        for l in range(y_lines):
            line = file.readline()[12:-3].split()
            y_axis += [ int(v) for v in line ]
        # Next lines are the matrix rows
        matrix = []
        for l in range(y_dimension):
            line = file.readline()[1:1+x_dimension]
            row = [ letter for letter in line ]
            matrix.append(row)
        # Check the final matrix size is as expected
        if len(matrix) != y_dimension:
            raise ValueError('Different number of rows than expected')
        sample_row = matrix[-1]
        if len(sample_row) != x_dimension:
            raise ValueError('Different number of columns than expected')
        # Return the output
        return { 'entities': entities, 'x_axis': x_axis, 'y_axis': y_axis, 'matrix': matrix }