from os import remove, rename
from os.path import exists
from shutil import copyfile
from subprocess import run, PIPE, Popen
from re import search

from model_workflow.utils.constants import GROMACS_EXECUTABLE, GREY_HEADER, COLOR_END
from model_workflow.utils.file import File
from model_workflow.utils.type_hints import *

from model_workflow.tools.fix_gromacs_masses import fix_gromacs_masses

# Set a function to call gromacs in a more confortable and standarized way:
# - Standard gromacs executable, which may be provided by the user
# - Gromacs mass fixes by using a custom atommass.dat file with extended atom names
# - Hidden unnecessary output logs and grey-colored necessary ones
# - Missing output checks
# Then return both output and error logs
def run_gromacs(command : str, user_input : Optional[str] = None,
    expected_output_filepath : Optional[str] = 'auto',
    show_output_logs : bool = False, show_error_logs : bool = False) -> Tuple[str, str]:

    # Run a fix for gromacs if not done before
    # Note that this is run always at the moment the code is read, no matter the command or calling origin
    fix_gromacs_masses()

    # In case we have user input we must open a process to then pipe it in the gromacs process
    if user_input:
        # The -e option allows the interpretation of '\', which is critial for the make_ndx command
        user_input_process = Popen([ "echo", "-e", *user_input.split() ], stdout=PIPE)

    # Set the gromacs process
    # Note that at this point the command is not yet run
    process = run([ GROMACS_EXECUTABLE, *command.split(), '-quiet' ],
        stdin = user_input_process.stdout if user_input else None,
        stdout = PIPE, stderr = PIPE if not show_error_logs else None)
    
    # If error is to be shown the color it in grey
    # This is usually used to see progress logs
    if show_error_logs: print(GREY_HEADER, end='\r')

    # Consume the gromacs process output thus running the command
    output_logs = process.stdout.decode()
    # Consume also error logs, but this will not run the command again
    error_logs = process.stderr.decode() if not show_error_logs else None

    # End the grey coloring
    if show_error_logs: print(COLOR_END, end='\r')

    # In case the expected output is set as 'auto' we must guess it from the command
    # Normally output comes after the '-o' option
    if expected_output_filepath == 'auto':
        command_splits = command.split()
        if '-o' not in command_splits: expected_output_filepath = None
        else:
            option_flag_index = command_splits.index('-o')
            expected_output_filepath = command_splits[option_flag_index + 1]

    # If an output file was expected then check it actually exists
    if expected_output_filepath and not exists(expected_output_filepath):
        # If we are missing the expetced output then report it
        print(output_logs)
        print(error_logs)
        # Recreate the exact command
        final_command = f'{GROMACS_EXECUTABLE} {command}'
        if user_input: final_command += f' (with user input "{user_input}")'
        raise SystemExit(f'Something went wrong with Gromacs while running "{command}"')
    
    # If all was good then show final logs but only if it was requested
    if show_output_logs: print(output_logs)

    # Return outputs
    return output_logs, error_logs


# Get the first frame from a trajectory
def get_first_frame (input_structure_filename : str, input_trajectory_filename : str, output_frame_filename : str):
    # Run Gromacs
    run_gromacs(f'trjconv -s {input_structure_filename} -f {input_trajectory_filename} \
                -o {output_frame_filename} -dump 0', user_input = 'System')

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
        run_gromacs(f'trjcat -f {" ".join(input_trajectory_filenames)} \
                    -o {single_trajectory_filename}')
    else:
        single_trajectory_filename = sample_trajectory
    # In case input and output formats are different we must convert the trajectory
    if input_trajectories_format != output_trajectory_format:
        run_gromacs(f'trjconv -f {single_trajectory_filename} \
                    -o {output_trajectory_filename}')
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
    run_gromacs(f'trjconv -f {input_trajectory_filename} -o {output_trajectory_filename} \
                -fr {auxiliar_ndx_filename}', user_input = 'System')

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
    run_gromacs(f'editconf -f {input_structure_file.path} -o {output_structure_file.path} \
                -n {filter_index_filename}', user_input = filter_selection_name)

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
    # Now run gromacs trjconv command in order to extract the desired frames
    run_gromacs(f'trjconv -s {input_structure_file.path} -f {input_trajectory_file.path} \
                -o {output_trajectory_file.path} -n {filter_index_filename}',
                user_input = filter_selection_name)

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
        print(logs)
        raise ValueError('Failed to mine Gromacs error logs')
    return int(system_atoms_match[1])

# Count TPR atoms
def get_tpr_atom_count (tpr_filepath : str) -> int:
    # Make sure the filepath is valid
    if not exists(tpr_filepath):
        raise ValueError('Trying to count atoms from a topology which does not exist')
    # Run Gromacs only to see the number of atoms in the TPR
    output_logs, error_logs = run_gromacs(f'convert-tpr -s {tpr_filepath}', user_input = "whatever")
    # Mine the number of atoms in the system from the logs
    atom_count = mine_system_atoms_count(error_logs)
    return atom_count

# Read a tpr file by converting it to ASCII
def get_tpr_content (tpr_filepath : str) -> Tuple[str, str]:
    # Read the tpr file making a 'dump'
    return run_gromacs(f'dump -s {tpr_filepath}')

# Regular expresion to mine atom charges
GROMACS_TPR_ATOM_CHARGES_REGEX = r"q=([0-9e+-. ]*),"

# Get tpr atom charges
# This works for the new tpr format (tested in 122)
def get_tpr_charges (tpr_filepath : str) -> List[float]:
    # Read the TPR
    tpr_content, tpr_error_logs = get_tpr_content(tpr_filepath)
    # Mine the atomic charges
    charges = []
    # Iterate tpr content lines
    for line in tpr_content.split('\n'):
        # Skip everything which is not atomic charges data
        if line[0:16] != '            atom': continue
        # Parse the line to get only charges
        match = search(GROMACS_TPR_ATOM_CHARGES_REGEX, line)
        if match: charges.append(float(match[1]))
    # If we successfully got atom charges then return them
    if len(charges) > 0: return charges
    # If there are no charges at the end then something went wrong
    print(tpr_error_logs)
    raise RuntimeError(f'Charges extraction from tpr file "{tpr_filepath}" has failed')

# Regular expresion to mine atom bonds
GROMACS_TPR_ATOM_BONDS_REGEX = r"^\s*([0-9]*) type=[0-9]* \((BONDS|CONSTR|CONNBONDS)\)\s*([0-9]*)\s*([0-9]*)$"
# Set a regular expression for SETTLE bonds, used for rigid waters
# ---- From the paper ------------------------------------------------------------------------------------
# The SETTLE can be applied to a four-point water model like TIP4P5 which has the fourth point with
# a certain charge and no mass if the force acting on the fourth point is distributed onto the other three
# points with masses in a reasonable manner.
# S. Miyamoto and P.A. Kollman, “SETTLE: An analytical version of the SHAKE and RATTLE algorithms for rigid
# water models,” J. Comp. Chem., 13 952–962 (1992)
# --------------------------------------------------------------------------------------------------------
# So it may happen that we encounter SETTLE with 4 atoms
# This is not yet supported, but at least we check if this is happening to raise an error when found
GROMACS_TPR_SETTLE_REGEX = r"^\s*([0-9]*) type=[0-9]* \(SETTLE\)\s*([0-9]*)\s*([0-9]*)\s*([0-9]*)\s*([0-9]*)$"

# Get tpr atom bonds
def get_tpr_bonds (tpr_filepath : str) -> List[ Tuple[int, int] ]:
    # Read the TPR
    tpr_content, tpr_error_logs = get_tpr_content(tpr_filepath)
    lines = tpr_content.split('\n')
    # Mine the atomic bonds
    bonds = []
    # Save the moment we already processed a bond of each class (e.g. bonds, constr, connbonds, etc.)
    processed_bond_classes = set()
    # Save already processed bond classes which have an index restart
    # This is a coomon problem in TPR bonding, we do not know the meaning but bonds after the index restart are always wrong
    restarted_bond_classes = set()
    # Iterate tpr content lines to find bonds
    for line in lines:
        # Parse the line to get only atomic bonds
        match = search(GROMACS_TPR_ATOM_BONDS_REGEX, line)
        # If this is not a line with bonding data then skip it
        if not match: continue
        # Mine bond index and class
        bond_index = int(match[1])
        bond_class = match[2]
        # If this class has a restart then skip it
        if bond_class in restarted_bond_classes: continue
        # If the index of the bond is 0 but we already processed bonds from this classes then this is a restart
        if bond_index == 0 and bond_class in processed_bond_classes:
            # If not, then add the class to the restarted classes
            restarted_bond_classes.add(bond_class)
            # And skip it
            continue
        # Add the class to the list
        processed_bond_classes.add(bond_class)
        # Get atom indices of bonded atoms
        bond = ( int(match[3]), int(match[4]) )
        bonds.append(bond)
    # Iterate bonds again to find water bonds with the SETTLE algorithm
    for line in lines:
        # Parse the line to get only atomic bonds
        match = search(GROMACS_TPR_SETTLE_REGEX, line)
        # If this is not a line with bonding data then skip it
        if not match: continue
        # If there is a 4rth atom we stop here
        if match[5]: raise ValueError('Found SETTLE bonds with 4 atoms. This is not yet supported')
        # Mine the settle
        settle = ( int(match[2]), int(match[3]), int(match[4]) )
        # Set actual bonds from the settles
        # The oxygen is always declared first
        # DANI: De esto no estoy 100% seguro, pero Miłosz me dijo que probablemente es así
        # DANI: Si no se cumple, el test de coherent bonds enconrará un hidrógeno con 2 enlaces
        bonds += [ (settle[0], settle[1]), (settle[0], settle[2]) ]
    # If we successfully got atom bonds then return them
    if len(bonds) > 0: return bonds
    # If there are no bonds at the end then something went wrong
    print(tpr_error_logs)
    raise RuntimeError(f'Bonds extraction from tpr file "{tpr_filepath}" has failed')

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
    run_gromacs(f'convert-tpr -s {input_structure_file.path} -o {output_structure_file.path} \
                -n {filter_index_filename}', user_input = filter_selection_name)

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
    run_gromacs(f'trjcat -f {new_file.path} {current_file} -o {current_file.path}')

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

# Filter atoms in a pdb file
# This method conserves maximum resolution and chains
def pdb_filter (
    input_pdb_filepath : str,
    output_pdb_filepath : str,
    index_filepath : str,
    filter_group_name : str
):
    # Filter the PDB
    run_gromacs(f'editconf -f {input_pdb_filepath} -o {output_pdb_filepath} \
                -n {index_filepath}', user_input = filter_group_name)

# Filter atoms in a xtc file
# Note that here we do not hide the stderr
# This is because it shows the progress
# Instead we color the output grey
def xtc_filter(
    structure_filepath : str,
    input_trajectory_filepath : str,
    output_trajectory_filepath : str,
    index_filepath : str,
    filter_group_name : str
):
    # Filter the trajectory
    run_gromacs(f'trjconv -s {structure_filepath} -f {input_trajectory_filepath} \
        -o {output_trajectory_filepath} -n {index_filepath}', user_input = filter_group_name)

# Filter atoms in both a pdb and a xtc file
def tpr_filter(
    input_tpr_filepath : str,
    output_tpr_filepath : str,
    index_filepath : str,
    filter_group_name : str
):
    # Filter the topology
    run_gromacs(f'convert-tpr -s {input_tpr_filepath} -o {output_tpr_filepath} \
                -n {index_filepath}', user_input = filter_group_name)

# Create a .ndx file from a complex mask
# e.g. no water and no ions -> !"Water"&!"Ion"
# This will return the group name to be further used
def make_index (input_structure_file : 'File', output_index_file : 'File', mask : str) -> str:
    # Run Gromacs
    run_gromacs(f'make_ndx -f {input_structure_file.path} -o {output_index_file.path}',
                user_input = f'{mask} \nq')
    # The group name is automatically assigned by gromacs
    # It equals the mask but removing andy " symbol
    group_name = mask.replace('"','')
    return group_name