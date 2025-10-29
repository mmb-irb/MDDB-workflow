# Functions powered by VMD
# Humphrey, W., Dalke, A. and Schulten, K., "VMD - Visual Molecular Dynamics", J. Molec. Graphics, 1996, vol. 14, pp. 33-38. 
# http://www.ks.uiuc.edu/Research/vmd/

import os
from os.path import exists

from subprocess import run, PIPE, STDOUT

from mddb_workflow.utils.file import File
from mddb_workflow.utils.type_hints import *
from mddb_workflow.utils.auxiliar import warn

# Set characters to be escaped since they have a meaning in TCL
TCL_RESERVED_CHARACTERS = ['"','[',']']

def escape_tcl_selection (selection : str) -> str:
    """ Given a VMD atom selection string, escape TCL meaningful characters and return the escaped string. """
    escaped_selection = selection
    for character in TCL_RESERVED_CHARACTERS:
        escaped_selection = escaped_selection.replace(character, '\\' + character)
    return escaped_selection

# Set the script filename with all commands to be passed to vmd
commands_filename = '.commands.vmd'

# List all the vmd supported trajectory formats
vmd_supported_structure_formats = {'pdb', 'prmtop', 'psf', 'parm', 'gro'} # DANI: Esto lo he hecho rápido, hay muchas más
vmd_supported_trajectory_formats = {'mdcrd', 'crd', 'dcd', 'xtc', 'trr', 'nc', 'netcdf', 'cdf', 'pdb', 'rst7'}

# Set a vmd format translator
vmd_format = {
    'pdb': 'pdb',
    'prmtop': 'prmtop',
    'psf': 'psf',
    'gro': 'gro',
    'crd': 'crd',
    'mdcrd': 'crd',
    'dcd': 'dcd',
    'xtc': 'xtc',
    'trr': 'trr',
    'netcdf': 'netcdf',
    'nc': 'netcdf',
    'cdf': 'netcdf',
    'rst7': 'rst7'
}

# Given a vmd supported topology with no coordinates and a single frame file, generate a pdb file
def vmd_to_pdb (
    input_structure_filename : str,
    input_trajectory_filename : str,
    output_structure_filename : str):

    # Get the trajectory format according to the vmd dictionary
    input_trajectory_file = File(input_trajectory_filename)
    trajectory_format = input_trajectory_file.format
    vmd_trajectory_format = vmd_format[trajectory_format]

    # Prepare a script for VMD to run. This is Tcl language
    with open(commands_filename, "w") as file:
        # Load only the first frame of the trajectory
        file.write(f'animate read {vmd_trajectory_format} {input_trajectory_filename} end 1\n')
        # Select all atoms
        file.write('set all [atomselect top "all"]\n')
        # Write the current topology in 'pdb' format
        file.write('$all frame first\n')
        file.write(f'$all writepdb {output_structure_filename}\n')
        file.write('exit\n')

    # Run VMD
    vmd_process = run([
        "vmd",
        input_structure_filename,
        "-e",
        commands_filename,
        "-dispdev",
        "none"
    ], stdout=PIPE, stderr=PIPE)
    output_logs = vmd_process.stdout.decode()

    if not exists(output_structure_filename):
        print(output_logs)
        error_logs = vmd_process.stderr.decode()
        print(error_logs)
        raise SystemExit('Something went wrong with VMD')

    os.remove(commands_filename)
# Set function supported formats
vmd_to_pdb.format_sets = [
    {
        'inputs': {
            'input_structure_filename': {'psf', 'parm', 'prmtop'},
            'input_trajectory_filename': vmd_supported_trajectory_formats
        },
        'outputs': {
            'output_structure_filename': {'pdb'}
        }
    }
]

def chainer (
    input_pdb_filename : str,
    atom_selection : Optional[str] = None,
    chain_letter : Optional[str] = None,
    output_pdb_filename : Optional[str] = None
):
    """ This tool allows you to set the chain of all atoms in a selection.
    This is powered by VMD and thus the selection lenguage must be the VMD's.
    Arguments are as follows:

    1. Input pdb filename
    2. Atom selection (All atoms by defualt)
    3. Chain letter (May be the flag 'fragment', which is the default indeed)
    4. Output pdb filename (Input filename by default)

    WARNING: When no selection is passed, if only a part of a "fragment" is missing the chain then the whole fragment will be affected
    WARNING: VMD only handles fragments if there are less fragments than letters in the alphabet
    DEPRECATED: Use the structures chainer instead. """

    # If no atom selection is provided then all atoms are selected
    if not atom_selection:
        atom_selection = 'all'

    # Escape TCL meaningful characters
    escaped_atom_selection = escape_tcl_selection(atom_selection)

    # If no chain letter is provided then the flag 'fragment' is used
    if not chain_letter:
        chain_letter = 'fragment'

    # If no output filename is provided then use input filename as output filename
    if not output_pdb_filename:
        output_pdb_filename = input_pdb_filename

    # Check the file exists
    if not exists(input_pdb_filename):
        raise SystemExit('ERROR: The file does not exist')
       
    with open(commands_filename, "w") as file:
        # Select the specified atoms and set the specified chain
        file.write(f'set atoms [atomselect top "{escaped_atom_selection}"]\n')
        # In case chain letter is not a letter but the 'fragment' flag, asign chains by fragment
        # Fragments are atom groups which are not connected by any bond
        if chain_letter == 'fragment':
            # Get all different chain names
            file.write('set chains_sample [lsort -unique [${atoms} get chain]]\n')
            # Set letters in alphabetic order
            file.write('set letters "A B C D E F G H I J K L M N O P Q R S T U V W X Y Z"\n')
            # Get the number of fragments
            file.write('set fragment_number [llength [lsort -unique -integer [${atoms} get fragment]]]\n')
            # For each fragment, set the chain of all atoms which belong to this fragment alphabetically
            # e.g. fragment 0 -> chain A, fragment 1 -> chain B, ...
            file.write('for {set i 0} {$i <= $fragment_number} {incr i} {\n')
            file.write('	set fragment_atoms [atomselect top "fragment $i"]\n')
            file.write('	$fragment_atoms set chain [lindex $letters $i]\n')
            file.write('}\n')
            # Otherwise, set the specified chain
        else:
            file.write(f'$atoms set chain {chain_letter}\n')
        # Write the current topology in 'pdb' format
        file.write('set all [atomselect top "all"]\n')
        file.write('$all frame first\n')
        file.write(f'$all writepdb {output_pdb_filename}\n')
        file.write('exit\n')

    # Run VMD
    logs = run([
        "vmd",
        input_pdb_filename,
        "-e",
        commands_filename,
        "-dispdev",
        "none"
    ], stdout=PIPE, stderr=PIPE).stdout.decode()

    # If the expected output file was not generated then stop here and warn the user
    if not exists(output_pdb_filename):
        print(logs)
        raise SystemExit('Something went wrong with VMD')

    # Remove the vmd commands file
    os.remove(commands_filename)

def merge_and_convert_trajectories (
    input_structure_filename : Optional[str],
    input_trajectory_filenames : list[str],
    output_trajectory_filename : str
    ):
    """ Get vmd supported trajectories merged and converted to a different format.
    WARNING: Note that this process is not memory efficient so beware the size of trajectories to be converted.
    WARNING: The input structure filename may be None. """

    warn('You are using a not memory efficient tool. If the trajectory is too big your system may not hold it.')
    print('Note that we cannot display the progress of the conversion since we are using VMD')

    # Get the format to export coordinates
    output_trajectory_file = File(output_trajectory_filename)
    output_trajectory_format = output_trajectory_file.format

    # Although 'crd' and 'mdcrd' may be the same, VMD only recognizes 'crd' as exporting coordinates file type
    if output_trajectory_format == 'mdcrd':
        output_trajectory_format = 'crd'

    # Prepare a script for the VMD to automate the data parsing. This is Tcl lenguage
    # In addition, if chains are missing, this script asigns chains by fragment
    # Fragments are atom groups which are not connected by any bond
    with open(commands_filename, "w") as file:
        # Select all atoms
        file.write('set all [atomselect top "all"]\n')
        # Write the current trajectory in the specified format format
        file.write(f'animate write {output_trajectory_format} {output_trajectory_filename} waitfor all sel $all\n')
        file.write('exit\n')

    inputs = [ input_structure_filename, *input_trajectory_filenames ] if input_structure_filename else input_trajectory_filenames

    # Run VMD
    logs = run([
        "vmd",
        *inputs,
        "-e",
        commands_filename,
        "-dispdev",
        "none"
    ], stdout=PIPE, stderr=STDOUT).stdout.decode() # Redirect errors to the output in order to dont show them in console

    # If the expected output file was not generated then stop here and warn the user
    if not exists(output_trajectory_filename):
        print(logs)
        raise SystemExit('Something went wrong with VMD')

    # Remove the vmd commands file
    os.remove(commands_filename)

# Set function supported formats
merge_and_convert_trajectories.format_sets = [
    {
        'inputs': {
            'input_structure_filename': None,
            'input_trajectory_filenames': {'dcd', 'xtc', 'trr', 'nc'}
        },
        'outputs': {
            # NEVER FORGET: VMD cannot write to all formats it supports to read
            'output_trajectory_filename': {'mdcrd', 'crd', 'dcd', 'trr'}
        }
    },
    {
        'inputs': {
            'input_structure_filename': vmd_supported_structure_formats,
            'input_trajectory_filenames': {'mdcrd', 'crd', 'rst7'}
        },
        'outputs': {
            # NEVER FORGET: VMD cannot write to all formats it supports to read
            'output_trajectory_filename': {'mdcrd', 'crd', 'dcd', 'trr'}
        }
    }
]

def get_vmd_selection_atom_indices (input_structure_filename : str, selection : str) -> list[int]:
    """ Given an atom selection in vmd syntax, return the list of atom indices it corresponds to. """

    # Escape TCL meaningful characters
    escaped_selection = escape_tcl_selection(selection)
    
    # Prepare a script for VMD to run. This is Tcl language
    # The output of the script will be written to a txt file
    atom_indices_filename = '.vmd_output.txt'
    with open(commands_filename, "w") as file:
        # Select the specified atoms
        file.write(f'set selection [atomselect top "{escaped_selection}"]\n')
        # Save atom indices from the selection
        file.write('set indices [$selection list]\n')
        # Write atom indices to a file
        file.write(f'set indices_file [open {atom_indices_filename} w]\n')
        file.write('puts $indices_file $indices\n')
        file.write('exit\n')

    # Run VMD
    logs = run([
        "vmd",
        input_structure_filename,
        "-e",
        commands_filename,
        "-dispdev",
        "none"
    ], stdout=PIPE, stderr=PIPE).stdout.decode()

    # If the expected output file was not generated then stop here and warn the user
    if not exists(atom_indices_filename):
        print(logs)
        raise SystemExit('Something went wrong with VMD')

    # Read the VMD output
    with open(atom_indices_filename, 'r') as file:
        raw_atom_indices = file.read()

    # Parse the atom indices string to an array of integers
    atom_indices = [ int(i) for i in raw_atom_indices.split() ]
    
    # Remove trahs files
    trash_files = [ commands_filename, atom_indices_filename ]
    for trash_file in trash_files:
        os.remove(trash_file)

    return atom_indices

def get_covalent_bonds (structure_filename : str, selection : Optional['Selection'] = None) -> list[ list[int] ]:
    """ Set a function to retrieve all covalent (strong) bonds in a structure using VMD.
    You may provide an atom selection as well. """

    # Parse the selection to vmd
    vmd_selection = 'all'
    if selection:
        vmd_selection = selection.to_vmd()

    # Prepare a script for the VMD to automate the commands. This is Tcl lenguage
    output_bonds_file = '.bonds.txt'
    with open(commands_filename, "w") as file:
        # Select atoms
        file.write(f'set atoms [atomselect top "{vmd_selection}"]\n')
        # Save covalent bonds
        file.write('set bonds [$atoms getbonds]\n')
        # Write those bonds to a file
        file.write(f'set bondsfile [open {output_bonds_file} w]\n')
        file.write('puts $bondsfile $bonds\n')
        file.write('exit\n')
        
    # Run VMD
    logs = run([
        "vmd",
        structure_filename,
        "-e",
        commands_filename,
        "-dispdev",
        "none"
    ], stdout=PIPE, stderr=PIPE).stdout.decode()

    # If the output file is missing at this point then it means something went wrong
    if not exists(output_bonds_file):
        print(logs)
        raise SystemExit('Something went wrong with VMD')
    
    # Read the VMD output
    with open(output_bonds_file, 'r') as file:
        raw_bonds = file.read()

    # Remove vmd files since they are no longer usefull
    for f in [ commands_filename, output_bonds_file ]:
        os.remove(f)

    # Sometimes there is a breakline at the end of the raw bonds string and it must be removed
    # Add a space at the end of the string to make the parser get the last character
    raw_bonds = raw_bonds.replace('\n', '') + ' '
    
    # Parse the raw bonds string to a list of atom bonds (i.e. a list of lists of integers)
    # Raw bonds format is (for each atom in the selection):
    # '{index1, index2, index3 ...}' with the index of each connected atom
    # 'index' if there is only one connected atom
    # '{}' if there are no connected atoms
    bonds_per_atom = []
    last_atom_index = ''
    last_atom_bonds = []
    in_brackets = False
    for character in raw_bonds:
        if character == ' ':
            if len(last_atom_index) > 0:
                if in_brackets:
                    last_atom_bonds.append(int(last_atom_index))
                else:
                    bonds_per_atom.append([int(last_atom_index)])
                last_atom_index = ''
            continue
        if character == '{':
            in_brackets = True
            continue
        if character == '}':
            if last_atom_index == '':
                bonds_per_atom.append([])
                in_brackets = False
                continue
            last_atom_bonds.append(int(last_atom_index))
            last_atom_index = ''
            bonds_per_atom.append(last_atom_bonds)
            last_atom_bonds = []
            in_brackets = False
            continue
        last_atom_index += character
                
    return bonds_per_atom

def get_covalent_bonds_between (
    structure_filename : str,
    # Selections may be either selection instances or selection strings already in VMD format
    selection_1 : Union['Selection', str],
    selection_2 : Union['Selection', str]
) -> list[ list[int] ]:
    """ Set a function to retrieve covalent (strong) bonds between 2 atom selections. """
    
    # Parse selections (if not parsed yet)
    parsed_selection_1 = selection_1 if type(selection_1) == str else selection_1.to_vmd()
    parsed_selection_2 = selection_2 if type(selection_2) == str else selection_2.to_vmd()

    # Prepare a script for the VMD to automate the commands. This is Tcl lenguage
    output_index_1_file = '.index1.txt'
    output_index_2_file = '.index2.txt'
    output_bonds_file = '.bonds.ext'
    with open(commands_filename, "w") as file:
        # Select the specified atoms in selection 1
        file.write(f'set sel1 [atomselect top "{parsed_selection_1}"]\n')
        # Save all atom index in the selection
        file.write('set index1 [$sel1 list]\n')
        # Write those index to a file
        file.write(f'set indexfile1 [open {output_index_1_file} w]\n')
        file.write('puts $indexfile1 $index1\n')
        # Save all covalent bonds in the selection
        file.write('set bonds [$sel1 getbonds]\n')
        # Write those bonds to a file
        file.write(f'set bondsfile [open {output_bonds_file} w]\n')
        file.write('puts $bondsfile $bonds\n')
        # Select the specified atoms in selection 2
        file.write(f'set sel2 [atomselect top "{parsed_selection_2}"]\n')
        # Save all atom index in the selection
        file.write('set index2 [$sel2 list]\n')
        # Write those index to a file
        file.write(f'set indexfile2 [open {output_index_2_file} w]\n')
        file.write('puts $indexfile2 $index2\n')
        file.write('exit\n')
        
    # Run VMD
    logs = run([
        "vmd",
        structure_filename,
        "-e",
        commands_filename,
        "-dispdev",
        "none"
    ], stdout=PIPE, stderr=PIPE).stdout.decode()

    # If the expected output file was not generated then stop here and warn the user
    if not exists(output_index_1_file) or not exists(output_bonds_file) or not exists(output_index_2_file):
        print(logs)
        raise SystemExit('Something went wrong with VMD')
    
    # Read the VMD output
    with open(output_index_1_file, 'r') as file:
        raw_index_1 = file.read()
    with open(output_bonds_file, 'r') as file:
        raw_bonds = file.read()
    with open(output_index_2_file, 'r') as file:
        raw_index_2 = file.read()

    # Remove vmd files since they are no longer usefull
    for f in [ commands_filename, output_index_1_file, output_index_2_file, output_bonds_file ]:
        os.remove(f)

    # Sometimes there is a breakline at the end of the raw bonds string and it must be removed
    # Add a space at the end of the string to make the parser get the last character
    raw_bonds = raw_bonds.replace('\n', '') + ' '
    
    # Raw indexes is a string with all indexes separated by spaces
    index_1 = [ int(i) for i in raw_index_1.split() ]
    index_2 = set([ int(i) for i in raw_index_2.split() ])
    
    # Parse the raw bonds string to a list of atom bonds (i.e. a list of lists of integers)
    # Raw bonds format is (for each atom in the selection):
    # '{index1, index2, index3 ...}' with the index of each connected atom
    # 'index' if there is only one connected atom
    # '{}' if there are no connected atoms
    bonds_per_atom = []
    last_atom_index = ''
    last_atom_bonds = []
    in_brackets = False
    for character in raw_bonds:
        if character == ' ':
            if len(last_atom_index) > 0:
                if in_brackets:
                    last_atom_bonds.append(int(last_atom_index))
                else:
                    bonds_per_atom.append([int(last_atom_index)])
                last_atom_index = ''
            continue
        if character == '{':
            in_brackets = True
            continue
        if character == '}':
            if last_atom_index == '':
                bonds_per_atom.append([])
                in_brackets = False
                continue
            last_atom_bonds.append(int(last_atom_index))
            last_atom_index = ''
            bonds_per_atom.append(last_atom_bonds)
            last_atom_bonds = []
            in_brackets = False
            continue
        last_atom_index += character
        
    # At this point indexes and bonds from the first selection should match in number
    if len(index_1) != len(bonds_per_atom):
        raise ValueError(f'Indexes ({len(index_1)}) and atom bonds ({len(bonds_per_atom)}) do not match in number')
        
    # Now get all covalent bonds which include an index from the atom selection 2
    crossed_bonds = []
    for i, index in enumerate(index_1):
        bonds = bonds_per_atom[i]
        for bond in bonds:
            if bond in index_2:
                # DANI: A tuple may make more sense than a list to define a bond
                # DANI: However tuples are converted to lists when JSON serialized
                # DANI: And the interactions are saved to json, so a list keeps things more coherent
                crossed_bond = [index, bond]
                crossed_bonds.append(crossed_bond)
                
    return crossed_bonds

def get_interface_atom_indices (
    input_structure_filepath : str,
    input_trajectory_filepath : str,
    selection_1 : str,
    selection_2 : str,
    distance_cutoff : float,
) -> list[int]:
    """ Given two atom selections, find interface atoms and return their indices
    Interface atoms are those atoms closer than the cutoff in at least 1 frame along a trajectory
    Return also atom indices for the whole selections. """

    # Set the interface selections
    interface_selection_1 = (f'({selection_1}) and within {distance_cutoff} of ({selection_2})')
    interface_selection_2 = (f'({selection_2 }) and within {distance_cutoff} of ({selection_1})')
    
    # Set the output txt files for vmd to write the atom indices
    # Note that these output files are deleted at the end of this function
    selection_1_filename = '.selection_1.txt'
    selection_2_filename = '.selection_2.txt'
    interface_selection_1_filename = '.interface_selection_1.txt'
    interface_selection_2_filename = '.interface_selection_2.txt'
    interacting_frames_filename = '.iframes.txt'
    total_frames_filename = '.nframes.txt'

    # Prepare a script for VMD to run. This is Tcl language
    commands_filename = '.commands.vmd'
    with open(commands_filename, "w") as file:
        # -------------------------------------------
        # First get the whole selection atom indices
        # -------------------------------------------
        # Select the specified atoms
        file.write(f'set selection [atomselect top "{selection_1}"]\n')
        # Save atom indices from the selection
        file.write('set indices [$selection list]\n')
        # Write atom indices to a file
        file.write(f'set indices_file [open {selection_1_filename} w]\n')
        file.write('puts $indices_file $indices\n')
        # Select the specified atoms
        file.write(f'set selection [atomselect top "{selection_2}"]\n')
        # Save atom indices from the selection
        file.write('set indices [$selection list]\n')
        # Write atom indices to a file
        file.write(f'set indices_file [open {selection_2_filename} w]\n')
        file.write('puts $indices_file $indices\n')
        # -------------------------------------------
        # Now get the interface selection atom indices
        # Also count the number of frames where there is at least one interacting residue
        # -------------------------------------------
        # Capture indices for each frame in the trajectory
        file.write('set accumulated_interface1_atom_indices []\n')
        file.write('set accumulated_interface2_atom_indices []\n')
        file.write(f'set interface1 [atomselect top "{interface_selection_1}"]\n')
        file.write(f'set interface2 [atomselect top "{interface_selection_2}"]\n')
        # Capture the number of frames where the interaction happens
        file.write('set iframes 0\n')
        # Get the number of frames in the trajectory
        file.write('set nframes [molinfo top get numframes]\n')
        # Iterate over each frame
        # Note that we skip the first frame (i = 1, not i = 0) since it belongs to the structure
        file.write('for { set i 1 } { $i < $nframes } { incr i } {\n')
        # Update the selection in the current frame
        file.write('    $interface1 frame $i\n')
        file.write('    $interface1 update\n')
        # Add its atom indices to the acumulated atom indices
        file.write('    set interface1_atom_indices [$interface1 list]\n')
        file.write('    set accumulated_interface1_atom_indices [concat $accumulated_interface1_atom_indices $interface1_atom_indices ]\n')
        # Repeat with the selection 2
        file.write('    $interface2 frame $i\n')
        file.write('    $interface2 update\n')
        file.write('    set interface2_atom_indices [$interface2 list]\n')
        file.write('    set accumulated_interface2_atom_indices [concat $accumulated_interface2_atom_indices $interface2_atom_indices ]\n')
        # If there was at least one atom in one of the interactions then add one to the interaction frame count
        # Note that checking both interactions would be redundant so one is enough
        file.write('    if { [llength $interface1_atom_indices] > 0 } {\n')
        file.write('        incr iframes\n')
        file.write('    }\n')
        file.write('}\n')
        # Write the number of interacting frames and total frames to files
        file.write(f'set iframes_file [open {interacting_frames_filename} w]\n')
        file.write('puts $iframes_file $iframes\n')
        file.write(f'set nframes_file [open {total_frames_filename} w]\n')
        # Note that we must substract 1 from the total frames count since the first frame was skipped
        file.write('set cframes [expr $nframes - 1]\n')
        file.write('puts $nframes_file $cframes\n')
        # Remove duplicated indices
        # Use the -integer argument to do the correct sorting
        # Otherwise numbers are sorted as strings so you have 1, 10, 100, 2, etc.
        file.write('set unique_interface1_atom_indices [lsort -integer -unique $accumulated_interface1_atom_indices]\n')
        file.write('set unique_interface2_atom_indices [lsort -integer -unique $accumulated_interface2_atom_indices]\n')
        # Write indices to files
        file.write(f'set indices_file [open {interface_selection_1_filename} w]\n')
        file.write('puts $indices_file $unique_interface1_atom_indices\n')
        file.write(f'set indices_file [open {interface_selection_2_filename} w]\n')
        file.write('puts $indices_file $unique_interface2_atom_indices\n')
        file.write('exit\n')

    # Run VMD
    logs = run([
        "vmd",
        input_structure_filepath,
        input_trajectory_filepath,
        "-e",
        commands_filename,
        "-dispdev",
        "none"
    ], stdout=PIPE, stderr=PIPE).stdout.decode()

    # If any of the output files do not exist at this point then it means something went wrong with vmd
    expected_output_files = [
        selection_1_filename,
        selection_2_filename,
        interface_selection_1_filename,
        interface_selection_2_filename,
        interacting_frames_filename,
        total_frames_filename
    ]
    for output_file in expected_output_files:
        if not os.path.exists(output_file):
            print(logs)
            raise SystemExit('Something went wrong with VMD')
    
    # Set a function to read the VMD output and parse the atom indices string to an array of integers
    def process_vmd_output (output_filename : str) -> list[int]:
        with open(output_filename, 'r') as file:
            raw_atom_indices = file.read()
        return [ int(i) for i in raw_atom_indices.split() ]

    # Read the VMD output
    selection_1_atom_indices = process_vmd_output(selection_1_filename)
    selection_2_atom_indices = process_vmd_output(selection_2_filename)
    selection_1_interface_atom_indices = process_vmd_output(interface_selection_1_filename)
    selection_2_interface_atom_indices = process_vmd_output(interface_selection_2_filename)
    interacting_frames = process_vmd_output(interacting_frames_filename)[0]
    total_frames = process_vmd_output(total_frames_filename)[0]
    
    # Remove trash files
    trash_files = [ commands_filename ] + expected_output_files
    for trash_file in trash_files:
        os.remove(trash_file)

    # Return the results
    return {
        'selection_1_atom_indices': selection_1_atom_indices,
        'selection_2_atom_indices': selection_2_atom_indices,
        'selection_1_interface_atom_indices': selection_1_interface_atom_indices,
        'selection_2_interface_atom_indices': selection_2_interface_atom_indices,
        'interacting_frames': interacting_frames,
        'total_frames': total_frames
    }