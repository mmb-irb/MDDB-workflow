# This script is used to merge trajectory files and convert both topology and trajectory files formats
# Everything is carried by VMD software
# Humphrey, W., Dalke, A. and Schulten, K., "VMD - Visual Molecular Dynamics", J. Molec. Graphics, 1996, vol. 14, pp. 33-38. 
# http://www.ks.uiuc.edu/Research/vmd/

import os

from subprocess import run, PIPE, Popen

# Set he path to a script with all commands needed for vmd to parse the topology file
commands_filename = 'commands.vmd'

# input_topology_filename - The name string of the input topology file (path)
# Tested supported formats are .mae (WARNING: May not be fully compatible with the workflow), .psf and .pdb
# input_trajectory_filenames - The name string of the input trajectory filename(s) (path)
# Multiple files may be selected with bash syntax (e.g. *.dcd)
# Tested supported formats are .dcd and .xtc
# output_topology_filename - The name string of the output topology file (path)
# Tested supported format is .pdb
# output_trajectory_filename - The name string of the output trajectory file (path)
# Tested supported format is .trr
def vmd_processor (
    input_topology_filename : str,
    input_trajectory_filenames : list,
    output_topology_filename : str,
    output_trajectory_filename : str
    ) -> str:

    # Set an execptional output filename name for the vmd processed trajectory
    # It must be in 'trr' format since VMD does not support 'xtc' format
    processed_trajectory = 'md.trr'

    # Prepare a script for the VMD to automate the data parsing. This is Tcl lenguage
    # In addition, if chains are missing, this script asigns chains by fragment
    # Fragments are atom groups which are not connected by any bond
    with open(commands_filename, "w") as file:
        # Select all atoms
        file.write('set all [atomselect top "all"]\n')
        # Get all different chain names
        file.write('set chains_sample [lsort -unique [${all} get chain]]\n')
        # If there are only 'X' chains it means there are no chains at all
        # VMD asigns 'X' to missing chains by default
        file.write('if {[string compare $chains_sample X] == 0} {\n')
        # Set letters in alphabetic order
        file.write('	set letters "A B C D E F G H I J K L M N O P Q R S T U V W X Y Z"\n')
        # Get the number of fragments
        file.write('	set fragment_number [llength [lsort -unique -integer [${all} get fragment]]]\n')
        # For each fragment, set the chain of all atoms which belong to this fragment alphabetically
        # e.g. fragment 0 -> chain A, fragment 1 -> chain B, ...
        file.write('	for {set i 0} {$i <= $fragment_number} {incr i} {\n')
        file.write('		set fragment_atoms [atomselect top "fragment $i"]\n')
        file.write('		$fragment_atoms set chain [lindex $letters $i]\n')
        file.write('	}\n')
        file.write('}\n')
        # Write the current trajectory in 'trr' format
        if not os.path.exists(output_trajectory_filename) :
            file.write('animate write trr ' + processed_trajectory + ' waitfor all sel $all\n')
        # Write the current topology in 'pdb' format
        file.write('$all frame first\n')
        file.write('$all writepdb ' + output_topology_filename + '\n')
        file.write('exit\n')

    # Run VMD
    logs = run([
        "vmd",
        input_topology_filename,
        *input_trajectory_filenames,
        "-e",
        commands_filename,
        "-dispdev",
        "none"
    ], stdout=PIPE).stdout.decode()

    # Check the output files have been created
    # Otherwise print logs and throw an error
    if not ( os.path.exists(output_topology_filename) and ( os.path.exists(processed_trajectory) or os.path.exists(output_trajectory_filename) ) ) :
        print(logs)
        raise SystemExit("ERROR: Something was wrong with VMD")

    # Convert the output processed trajectory from 'trr' to 'xtc' format
    if not os.path.exists(output_trajectory_filename) and os.path.exists(processed_trajectory):
        p = Popen([
            "echo",
            "System",
        ], stdout=PIPE)
        logs = run([
            "gmx",
            "trjconv",
            "-s",
            output_topology_filename,
            "-f",
            processed_trajectory,
            '-o',
            output_trajectory_filename,
            '-quiet'
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()

    # Remove the trajectory in 'trr' format to free up space
    if os.path.exists(processed_trajectory):
        logs = run([
            "rm",
            processed_trajectory,
        ], stdout=PIPE).stdout.decode()

    # Return VMD logs
    return logs

# An abbreviated version of the function above
# Add chains to a given topology which is missing chains
def vmd_chainer (
    input_topology_filename : str,
    output_topology_filename : str,
    ) -> str:

    # Prepare a script for the VMD to automate the data parsing. This is Tcl lenguage
    # In addition, if chains are missing, this script asigns chains by fragment
    # Fragments are atom groups which are not connected by any bond
    with open(commands_filename, "w") as file:
        # Select all atoms
        file.write('set all [atomselect top "all"]\n')
        # Get all different chain names
        file.write('set chains_sample [lsort -unique [${all} get chain]]\n')
        # If there are only 'X' chains it means there are no chains at all
        # VMD asigns 'X' to missing chains by default
        file.write('if {[string compare $chains_sample X] == 0} {\n')
        # Set letters in alphabetic order
        file.write('	set letters "A B C D E F G H I J K L M N O P Q R S T U V W X Y Z"\n')
        # Get the number of fragments
        file.write('	set fragment_number [llength [lsort -unique -integer [${all} get fragment]]]\n')
        # For each fragment, set the chain of all atoms which belong to this fragment alphabetically
        # e.g. fragment 0 -> chain A, fragment 1 -> chain B, ...
        file.write('	for {set i 0} {$i <= $fragment_number} {incr i} {\n')
        file.write('		set fragment_atoms [atomselect top "fragment $i"]\n')
        file.write('		$fragment_atoms set chain [lindex $letters $i]\n')
        file.write('	}\n')
        file.write('}\n')
        # Write the current topology in 'pdb' format
        file.write('$all frame first\n')
        file.write('$all writepdb ' + output_topology_filename + '\n')
        file.write('exit\n')

    # Run VMD
    logs = run([
        "vmd",
        input_topology_filename,
        "-e",
        commands_filename,
        "-dispdev",
        "none"
    ], stdout=PIPE).stdout.decode()

    # Check the output files have been created
    # Otherwise print logs and throw an error
    if not os.path.exists(output_topology_filename) :
        print(logs)
        raise SystemExit("ERROR: Something was wrong with VMD")

    # Remove the commands file
    run([
        "rm",
        commands_filename,
    ], stdout=PIPE).stdout.decode()

    # Return VMD logs
    return logs

# Given a psf topology and a single frame, generate a pdb file
def psf_to_pdb (
    input_topology_filename : str,
    input_frame_filename : str,
    output_filename : str):

    # Prepare a script for the VMD to automate the data parsing. This is Tcl lenguage
    # In addition, if chains are missing, this script asigns chains by fragment
    # Fragments are atom groups which are not connected by any bond
    with open(commands_filename, "w") as file:
        # Select all atoms
        file.write('set all [atomselect top "all"]\n')
        # Write the current topology in 'pdb' format
        file.write('$all frame first\n')
        file.write('$all writepdb ' + output_filename + '\n')
        file.write('exit\n')

    # Run VMD
    logs = run([
        "vmd",
        input_topology_filename,
        input_frame_filename,
        "-e",
        commands_filename,
        "-dispdev",
        "none"
    ], stdout=PIPE).stdout.decode()