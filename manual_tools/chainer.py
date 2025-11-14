import sys
import os
from subprocess import run, PIPE

from mddb_workflow.utils.auxiliar import InputError

# DANI: This is deprecared

# This tool allows you to set the chain of all atoms in a selection
# This is powered by VMD and thus the selection lenguage must be the VMD's
# Arguments are as follows:
# 1 - Input pdb filename
# 2 - Atom selection
# 3 - Chain letter
# 4 - Output pdb filename (Optional. Input filename by default)
# e.g. python chainer.py example.pdb 'resname POPS' 'M'

# This arguments are passed when the script is called
# The input pdb filename
input_pdb_filename = sys.argv[1]
# The atom selection
atom_selection = sys.argv[2]
# The new chain letter
chain_letter = sys.argv[3]
# The output pdb filename
if len(sys.argv) > 4:
    output_pdb_filename = sys.argv[4]
else:
    output_pdb_filename = input_pdb_filename

# Check the file exists
if not os.path.exists(input_pdb_filename):
    raise InputError(f'The file {input_pdb_filename} does not exist')

# Set he path to a script with all commands needed for vmd to parse the topology file
commands_filename = '.commands.vmd'

# Prepare a script for the VMD to automate the data parsing. This is Tcl lenguage
# In addition, if chains are missing, this script asigns chains by fragment
# Fragments are atom groups which are not connected by any bond
with open(commands_filename, "w") as file:
    # Select the specified atoms and set the specified chain
    file.write('set atoms [atomselect top "' + atom_selection + '"]\n')
    file.write('$atoms set chain ' + chain_letter + '\n')
    # Write the current topology in 'pdb' format
    file.write('set all [atomselect top "all"]\n')
    file.write('$all frame first\n')
    file.write('$all writepdb ' + output_pdb_filename + '\n')
    file.write('exit\n')

# Run VMD
logs = run([
    "vmd",
    input_pdb_filename,
    "-e",
    commands_filename,
    "-dispdev",
    "none"
], stdout=PIPE).stdout.decode()

# Remove the vmd commands file
os.remove(commands_filename)

