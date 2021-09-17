import sys
import os
from subprocess import run, PIPE, Popen

import prody

# This tool allows you to set a gromax indexes file (.ndx) including an atom selection from prody selection syntax
# If the append argument is passed the selections will be appended to the existing output ndx file
# If the append argument is passed and the output ndx does not exist then it is created as the default Gromacs selections
# Arguments are as follows:
# 1 - Input pdb filename
#     e.g. 'example.pdb'
# 2 - Selection
#     e.g. 'chain E'
# 3 - Group name (Optional. 'custom' by default)
#     e.g. 'example selection'
# 4 - Output ndx filename (Optional. 'index.ndx' by default)
#     e.g. 'example.ndx'
# 5 - Append (Optional. False by default)
# e.g. python selection_to_ndx.py example.pdb 'chain A' 'custom selection' 'custom.ndx' True

# This arguments are passed when the script is called
# The input pdb filename
pdb_filename : str = sys.argv[1]
# The selection in prody syntax
selection : str = sys.argv[2]
# The atom group name for gromacs
if len(sys.argv) > 3:
    group_name : str = sys.argv[3]
else:
    group_name = 'custom'
# The output .ndx filename
if len(sys.argv) > 4:
    output_filename : str = sys.argv[4]
else:
    output_filename = 'index.ndx'
# The output .ndx filename
if len(sys.argv) > 5:
    append : bool = sys.argv[5]
else:
    append = False

# Check the file exists
if not os.path.exists(pdb_filename):
    raise SystemExit('ERROR: The file does not exist')

# Parse hte input topology
pdb = prody.parsePDB(pdb_filename)
# In case the append is true, create the default ndx file to append the new selection
if append and not os.path.exists(output_filename):
    p = Popen([
        "echo",
        'q',
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "make_ndx",
        "-f",
        pdb_filename,
        '-o',
        output_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()
# Add each selection to the ndx file
with open(output_filename, "a" if append else "w") as file:
    selected_pdb = pdb.select(selection)
    if not selected_pdb:
        raise SystemExit('WARNING: The fit center selection "' + selection + '" matches no atoms')
    print(selection + ' -> ' + str(selected_pdb.numAtoms()) + ' atoms selected')
    selected_atoms = list(selected_pdb.iterAtoms())
    selected_atom_indexes = [ atom.getIndex() for atom in selected_atoms ]
    # Add a header with the name for each group
    content = '[ ' + group_name + ' ]\n'
    count = 0
    for index in selected_atom_indexes:
        # Add a breakline each 15 indices
        count += 1
        if count == 15:
            content += '\n'
            count = 0
        # Add a space between indices
        # Atom indices go from 0 to n-1
        # Add +1 to the index since gromacs counts from 1 to n
        content += str(index + 1) + ' '
    file.write(content + '\n')