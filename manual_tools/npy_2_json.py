from sys import argv
from os.path import exists
from numpy import load
from json import dump

# Get input filename
input_filename = argv[1]
# Get output filename
output_filename = argv[2]

# Convert a '.npy' file to a '.json' file
def npy_2_json (input_filename : str, output_filename : str):
    if not exists(input_filename):
        raise SystemExit('Input file "' + input_filename + '" does not exist')
    content = list(load(input_filename))
    with open(output_filename, 'w') as file:
        dump(content, file, indent=4)

npy_2_json(input_filename, output_filename)
