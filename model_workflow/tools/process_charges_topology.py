import os

def is_raw (filename : str) -> bool:
    return filename == 'charges.txt'

def is_prmtop (filename : str) -> bool:
    return filename[-7:] == '.prmtop'

def is_top (filename : str) -> bool:
    return filename[-4:] == '.top'

def is_psf (filename : str) -> bool:
    return filename[-4:] == '.psf'

# This script is used to get standarized charges topology
# In case charges are in a 'charges.txt' this process does nothing
# In case charges are in .prmtop, .top or .psf file we rename the file
def process_charges_topology (input_charges_filename : str) -> str:

    if not input_charges_filename:
        return None

    # If it is a 'charges.txt' stop here
    if is_raw(input_charges_filename):
        return input_charges_filename

    # In other cases rename the topology and return the standard filename
    if is_prmtop(input_charges_filename):
        standard_filename = 'topology.prmtop'
    elif is_top(input_charges_filename):
        standard_filename = 'topology.top'
    elif is_psf(input_charges_filename):
        standard_filename = 'topology.psf'
    
    if not os.path.exists(standard_filename):
        os.rename(input_charges_filename, standard_filename)

    return standard_filename

# Find out if there is any of the standard topology filenames in the current directory
# In that case, return the first standard filename found
# Raw energies have priority
def find_energies_filename():

    # Set the standard filenames
    standard_filenames = [
        'charges.txt',
        'topology.prmtop',
        'topology.top',
        'topology.psf',
    ]

    # return the first existing standard filename
    for standard_filename in standard_filenames:
        if os.path.exists(standard_filename):
            return standard_filename

    return None