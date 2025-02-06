from sys import argv

# Make sure the target filename was passed and there are no additional arguments
if len(argv) < 2:
    raise SystemExit('This script expects arguments: the target filepaths\n'
        'e.g. python get_energies_sum.py target_1.cmip.pdb target_2.cmip.pdb')

# Get target filepaths from user arguments
target_filepaths = argv[1:]

# Set a function to round numbers to cents
def rount_to_cents (value : float) -> float:
    return round(value * 100) / 100

# Mine cmip output pdb files, sum their energies and log the result
def mine_cmip_output (target_filepath : str):
    atom_energies = []
    with open(target_filepath, 'r') as file:
        lines = list(file)
        # Mine energies line by line (i.e. atom by atom)
        for line in lines:
            # Numbers may become larger than expected when they are close to 0 since they are expressed in E notation
            # e.g. 0.1075E-04
            # For this reason we must mine these last values relying on whitespaces between them
            line_splits = line.split()
            vdw = float(line_splits[-3])
            es = float(line_splits[-2])
            both = float(line_splits[-1])
            # Values greater than 100 are represented as 0
            # This step is performed to filter 'infinity' values
            energies = (vdw, es, both)
            if both > 100:
                energies = (0, 0, 0)
            # Add current atom energy values to the atom energies list
            atom_energies.append(energies)
    # Sum energies
    total_vdw = total_es = total_both = 0
    for energies in atom_energies:
        total_vdw += energies[0]
        total_es += energies[1]
        total_both += energies[2]
    # Round to cents for a prettier display
    total_vdw = rount_to_cents(total_vdw)
    total_es = rount_to_cents(total_es)
    total_both = rount_to_cents(total_both)
    print(f' Total energies in {target_filepath}: vmd {total_vdw}, es {total_es}, both {total_both}')

# Mine and log energies for every input file
for target_filepath in target_filepaths:
    try:
        mine_cmip_output(target_filepath)
    except:
        raise SystemExit(f'Failed to mine {target_filepath}. Are you sure this is a CMIP output file?')
