from sys import argv

target_filepath = argv[1]

# Mine cmip output pdb files
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
print(f' Total energies: vmd {total_vdw}, es {total_es}, both {total_both}')
