import pytraj as pt

from pathlib import Path

# Given a topology which includes charges
# Extract those charges and save them in a list to be returned
# Supported formats (tested): prmtop, top, psf (standard psf, not from DESRES)
# Supported formats (not tested): mol2, cif, sdf
# Non supported formats: mae, tpr, pdb (has no charges)
def get_topology_charges (input_topology_filename : str):
    topology = pt.load_topology(filename=input_topology_filename)
    # WARNING: We must convert this numpy ndarray to a normal list
    # Otherwise the search by index is extremly ineficient
    topology_charges = list(topology.charge)
    return topology_charges

# Given a raw file with listed charges
# Extract those charges and save them in a list to be returned
def get_raw_charges (input_charges_filename : str):
    charges = []
    with open(input_charges_filename, 'r') as file:
        lines = file.readlines()
        for line in lines:
            charges.append(float(line))
    return charges

# Set the path to the source res.lib file
resources = str(Path(__file__).parent.parent / "utils" / "resources")
reslib_source = resources + '/res.lib'

# Given a topology which does NOT include charges
# Find which atoms are already describen in the source res.lib and write them in a new local res.lib
def get_reslib_charges (
    input_topology_filename : str,
    output_reslib_filename : str,
    source_reslib_filename : str = reslib_source):
    
    topology = pt.load_topology(filename=input_topology_filename)
    # Get the source reslib atoms
    source_atoms = get_reslib_atoms(source_reslib_filename)
    # Trak if some atom is missing
    missing_atoms = False
    # Save the name of already found unknown residues so we can skip them when repeated
    already_parsed_atoms = []
    with open(output_reslib_filename, 'w') as file:
        atoms = list(topology.atoms)
        for a, atom in enumerate(atoms):
            residue = atom.resname
            name = atom.name
            # Skip this atom if we already found it
            if next((at for at in already_parsed_atoms if at['residue'] == residue and at['name'] == name), None):
                continue
            # If the atom is found in the source res lib copy the current reslib line in our new reslib
            source_atom = next((at for at in source_atoms if at['residue'] == residue and at['name'] == name), None)
            if source_atom:
                already_parsed_atoms.append(source_atom)
                file.write(source_atom['line'])
                continue
            missing_atoms = True
            print('res.lib missing atom ->  ' + residue + '  ' + name)
    return missing_atoms

# Save all atoms in res.lib
# The residue and atom names are enought to check if we have them in the topology
def get_reslib_atoms(reslib_filename : str) -> list:
    atoms = []
    with open(reslib_filename, 'r') as file:
        lines = file.readlines()
        for line in lines:
            # Skip comment lines
            if line[0] == '#':
                continue
            atom = { 'line': line }
            atom['residue'] = line[0:4].strip()
            atom['name'] = line[5:9].strip()
            atoms.append(atom)
    return atoms