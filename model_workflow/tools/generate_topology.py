import json

# Generate a json file to be uploaded to the database with topology data
def generate_topology (
    structure,
    charges : list,
    residues_map : dict,
    output_topology_filename : str
):
    print(structure)
    # The structure will be a bunch of arrays
    # Atom data
    atom_names = []
    atom_elements = []
    for atom in structure.atoms:
        atom_names.append(atom.name)
        if atom.element:
            atom_elements.append(atom.element)
    # Residue data
    residue_names = []
    residue_numbers = []
    residue_icodes = []
    residue_atom_indices = []
    for residue in structure.residues:
        residue_names.append(residue.name)
        residue_numbers.append(residue.number)
        if residue.icode:
            residue_icodes.append(residue.icode)
        residue_atom_indices.append(residue.atom_indices)
    # Chain data
    chain_names = []
    chain_residue_indices = []
    for chain in structure.chains:
        chain_names.append(chain.name)
        chain_residue_indices.append(chain.residue_indices)

    # Setup the final output
    topology = {
        'atom_names': atom_names,
        'atom_elements': atom_elements,
        'atom_charges': charges,
        'residue_names': residue_names,
        'residue_numbers': residue_numbers,
        'residue_icodes': residue_icodes,
        'residue_atom_indices': residue_atom_indices,
        'chain_names': chain_names,
        'chain_residue_indices': chain_residue_indices,
        **residues_map,
    }
    with open(output_topology_filename, 'w') as file:
        json.dump(topology, file)