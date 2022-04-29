import json

# Generate a json file to be uploaded to the database with topology data
def generate_topology (
    structure,
    charges : list,
    residues_map : dict,
    output_topology_filename : str
):
    # The structure will be a bunch of arrays
    # Atom data
    atom_names = []
    atom_elements = []
    atom_residue_indices = []
    for atom in structure.atoms:
        atom_names.append(atom.name)
        if atom.element:
            atom_elements.append(atom.element)
        atom_residue_indices.append(atom.residue.index)
        
    # Residue data
    residue_names = []
    residue_numbers = []
    residue_icodes = []
    residue_chain_indices = []
    for residue in structure.residues:
        residue_names.append(residue.name)
        residue_numbers.append(residue.number)
        if residue.icode:
            residue_icodes.append(residue.icode)
        residue_chain_indices.append(residue.chain.index)
    # Chain data
    chain_names = []
    for chain in structure.chains:
        chain_names.append(chain.name)

    # Setup the final output
    topology = {
        'atom_names': atom_names,
        'atom_elements': atom_elements,
        'atom_charges': charges,
        'atom_residue_indices': atom_residue_indices,
        'residue_names': residue_names,
        'residue_numbers': residue_numbers,
        'residue_icodes': residue_icodes,
        'residue_chain_indices': residue_chain_indices,
        'chain_names': chain_names,
        **residues_map,
    }
    with open(output_topology_filename, 'w') as file:
        json.dump(topology, file)