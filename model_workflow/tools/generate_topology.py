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
    structure_atoms = structure.atoms
    atom_count = len(structure_atoms)
    atom_names = [ None ] * atom_count
    atom_elements = [ None ] * atom_count
    atom_residue_indices = [ None ] * atom_count
    for index, atom in enumerate(structure_atoms):
        atom_names[index] = atom.name
        atom_elements[index] = atom.element
        atom_residue_indices[index] = atom.residue.index
        
    # Residue data
    structure_residues = structure.residues
    residue_count = len(structure_residues)
    residue_names = [ None ] * residue_count
    residue_numbers = [ None ] * residue_count
    residue_icodes = [ None ] * residue_count
    residue_chain_indices = [ None ] * residue_count
    for index, residue in enumerate(structure_residues):
        residue_names[index] = residue.name
        residue_numbers[index] = residue.number
        if residue.icode:
            residue_icodes[index] = residue.icode
        residue_chain_indices[index] = residue.chain.index

    # In case there are not icodes at all set the icodes list empty
    if all(icode == None for icode in residue_icodes):
        residue_icodes = []

    # Chain data
    structure_chains = structure.chains
    chain_count = len(structure_chains)
    chain_names = [ None ] * chain_count
    for index, chain in enumerate(structure_chains):
        chain_names[index] = chain.name

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