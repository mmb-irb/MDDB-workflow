import json
from typing import List

from model_workflow.utils.auxiliar import save_json
from model_workflow.utils.constants import NO_REFERABLE_FLAG, NOT_FOUND_FLAG

# Generate a json file to be uploaded to the database (mongo) with topology data
def generate_topology (
    structure : 'Structure',
    charges : List[int],
    protein_map : dict,
    ligand_map : List[dict],
    pbc_residues : List[int],
    output_topology_filepath : str
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

    # Set the atom bonds
    # In order to make it more standard sort atom bonds by their indices
    atom_bonds = [ sorted(atom_indices) for atom_indices in structure.bonds ]

    # Residue data
    structure_residues = structure.residues
    residue_count = len(structure_residues)
    residue_names = [ None ] * residue_count
    residue_numbers = [ None ] * residue_count
    # Icodes are saved as a dictionary since usually only a few residues have icodes (or none)
    # Resiude ids are used as keys and, when loaded to mongo, they will become strings
    # Saving icodes as an array would be unefficient since it will result in an array filled with nulls
    residue_icodes = {}
    residue_chain_indices = [ None ] * residue_count
    for index, residue in enumerate(structure_residues):
        residue_names[index] = residue.name
        residue_numbers[index] = residue.number
        if residue.icode:
            residue_icodes[index] = residue.icode
        residue_chain_indices[index] = residue.chain.index

    # In case there are not icodes at all set the icodes as None (i.e. null for mongo)
    if len(list(residue_icodes.keys())) == 0:
        residue_icodes = None

    # Chain data
    structure_chains = structure.chains
    chain_count = len(structure_chains)
    chain_names = [ None ] * chain_count
    for index, chain in enumerate(structure_chains):
        chain_names[index] = chain.name

    # Check we have charges and, if not, set charges as None (i.e. null for mongo)
    atom_charges = charges if charges and len(charges) > 0 else None
    if not atom_charges:
        print('WARNING: Topology is missing atom charges')

    # Reformat mapping data to the topology system

    # Add the reference type to each reference object
    for data in protein_map:
        data['type'] = 'protein'
    for data in ligand_map:
        data['type'] = 'ligand'

    # Get the count of residues from the structure
    residues_count = len(structure.residues)

    # Now format data
    reference_ids = []
    reference_types = []
    residue_reference_indices = [ None ] * residues_count
    residue_reference_numbers = [ None ] * residues_count

    for data in protein_map + ligand_map:
        match = data['match']
        # Get the reference index
        # Note that several matches may belong to the same reference and thus have the same index
        reference = match['ref']
        # If reference is missing at this point then it means we failed to find a matching reference
        if reference == None:
            continue
        # If we have the "no referable" flag
        if reference == NO_REFERABLE_FLAG:
            if NO_REFERABLE_FLAG not in reference_ids:
                reference_ids.append(NO_REFERABLE_FLAG)
            reference_index = reference_ids.index(NO_REFERABLE_FLAG)
            for residue_index in data['residue_indices']:
                residue_reference_indices[residue_index] = reference_index
            continue
        # If we have the "not found" flag
        if reference == NOT_FOUND_FLAG:
            if NOT_FOUND_FLAG not in reference_ids:
                reference_ids.append(NOT_FOUND_FLAG)
            reference_index = reference_ids.index(NOT_FOUND_FLAG)
            for residue_index in data['residue_indices']:
                residue_reference_indices[residue_index] = reference_index
            continue
        # Get the reference type
        reference_type = data['type']
        reference_types.append(reference_type)
        # Get the reference id
        reference_id = None
        if reference_type == 'protein':
            reference_id = reference['uniprot']
        elif reference_type == 'ligand':
            reference_id = reference['pubchem']
        else:
            raise ValueError('Not supported type ' + reference_type)
        # If we have a regular reference id (i.e. not a no referable / not found flag)
        if reference_id not in reference_ids:
            reference_ids.append(reference_id)
        reference_index = reference_ids.index(reference_id)
        # Set the topology reference number and index for each residue
        # Note that ligands do not have any residue reference numbering
        if reference_type == 'protein':
            for residue_index, residue_number in zip(data['residue_indices'], match['map']):
                if residue_number == None:
                    continue
                residue_reference_indices[residue_index] = reference_index
                residue_reference_numbers[residue_index] = residue_number
        if reference_type == 'ligand':
            for residue_index in data['residue_indices']:
                residue_reference_indices[residue_index] = reference_index
    # If there are not references at the end then set all fields as None, in order to save space
    if len(reference_ids) == 0:
        reference_ids = None
        residue_reference_indices = None
        residue_reference_numbers = None

    # Setup the final output
    topology = {
        'atom_names': atom_names,
        'atom_elements': atom_elements,
        'atom_charges': atom_charges,
        'atom_residue_indices': atom_residue_indices,
        'atom_bonds': atom_bonds,
        'residue_names': residue_names,
        'residue_numbers': residue_numbers,
        'residue_icodes': residue_icodes,
        'residue_chain_indices': residue_chain_indices,
        'chain_names': chain_names,
        # Residues map
        'references': reference_ids,
        'reference_types': reference_types,
        'residue_reference_indices': residue_reference_indices,
        'residue_reference_numbers': residue_reference_numbers,
        # Save also the pbc residues here
        'pbc_residues': pbc_residues
    }
    save_json(topology, output_topology_filepath)