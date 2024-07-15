from typing import List
from model_workflow.utils.constants import NO_REFERABLE_FLAG, NOT_FOUND_FLAG

# Build the residue map from both proteins and ligands maps
# This is formatted as both the standard topology and metadata generators expect them
def generate_residue_mapping(
    protein_map : List[dict],
    ligand_map : List[dict],
    structure : 'Structure',
) -> dict:

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
            reference_types.append(reference_type)
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
        reference_types = None
        residue_reference_indices = None
        residue_reference_numbers = None

    residue_map = {
        'references': reference_ids,
        'reference_types': reference_types,
        'residue_reference_indices': residue_reference_indices,
        'residue_reference_numbers': residue_reference_numbers
    }

    return residue_map