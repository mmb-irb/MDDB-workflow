import os
import json
import urllib.request

from typing import List, Tuple, Optional, Union

from pathlib import Path

from model_workflow.tools.residues_library import residue_name_2_letter

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo
from Bio.Blast import NCBIWWW
import xmltodict

# from Bio import Align
# from Bio.Align import substitution_matrices

# # Set the aligner
# aligner = Align.PairwiseAligner(mode='local')
# aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
# aligner.open_gap_score = -10
# aligner.extend_gap_score = -0.5

# Set the name of the json file which will store the used reference objects
references_filename = 'references.json'

# Set a flag to represent a synthetic construct reference (i.e. no reference although it is protein)
synthetic_construct_flag = 'sc'

# Map the structure aminoacids sequences against the standard reference sequences
# References are uniprot accession ids and they are optional
# For each reference, align the reference sequence with the topology sequence
# Chains which do not match any reference sequence will be blasted
# Note that an internet connection is required both to retireve the uniprot reference sequence and to do the blast
# NEVER FORGET: This system relies on the fact that topology chains are not repeated
def generate_map_online (
    structure : 'Structure',
    forced_references : Optional[ Union[list,dict] ] = None,
    pdb_ids : List[str] = []
) -> dict:
    # Check if the forced references are strict (i.e. reference per chain, as a dictionary) or flexible (list of references)
    strict_references = type(forced_references) == dict
    # Check the synthetic construct flag not to be passed when references are not strict
    if not strict_references and synthetic_construct_flag in forced_references:
        raise SystemExit('WRONG INPUT: The "synthetic construct" flag cannot be passed in a list. You must use a chain keys dictionary (e.g. {"A":"sc"})')
    # Store all the references which are got through this process
    # Note that not all references may be used at the end
    references = {}
    # Given a uniprot accession, get the reference object
    # Try first asking to the MDposit database in case the reference exists already
    # If not, retrieve UniProt data and build the reference object
    # Return also a boolean to set if the reference already existed (True) or not (False)
    def get_reference (uniprot_accession : str) -> Tuple[dict, bool]:
        reference = references.get(uniprot_accession, None)
        if reference:
            return reference, True
        reference = get_mdposit_reference(uniprot_accession)
        if reference:
            return reference, True
        reference = get_uniprot_reference(uniprot_accession)
        return reference, False
    # Import local references, in case the references json file already exists
    if os.path.exists(references_filename):
        references = import_references()
    # Get the structure chain sequences
    structure_sequences = get_chain_sequences(structure)
    # Find out which chains are protein
    protein_sequences = []
    for structure_sequence in structure_sequences:
        sequence = structure_sequence['sequence']
        if next((letter for letter in sequence if letter != 'X'), None):
            structure_sequence['match'] = { 'ref': None, 'map': None, 'score': 0 }
            protein_sequences.append(structure_sequence)
    # For each input forced reference, get the reference sequence
    reference_sequences = {}
    if forced_references:
        forced_uniprot_ids = list(forced_references.values()) if strict_references else forced_references
        for uniprot_id in forced_uniprot_ids:
            # If instead of a uniprot id there is a 'synthetic construct' flag
            # A synthetic construct does not have a reference sequence by definition so we skip this process
            if uniprot_id == synthetic_construct_flag:
                continue
            # If reference is already in the list (i.e. it has been imported) then skip this process
            reference = references.get(uniprot_id, None)
            if reference:
                reference_sequences[uniprot_id] = reference['sequence']
                continue
            reference, already_loaded = get_reference(uniprot_id)
            reference_sequences[uniprot_id] = reference['sequence']
            # Save the current whole reference object for later
            references[reference['uniprot']] = reference
    # Save already tried alignments to not repeat the alignment further
    tried_alignments = { structure_sequence['name']: [] for structure_sequence in protein_sequences }
    # Try to match all protein sequences with the available reference sequences
    # In case of match, objects in the 'protein_sequences' list are modified by adding the result
    # Finally, return True if all protein sequences were matched with the available reference sequences or False if not
    def match_sequences () -> bool:
        # Track each chain-reference alignment match and keep the score of successful alignments
        # Now for each structure sequence, align all reference sequences and keep the best alignment (if it meets the minimum)
        for structure_sequence in protein_sequences:
            chain = structure_sequence['name']
            chain_tried_alignments = tried_alignments[chain]
            # In case references are forced per chain check if there is a reference for this chain and match according to this
            if strict_references:
                # Get the forced specific chain for this sequence, if any
                forced_reference = forced_references.get(chain, None)
                if forced_reference:
                    # If the chain has a specific forced reference then we must align it just once
                    # Skip this process in further matches
                    if structure_sequence['match']['ref']:
                        continue
                    # In case the forced reference is the synthetic construct flag
                    # Thus it has no reference sequence and we must not try to match it
                    # Actually, any match would be accidental and not correct
                    if forced_reference == synthetic_construct_flag:
                        structure_sequence['match'] = { 'ref': synthetic_construct_flag }
                        continue
                    # Get the forced reference sequence and align it to the chain sequence in order to build the map
                    reference_sequence = reference_sequences[forced_reference]
                    print(' Aligning chain ' + chain + ' with ' + forced_reference + ' reference sequence')
                    align_results = align(reference_sequence, structure_sequence['sequence'])
                    # The align must match or we stop here and warn the user
                    if not align_results:
                        raise SystemExit('Forced reference ' + chain + ' -> ' + forced_reference + ' does not match in sequence')
                    sequence_map, align_score = align_results
                    reference = references[forced_reference]
                    structure_sequence['match'] = { 'ref': reference, 'map': sequence_map, 'score': align_score }
                    continue
            for uniprot_id, reference_sequence in reference_sequences.items():
                # If this alignment has been tried already then skip it
                if uniprot_id in chain_tried_alignments:
                    continue
                # Align the structure sequence with the reference sequence
                print(' Aligning chain ' + chain + ' with ' + uniprot_id + ' reference sequence')
                align_results = align(reference_sequence, structure_sequence['sequence'])
                tried_alignments[chain].append(uniprot_id) # Save the alignment try, no matter if it works or not
                if not align_results:
                    continue
                # In case we have a valid alignment, check the alignment score is better than the current reference score (if any)
                sequence_map, align_score = align_results
                current_reference = structure_sequence['match']
                if current_reference['score'] > align_score:
                    continue
                reference = references[uniprot_id]
                # If the alignment is better then we impose the new reference
                structure_sequence['match'] = { 'ref': reference, 'map': sequence_map, 'score': align_score }
        # Sum up the current matching
        print(' Reference summary:')
        for structure_sequence in structure_sequences:
            name = structure_sequence['name']
            match = structure_sequence.get('match', None)
            if not match:
                print('   ' + name + ' -> Not protein')
                continue
            reference = structure_sequence['match'].get('ref', None)
            if not reference:
                print('   ' + name + ' -> ¿?')
                continue
            if reference == synthetic_construct_flag:
                print('   ' + name + ' -> Synthetic construct')
                continue
            uniprot_id = reference['uniprot']
            print('   ' + name + ' -> ' + uniprot_id)
        # Finally, return True if all protein sequences were matched with the available reference sequences or False if not
        return all([ structure_sequence['match']['ref'] for structure_sequence in protein_sequences ])
    # If we have every protein chain matched with a reference already then we stop here
    if match_sequences():
        export_references(protein_sequences)
        return format_topology_data(structure, protein_sequences)
    # If there are still any chain which is not matched with a reference then we need more references
    # To get them, retrieve all uniprot codes associated to the pdb codes, if any
    for pdb_id in pdb_ids:
        # Ask PDB
        uniprot_ids = pdb_to_uniprot(pdb_id)
        for uniprot_id in uniprot_ids:
            # Build a new reference from the resulting uniprot
            reference, already_loaded = get_reference(uniprot_id)
            if reference == None:
                continue
            reference_sequences[reference['uniprot']] = reference['sequence']
            # Save the current whole reference object for later
            references[reference['uniprot']] = reference
    # If we have every protein chain matched with a reference already then we stop here
    if match_sequences():
        export_references(protein_sequences)
        return format_topology_data(structure, protein_sequences)
    # If there are still any chain which is not matched with a reference then we need more references
    # To get them, we run a blast with each orphan chain sequence
    for structure_sequence in protein_sequences:
        # Skip already references chains
        if structure_sequence['match']['ref']:
            continue
        # Run the blast
        sequence = structure_sequence['sequence']
        uniprot_id = blast(sequence)
        # Build a new reference from the resulting uniprot
        reference, already_loaded = get_reference(uniprot_id)
        reference_sequences[reference['uniprot']] = reference['sequence']
        # Save the current whole reference object for later
        references[reference['uniprot']] = reference
        # If we have every protein chain matched with a reference already then we stop here
        if match_sequences():
            export_references(protein_sequences)
            return format_topology_data(structure, protein_sequences)
    raise RuntimeError('The BLAST failed to find a matching reference sequence for at least one protein sequence')

# Export reference objects data to a json file
# This file is used by the loader to load new references to the database
# Note that all references are saved to this json file, even those which are already in the database
# It is the loader who is the responsible to check which references must be loaded and which ones are loaded already
# Note that mapping data (i.e. which residue belongs to each reference) is not saved
def export_references (mapping_data : list):
    final_references = []
    final_uniprots = []
    for data in mapping_data:
        match = data['match']
        ref = match['ref']
        if ref == synthetic_construct_flag:
            continue
        uniprot = ref['uniprot']
        if uniprot in final_uniprots:
            continue
        final_references.append(ref)
        final_uniprots.append(uniprot)
    with open(references_filename, 'w') as file:
        json.dump(final_references, file, indent=4)

# Import reference json file so we do not have to rerun this process
def import_references () -> list:
    print(' Importing references from ' + references_filename)
    # Read the file
    with open(references_filename, 'r') as file:
        file_content = json.load(file)
    # Format data as the process expects to find it
    references = {}
    for reference in file_content:
        uniprot = reference['uniprot']
        references[uniprot] = reference
    return references

# Reformat mapping data to the topology system (introduced later)
def format_topology_data (structure : 'Structure', mapping_data : list) -> dict:
    # Get the count of residues from the structure
    residues_count = len(structure.residues)
    # Now format data
    reference_ids = []
    residue_reference_indices = [ None ] * residues_count
    residue_reference_numbers = [ None ] * residues_count
    for data in mapping_data:
        match = data['match']
        # Get the reference index
        # Note that several matches may belong to the same reference and thus have the same index
        reference = match['ref']
        # If we have the synthetic construct flag
        if reference == synthetic_construct_flag:
            if synthetic_construct_flag not in reference_ids:
                reference_ids.append(synthetic_construct_flag)
            reference_index = reference_ids.index(synthetic_construct_flag)
            for residue_index in data['residue_indices']:
                residue_reference_indices[residue_index] = reference_index
            continue
        # If we have a regular uniprot id
        uniprot_id = reference['uniprot']
        if uniprot_id not in reference_ids:
            reference_ids.append(reference['uniprot'])
        reference_index = reference_ids.index(uniprot_id)
        # Set the topology reference number and index for each residue
        for residue_index, residue_number in zip(data['residue_indices'], match['map']):
            if residue_number == None:
                continue
            residue_reference_indices[residue_index] = reference_index
            residue_reference_numbers[residue_index] = residue_number
    # If there are not references at the end then set all fields as None, in order to save space
    if len(reference_ids) == 0:
        reference_ids = None
        residue_reference_indices = None
        residue_reference_numbers = None
    # Return the 3 topology fields as they are in the database
    residues_map = {
        'references': reference_ids,
        'residue_reference_indices': residue_reference_indices,
        'residue_reference_numbers': residue_reference_numbers,
    }
    return residues_map

# Align a reference aminoacid sequence with each chain sequence in a topology
# NEVER FORGET: This system relies on the fact that topology chains are not repeated
def map_sequence (ref_sequence : str, structure : 'Structure') -> list:
    sequences = get_chain_sequences(structure)
    mapping = []
    for s in sequences:
        sequence = sequences[s]
        sequence_map = align(ref_sequence, sequence)
        mapping += sequence_map
    return mapping

# Get each chain name and aminoacids sequence in a topology
# Output format example: [ { 'sequence': 'VNLTT', 'indices': [1, 2, 3, 4, 5] }, ... ]
def get_chain_sequences (structure : 'Structure') -> list:
    sequences = []
    chains = structure.chains
    for chain in chains:
        name = chain.name
        sequence = ''
        residue_indices = []
        for residue in chain.residues:
            letter = residue_name_2_letter(residue.name, 'aminoacids')
            sequence += letter
            residue_indices.append(residue.index)
        # Save sequences by chain name (i.e. chain id or chain letter)
        sequence_object = { 'name': name, 'sequence': sequence, 'residue_indices': residue_indices }
        sequences.append(sequence_object)
    return sequences

# Align two aminoacid sequences
# Return a list with the reference residue indexes (values)
# which match each new sequence residues indexes (indexes)
# Return also the score of the alignment
# Return None when there is not valid alignment at all
# Set verbose = True to see a visual summary of the sequence alignments in the logs
def align (ref_sequence : str, new_sequence : str, verbose : bool = False) -> Optional[ Tuple[list, float] ]:

    #print('- REFERENCE\n' + ref_sequence + '\n- NEW\n' + new_sequence)

    # If the new sequence is all 'X' stop here, since this would make the alignment infinite
    # Then an array filled with None is returned
    if all([ letter == 'X' for letter in new_sequence ]):
        return None

    # Return the new sequence as best aligned as possible with the reference sequence
    alignments = pairwise2.align.localds(ref_sequence, new_sequence, MatrixInfo.blosum62, -10, -0.5)
    # DANI: Habría que hacerlo de esta otra forma según el deprecation warning (arriba hay más código)
    # DANI: El problema es que el output lo tiene todo menos la sequencia en formato alienada
    # DANI: i.e. formato '----VNLTT', que es justo el que necesito
    #alignments = aligner.align(ref_sequence, new_sequence)

    # In case there are no alignments it means the current chain has nothing to do with this reference
    # Then an array filled with None is returned
    if len(alignments) == 0:
        return None

    # Several alignments may be returned, specially when it is a difficult or impossible alignment
    # Output format example: '----VNLTT'
    best_alignment = alignments[0]
    aligned_sequence = best_alignment[1]
    score = alignments[0][2]
    # WARNING: Do not use 'aligned_sequence' length here since it has the total sequence length
    normalized_score = score / len(new_sequence)
    if verbose:
        print(format_alignment(*alignments[0]))

    # If the normalized score does not reaches the minimum we consider the alignment is not valid
    # It may happen when the reference goes for a specific chain but we must map all chains
    # This 1 has been found experimentally
    # Non maching sequence may return a 0.1-0.3 normalized score
    # Matching sequence may return >4 normalized score
    if normalized_score < 1:
        print('    Not valid alignment')
        return None

    # Tell the user about the success
    beautiful_normalized_score = round(normalized_score * 100) / 100
    print('    Valid alignment -> Normalized score = ' + str(beautiful_normalized_score))

    # Match each residue
    aligned_mapping = []
    aligned_index = 0
    for l, letter in enumerate(aligned_sequence):
        # Guions are skipped
        if letter == '-':
            continue
        # Get the current residue of the new sequence
        equivalent_letter = new_sequence[aligned_index]
        if not letter == equivalent_letter:
            raise SystemExit('Something was wrong :S')
        # 'X' residues cannot be mapped since reference sequences should never have any 'X'
        if letter == 'X':
            aligned_mapping.append(None)
            aligned_index += 1
            continue
        # Otherwise add the equivalent aligned index to the mapping
        # WARNING: Add +1 since uniprot residue counts start at 1, not 0
        aligned_mapping.append(l + 1)
        aligned_index += 1

    return aligned_mapping, normalized_score

# Given an aminoacids sequence, return a list of uniprot ids
# Note that we are blasting against UniProtKB / Swiss-Prot so results will always be valid UniProt accessions
# WARNING: This always means results will correspond to curated entries only
#   If your sequence is from an exotic organism the result may be not from it but from other more studied organism
def blast (sequence : str) -> List[str]:
    print('Throwing blast...')
    result = NCBIWWW.qblast(
        program = "blastp",
        database = "swissprot", # UniProtKB / Swiss-Prot
        sequence = sequence,
    )
    parsed_result = xmltodict.parse(result.read())
    hits = parsed_result['BlastOutput']['BlastOutput_iterations']['Iteration']['Iteration_hits']['Hit']
    # Get the first result only
    result = hits[0]
    # Return the accession
    # DANI: Si algun día tienes problemas porque te falta el '.1' al final del accession puedes sacarlo de Hit_id
    accession = result['Hit_accession']
    print('Result: ' + accession)
    print(result['Hit_def'])
    return accession

# Given a uniprot accession, use the MDposit API to request its data in case it is already in the database
def get_mdposit_reference (uniprot_accession : str) -> Optional[dict]:
    # Request MDposit
    #request_url = 'https://mdposit-dev.bsc.es/api/rest/v1/references/' + uniprot_accession
    # DANI: De momento la query hay que hacerla a cv19 para ahorrarnos problemas con los certificados
    # DANI: Las references no se filtran por collections así que las tendremos todas aunque estos en cv19
    request_url = 'https://bioexcel-cv19-dev.bsc.es/api/rest/v1/references/' + uniprot_accession
    try:
        with urllib.request.urlopen(request_url) as response:
            parsed_response = json.loads(response.read().decode("utf-8"))
    # If the accession is not found in the database then we stop here
    except urllib.error.HTTPError as error:
        # If the uniprot accession is not yet in the MDposit references then return None
        if error.code == 404:
            return None
        else:
            print('Error when requesting ' + request_url)
            raise ValueError('Something went wrong with the MDposit request (error ' + str(error.code) + ')')
    return parsed_response

# Given a uniprot accession, use the uniprot API to request its data and then mine what is needed for the database
def get_uniprot_reference (uniprot_accession : str) -> dict:
    # Request Uniprot
    request_url = 'https://www.ebi.ac.uk/proteins/api/proteins/' + uniprot_accession
    parsed_response = None
    try:
        with urllib.request.urlopen(request_url) as response:
            parsed_response = json.loads(response.read().decode("utf-8"))
    # If the accession is not found in UniProt then the id is not valid
    except urllib.error.HTTPError as error:
        print('Error when requesting ' + request_url)
        raise ValueError('Something went wrong with the Uniprot request (error ' + str(error.code) + ')')
    # If we have not a response at this point then it may mean we are trying to access an obsolete entry (e.g. P01607)
    if parsed_response == None:
        print('WARNING: Cannot find UniProt entry for accession ' + uniprot_accession)
        return None
    # Get the full protein name
    protein_data = parsed_response['protein']
    protein_name_data = protein_data.get('recommendedName', None)
    # DANI: It is possible that the 'recommendedName' is missing if it is not a reviewed UniProt entry
    if not protein_name_data:
        print('WARNING: The UniProt accession ' + uniprot_accession + ' is missing the recommended name. You should consider changing the reference.')
        protein_name_data = protein_data.get('submittedName', None)[0]
    if not protein_name_data:
        raise ValueError('Unexpected structure in UniProt response for accession ' + uniprot_accession)
    protein_name = protein_name_data['fullName']['value']
    # Get the gene names as a single string
    gene_names = []
    # WARNING: Some uniprot entries are missing gene names (e.g. P00718)
    genes = parsed_response.get('gene', [])
    for gene in genes:
        gene_name = gene.get('name', None)
        if not gene_name:
            gene_name = gene.get('orfNames', [])[0]
        if not gene_name:
            raise ValueError('The uniprot response for ' + uniprot_accession + ' has an unexpected format')
        gene_names.append(gene_name['value'])
    gene_names = ', '.join(gene_names) if len(gene_names) > 0 else None
    # Get the organism name
    organism = parsed_response['organism']['names'][0]['value']
    # Get the aminoacids sequence
    sequence = parsed_response['sequence']['sequence']
    # Get interesting regions to be highlighted in the client
    domains = []
    for feature in parsed_response['features']:
        if feature['type'] != "CHAIN":
            continue
        name = feature['description']
        comments = [ comment for comment in parsed_response['comments'] if name == comment.get('molecule', None) ]
        comment_text = [ comment['text'][0]['value'] for comment in comments if comment.get('text', False) ]
        description = '\n\n'.join(comment_text)
        domains.append({
            'name': name,
            'description': description,
            # Set the representations to be configured in the client viewer to show this domain
            'representations':[{
                'name': name,
                'selection': feature['begin'] + '-' + feature['end']
            }]
        })
    return {
        'name': protein_name,
        'gene': gene_names,
        'organism': organism,
        'uniprot': uniprot_accession,
        'sequence': sequence,
        'domains': domains
    }

# Given a pdb Id, get its uniprot id
# e.g. 6VW1 -> Q9BYF1, P0DTC2, P59594
def pdb_to_uniprot (pdb_id : str) -> List[str]:
    # Request the MMB service to retrieve pdb data
    request_url = 'https://mmb.irbbarcelona.org/api/pdb/' + pdb_id + '/entry'
    try:
        with urllib.request.urlopen(request_url) as response:
            parsed_response = json.loads(response.read().decode("utf-8"))
    # If the accession is not found in the PDB then we can stop here
    except urllib.error.HTTPError as error:
        if error.code == 404:
            return None
        else:
            raise ValueError('Something went wrong with the PDB request: ' + request_url)
    # Get the uniprot accessions
    uniprot_ids = [ uniprot['_id'] for uniprot in parsed_response['uniprotRefs'] ]
    print(' References for PDB code ' + pdb_id + ': ' + ', '.join(uniprot_ids))
    return uniprot_ids