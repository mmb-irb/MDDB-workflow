import json
import time
import re
import urllib.request
import requests

from typing import List, Tuple, Optional

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

# Map the structure aminoacids sequences against the standard reference sequences
# References are uniprot accession ids and they are optional
# For each reference, align the reference sequence with the topology sequence
# Chains which do not match any reference sequence will be blasted
# Note that an internet connection is required both to retireve the uniprot reference sequence and to do the blast
# NEVER FORGET: This system relies on the fact that topology chains are not repeated
def generate_map_online (structure : 'Structure', forced_references : List[str] = []) -> dict:
    # Store all the references which are got through this process
    # Note that not all references may be used at the end
    references = {}
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
        for uniprot_id in forced_references:
            reference, already_loaded = get_reference(uniprot_id)
            reference_sequences[reference['uniprot']] = reference['sequence']
            # Save the current whole reference object for later
            references[reference['uniprot']] = reference
    # Try to match all protein sequences with the available reference sequences
    # In case of match, objects in the 'protein_sequences' list are modified by adding the result
    # Finally, return True if all protein sequences were matched with the available reference sequences or False if not
    def match_sequences () -> bool:
        # Track each chain-reference alignment match and keep the score of successful alignments
        # Now for each structure sequence, align all reference sequences and keep the best alignment (if it meets the minimum)
        for structure_sequence in protein_sequences:
            for uniprot_id, reference_sequence in reference_sequences.items():
                # Align the structure sequence with the reference sequence
                align_results = align(reference_sequence, structure_sequence['sequence'])
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
        # Finally, return True if all protein sequences were matched with the available reference sequences or False if not
        return all([ structure_sequence['match']['ref'] for structure_sequence in protein_sequences ])
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
        uniprot_ids = blast(sequence)
        # Build a new reference from each resulting uniprot
        for uniprot_id in uniprot_ids:
            reference, already_loaded = get_reference(uniprot_id)
            reference_sequences[reference['uniprot']] = reference['sequence']
            # Save the current whole reference object for later
            references[reference['uniprot']] = reference
    # If we have every protein chain matched with a reference already then we stop here
    if match_sequences():
        export_references(protein_sequences)
        return format_topology_data(structure, protein_sequences)
    raise RuntimeError('The BLAST failed to find a matching reference sequence for at least one protein sequence')

# Try to match all protein sequences with the available reference sequences
# In case of match, objects in the 'protein_sequences' list are modified by adding the result
# Finally, return True if all protein sequences were matched with the available reference sequences or False if not
def match_sequences (protein_sequences : list, reference_sequences : dict) -> bool:
    # Track each chain-reference alignment match and keep the score of successful alignments
    # Now for each structure sequence, align all reference sequences and keep the best alignment (if it meets the minimum)
    for structure_sequence in protein_sequences:
        for uniprot_id, reference_sequence in reference_sequences.items():
            # Align the structure sequence with the reference sequence
            align_results = align(reference_sequence, structure_sequence['sequence'])
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
    # Finally, return True if all protein sequences were matched with the available reference sequences or False if not
    return all([ structure_sequence['match']['ref'] for structure_sequence in protein_sequences ])

# Export reference objects data to a json file
# This file is used by the loader to load new references to the database
# Note that all references are saved to this json file, even those which are already in the database
# It is the loader who is the responsible to check which references must be loaded and which ones are loaded already
def export_references (mapping_data : list):
    final_references = []
    final_uniprots = []
    for data in mapping_data:
        match = data['match']
        ref = match['ref']
        uniprot = ref['uniprot']
        if uniprot in final_uniprots:
            continue
        final_references.append(ref)
        final_uniprots.append(uniprot)
    with open(references_filename, 'w') as file:
        json.dump(final_references, file, indent=4)

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
        uniprot_id = reference['uniprot']
        if uniprot_id not in reference_ids:
            reference_ids.append(reference['uniprot'])
        reference_index = reference_ids.index(uniprot_id)
        residue_indices = data['residue_indices']
        for r, residue_number in enumerate(match['map']):
            if residue_number == None:
                continue
            residue_index = residue_indices[r]
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
        sequence_object = { 'sequence': sequence, 'residue_indices': residue_indices }
        sequences.append(sequence_object)
    return sequences

# Align two aminoacid sequences
# Return a list with the reference residue indexes (values)
# which match each new sequence residues indexes (indexes)
# Return also the score of the alignment
# Return None when there is not valid alignment at all
def align (ref_sequence : str, new_sequence : str) -> Optional[ Tuple[list, float] ]:

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
    print(format_alignment(*alignments[0]))
    score = alignments[0][2]
    # WARNING: Do not use 'aligned_sequence' length here since it has the total sequence length
    normalized_score = score / len(new_sequence)
    print('Normalized score: ' + str(normalized_score))

    # If the normalized score does not reaches the minimum we consider the alignment is not valid
    # It may happen when the reference goes for a specific chain but we must map all chains
    # This 1 has been found experimentally
    # Non maching sequence may return a 0.1-0.3 normalized score
    # Matching sequence may return >4 normalized score
    if normalized_score < 1:
        print('Not valid alignment')
        return None

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

# Given an EMBL accession, return the corresponding UniProt accession
# e.g. BAH02663.1 -> B6ZGN7
# None is return when the embl accession has not a UniProt id
def embl_to_uniprot (embl_id : str) -> str:
    # Request the EMBL API to retrieve the id associated data
    request_url = 'https://www.ebi.ac.uk/ena/browser/api/embl/' + embl_id + '?lineLimit=0&annotationOnly=true'
    try:
        with urllib.request.urlopen(request_url) as response:
            parsed_response = response.read().decode("utf-8")
    # If the accession is not found in the EMBL then we can stop here
    except urllib.error.HTTPError as error:
        if error.code == 400 or error.code == 404:
            return None
        else:
            raise ValueError('Something went wrong with the ENA request: ' + request_url)
    # Mine the uniprot accession
    result = re.search('db_xref="UniProtKB/TrEMBL:([0-9A-Z]*)"', parsed_response)
    if result == None:
        return None
    return result[1]

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
    return uniprot_ids

# Given an accession from the blast, try to find and equivalent uniprot accession
# Note that a single accession may be related to several uniprot ids
# In those cases we return all uniprot ids
# e.g. 6VW1 -> Q9BYF1, P0DTC2, P59594
# e.g. 5cnn -> P00533
# e.g. BAH02663.1 -> B6ZGN7
def accession_to_uniprot (accession : str) -> List[str]:
    print('Searching UniProt ID for ' + accession)
    # Try with UniProt ID mapping for GenBank accessions
    uniprot_id = genbank_to_uniprot(accession)
    if uniprot_id:
        print('    Got it from GenBank -> ' + uniprot_id)
        return [ uniprot_id ]
    # Try with ENA
    uniprot_id = embl_to_uniprot(accession)
    if uniprot_id:
        print('    Got it from ENA -> ' + uniprot_id)
        return [ uniprot_id ]
    # Try with PDB
    uniprot_ids = pdb_to_uniprot(accession)
    if uniprot_ids:
        print('    Got it from PDB -> ' + ', '.join(uniprot_ids))
        return uniprot_ids
    print('    Not found')
    return []

# Given an aminoacids sequence, return a list of uniprot ids
# At least 1 of the uniprot ids will belong to a very similar sequence
# This happens because the blast results return not only UniProt IDs but also other database accessions
# In those cases we find the equivalent (or most related) uniprot id
# In some databases a single accession may correspond to several uniprot accessions (e.g. PDB)
def blast (sequence : str) -> List[str]:
    print('Throwing blast...')
    result = NCBIWWW.qblast(
        program = "blastp",
        database = "nr",
        sequence = sequence,
    )
    parsed_result = xmltodict.parse(result.read())
    hits = parsed_result['BlastOutput']['BlastOutput_iterations']['Iteration']['Iteration_hits']['Hit']
    # Accessions for each sequence in the blast results are for the EMBL database
    # DANI: Si algun día tienes problemas porque te falta el '.1' al final del accession puedes sacarlo de Hit_id
    accessions = [ hit['Hit_accession'] for hit in hits ]
    # Return the first UniProt accesssion
    # Note that not all accessions have an equivalent UniProt accession (e.g. XP_032963186)
    for accession in accessions:
        uniprot_ids = accession_to_uniprot(accession)
        if uniprot_ids:
            return uniprot_ids
    raise ValueError('None of the blast result accessions corresponds to a UniProt reference sequence')

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
        if error.code == 404:
            return None
        else:
            raise ValueError('Something went wrong with the MDposit request: ' + request_url)
    return parsed_response

# Given a uniprot accession, use the uniprot API to request its data and then mine what is needed for the database
def get_uniprot_reference (uniprot_accession : str) -> dict:
    # Request Uniprot
    request_url = 'https://www.ebi.ac.uk/proteins/api/proteins/' + uniprot_accession
    try:
        with urllib.request.urlopen(request_url) as response:
            parsed_response = json.loads(response.read().decode("utf-8"))
    # If the accession is not found in UniProt then the id is not valid
    except urllib.error.HTTPError as error:
        if error.code == 400:
            raise ValueError('Something went wrong with the Uniprot request: ' + request_url)
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
    gene_names = ', '.join([ gene['name']['value'] for gene in parsed_response['gene'] ])
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

# Given a uniprot accession, get the reference object
# Try first asking to the MDposit database in case the reference exists already
# If not, retrieve UniProt data and build the reference object
# Return also a boolean to set if the reference already existed (True) or not (False)
def get_reference (uniprot_accession : str) -> Tuple[dict, bool]:
    reference = get_mdposit_reference(uniprot_accession)
    if reference:
        return reference, True
    return get_uniprot_reference(uniprot_accession), False

# -----------------------------------------------------------------------------------------------
# The following code has been copied from https://www.uniprot.org/help/id_mapping (10/10/2022)

# Set the UniProt API URL
API_URL = "https://rest.uniprot.org"
# Set the seconds to wait between every retry
POLLING_INTERVAL = 3

# Check if the job has already a result
def check_response (response):
    try:
        response.raise_for_status()
    except requests.HTTPError:
        print(response.json())
        raise

# Send the job
def submit_id_mapping (from_db, to_db, ids):
    request = requests.post(
        f"{API_URL}/idmapping/run",
        data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
    )
    check_response(request)
    return request.json()["jobId"]

# Set the session, which is used firther by 2 functions
session = requests.Session()

# Check repeatedly until results are returned
def check_id_mapping_results_ready (job_id):
    while True:
        request = session.get(f"{API_URL}/idmapping/status/{job_id}")
        check_response(request)
        j = request.json()
        if "jobStatus" in j:
            if j["jobStatus"] == "RUNNING":
                #print(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(j["jobStatus"])
        else:
            return True

# Find the link to the results
def get_id_mapping_results_link (job_id):
    url = f"{API_URL}/idmapping/details/{job_id}"
    request = session.get(url)
    check_response(request)
    return request.json()["redirectURL"]

# -----------------------------------------------------------------------------------------------

# Mine the uniprot id from the results
def mine_uniprot (results_url : str) -> Optional[str]:
    with urllib.request.urlopen(results_url) as response:
        parsed_response = json.loads(response.read().decode("utf-8"))
    results = parsed_response['results']
    # Return None if there are not results
    if len(results) == 0:
        return None
    # Get the first result accesion
    first_result = results[0]
    uniprot_id = first_result['to']['primaryAccession']
    return uniprot_id
    

# Main function: Given a GenBank protein accession find its corresponding UniProt id
def genbank_to_uniprot (genbank_accession : str) -> str:
    # Send the job and store its job id
    job_id = submit_id_mapping(from_db="EMBL-GenBank-DDBJ_CDS", to_db="UniProtKB", ids=[genbank_accession])
    # Wait until we have results
    if check_id_mapping_results_ready(job_id):
        # Find the link to the results
        results_url = get_id_mapping_results_link(job_id)
        # Mine the uniprot id from the results
        return mine_uniprot(results_url)
    raise SystemExit('Something went wrong with the results')