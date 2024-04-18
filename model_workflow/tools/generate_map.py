import os
import sys
import json
import urllib.request

from typing import List, Tuple, Optional, Union

from pathlib import Path

from model_workflow.tools.residues_library import residue_name_2_letter
from model_workflow.utils.auxiliar import load_json, save_json
from model_workflow.utils.constants import OUTPUT_REFERENCES_FILENAME, REFERENCE_SEQUENCE_FLAG, NO_REFERABLE_FLAG, NOT_FOUND_FLAG

import xmltodict

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Blast import NCBIWWW

# Bio.SubsMat library shows a huge deprecation warning
# This is disturbing so we supress its output while importing it

# Save current stderr to further restore it
stderr_backup = sys.stderr
# Suppress stderr
sys.stderr = None
# Import the function
from Bio.SubsMat import MatrixInfo
# Restore stderr
sys.stderr = stderr_backup

# from Bio import Align
# from Bio.Align import substitution_matrices

# # Set the aligner
# aligner = Align.PairwiseAligner(mode='local')
# aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
# aligner.open_gap_score = -10
# aligner.extend_gap_score = -0.5

# Map the structure aminoacids sequences against the standard reference sequences
# References are uniprot accession ids and they are optional
# For each reference, align the reference sequence with the topology sequence
# Chains which do not match any reference sequence will be blasted
# Note that an internet connection is required both to retireve the uniprot reference sequence and to do the blast
# NEVER FORGET: This system relies on the fact that topology chains are not repeated
def generate_protein_mapping (
    structure : 'Structure',
    register : dict,
    mercy : List[str] = [],
    forced_references : Union[list,dict] = [],
    pdb_ids : List[str] = [],
) -> dict:
    # Remove previous warnings, if any
    register.remove_warnings(REFERENCE_SEQUENCE_FLAG)
    # Forced references must be list or dict
    # If it is none then we set it as an empty list
    if forced_references == None:
        forced_references = []
    # Check if the forced references are strict (i.e. reference per chain, as a dictionary) or flexible (list of references)
    strict_references = type(forced_references) == dict
    # Check the "no referable" flag not to be passed when references are not strict
    if not strict_references and NO_REFERABLE_FLAG in forced_references:
        raise SystemExit('WRONG INPUT: The "no referable" flag cannot be passed in a list. You must use a chain keys dictionary (e.g. {"A":"' + NO_REFERABLE_FLAG + '"})')
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
    imported_references = None
    if os.path.exists(OUTPUT_REFERENCES_FILENAME):
        imported_references = import_references()
        # Append the imported references to the overall references pool
        for k,v in imported_references.items():
            references[k] = v
    # Get the structure chain sequences
    parsed_chains = get_parsed_chains(structure)
    # Find out which chains are protein
    protein_parsed_chains = []
    for chain_data in parsed_chains:
        sequence = chain_data['sequence']
        if next((letter for letter in sequence if letter != 'X'), None):
            chain_data['match'] = { 'ref': None, 'map': None, 'score': 0 }
            protein_parsed_chains.append(chain_data)
    # If there are no protein sequences then there is no need to map anything
    if len(protein_parsed_chains) == 0:
        print(' There are no protein sequences')
        return protein_parsed_chains
    # For each input forced reference, get the reference sequence
    reference_sequences = {}
    # Save already tried alignments to not repeat the alignment further
    tried_alignments = { chain_data['name']: [] for chain_data in protein_parsed_chains }
    # Set a function to try to match all protein sequences with the available reference sequences
    # In case of match, objects in the 'protein_parsed_chains' list are modified by adding the result
    # Finally, return True if all protein sequences were matched with the available reference sequences or False if not
    def match_sequences () -> bool:
        # Track each chain-reference alignment match and keep the score of successful alignments
        # Now for each structure sequence, align all reference sequences and keep the best alignment (if it meets the minimum)
        for chain_data in protein_parsed_chains:
            chain = chain_data['name']
            chain_tried_alignments = tried_alignments[chain]
            # In case references are forced per chain check if there is a reference for this chain and match according to this
            if strict_references:
                # Get the forced specific chain for this sequence, if any
                forced_reference = forced_references.get(chain, None)
                if forced_reference:
                    # If the chain has a specific forced reference then we must align it just once
                    # Skip this process in further matches
                    if chain_data['match']['ref']:
                        continue
                    # In case the forced reference is the "no referable" flag
                    # Thus it has no reference sequence and we must not try to match it
                    # Actually, any match would be accidental and not correct
                    if forced_reference == NO_REFERABLE_FLAG:
                        chain_data['match'] = { 'ref': NO_REFERABLE_FLAG }
                        continue
                    # In case the forced reference is the "not found" flag
                    # This should not happend but we may be using references as forced references, so just in case
                    if forced_reference == NOT_FOUND_FLAG:
                        chain_data['match'] = { 'ref': NOT_FOUND_FLAG }
                        continue
                    # Get the forced reference sequence and align it to the chain sequence in order to build the map
                    reference_sequence = reference_sequences[forced_reference]
                    print(' Aligning chain ' + chain + ' with ' + forced_reference + ' reference sequence')
                    align_results = align(reference_sequence, chain_data['sequence'])
                    # The align must match or we stop here and warn the user
                    if not align_results:
                        raise SystemExit('Forced reference ' + chain + ' -> ' + forced_reference + ' does not match in sequence')
                    sequence_map, align_score = align_results
                    reference = references[forced_reference]
                    chain_data['match'] = { 'ref': reference, 'map': sequence_map, 'score': align_score }
                    continue
            # If the chain has already a match and this match is among the forced references then stop here
            # Forced references have priority and this avoids having a match with a not forced reference further
            if chain_data['match']['ref'] and chain_data['match']['ref']['uniprot'] in forced_references:
                continue
            # Iterate over the deifferent available reference sequences
            for uniprot_id, reference_sequence in reference_sequences.items():
                # If this alignment has been tried already then skip it
                if uniprot_id in chain_tried_alignments:
                    continue
                # Align the structure sequence with the reference sequence
                # NEVER FORGET: This system relies on the fact that topology chains are not repeated
                print(' Aligning chain ' + chain + ' with ' + uniprot_id + ' reference sequence')
                align_results = align(reference_sequence, chain_data['sequence'])
                tried_alignments[chain].append(uniprot_id) # Save the alignment try, no matter if it works or not
                if not align_results:
                    continue
                # In case we have a valid alignment, check the alignment score is better than the current reference score (if any)
                sequence_map, align_score = align_results
                current_reference = chain_data['match']
                if current_reference['score'] > align_score:
                    continue
                reference = references[uniprot_id]
                # If the alignment is better then we impose the new reference
                chain_data['match'] = { 'ref': reference, 'map': sequence_map, 'score': align_score }
        # Sum up the current matching
        print(' Reference summary:')
        for chain_data in parsed_chains:
            name = chain_data['name']
            match = chain_data.get('match', None)
            if not match:
                print('   ' + name + ' -> Not protein')
                continue
            reference = chain_data['match'].get('ref', None)
            if not reference:
                print('   ' + name + ' -> ¿?')
                continue
            if reference == NO_REFERABLE_FLAG:
                print('   ' + name + ' -> No referable')
                continue
            if reference == NOT_FOUND_FLAG:
                print('   ' + name + ' -> Not found')
                continue
            uniprot_id = reference['uniprot']
            print('   ' + name + ' -> ' + uniprot_id)
        # Export already matched references
        export_references(protein_parsed_chains)
        # Finally, return True if all protein sequences were matched with the available reference sequences or False if not
        allright = all([ chain_data['match']['ref'] for chain_data in protein_parsed_chains ])            
        return allright
    # --- End of match_sequences function --------------------------------------------------------------------------------
    # First use the forced references for the matching
    if forced_references:
        forced_uniprot_ids = list(forced_references.values()) if strict_references else forced_references
        for uniprot_id in forced_uniprot_ids:
            # If instead of a uniprot id there is a 'no referable' flag then we skip this process
            if uniprot_id == NO_REFERABLE_FLAG:
                continue
            # If instead of a uniprot id there is a 'not found' flag then we skip this process
            # This should not happend but we may be using references as forced references, so just in case
            if uniprot_id == NOT_FOUND_FLAG:
                continue
            # If reference is already in the list (i.e. it has been imported) then skip this process
            reference = references.get(uniprot_id, None)
            if reference:
                reference_sequences[uniprot_id] = reference['sequence']
                continue
            # Find the reference data for the given uniprot id
            reference, already_loaded = get_reference(uniprot_id)
            reference_sequences[uniprot_id] = reference['sequence']
            # Save the current whole reference object for later
            references[reference['uniprot']] = reference
        # Now that we have all forced references data perform the matching
        # If we have every protein chain matched with a reference then we stop here
        print(' Using only forced references from the inputs file')
        if match_sequences():
            return protein_parsed_chains
    # Now add the imported references to reference sequences. Thus now they will be 'matchable'
    # Thus now they will be 'matchable', so try to match sequences again in case any of the imported references has not been tried
    # It was not done before since we want forced references to have priority
    if imported_references:
        need_rematch = False
        for uniprot_id, reference in imported_references.items():
            # If the imported reference has been aligned already (i.e. it was a forced reference)
            if uniprot_id in reference_sequences:
                continue
            # Otherwise, include it
            need_rematch = True
            reference_sequences[uniprot_id] = reference['sequence']
        # If there was at least one imported reference missing then rerun the matching
        if need_rematch:
            print(' Using also references imported from references.json')
            if match_sequences():
                return protein_parsed_chains
    # If there are still any chain which is not matched with a reference then we need more references
    # To get them, retrieve all uniprot codes associated to the pdb codes, if any
    if pdb_ids and len(pdb_ids) > 0:
        for pdb_id in pdb_ids:
            # Ask PDB
            uniprot_ids = pdb_to_uniprot(pdb_id)
            if uniprot_ids == None:
                continue 
            for uniprot_id in uniprot_ids:
                # Build a new reference from the resulting uniprot
                reference, already_loaded = get_reference(uniprot_id)
                if reference == None:
                    continue
                reference_sequences[reference['uniprot']] = reference['sequence']
                # Save the current whole reference object for later
                references[reference['uniprot']] = reference
        # If we have every protein chain matched with a reference already then we stop here
        print(' Using also references related to PDB ids from the inputs file')
        if match_sequences():
            return protein_parsed_chains
    # If there are still any chain which is not matched with a reference then we need more references
    # To get them, we run a blast with each orphan chain sequence
    for chain_data in protein_parsed_chains:
        # Skip already references chains
        if chain_data['match']['ref']:
            continue
        # Run the blast
        sequence = chain_data['sequence']
        uniprot_id = blast(sequence)
        if not uniprot_id:
            chain_data['match'] = { 'ref': NOT_FOUND_FLAG }
            continue
        # Build a new reference from the resulting uniprot
        reference, already_loaded = get_reference(uniprot_id)
        reference_sequences[reference['uniprot']] = reference['sequence']
        # Save the current whole reference object for later
        references[reference['uniprot']] = reference
        # If we have every protein chain matched with a reference already then we stop here
        print(' Using also references from blast')
        if match_sequences():
            return protein_parsed_chains
    # At this point we should have macthed all sequences
    # If not, kill the process unless mercy was given
    must_be_killed = REFERENCE_SEQUENCE_FLAG not in mercy
    if must_be_killed:
        raise SystemExit('BLAST failed to find a matching reference sequence for at least one protein sequence')
    print('WARNING: BLAST failed to find a matching reference sequence for at least one protein sequence')
    register.add_warning(REFERENCE_SEQUENCE_FLAG, 'There is at least one protein region which is not mapped to any reference sequence')
    return protein_parsed_chains

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
        if not ref or ref == NO_REFERABLE_FLAG or ref == NOT_FOUND_FLAG:
            continue
        uniprot = ref['uniprot']
        if uniprot in final_uniprots:
            continue
        final_references.append(ref)
        final_uniprots.append(uniprot)
    # If there are no references to exported then do not egenrate the json file
    if len(final_references) == 0:
        return
    # Write references to a json file
    save_json(final_references, OUTPUT_REFERENCES_FILENAME, indent = 4)

# Import reference json file so we do not have to rerun this process
def import_references () -> list:
    print(' Importing references from ' + OUTPUT_REFERENCES_FILENAME)
    # Read the file
    file_content = load_json(OUTPUT_REFERENCES_FILENAME)
    # Format data as the process expects to find it
    references = {}
    for reference in file_content:
        uniprot = reference['uniprot']
        references[uniprot] = reference
    return references

# Get each chain name and aminoacids sequence in a topology
# Output format example: [ { 'sequence': 'VNLTT', 'indices': [1, 2, 3, 4, 5] }, ... ]
def get_parsed_chains (structure : 'Structure') -> list:
    parsed_chains = []
    chains = structure.chains
    for chain in chains:
        name = chain.name
        sequence = ''
        residue_indices = []
        for residue in chain.residues:
            letter = residue_name_2_letter(residue.name, 'aminoacids')
            sequence += letter
            residue_indices.append(residue.index)
        sequence_object = { 'name': name, 'sequence': sequence, 'residue_indices': residue_indices }
        parsed_chains.append(sequence_object)
    return parsed_chains

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
def blast (sequence : str) -> Optional[str]:
    print('Throwing blast...')
    result = NCBIWWW.qblast(
        program = "blastp",
        database = "swissprot", # UniProtKB / Swiss-Prot
        sequence = sequence,
    )
    parsed_result = xmltodict.parse(result.read())
    hits = parsed_result['BlastOutput']['BlastOutput_iterations']['Iteration']['Iteration_hits']
    # When there is no result return None
    # Note that this is possible although hardly unprobable
    if not hits:
        return None
    # Get the first result only
    result = hits['Hit'][0]
    # Return the accession
    # DANI: Si algun día tienes problemas porque te falta el '.1' al final del accession puedes sacarlo de Hit_id
    accession = result['Hit_accession']
    print('Result: ' + accession)
    return accession

# Given a uniprot accession, use the MDposit API to request its data in case it is already in the database
def get_mdposit_reference (uniprot_accession : str) -> Optional[dict]:
    # Request MDposit
    request_url = 'https://mdposit-dev.mddbr.eu/api/rest/v1/references/' + uniprot_accession
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
    # This error may occur if there is no internet connection
    except urllib.error.URLError as error:
        print('Error when requesting ' + request_url)
        raise ValueError('Something went wrong with the MDposit request')
    return parsed_response

# Given a uniprot accession, use the uniprot API to request its data and then mine what is needed for the database
def get_uniprot_reference (uniprot_accession : str) -> Optional[dict]:
    # Request Uniprot
    request_url = 'https://www.ebi.ac.uk/proteins/api/proteins/' + uniprot_accession
    parsed_response = None
    try:
        with urllib.request.urlopen(request_url) as response:
            parsed_response = json.loads(response.read().decode("utf-8"))
    # If the accession is not found in UniProt then the id is not valid
    except urllib.error.HTTPError as error:
        if error.code == 404:
            print('WARNING: Cannot find UniProt entry for accession ' + uniprot_accession)
            return None
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
    # WARNING: Some uniprot entries have names in non-canonical keys (e.g. orfNames, olnNames)
    #       These names are not shown event in the uniprot web page so we also skip them
    genes = parsed_response.get('gene', [])
    for gene in genes:
        gene_name = gene.get('name', None)
        if not gene_name:
            continue
        gene_names.append(gene_name['value'])
    gene_names = ', '.join(gene_names) if len(gene_names) > 0 else None
    # Get the organism name
    organism = parsed_response['organism']['names'][0]['value']
    # Get the aminoacids sequence
    sequence = parsed_response['sequence']['sequence']
    # Get interesting regions to be highlighted in the client
    domains = []
    # WARNING: Some uniprot entries are missing features (e.g. O27908)
    features = parsed_response.get('features', [])
    # Get comments data which may be useful to further build domain descriptions
    comments_data = parsed_response.get('comments', None)
    # There are many types of features but in this case we will focus on domanins and similar
    target_types = ['CHAIN', 'REGION', 'DOMAIN', 'MOTIF', 'SITE']
    for feature in features:
        # Skip features of other types
        if feature['type'] not in target_types:
            continue
        # Get the domain name
        name = feature['description']
        # Build the domain description from the coments data
        description = None
        if comments_data:
            comments = [ comment for comment in parsed_response['comments'] if name == comment.get('molecule', None) ]
            comment_text = [ comment['text'][0]['value'] for comment in comments if comment.get('text', False) ]
            description = '\n\n'.join(comment_text)
        # Set the domain selection
        # The domain 'begin' and 'end' values may include non-numeric symbols such as '~', '>' or '<'
        # These values are usually ignored or replaced by '?' in the UniProt web client
        # There may be not numeric value at all (e.g. Q99497)
        # In these cases uniprot shows the domain but it does not allow to use its functionallities 
        # e.g. you can not blast using the domain sequence
        # It makes not sense having a domain with unkown range to me so we skip these domains
        begin = feature['begin'].replace('~', '').replace('>', '').replace('<', '')
        end = feature['end'].replace('~', '').replace('>', '').replace('<', '')
        if begin == '' or end == '':
            continue
        selection = begin + '-' + end
        # If we already have a domain with the same name then join both domains
        # For instance, you may have several repetitions of the 'Disordered' region
        already_existing_domain = next((domain for domain in domains if domain['name'] == name), None)
        if already_existing_domain:
            already_existing_domain['selection'] += ', ' + selection
            continue
        # Otherwise, create a new domain
        domain = {
            'name': name,
            'selection': selection
        }
        # Add a description only if we succesfully mined it
        # Note that only features tagged as CHAIN have comments (not sure about this)
        if description:
            domain['description'] = description
        # Add the new domain to the list
        domains.append(domain)
    # Mine protein functions from Gene Ontology references
    # Get database references
    db_references = parsed_response['dbReferences']
    # Get references from Gene Ontology (GO) only
    go_references = [ ref for ref in db_references if ref['type'] == 'GO' ]
    # A Gene Ontology entry may be one of three types:
    # Cellular Component (C), Molecular Function (F) and Biological Process (P)
    # In this case we are interested in protein function only so we will keep those with the 'F' header
    functions = [ ref['properties']['term'][2:] for ref in go_references if ref['properties']['term'][0:2] == 'F:' ]
    # Set the final reference to be uploaded to the database
    return {
        'name': protein_name,
        'gene': gene_names,
        'organism': organism,
        'uniprot': uniprot_accession,
        'sequence': sequence,
        'domains': domains,
        'functions': functions
    }

# Given a pdb Id, get its uniprot id
# e.g. 6VW1 -> Q9BYF1, P0DTC2, P59594
def pdb_to_uniprot (pdb_id : str) -> Optional[ List[str] ]:
    # Request the MMB service to retrieve pdb data
    request_url = 'http://mdb-login.bsc.es/api/pdb/' + pdb_id + '/entry'
    try:
        with urllib.request.urlopen(request_url) as response:
            parsed_response = json.loads(response.read().decode("utf-8"))
    # If the accession is not found in the PDB then we can stop here
    except urllib.error.HTTPError as error:
        if error.code == 404:
            print(' PDB code ' + pdb_id + ' not found')
            return None
        else:
            raise ValueError('Something went wrong with the PDB request: ' + request_url)
    # Get the uniprot accessions
    uniprot_ids = [ uniprot['_id'] for uniprot in parsed_response['uniprotRefs'] ]
    print(' References for PDB code ' + pdb_id + ': ' + ', '.join(uniprot_ids))
    return uniprot_ids

# This function is used by the generate_metadata script
# 1. Get structure sequences
# 2. Calculate which reference domains are covered by the previous sequence
# 3. In case it is a covid spike, align the previous sequence against all saved variants (they are in 'utils')
from model_workflow.resources.covid_variants import get_variants
def get_sequence_metadata (structure : 'Structure', residue_map : dict) -> dict:
    # First mine the sequences from the structure
    # Set a different sequence for each chain
    sequences = [ chain.get_sequence() for chain in structure.chains ]
    # Get values from the residue map
    # Get protein references from the residues map
    reference_ids = []
    for ref, ref_type in zip(residue_map['references'], residue_map['reference_types']):
        if ref_type == 'protein':
            reference_ids.append(ref)
    residue_reference_numbers = residue_map['residue_reference_numbers']
    residue_reference_indices = residue_map['residue_reference_indices']
    # Load references data, which should already be save to the references data file
    references_data = import_references() if os.path.exists(OUTPUT_REFERENCES_FILENAME) else {}
    # In case we have the SARS-CoV-2 spike among the references, check also which is the variant it belongs to
    # DANI: Esto es un parche. En un futuro buscaremos una manera mejor de comprovar variantes en cualquier contexto
    variant = None
    covid_spike_reference_id = 'P0DTC2'
    if covid_spike_reference_id in reference_ids:
        # Load the reference variant sequences
        covid_spike_variants = get_variants()
        # Load the reference data and find a chain which belongs to the spike
        # We consider having only one variant, so all chains should return the same result
        reference_data = references_data[covid_spike_reference_id]
        covid_spike_reference_index = reference_ids.index(covid_spike_reference_id)
        sample_residue_index = next((residue_index for residue_index, reference_index in enumerate(residue_reference_indices) if reference_index == covid_spike_reference_index), None)
        if sample_residue_index == None:
            raise SystemExit('Failed to find residue belonging to SARS-CoV-2 variant') 
        sample_chain_sequence = structure.residues[sample_residue_index].chain.get_sequence()
        # Align the sample chain sequence against all variants to find the best match
        highest_score = 0
        for variant_name, variant_sequence in covid_spike_variants.items():
            align_results = align(variant_sequence, sample_chain_sequence)
            if not align_results:
                continue
            mapping, score = align_results
            if score > highest_score:
                highest_score = score
                variant = variant_name
        # At this point there should be a result
        if not variant:
            raise SystemExit('Something went wrong trying to find the SARS-CoV-2 variant')
    # Set which reference domains are present in the structure
    domains = []
    for reference_index, reference_id in enumerate(reference_ids):
        # If this is not referable then there are no domains to mine
        if reference_id == NO_REFERABLE_FLAG:
            continue
        # If the reference was not found then there are no domains to mine
        if reference_id == NOT_FOUND_FLAG:
            continue
        # Get residue numbers of residues which belong to the current reference in th residue map
        residue_numbers = []
        for residue_index, residue_reference_index in enumerate(residue_reference_indices):
            if residue_reference_index == reference_index:
                residue_numbers.append(residue_reference_numbers[residue_index])
        # Set a function to find if any of the residue numbers is within a given range
        def in_range (start : int, end : int) -> bool:
            for residue_number in residue_numbers:
                if residue_number >= start and residue_number <= end:
                    return True
            return False
        # Get the reference data for the current reference uniprot id
        reference_data = references_data.get(reference_id, None)
        if not reference_data:
            raise SystemExit(reference_id + ' is not in the references data file')
        # Iterate over data domains
        for domain in reference_data['domains']:
            selection = domain['selection']
            # Iterate over the residue ranges in the selection field
            residue_ranges = selection.split(', ')
            for residue_range in residue_ranges:
                start, end = [ int(num) for num in residue_range.split('-') ]
                # If this range includes any resiudes in the residue map then the domain is included in the list
                if in_range(start, end):
                    domains.append(domain['name'])
                    break
    # Return the sequence matadata
    return {
        'sequences': sequences,
        'domains': domains,
        'cv19_variant': variant
    }
