
import json
import time

from urllib.request import urlopen
from urllib.parse import urlencode
from urllib.error import HTTPError
from urllib import request

from model_workflow.tools.residues_library import residue_name_2_letter
from model_workflow.utils.auxiliar import load_json, save_json, RemoteServiceError
from model_workflow.utils.type_hints import *

# Set analysis version
CHAINS_VERSION = '0.1'

# Get the sequence and name of the chain in the structure and request the InterProScan 
def request_interpsocan (sequence : str) -> str:
    # Set the request URL
    request_url = 'https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run'
    # Set the POST data
    data = urlencode({
        'email': 'daniel.beltran@irbbarcelona.org',
        'title': f'Chain X',
        'sequence': f'>chain X\n{sequence}'
    }).encode()
    parsed_response = None
    try:
        with urlopen(request_url, data=data) as response:
            parsed_response = response.read().decode("utf-8")
    except HTTPError as error:
        print(error.read().decode())
        if error.code == 404:
            print(f' Not found')
            return None
        else:
            raise ValueError('Something went wrong with the InterProScan request: ' + request_url)
    return parsed_response

# Check the status of the InterProScan job
def check_interproscan_status (jobid : str) -> str:
    # Set the request URL
    request_url = f'https://www.ebi.ac.uk/Tools/services/rest/iprscan5/status/{jobid}'   
    parsed_response = None
    try:
        with urlopen(request_url) as response:
            parsed_response = response.read().decode("utf-8")
    except HTTPError as error:
        print(error.read().decode())
        if error.code == 404:
            print(f' Not found')
            return None
        else:
            raise ValueError('Something went wrong with the InterProScan status request: ' + request_url)
    return parsed_response

# Obtain the result of the InterProScan job in json format
def check_interproscan_result (jobid : str) -> dict:
    # Set the request URL
    request_url = f'https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/{jobid}/json'   
    parsed_response = None
    try:
        with urlopen(request_url) as response:
            parsed_response = json.loads(response.read().decode("utf-8"))
    except HTTPError as error:
        print(error.read().decode())
        if error.code == 404:
            print(f' Not found')
            return None
        elif error.code == 503:
            raise SystemExit('InterProScan Service unavailable. Please try again later.')
        else:
            raise ValueError('Something went wrong with the InterProScan results request: ' + request_url)
    return parsed_response

# Check if the required sequence is already in the MDDB database
def check_sequence_in_mddb (sequence : str, database_url : str) -> dict:
    # Request PubChem
    request_url = f'{database_url}rest/v1/references/chains/{sequence}'
    parsed_response = None
    try:
        with urlopen(request_url) as response:
            #parsed_response = json.loads(response.read().decode("windows-1252"))
            parsed_response = json.loads(response.read().decode("utf-8", errors='ignore'))
    # If the accession is not found in PubChem then the id is not valid
    # This may happen with pubchem ids of non-discrete compounds (e.g. 483927498)
    except HTTPError as error:
        if error.code == 404:
            return None
        elif error.code == 503:
            print('MDDB Service unavailable. Please try again later.')
            return None
        print('Error when requesting ', request_url)
        return None
    return parsed_response


# Get the parsed chains from the structure
def get_protein_parsed_chains (structure : 'Structure') -> list:
    parsed_chains = []
    chains = structure.chains
    # Iterate over the chains in the structure
    for chain in chains:
        # Get the name of the chain
        name = chain.name
        sequence = ''
        # Iterate over the residues in the chain
        for residue in chain.residues:
            # Translate the residue letter to his equivalent in the aminoacids library
            letter = residue_name_2_letter(residue.name, 'aminoacids')
            sequence += letter
        # If all residues are 'X' then it means this is not a protein
        if all(letter == 'X' for letter in sequence):
            continue
        # Create a dictionary with the chain name, sequence and residue indices that be returned
        sequence_object = { 'name': name, 'sequence': sequence }
        parsed_chains.append(sequence_object)
    return parsed_chains

# Set the expected ligand data fields
CHAIN_DATA_FIELDS = set(['sequence', 'interproscan'])

# Import the chains data from a file if exists
def import_chains (chains_references_file : 'File') -> dict:
    # Read the file
    imported_chains = load_json(chains_references_file.path)
    # Format data as the process expects to find it
    for imported_chain in imported_chains:
        for expected_field in CHAIN_DATA_FIELDS:
            if expected_field not in imported_chain:
                imported_chain[expected_field] = None
    return imported_chains

# Define the main function that will be called from the main script
# This function will get the parsed chains from the structure and request the InterProScan service
# to obtain the data for each chain
def generate_chain_references (
    structure : 'Structure',
    chains_references_file : 'File',
    database_url : str,
) -> dict:
    
    print('-> Getting protein chains data')

    # Obtain the protein parsed chains from the structure
    protein_parsed_chains = get_protein_parsed_chains(structure)

    # Get unique sequences
    protein_sequences = set([ chain['sequence'] for chain in protein_parsed_chains ])

    print(f' Found {len(protein_parsed_chains)} protein chains with {len(protein_sequences)} unique sequences')
        
    # Save data from all chains to be saved in a file
    chains_data = []
    # Load the chains file if exists already
    if chains_references_file.exists:
        chains_data += import_chains(chains_references_file)

    # Save the jobids of every call to InterProScan
    interproscan_jobids = {}

    # Iterate protein sequences
    for sequence in protein_sequences:
        # Check if the chain data already exists in the chains file
        chain_data = next((data for data in chains_data if data['sequence'] == sequence), None)
        # If we have no previous chain data then check if the sequence is already in the MDDB database
        if chain_data == None:
            chain_data = check_sequence_in_mddb(sequence, database_url)
            if chain_data is not None:
                chains_data.append(chain_data)
                # Save the chains data at this point
                # This may seem redundant since data will be not loaded further in the database
                # However, saving the chain in the local backup file is useufl to further run without internet connection
                save_json(chains_data, chains_references_file.path)
        # If we still have no chain data then create a new chain data dict
        # Set an object with the results of every call to InterProScan
        if chain_data == None:
            chain_data = {
                'sequence': sequence,
                'interproscan': None
            }
            chains_data.append(chain_data)
        # If chain data is missing any analysis then send a job
        # Request the InterProScan service
        # Keep the returned job ids to check the status and get the results later
        if chain_data['interproscan'] == None:
            interproscan_jobid = request_interpsocan(sequence)
            interproscan_jobids[sequence] = interproscan_jobid

    # Get the pending interpsocan jobids 
    pending_jobids = list(interproscan_jobids.values())

    # If we already have the results of all the chains then we can skip the next steps
    if len(pending_jobids) == 0:
        print(' All reference chains are already in the backup file')
        return chains_data

    # Iterate over the jobids to check the status and get the results
    # If the status is 'FINISHED' then we can get the results and eliminate the jobid from the list 
    # until there are no more jobids in either list
    while len(pending_jobids) >= 1:
        time.sleep(3)  # Wait for 3 seconds
        print(f' We are still waiting for {len(pending_jobids)} jobs to finish', end='\r')
        for sequence, interproscan_jobid in interproscan_jobids.items():
            # If the jobid is already processed then skip it
            if interproscan_jobid not in pending_jobids:
                continue
            status = check_interproscan_status(interproscan_jobid)
            # We only know four possible status for interproscan jobs, but it could be more
            if status == 'RUNNING' or status == 'PENDING' or status == 'QUEUED':
                continue
            # If the status is something that we don´t know then we raise an error in order to solucionate this problem
            if status != 'FINISHED':
                raise ValueError('Something went wrong with the InterProScan job: ' + interproscan_jobid)
            # Retrive the results from InterProScan
            interproscan_result = check_interproscan_result(interproscan_jobid)
            # Get corresponding chain data and add the InterProScan results
            chain_data = next(data for data in chains_data if data['sequence'] == sequence)
            chain_data['version'] = CHAINS_VERSION
            chain_data['interproscan'] = interproscan_result
            # Remove version and pathways so Mongo don't get confused when they change
            del chain_data['interproscan']['interproscan-version']
            # RUBEN: creo que results siempre tiene un solo elemento, pero por si acaso iteramos
            for result in chain_data['interproscan']['results']:
                for match in result['matches']:
                    if match['signature']['entry'] is not None:
                        del match['signature']['entry']['pathwayXRefs']
            # Remove the jobid from the queue list
            pending_jobids.remove(interproscan_jobid)
            # Save the result
            save_json(chains_data, chains_references_file.path)
    
    print(' Protein chains data obtained              ')
