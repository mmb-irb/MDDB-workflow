import os
import sys
import json
import urllib.request
from rdkit import Chem
from rdkit.Chem import AllChem
from mordred import Calculator, descriptors
from typing import List, Tuple, Optional, Union

from pathlib import Path
from model_workflow.tools.residues_library import residue_name_2_letter
from model_workflow.utils.auxiliar import load_json, save_json
from urllib.request import Request, urlopen
from model_workflow.utils.structures import Structure
from model_workflow.utils.auxiliar import load_json, save_json
from urllib.request import Request, urlopen
from urllib.parse import urlencode
from urllib.error import HTTPError, URLError
from urllib import request, parse
import json
import time

# Get the sequence and name of the chain in the structure and request the InterProScan 
def request_interpsocan (sequence : str, chain_name : str) -> str:
    # Set the request URL
    request_url = 'https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run'
    # Set the POST data
    data = urlencode({
        'email': 'daniel.beltran@irbbarcelona.org',
        'title': f'Chain {chain_name}',
        'sequence': f'>chain {chain_name}\n{sequence}'
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
def check_interproscan_result (jobid : str) -> str:
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
        else:
            raise ValueError('Something went wrong with the InterProScan results request: ' + request_url)
    return parsed_response

# Get the sequence and name of the chain in the structure and request the HMMER services
def request_hmmer (sequence : str, chain_name : str) -> str:
    # Set the request URL
    request_url = 'https://www.ebi.ac.uk/Tools/hmmer/search/phmmer'
    # Set the POST data
    data = urlencode({
        'seqdb': 'pdb',
        'seq': f'>chain {chain_name}\n{sequence}'
    }).encode()
    # Set the headers (they are needed in hmmer services)
    headers = {
    'Accept': 'application/json'
    }
    req = request.Request(request_url, data=data, headers=headers)
    parsed_response = None
    try:
        with request.urlopen(req) as response:
            parsed_response = json.loads(response.read().decode("utf-8"))
    except HTTPError as error:
        print(error.read().decode())
        if error.code == 404:
            print('Not found')
            parsed_response = None
        else:
            raise ValueError('Something went wrong with the request: ' + request_url)
    return parsed_response

# Check the status of the HMMER job
def check_hmmer_result (jobid : str) -> str:
    # Set the request URL
    request_url = f'https://www.ebi.ac.uk/Tools/hmmer/results/{jobid}.1'   
    parsed_response = None
    headers = {
    'Accept': 'application/json'
    }
    req = request.Request(request_url, headers=headers)
    try:
        with request.urlopen(req) as response:
            parsed_response = json.loads(response.read().decode("utf-8"))
    except HTTPError as error:
        print(error.read().decode())
        if error.code == 404:
            print('Not found')
            parsed_response = None
        else:
            raise ValueError('Something went wrong with the request: ' + request_url)
    return parsed_response

# Get the parsed chains from the structure
def get_parsed_chains (structure : 'Structure') -> list:
    parsed_chains = []
    chains = structure.chains
    # Iterate over the chains in the structure
    for chain in chains:
        # Get the name of the chain
        name = chain.name
        sequence = ''
        residue_indices = []
        # Iterate over the residues in the chain
        for residue in chain.residues:
            # Translate the residue letter to his equivalent in the aminoacids library
            letter = residue_name_2_letter(residue.name, 'aminoacids')
            sequence += letter
            residue_indices.append(residue.index)
        # Create a dictionary with the chain name, sequence and residue indices that be returned
        sequence_object = { 'name': name, 'sequence': sequence, 'residue_indices': residue_indices }
        parsed_chains.append(sequence_object)
    return parsed_chains

# Set the expected ligand data fields
CHAIN_DATA_FIELDS = set(['sequence', 'interproscan', 'hmmer'])

# Import the chains data from a file if exists
def import_chains(chains_references_filepath : str) -> dict:
    # Read the file
    imported_chains = load_json(chains_references_filepath)
    # Format data as the process expects to find it
    for imported_chain in imported_chains:
        for expected_field in CHAIN_DATA_FIELDS:
            if expected_field not in imported_chain:
                imported_chain[expected_field] = None
    return imported_chains

# Define the main function that will be called from the main script
# This function will get the parsed chains from the structure and request the InterProScan and HMMER services 
# to obtain the data for each chain
def chains_data (
        structure : Structure,
        chains_references_filepath : str,
        ) -> dict:
    
    print('-> Getting protein chains data') 

    # Obtain the parsed chains from the structure
    parsed_chains = get_parsed_chains(structure)
    protein_parsed_chains = []
    # Iterate over the parsed chains to get the protein sequences
    for chain_data in parsed_chains:
        sequence = chain_data['sequence']
        # We have to check if the sequence is a protein sequence
        # If the sequence has any letter different from 'X' then it is a protein sequence
        if next((letter for letter in sequence if letter != 'X'), None):
            chain_data['match'] = { 'ref': None, 'map': None, 'score': 0 }
            protein_parsed_chains.append(chain_data)
    # If there are no protein sequences then there is no need to analyce anything
    if len(protein_parsed_chains) == 0:
        print(' There are no protein sequences')
        return protein_parsed_chains
    
    # Save data from all chains to be saved in a file
    json_chains_data = []
    # Load the chains file if exists already
    if os.path.exists(chains_references_filepath):
        json_chains_data += import_chains(chains_references_filepath)

    # Iterate over the protein parsed chains to request the InterProScan and HMMER services
    chains_data_list = []
    interproscan_jobids = []
    hmmer_jobids = []
    for chain in protein_parsed_chains:
        # Get the chain name and sequence
        chain_name = chain['name']
        sequence = chain['sequence']

        # Check if the chain data already exists in the chains file
        chain_exists = False
        for chain_data in json_chains_data:
            # If the sequence is the same then the chain data already exists
            # We can skip the request to the services and use the data from the file
            if chain['sequence'] == chain_data['sequence']:
                chains_data_list.append(chain_data)
                chain_exists = True
                break
        # If the chain data already exists then we need to remove it from the protein parsed chains list   
        if chain_exists:
            protein_parsed_chains.remove(chain)
            continue
        
        # If the chain data does not exist then we need to request the services
        # Request the InterProScan and HMMER services
        interproscan_jobid = request_interpsocan(sequence, chain_name)
        hmmer_jobid = request_hmmer(sequence, chain_name)
        # Add the jobids to the corresponding lists
        interproscan_jobids.append(interproscan_jobid)
        hmmer_jobids.append(hmmer_jobid)

    # Save the final results to each list
    interproscan_results = []
    hmmer_results = []
    # Iterate over the jobids to check the status and get the results
    # If the status is 'FINISHED' then we can get the results and eliminate the jobid from the list 
    # until there are no more jobids in either list
    while len(interproscan_jobids) >= 1 or len(hmmer_jobids) >= 1:
        time.sleep(5)  # Wait for 5 seconds
        for interproscan_jobid in interproscan_jobids:
            status = check_interproscan_status(interproscan_jobid)
            # We only know four possible status for interproscan jobs, but it could be more
            if status == 'RUNNING' or status == 'PENDING' or status == 'QUEUED':
                continue
            else:
                # If the status is 'FINISHED' then we can get the results and delete the jobid from the list
                if status == 'FINISHED':
                    interproscan_result = check_interproscan_result(interproscan_jobid)
                    interproscan_results.append(interproscan_result)
                    interproscan_jobids.remove(interproscan_jobid)
                # If the status is something that we don´t know then we raise an error in order to solucionate this problem
                else:
                    raise ValueError('Something went wrong with the InterProScan job: ' + interproscan_jobid)
        
        # AGUS: HMMER devuelve directamente el resultado, no hay ningún status
        # AGUS: La única forma que se me ocurrió para comprobar si estaba acabado el job es con el campo 'results'
        # AGUS: pero si por lo que fuera fallara el sistema de EBI, no sé que devuelve (no me ha pasado todavía). Podría pasar en un futuro
        for hmmer_jobid in hmmer_jobids:
            if 'results' in hmmer_jobid:
                hmmer_result = check_hmmer_result(hmmer_jobid['results']['uuid'])
                hmmer_results.append(hmmer_result)
                hmmer_jobids.remove(hmmer_jobid)

    # AGUS: según el orden en el que acaben los jobs, los resultados se pueden mezclar
    # AGUS: entre las listas de interproscan_results y hmmer_results,  SOLUCIONAR ESTO
    for i in range(len(protein_parsed_chains)):
        chain_result = {
            'sequence': protein_parsed_chains[i]['sequence'],
            'interproscan': interproscan_results[i],
            'hmmer': hmmer_results[i]
        }
        chains_data_list.append(chain_result)
    
    with open(chains_references_filepath, 'w') as file:
        json.dump(chains_data_list, file, indent=4)
    
    return print('-> Protein chains data obtained')
