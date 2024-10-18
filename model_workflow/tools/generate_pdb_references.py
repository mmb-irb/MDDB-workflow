import json
import urllib.request
from typing import List

from model_workflow.utils.auxiliar import load_json, save_json

# Generate a json file including all PDB references and return these references to the workflow for further usage
def generate_pdb_references (pdb_ids : List, pdb_references_file : 'File') -> List[dict]:
    # If we already have PDB references then load them
    previous_pdb_references = {}
    if pdb_references_file.exists:
        previous_content = load_json(pdb_references_file.path)
        for pdb_reference in previous_content:
            pdb_id = pdb_reference['id']
            previous_pdb_references[pdb_id] = pdb_reference
    # Mine PDB data for every PDB id
    pdb_references = []
    for pdb_id in pdb_ids:
        # Find if we already have this PDB id among hte previous references
        pdb_data = previous_pdb_references.get(pdb_id, None)
        if pdb_data:
            pdb_references.append(pdb_data)
            continue
        # Otherwise download and mine the PDB data
        pdb_data = mine_pdb_data(pdb_id)
        pdb_references.append(pdb_data)
    # Write references to a json file
    save_json(pdb_references, pdb_references_file.path, indent = 4)
    # Return the PDB references
    return pdb_references

# Mine PDB data
def mine_pdb_data (pdb_id : str) -> dict:
    # Request the MMB service to retrieve pdb data
    # request_url = f'http://mdb-login.bsc.es/api/pdb/{pdb_id}/entry'
    request_url = f'https://mmb.irbbarcelona.org/api/pdb/{pdb_id}/entry'
    try:
        with urllib.request.urlopen(request_url) as response:
            parsed_response = json.loads(response.read().decode("utf-8"))
    # If the accession is not found in the PDB then we can stop here
    except urllib.error.HTTPError as error:
        if error.code == 404:
            print(f' PDB code {pdb_id} not found')
            return None
        else:
            raise ValueError('Something went wrong with the PDB request: {request_url}')
    # Mine the PDB data we need
    pdb_data = { 'id': pdb_id }
    pdb_data['title'] = parsed_response['compound']
    pdb_data['class'] = parsed_response['header']
    pdb_data['authors'] = parsed_response['authors']
    pdb_data['date'] = parsed_response['ascDate']
    pdb_data['organism'] = parsed_response['sources']
    pdb_data['method'] = parsed_response['expType']
    pdb_data['resolution'] = parsed_response['resol']
    pdb_data['date'] = parsed_response['ascDate']
    chain_uniprots = {}
    for chain in parsed_response['chains']:
        chain_type = chain['type']
        # Mine protein chains only
        if chain_type != 'protein':
            continue
        # Get the chain letter
        chain_id = chain['_id']
        letter = chain_id[-1]
        # Get the uniprot id associated to this chain
        hit = chain['swpHit']
        uniprot_id = hit['idHit']
        chain_uniprots[letter] = uniprot_id
    pdb_data['chain_uniprots'] = chain_uniprots
    return pdb_data