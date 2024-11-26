import json
import urllib.request
from typing import List

from model_workflow.utils.auxiliar import load_json, save_json, RemoteServiceError

# Generate a json file including all PDB references and return these references to the workflow for further usage
def generate_pdb_references (pdb_ids : List, pdb_references_file : 'File') -> List[dict]:
    print('-> Generating PDB references file')
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

# Set service URLs to be requested
pdb_data_services = {
    'IRB': f'https://mmb.irbbarcelona.org',
    'BSC': f'https://mdb-login.bsc.es'
}
# Get PDB data from a remote service
def get_pdb_data (pdb_id : str, service = 'IRB') -> dict:
    # Set the request URL
    service_url = pdb_data_services[service]
    request_url = f'{service_url}/api/pdb/{pdb_id}/entry'
    # Send the request and save the response
    print(f'Requesting {request_url} (...)', end='\r')
    parsed_response = None
    try:
        with urllib.request.urlopen(request_url) as response:
            parsed_response = json.loads(response.read().decode("utf-8"))
    # Handle HTTP errors
    except urllib.error.HTTPError as error:
        print(f'There was a problem when requesting {request_url}')
        print(f'  Error code: {error.code}')
        # If the PDB id is not found then we can stop here
        if error.code == 404:
            print(f'  PDB code {pdb_id} not found')
            return None
        # If the API is not responding try another service
        elif error.code == 500 or error.code == 502 or error.code == 503:
            print(f'  {service} API to retrieve PDB data may be out of service')
            # Before we surrender we try with the other available service
            if service == 'IRB':
                print('  Retrying with a different service')
                return get_pdb_data(pdb_id, service='BSC')
            # If we already tried with the other service then surrender
            raise RemoteServiceError('All APIs to retrieve PDB data may be out of service')
        # If the error is not known then stop here
        else:
            raise RemoteServiceError(f'Something went wrong with the PDB data request')
    # Handle URL errors
    except urllib.error.URLError as error:
        print(f'There was a problem when requesting {request_url}')
        print(f'  Error reason: {error.reason}')
        # These errors are not our fault, but the service is unsatable
        print(f'  {service} API to retrieve PDB data may be out of service')
        # Try with a different service
        if service == 'IRB':
            print('  Retrying with a different service')
            return get_pdb_data(pdb_id, service='BSC')
        raise RemoteServiceError(f'Something went wrong with the PDB data request')
    # Return the response
    print(f'Successfully requested {request_url}')
    return parsed_response

# Mine PDB data
def mine_pdb_data (pdb_id : str) -> dict:
    # Download the PDB data
    parsed_response = get_pdb_data(pdb_id)
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
        if hit == None:
            continue
        uniprot_id = hit['idHit']
        chain_uniprots[letter] = uniprot_id
    pdb_data['chain_uniprots'] = chain_uniprots
    return pdb_data