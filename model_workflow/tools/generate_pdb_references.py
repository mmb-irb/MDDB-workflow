import json
import urllib.request
import urllib.parse
from model_workflow.utils.type_hints import *

from model_workflow.utils.auxiliar import RemoteServiceError, load_json, save_json, request_pdb_data

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
        pdb_reference = previous_pdb_references.get(pdb_id, None)
        if pdb_reference:
            pdb_references.append(pdb_reference)
            continue
        # Otherwise download and mine the PDB data
        pdb_reference = get_pdb_reference(pdb_id)
        pdb_references.append(pdb_reference)
    # Write references to a json file
    save_json(pdb_references, pdb_references_file.path, indent = 4)
    # Return the PDB references
    return pdb_references

# Download PDB data from the PDB API
def get_pdb_reference (pdb_id : str) -> dict:
    # Set the request query
    query = '''query ($id: String!) {
        entry(entry_id: $id) {
            rcsb_id
            struct { title }
            struct_keywords { pdbx_keywords }
            refine { pdbx_refine_id ls_d_res_high }
            rcsb_accession_info { initial_release_date }
            audit_author { name }
            polymer_entities {
                rcsb_polymer_entity_container_identifiers { asym_ids uniprot_ids }
                rcsb_entity_source_organism { scientific_name }
            }
        }
    }'''
    # Request PDB data
    parsed_response = request_pdb_data(pdb_id, query)
    # Mine data
    pdb_data = {}
    pdb_data['id'] = parsed_response['rcsb_id']
    pdb_data['title'] = parsed_response['struct']['title']
    pdb_data['class'] = parsed_response['struct_keywords']['pdbx_keywords']
    pdb_data['authors'] = [ author['name'] for author in parsed_response['audit_author'] ]
    pdb_data['date'] = parsed_response['rcsb_accession_info']['initial_release_date']
    pdb_data['method'] = parsed_response['refine'][0]['pdbx_refine_id']
    pdb_data['resolution'] = parsed_response['refine'][0]['ls_d_res_high']
    chain_uniprots = {}
    organisms = []
    for polymer in parsed_response['polymer_entities']:
        identifier = polymer['rcsb_polymer_entity_container_identifiers']
        # Get the organisms
        organisms += [ organism['scientific_name'] for organism in polymer['rcsb_entity_source_organism'] ]
        # Get the uniprot
        uniprots = identifier.get('uniprot_ids', None)
        if not uniprots: continue
        if len(uniprots) > 1: raise ValueError('Multiple PDB ids in the same chain?')
        uniprot_id = uniprots[0]
        chains = identifier['asym_ids']
        for chain in chains:
            chain_uniprots[chain] = uniprot_id
    pdb_data['chain_uniprots'] = chain_uniprots
    pdb_data['organisms'] = organisms
    return pdb_data

# Set service URLs to be requested
pdb_data_services = {
    'IRB': f'https://mmb.irbbarcelona.org',
    'BSC': f'https://mdb-login.bsc.es'
}
# Download PDB data from a remote service
# DEPRECATED: our custom PDB API may return UniProt ids not matching the ones in the PDB
# e.g. 1AK4 and 4Z80
def DEPRECATED_download_pdb_data (pdb_id : str, service = 'IRB') -> dict:
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
            print(f'  PDB id {pdb_id} not found')
            return None
        # If the API is not responding try another service
        elif error.code == 500 or error.code == 502 or error.code == 503:
            print(f'  {service} API to retrieve PDB data may be out of service')
            # Before we surrender we try with the other available service
            if service == 'IRB':
                print('  Retrying with a different service')
                return DEPRECATED_download_pdb_data(pdb_id, service='BSC')
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
            return DEPRECATED_download_pdb_data(pdb_id, service='BSC')
        raise RemoteServiceError(f'Something went wrong with the PDB data request')
    # Return the response
    print(f'Successfully requested {request_url}')
    return parsed_response