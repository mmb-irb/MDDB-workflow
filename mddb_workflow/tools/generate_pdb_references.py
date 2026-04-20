import os
import json
import urllib.request

from mddb_workflow.utils.auxiliar import RemoteServiceError, load_json, save_json, request_pdb_data, round_to_thousandths, warn
from mddb_workflow.utils.structures import Structure
from mddb_workflow.utils.gmx_spells import run_gromacs
from mddb_workflow.utils.type_hints import *

from mddb_workflow.tools.xvg_parse import xvg_parse
from mddb_workflow.tools.generate_map import normalize_protein_sequence, get_uniprot_reference, align


def prepare_pdb_references (pdb_ids : list[str], pdb_references_file : 'File'):
    """Prepare the PDB references json file to be uploaded to the database."""
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
    # DANI: There is no need to return PDB references since it is not used in the workflow
    #return pdb_references

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
            exptl { method }
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
    # pdbx_refine_id not on every PDB like in 1FI3
    #pdb_data['method'] = parsed_response['refine'][0]['pdbx_refine_id']
    pdb_data['method'] = parsed_response['exptl'][0]['method']
    # ls_d_res_high not on every PDB like in 1FI3
    # pdb_data['resolution'] = parsed_response['refine'][0]['ls_d_res_high']
    chain_uniprots = {}
    organisms = []
    for polymer in parsed_response['polymer_entities']:
        identifier = polymer['rcsb_polymer_entity_container_identifiers']
        # Get all chains with this configuration
        chains = identifier['asym_ids']
        # Get the organisms
        organism_entries = polymer['rcsb_entity_source_organism']
        if organism_entries != None:
            organisms += [ organism['scientific_name'] for organism in organism_entries ]
        # Get the uniprot
        uniprots = identifier.get('uniprot_ids', None)
        if not uniprots: continue
        if len(uniprots) > 1:
            warn(f'PDB {pdb_id} has multiple uniprots: {uniprots}. Using only the first one')
        uniprot_id = uniprots[0]
        # Get all chains with this configuration
        chains = identifier['asym_ids']
        for chain in chains:
            chain_uniprots[chain] = uniprot_id
    pdb_data['organisms'] = list(set(organisms))
    pdb_data['chain_uniprots'] = chain_uniprots
    # Download the structure and calculate the solvent accessible surface in the PDB structure as reference
    # Do it for protein chains only and get also their aminoacid sequences
    chain_sas, chain_seq = calculate_pdb_chain_sas(pdb_id)
    # Find the UniProt reference numeration for every residue in the sequence
    chain_residues_uniprot_numeration = {}
    # Iterate protein chains
    for chain, uniprot_id in chain_uniprots.items():
        # Get the chain sequence
        sequence = chain_seq[chain]
        # Get the uniprot reference, including the reference sequence
        uniprot_reference = get_uniprot_reference(uniprot_id)
        reference_sequence = uniprot_reference['sequence']
        # Align the PDB sequence against the uniprot reference sequence to know the reference residue numeration
        residue_numeration = None
        align_match = align(reference_sequence, sequence) if sequence else None
        if not align_match:
            warn(f'PDB {pdb_id} has a polymer tagged with the UniProt {uniprot_id} but its reference sequence does not match.')
            chain_residues_uniprot_numeration[chain] = 'Missmatched'
            continue
        residue_numeration, score = align_match
        chain_residues_uniprot_numeration[chain] = residue_numeration
    # Finally assign these values to the reference data
    pdb_data['chain_seq'] = chain_seq
    pdb_data['chain_resnum'] = chain_residues_uniprot_numeration
    pdb_data['chain_sas'] = chain_sas
    # Sort all dictionaries
    # Otherwise the loader will complain about a change every time since the order may change
    for key, value in pdb_data.items():
        if type(value) == dict:
            pdb_data[key] = dict(sorted(value.items()))
    # Set a version number so new references can replace old references
    pdb_data['version'] = '0.0.1'
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

# Set a residual output filepath
RESIDUAL_AREA_FILENAME = '.residual_area.xvg'
def calculate_pdb_chain_sas (pdb_id : str) -> dict:
    """Calculate the solvent accessible surface in the PDB structure as reference."""
    # Download the structure
    # Download the mmcif instead of the PDB to avoid having problems
    # Some PDB entries may be not available in PDB format anymore
    auxiliar_mmcif_filepath = f'.{pdb_id}.cif'
    download_mmcif(pdb_id, auxiliar_mmcif_filepath)
    # Load the structure
    structure = Structure.from_mmcif_file(auxiliar_mmcif_filepath)
    # Now export it to PDB so gromacs can read it
    auxiliar_pdb_filepath = f'.{pdb_id}.pdb'
    structure.generate_pdb_file(auxiliar_pdb_filepath)
    # Target only protein chains
    protein_selection = structure.select_protein()
    protein_chains = structure.get_selection_chains(protein_selection)
    # Keep the SAS of every chain
    pdb_chains_sas = {}
    # Iterate PDB chains
    for chain in protein_chains:
        # Run a SAS analysis only for this specific chain
        # Convert the chain selection to a ndx file
        chain_selection = chain.get_selection()
        sasa_selection_name = f'chain_{chain.name}'
        sasa_ndx_selection = chain_selection.to_ndx(sasa_selection_name)
        ndx_filename = f'.indices.ndx'
        with open(ndx_filename, 'w') as file:
            file.write(sasa_ndx_selection)
        # Set the main output filepath
        current_chain_sasa = f'sasa_{chain.name}.xvg'
        # If not defined, then by default it creates a 'area.xvg' file where we are running the workflow
        run_gromacs(f'sasa -s {auxiliar_pdb_filepath} -oa {current_chain_sasa} -o {RESIDUAL_AREA_FILENAME}\
            -n {ndx_filename} -surface {sasa_selection_name}', expected_output_filepath = current_chain_sasa)
        # Mine the sasa results (.xvg file)
        # Areas from excluded atoms are not recorded in the xvg file
        sasa = xvg_parse(current_chain_sasa, ['n', 'area', 'sd'])
        # Restructure data by adding all atoms sas per residue
        atom_numbers = sasa['n']
        atom_areas = sasa['area']
        sas_per_residues = [0.0] * len(structure.residues)
        for atom_number, atom_area in zip(atom_numbers, atom_areas):
            atom_index = atom_number - 1
            atom = structure.atoms[atom_index]
            residue_index = atom.residue_index
            sas_per_residues[residue_index] += atom_area
        # Round to thousands
        sas_per_residues = [ round_to_thousandths(sas) for sas in sas_per_residues ]
        pdb_chains_sas[chain.name] = sas_per_residues
        # Remove files which are no longer required
        os.remove(ndx_filename)
        os.remove(current_chain_sasa)
        os.remove(RESIDUAL_AREA_FILENAME)
    # Keep also the chain sequence
    # Note that the sequence coming from the GraphQL is the reference sequence, including also the endings
    # And only god knows how to ask GraphQL to get the actual present sequence in the PDB structure
    # In order to get the structure we read it directly from the parsed structure
    pdb_chains_seq = {}
    for chain in protein_chains:
        pdb_chains_seq[chain.name] = chain.get_sequence()
    # Remove the auxiliar file
    os.remove(auxiliar_pdb_filepath)
    os.remove(auxiliar_mmcif_filepath)
    return pdb_chains_sas, pdb_chains_seq
    
def download_mmcif (pdb_id : str, output_filepath : str):
    """Download a mmCIF file corresponding to a PDB entry."""
    request_url = f'https://files.rcsb.org/download/{pdb_id}.cif'
    try:
        parsed_response = None
        with urllib.request.urlopen(request_url) as response:
            parsed_response = response.read().decode("utf-8")
        with open(output_filepath, 'w') as file:
            file.write(parsed_response)
    except Exception as error:
        raise Exception(f'Something went wrong when downloading {pdb_id} structure: {request_url} with error: {error}')