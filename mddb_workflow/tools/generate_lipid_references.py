import requests
from functools import lru_cache
from mddb_workflow.utils.auxiliar import save_json
from mddb_workflow.utils.auxiliar import warn
from mddb_workflow.utils.type_hints import *

def generate_lipid_references(inchikeys: dict,
                              output_filepath : str,
                              ) -> list[dict]:
    """Generate the lipid references.

    Returns:
        dict: A list of dictionaries containing lipid references and
        the residues indices. For example:
            [{'name': 'CHL1',
              'residue_indices': [935, 936, 937, ...],
              'fragments': [],
              'match': {'ref': {'inchikey': 'HVYWMOMLDIMFJA-DPAQBDIFSA-N'}}}, ...]
    """
    # Patch case where there no internet
    try:
        # This would return a ConnectionError
        is_in_swisslipids('test')
    except Exception as e:
        # Then we map the lipids/membrane
        warn(f'There was a problem connecting to the SwissLipids database: {e}')
        return None

    lipid_references = []
    lipid_map = []
    for inchikey, res_data in inchikeys.items():
        SL_data = is_in_swisslipids(inchikey)
        LM_data = is_in_LIPID_MAPS(inchikey)
        # If we dont find it, we try without stereochemistry
        if not SL_data:
            SL_data = is_in_swisslipids(inchikey, only_first_layer=True)
        if not LM_data:
            LM_data = is_in_LIPID_MAPS(inchikey, only_first_layer=True)

        # We don't use lipid data for now, if we have it it is present in LIPID MAPS
        if SL_data or LM_data:
            lipid_references.append({'inchikey': inchikey,
                                     'inchi': res_data['inchi'],
                                     'swisslipids': SL_data,
                                     'lipidmaps': LM_data,
                                     })
            # Format needed for generate_residue_mapping
            lipid_map.append({ 
                'name': list(res_data['resname'])[0], 
                'residue_indices': list(map(int, res_data['resindices'])),
                'fragments': res_data['fragments'], 
                'match': { 
                    'ref': { 'inchikey': inchikey } } 
                    })
 
            # QUALITY CHECKS
            cls = res_data['classification']
            # If the residue is a lipid, we check if it is classified as fatty/steroid
            if all('fatty' not in classes for classes in cls) and \
                all('steroid' not in classes for classes in cls):
                warn(f'The residue {str(res_data["resname"])} is classified as {cls}, but the InChIKey "{inchikey}" is a lipid.')

        else:
            # If the InChIKey is not in SwissLipids or LIPID MAPS, we check if it is classified as fatty
            if any('fatty' in classes for classes in res_data['classification']):
                warn(f'The InChIKey {inchikey} of {str(res_data["resname"])} is '
                     f'classified as fatty but is not a lipid.\n'
                     f'Resindices: {str(res_data["resindices"])}')

    save_json(lipid_references, output_filepath)
    return lipid_map

@lru_cache(maxsize=None)
def is_in_LIPID_MAPS(inchikey, only_first_layer=False) -> dict:
    """Search the InChi keys in LIPID MAPS"""
    headers = {'accept': 'json'}
    # https://www.lipidmaps.org/resources/rest
    # Output item = physchem, is the only one that returns data for the inchi key
    # for only the two first layers (main and atom connection)
    # To see InChiKey layers: 
    # https://www.inchi-trust.org/about-the-inchi-standard/
    # Or https://www.rhea-db.org/help/inchi-inchikey#What_is_an_InChIKey_
    key = inchikey[:14] if only_first_layer else inchikey
    url = f"https://www.lipidmaps.org/rest/compound/inchi_key/{key}/all"
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        js = response.json()
        if js != []:
           return js      
        else:
            return False  
    else:
        print(f"Error for {inchikey}: {response.status_code}")

@lru_cache(maxsize=None)
def is_in_swisslipids(inchikey, only_first_layer=False, 
                       protonation=True) -> dict:
    """Search the InChi keys in SwissLipids. Documentation: https://www.swisslipids.org/#/api"""
    key = inchikey[:14] if only_first_layer else inchikey
    headers = {'accept': 'json'}
    url = f"https://www.swisslipids.org/api/index.php/search?term={key}"
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        data = response.json()[0]
        detailed_data = get_swisslipids_info(data['entity_id'])
        data['synonyms'] = detailed_data.get('synonyms', [])
        return data
    else:
        return False

def get_swisslipids_info(entity_id) -> dict:
    """Get information about a SwissLipids entry."""
    headers = {'accept': 'json'}
    url = f"https://www.swisslipids.org/api/index.php/entity/{entity_id}"
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        return response.json()
    else:
        return False
