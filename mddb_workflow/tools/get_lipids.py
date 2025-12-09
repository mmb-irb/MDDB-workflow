import requests
from functools import lru_cache
from mddb_workflow.utils.auxiliar import warn
from mddb_workflow.utils.type_hints import *


def generate_lipid_references(inchikeys: dict[str, 'InChIKeyData']) -> dict[str, dict]:
    """Add lipid-specific database information to InChIKeyData objects.

    This function queries SwissLipids and LIPID MAPS databases for each InChI key
    and adds the results directly to the InChIKeyData objects. It also performs
    quality checks on lipid classifications.

    Args:
        inchikeys: Dictionary mapping InChI keys to InChIKeyData objects (modified in-place).

    Returns:
        set: Set of InChI keys that were identified as lipids.

    """
    # Check internet connection
    try:
        is_in_swisslipids('test')
    except Exception as e:
        warn(f'There was a problem connecting to the SwissLipids database: {e}')
        return {}
    lipid_references = {}

    for inchikey, mol_data in inchikeys.items():
        # If we dont find it, we try without stereochemistry (only connection layer)
        SL_data = is_in_swisslipids(inchikey) or is_in_swisslipids(inchikey, only_first_layer=True)
        LM_data = is_in_LIPID_MAPS(inchikey) or is_in_LIPID_MAPS(inchikey, only_first_layer=True)

        # Add lipid database data to InChIKeyData
        if SL_data or LM_data:
            lipid_references[inchikey] = {'swisslipids': SL_data, 'lipidmaps': LM_data}

            # QUALITY CHECKS
            clasi = mol_data.classification
            # If the residue is a lipid, we check if it is classified as fatty/steroid
            if all('fatty' not in classes for classes in clasi) and \
                all('steroid' not in classes for classes in clasi):
                warn(f'The {mol_data.moltype} {mol_data.molname} is classified as {clasi}, '
                     f'but the InChIKey "{inchikey}" may be from a lipid.')
        else:
            # If the InChIKey is not in SwissLipids or LIPID MAPS, check classification
            if any('fatty' in classes for classes in mol_data.classification):
                warn(f'The InChIKey {inchikey} of {mol_data.molname} is '
                     f'classified as fatty but is not a lipid.\n'
                     f'Resindices: {str(mol_data.resindices)}')

    return lipid_references


@lru_cache(maxsize=None)
def is_in_LIPID_MAPS(inchikey : str, only_first_layer : bool = False) -> dict:
    """Search the InChI keys in LIPID MAPS."""
    if not inchikey: raise RuntimeError('Empty inchikey')
    headers = {'accept': 'json'}
    # https://www.lipidmaps.org/resources/rest
    # Output item = physchem, is the only one that returns data for the inchi key
    # for only the two first layers (main and atom connection)
    # To see InChIKey layers:
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


def get_swisslipids_info(entity_id) -> dict:
    """Get information about a SwissLipids entry."""
    headers = {'accept': 'json'}
    url = f"https://www.swisslipids.org/api/index.php/entity/{entity_id}"
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        return response.json()
    else:
        return False


@lru_cache(maxsize=None)
def is_in_swisslipids(inchikey, only_first_layer=False) -> dict:
    """Search the InChI keys in SwissLipids.
    Documentation: https://www.swisslipids.org/#/api.
    """
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
