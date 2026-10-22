import requests
from functools import lru_cache
from mddb_workflow.utils.auxiliar import warn, retry_request, RequestError, RateLimitError, get_retry_after, InputError
from mddb_workflow.utils.constants import LIPID_REFERENCE_VERSION, LIPIDS_MATCH_FLAG
from mddb_workflow.tools.get_inchi_keys import get_database_inchikey_reference
from mddb_workflow.utils.selections import Selection
from mddb_workflow.utils.type_hints import *


def generate_lipid_references(
    inchikeys: dict[str, 'InChIKeyData'],
    cache: 'Cache',
    input_forced_lipids: Optional[list[str]] = None,
    mercy: list[str] = [],
    database: Optional['Database'] = None,
) -> dict[str, dict]:
    """Add lipid-specific database information to InChIKeyData objects.

    This function queries SwissLipids and LIPID MAPS databases for each InChI key
    and adds the results directly to the InChIKeyData objects. It also performs
    quality checks on lipid classifications. InChI keys already in the MDDB
    database reuse its lipid data instead.

    Residues with names in the forced lipids are treated as lipids although they are not in the databases.
    Residues classified as fatty which are not in the databases nor forced make the workflow fail,
    unless the lipids mercy flag is passed.

    Args:
        inchikeys: Dictionary mapping InChI keys to InChIKeyData objects (modified in-place).
        cache: Cache to store the MDDB database requests.
        input_forced_lipids: Residue names to be treated as lipids.
        mercy: Failures which do not kill the process.
        database: MDDB database handler to reuse existing references.

    Returns:
        dict: Lipid references by InChI key.

    """
    lipid_references = {}
    forced_lipids = set(input_forced_lipids or [])
    not_matched_fatty = []
    # Check internet connection only once we need to request the lipid databases
    is_connected = None

    for inchikey, mol_data in inchikeys.items():
        clasi = mol_data.classification
        if len(clasi) == 1 and 'ion' in clasi:
            # Skip ions
            continue
        # Residues forced by the user are lipids with no database references
        if set(mol_data.resnames) & forced_lipids:
            lipid_references[inchikey] = {
                'swisslipids': None,
                'lipidmaps': None,
                'forced': True,
                'version': LIPID_REFERENCE_VERSION
            }
            continue
        # Check if the lipid data is already in the MDDB database
        database_reference = get_database_inchikey_reference(inchikey, cache, database)
        if database_reference:
            SL_data = database_reference.get('swisslipids')
            LM_data = database_reference.get('lipidmaps')
        else:
            if is_connected is None:
                try:
                    is_in_swisslipids('test')
                    is_connected = True
                except Exception as e:
                    warn(f'There was a problem connecting to the SwissLipids database: {e}')
                    is_connected = False
            if not is_connected: continue
            # If we dont find it, we try without stereochemistry (only connection layer)
            SL_data = is_in_swisslipids(inchikey) or is_in_swisslipids(inchikey, only_first_layer=True)
            LM_data = is_in_LIPID_MAPS(inchikey) or is_in_LIPID_MAPS(inchikey, only_first_layer=True)

        # Add lipid database data to InChIKeyData
        if SL_data or LM_data:
            lipid_references[inchikey] = {
                'swisslipids': SL_data,
                'lipidmaps': LM_data,
                'version': LIPID_REFERENCE_VERSION
            }

            # QUALITY CHECKS
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
                not_matched_fatty.append((mol_data.molname, inchikey, mol_data.resindices))

    if not_matched_fatty and LIPIDS_MATCH_FLAG not in mercy:
        raise InputError('Some residues are classified as fatty but they are not found in SwissLipids or LIPID MAPS:\n\t'
            + '\n\t'.join(f'{name} - {key} (residues {resindices})' for name, key, resindices in not_matched_fatty)
            + '\n\tAdd their residue names to the "forced_lipids" input field to treat them as lipids'
            + f' or use the flag "-m {LIPIDS_MATCH_FLAG}" to bypass this check.')

    return lipid_references


def get_lipids_selection(
    structure: 'Structure',
    inchikey_map: list[dict],
    lipid_references: dict[str, dict],
    cg_selection: Optional['Selection'] = None,
) -> 'Selection':
    """Select lipids using lipid references (InChI keys found in SwissLipids/LIPID MAPS) for all-atom regions.
    Coarse grained regions have no InChI keys, so here we rely on the residue classification.
    """
    aa_lipids_selection = Selection()
    for inchikey_data in inchikey_map:
        if inchikey_data['generated_inchikey'] in lipid_references:
            aa_lipids_selection += structure.select_residue_indices(inchikey_data['residue_indices'])
    cg_lipids_selection = structure.select_lipids() & cg_selection
    return (aa_lipids_selection - cg_selection) + cg_lipids_selection


def check_lipid_response(response: 'requests.Response', service: str, key: str):
    """Raise an error if the lipid database response is not successful."""
    if response.status_code == 200: return
    if response.status_code == 429:
        raise RateLimitError(f'Too many requests to {service} for {key}',
            retry_after=get_retry_after(response.headers))
    raise RequestError(f'Something went wrong with the {service} request for {key} '
        f'(error {response.status_code}): {response.reason}')


@retry_request
@lru_cache(maxsize=None)
def is_in_LIPID_MAPS(inchikey: str, only_first_layer: bool = False) -> dict:
    """Search the InChI keys in LIPID MAPS."""
    if not inchikey: raise RuntimeError('Empty inchikey')
    # https://www.lipidmaps.org/resources/rest
    # To see InChIKey layers:
    # https://www.inchi-trust.org/about-the-inchi-standard/
    # Or https://www.rhea-db.org/help/inchi-inchikey#What_is_an_InChIKey_
    # Output item = all only matches full InChI keys
    # Output item = physchem is the only one that returns data for the first block of the key (atom connections)
    context = 'physchem' if only_first_layer else 'all'
    key = inchikey[:14] if only_first_layer else inchikey
    url = f"https://www.lipidmaps.org/rest/compound/inchi_key/{key}/{context}/json"
    headers = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64)'}
    response = requests.get(url, headers=headers)
    check_lipid_response(response, 'LIPID MAPS', key)
    js = response.json()
    # LIPID MAPS answers an empty list when the InChI key is not found
    if js == []: return False
    if not only_first_layer: return js
    # Several matches come as {'Row1': {...}, 'Row2': {...}} and a single match as the entry itself
    rows = list(js.values()) if 'Row1' in js else [js]
    # Prefer the entry with the same stereo layer (first two blocks), so only the protonation changes
    row = next((row for row in rows if row['inchi_key'][:25] == inchikey[:25]), rows[0])
    # Get the whole entry data with its full InChI key
    return is_in_LIPID_MAPS(row['inchi_key'])


def get_swisslipids_info(entity_id) -> dict:
    """Get information about a SwissLipids entry."""
    headers = {'accept': 'json'}
    url = f"https://www.swisslipids.org/api/index.php/entity/{entity_id}"
    response = requests.get(url, headers=headers)
    check_lipid_response(response, 'SwissLipids', entity_id)
    return response.json()


@retry_request
@lru_cache(maxsize=None)
def is_in_swisslipids(inchikey, only_first_layer=False) -> dict:
    """Search the InChI keys in SwissLipids.
    Documentation: https://www.swisslipids.org/#/api.
    """
    key = inchikey[:14] if only_first_layer else inchikey
    headers = {'accept': 'json'}
    url = f"https://www.swisslipids.org/api/index.php/search?term={key}"
    response = requests.get(url, headers=headers)
    # SwissLipids answers 404 when the InChI key is not found
    if response.status_code == 404: return None
    check_lipid_response(response, 'SwissLipids', key)
    data = response.json()[0]
    detailed_data = get_swisslipids_info(data['entity_id'])
    data['synonyms'] = detailed_data.get('synonyms', [])
    return data
