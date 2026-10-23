import requests
from functools import lru_cache
from mddb_workflow.utils.auxiliar import warn, retry_request, RequestError, RateLimitError, get_retry_after, InputError
from mddb_workflow.utils.constants import LIPID_REFERENCE_VERSION, INCOMPLETE_REFERENCE_VERSION, LIPIDS_MATCH_FLAG
from mddb_workflow.tools.get_inchi_keys import get_database_inchikey_reference
from mddb_workflow.tools.lipid_classification import classify_lipid, LIPID_CATEGORIES
from mddb_workflow.tools.membrane_mapping import MIN_LEAFLET_LIPIDS
from mddb_workflow.utils.selections import Selection
from mddb_workflow.utils.type_hints import *


def generate_lipid_references(
    inchikeys: dict[str, 'InChIKeyData'],
    cache: 'Cache',
    input_forced_lipids: Optional[list[str]] = None,
    mercy: list[str] = [],
    database: Optional['Database'] = None,
) -> dict[str, dict]:
    """Find the lipids among the InChI keys and get their references in lipid databases.

    Lipids are recognised from their structure with SMARTS of the LIPID MAPS categories (see lipid_classification).
    Then SwissLipids and LIPID MAPS are queried only for these lipids to get their references.
    Lipids which are not in these databases (e.g. synthetic ionizable lipids) have no references.
    InChI keys already in the MDDB database reuse its lipid data instead.
    When SwissLipids is not available lipids have no references either, but they get an incomplete
    reference version so their InChI key references are not reused from the MDDB database and they are remade later.

    Residues with names in the forced lipids are treated as lipids although their structure is not recognised.
    Residues classified as fatty whose structure is not recognised as a lipid are not lipids.
    When there are enough of them to be part of a membrane, a missed lipid would break the membrane mapping
    so the workflow fails, unless the lipids mercy flag is passed. Otherwise they are just warned.

    Args:
        inchikeys: Dictionary mapping InChI keys to InChIKeyData objects.
        cache: Cache to store the MDDB database requests.
        input_forced_lipids: Residue names to be treated as lipids.
        mercy: Failures which do not kill the process.
        database: MDDB database handler to reuse existing references.

    Returns:
        dict: Lipid references by InChI key with their LIPID MAPS category and database data.

    """
    lipid_references = {}
    forced_lipids = set(input_forced_lipids or [])
    not_matched_fatty = []

    for inchikey, mol_data in inchikeys.items():
        clasi = mol_data.classification
        if len(clasi) == 1 and 'ion' in clasi:
            # Skip ions
            continue
        category = classify_lipid(mol_data.inchi)
        # Residues forced by the user are lipids with no database references
        if set(mol_data.resnames) & forced_lipids:
            lipid_references[inchikey] = {
                'category': category,
                'swisslipids': None,
                'lipidmaps': None,
                'forced': True,
                'version': LIPID_REFERENCE_VERSION
            }
            continue
        if not category:
            # If the structure is not a lipid, check the residue classification
            if any('fatty' in classes for classes in clasi):
                molecules = len(mol_data.resindices) // mol_data.frag_len
                warn(f'The InChI key {inchikey} of {mol_data.molname} is classified as fatty '
                     f'but its structure is not a lipid ({molecules} molecules).\n'
                     f'Resindices: {str(mol_data.resindices)}')
                if molecules >= MIN_LEAFLET_LIPIDS:
                    not_matched_fatty.append((mol_data.molname, inchikey, molecules))
            continue
        # Check if the lipid data is already in the MDDB database
        database_reference = get_database_inchikey_reference(inchikey, cache, database)
        SL_data = LM_data = None
        version = LIPID_REFERENCE_VERSION
        if database_reference:
            SL_data = database_reference.get('swisslipids')
            LM_data = database_reference.get('lipidmaps')
        # Check internet connection only once we need to request the lipid databases
        elif is_swisslipids_available():
            # If we dont find it, we try without stereochemistry (only connection layer)
            SL_data = is_in_swisslipids(inchikey) or is_in_swisslipids(inchikey, only_first_layer=True)
            LM_data = is_in_LIPID_MAPS(inchikey) or is_in_LIPID_MAPS(inchikey, only_first_layer=True)
        else:
            version = INCOMPLETE_REFERENCE_VERSION
        lipid_references[inchikey] = {
            'category': category,
            'swisslipids': SL_data or None,
            'lipidmaps': LM_data or None,
            'version': version
        }

        # QUALITY CHECKS
        # If the residue is a lipid, we check if it is classified as fatty/steroid
        # InChI keys with no residues (e.g. when remaking database references) have no classification
        if clasi and all('fatty' not in classes for classes in clasi) and \
            all('steroid' not in classes for classes in clasi):
            warn(f'The {mol_data.moltype} {mol_data.molname} is classified as {clasi}, '
                 f'but the InChIKey "{inchikey}" is from a lipid ({LIPID_CATEGORIES[category]}).')

    if not_matched_fatty and LIPIDS_MATCH_FLAG not in mercy:
        raise InputError('Some residues are classified as fatty but their structure is not recognised as a lipid:\n\t'
            + '\n\t'.join(f'{name} - {key} ({molecules} molecules)' for name, key, molecules in not_matched_fatty)
            + '\n\tAdd their residue names to the "forced_lipids" input field to treat them as lipids'
            + f' or use the flag "-m {LIPIDS_MATCH_FLAG}" to bypass this check.')

    return lipid_references


def get_lipids_selection(
    structure: 'Structure',
    inchikey_map: list[dict],
    lipid_references: dict[str, dict],
    cg_selection: Optional['Selection'] = None,
) -> 'Selection':
    """Select lipids using lipid references (InChI keys recognised as lipids) for all-atom regions.
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


def is_swisslipids_available() -> bool:
    """Check if SwissLipids answers.
    Once it answers the test request is cached, so we do not request SwissLipids again for every lipid.
    Failures are not cached, since they may be transient (e.g. during a long references update).
    """
    try:
        is_in_swisslipids('test')
        return True
    except Exception as e:
        warn(f'There was a problem connecting to the SwissLipids database: {e}')
        return False


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
