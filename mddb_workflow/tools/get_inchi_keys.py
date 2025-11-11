import requests
import MDAnalysis
import multiprocessing
from rdkit import Chem
from functools import lru_cache
from dataclasses import dataclass, field
from mddb_workflow.utils.structures import Structure
from mddb_workflow.utils.auxiliar import save_json, warn
from mddb_workflow.utils.type_hints import *


@dataclass
class InChIKeyData:
    """Data structure for InChI key information."""
    inchi: str
    resindices: list[int] = field(default_factory=list)
    fragments: list[list[int]] = field(default_factory=list)
    resname: set[str] = field(default_factory=set)
    classification: set = field(default_factory=set)
    frag_len: int = 1


def residue_to_inchi(task: tuple['MDAnalysis.AtomGroup', int]) -> tuple[str, str, int]:
    """Process a single residue to get its InChI key and related information."""
    resatoms, resindices = task
    # Convert to RDKIT and get InChI data
    res_RD = resatoms.convert_to.rdkit()
    # Calculate InChI key and string
    inchikey = Chem.MolToInchiKey(res_RD)
    # rdinchi.MolToInchi so it doesnt print the warnings
    inchi, retcode, message, logs, aux = Chem.rdinchi.MolToInchi(res_RD)
    return (inchikey, inchi, resindices)


@lru_cache(maxsize=None)
def is_in_LIPID_MAPS(inchikey, only_first_layer=False) -> dict:
    """Search the InChI keys in LIPID MAPS."""
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


def add_lipid_references(key_2_data: dict[str, InChIKeyData]) -> set[str]:
    """Add lipid-specific database information to InChIKeyData objects.

    This function queries SwissLipids and LIPID MAPS databases for each InChI key
    and adds the results directly to the InChIKeyData objects. It also performs
    quality checks on lipid classifications.

    Args:
        key_2_data: Dictionary mapping InChI keys to InChIKeyData objects (modified in-place).

    Returns:
        set: Set of InChI keys that were identified as lipids.

    """
    # Check internet connection
    try:
        is_in_swisslipids('test')
    except Exception as e:
        warn(f'There was a problem connecting to the SwissLipids database: {e}')
        return set()

    lipid_inchikeys = set()

    for inchikey, res_data in key_2_data.items():
        # If we dont find it, we try without stereochemistry
        SL_data = is_in_swisslipids(inchikey) or is_in_swisslipids(inchikey, only_first_layer=True)
        LM_data = is_in_LIPID_MAPS(inchikey) or is_in_LIPID_MAPS(inchikey, only_first_layer=True)

        # Add lipid database data to InChIKeyData
        if SL_data or LM_data:
            res_data.swisslipids = SL_data
            res_data.lipidmaps = LM_data
            lipid_inchikeys.add(inchikey)

            # QUALITY CHECKS
            clasi = res_data.classification
            # If the residue is a lipid, we check if it is classified as fatty/steroid
            if all('fatty' not in classes for classes in clasi) and \
                all('steroid' not in classes for classes in clasi):
                warn(f'The residue {str(res_data.resname)} is classified as {clasi}, '
                     f'but the InChIKey "{inchikey}" is a lipid.')
        else:
            # If the InChIKey is not in SwissLipids or LIPID MAPS, check classification
            if any('fatty' in classes for classes in res_data.classification):
                warn(f'The InChIKey {inchikey} of {str(res_data.resname)} is '
                     f'classified as fatty but is not a lipid.\n'
                     f'Resindices: {str(res_data.resindices)}')

    return lipid_inchikeys


def get_inchikeys(
    universe: 'MDAnalysis.Universe',
    structure: 'Structure',
    output_filepath: str,
) -> dict:
    """Generate a dictionary mapping InChI keys to residue information for non-standard residues.

    This function uses MDAnalysis to parse the input structure and topology files and identifies
    residues that are not classified as 'ion', 'solvent', 'nucleic', or 'protein'. For each
    identified residue, it converts the structure to RDKit format to obtain the InChI key
    and InChI string. The resulting data is stored in dictionaries to map InChI keys to residue
    details and residue names to InChI keys. PDB coordinates are necesary to distinguish stereoisomers.

    Args:
        universe (Universe): The MDAnalysis Universe object containing the structure and topology.
        structure (Structure): The Structure object containing residues.
        output_filepath (str): Path to save the output JSON file containing InChI key references.

    Returns:
        dict: A dictionary where keys are InChI keys and values are dictionaries containing:
            - 'inchi' (str): The InChI string for the residue
            - 'resindices' (list): A list of all residue indices with this InChI key
            - 'fragments' (list[list]): Lists of residue indices that are connected as a single group
            - 'frag_len' (int): Length of the fragments. 1 if no fragments are present.
            - 'resname' (set): Set of residue names associated with this InChI key
            - 'classification' (set): Set of residue classifications for this InChI key

    Notes:
        The function also performs consistency checks, warning if multiple residue names
        map to the same InChI key or if multiple InChI keys map to the same residue name,
        which can indicate mismatched residue definitions or stereoisomers.

    """
    try:
        universe.universe.atoms.charges
    except Exception:
        warn('Topology file does not have charges, InChI keys may be unreliable.')

    # 1) Prepare residue data for parallel processing
    # First group residues that are bonded together
    tasks = []
    residues = structure.residues

    # Fragment = residues that are bonded together
    fragments = universe.atoms.fragments
    for i, fragment in enumerate(fragments):
        resindices = fragment.residues.resindices.tolist()

        # Continue to the next fragment if any of its
        # residues are of a disallowed classification
        classes = {residues[resindex].classification for resindex in resindices}
        if (classes.intersection({'ion', 'solvent'}) or
            (classes.intersection({'dna', 'rna', 'protein'}) and len(resindices) > 1)):
            continue

        # Select residues atoms with MDAnalysis
        resatoms = universe.residues[resindices].atoms
        if 'Cg' in resatoms.types:
            # Skip coarse grain residues
            continue
        # If you pass a residue selection to a parallel worker, you a passing a whole MDAnalysis
        # universe, slowing the process down because you have to pickle the object
        # To avoid this we create
        resatoms = MDAnalysis.Merge(resatoms).universe.atoms
        # Convert to RDKit and get InChI data
        tasks.append((resatoms, resindices))

    results = []
    # Execute tasks in parallel
    with multiprocessing.Pool() as pool:
        results = pool.map(residue_to_inchi, tasks)

    # 2) Process results and build dictionaries
    key_2_data: dict[str, InChIKeyData] = {}  # To see if different name for same residue
    name_2_key = {}  # To see if different residues under same name
    for (inchikey, inchi, resindices) in results:
        # Get or create the entry for this InChI key
        data = key_2_data.setdefault(inchikey, InChIKeyData(inchi=inchi))

        # Add residue index to the list
        data.resindices.extend(resindices)
        # Add residue name to the list. For multi residues we join the names
        resname = '-'.join(sorted([residues[index].name for index in resindices]))
        data.resname.add(resname)
        # Add residue class to the list
        if len(resindices) > 1:
            classes = tuple(set([residues[index].classification for index in resindices]))
            data.classification.add(classes)
            # Glucolipids saved the groups of residues the form a 'fragment' to solve a
            # problem with FATSLiM later (ex: A01IR, A01J5)
            data.fragments.append(list(map(int, resindices)))
        else:
            data.classification.add(residues[resindices[0]].classification)

        # Incorrect residue name, estereoisomers, lose of atoms...
        name_2_key.setdefault(resname, []).append(inchikey)

    # 3) Check data coherence
    for inchikey, data in key_2_data.items():
        # Check if there are multiple names for the same InChI key
        if len(data.resname) > 1:
            warn('Same residue with different names:\n'
                f'{inchikey} -> {str(data.resname)}')
        # Check if there are multiple classifications for the same InChI key
        if len(data.classification) > 1:
            warn('Same residue with different classifications:\n'
                 f'{inchikey} + -> {str(data.classification)} for names {str(data.resname)}')
        # Check if there are multiple fragments length for the same InChI key
        if len(data.fragments) == 0:
            data.frag_len = 1
        else:
            frag_len = len(set([len(fragment) for fragment in data.fragments]))
            assert frag_len == 1, \
                f'Fragments of different lengths for InChI key {inchikey}: {str(frag_len)}'
            data.frag_len = frag_len

    # Check if there are multiple InChI keys for the same name
    for name, inchikeys in name_2_key.items():
        inchikeys = list(set(inchikeys))
        if len(inchikeys) < 2: continue
        counts = {}
        # Count the number of fragments
        for key in inchikeys:
            # If there are not fragments, we use the number of residues
            counts[key] = (key_2_data[key].frag_len
                           if key_2_data[key].frag_len > 1
                           else len(key_2_data[key].resindices))
        # Format the counts for printing
        key_counts = '\n'.join([f'\t{key}: {counts[key]: >4}' for key in inchikeys])
        warn(f'The fragment {name} has more than one InChi key:\n'
                f'{key_counts}')

    # Add lipid-specific data to key_2_data
    lipid_inchikeys = add_lipid_references(key_2_data)

    # Build InChI key references and map from enriched data
    inchikey_references = []
    inchikey_map = {}
    for inchikey, res_data in key_2_data.items():
        inchikey_references.append({
            'inchikey': inchikey,
            'inchi': res_data.inchi,
            'swisslipids': getattr(res_data, 'swisslipids', None),
            'lipidmaps': getattr(res_data, 'lipidmaps', None),
        })
        inchikey_map[inchikey] = {
            'name': list(res_data.resname)[0],
            'inchi': res_data.inchi,
            'residue_indices': list(map(int, res_data.resindices)),
            'fragments': res_data.fragments,
            'is_lipid': inchikey in lipid_inchikeys,
            'match': {
                'ref': {'inchikey': inchikey}
            }
        }

    save_json(inchikey_references, output_filepath)
    return inchikey_map


def generate_lipid_mapping(inchikeys: dict) -> dict:
    """Generate a mapping of lipid InChI keys."""
    return {key: value for key, value in inchikeys.items() if value.get('is_lipid') is True}
