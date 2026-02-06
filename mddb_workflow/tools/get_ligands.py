import re
import json

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from mordred import Calculator, descriptors
from mddb_workflow.utils.constants import LIGANDS_MATCH_FLAG, PDB_TO_PUBCHEM
from mddb_workflow.utils.auxiliar import InputError, request_pdb_data, warn, retry_request, handle_http_request
from mddb_workflow.utils.type_hints import *


from functools import lru_cache
import requests
from scipy.optimize import linear_sum_assignment
import numpy as np

# Set the expected ligand data fields
LIGAND_DATA_FIELDS = set(['name', 'pubchem', 'drugbank', 'chembl', 'smiles', 'formula', 'morgan', 'mordred', 'pdbid', 'inchikey'])
MINIMUM_TANIMOTO_THRESHOLD = 0.3


@retry_request
def _make_get_request(url: str) -> requests.Response:
    return requests.get(url, timeout=30)


@retry_request
def _make_post_request(url, payload):
    return requests.post(url, data=payload, timeout=30)


def record_pubchem_match(
    inchikey: str,
    pubchem_id: str,
    match_mol: Chem.Mol,
    reference_fg: dict,
    ligand_references: dict,
    reason: str,
    threshold: float = MINIMUM_TANIMOTO_THRESHOLD,
) -> None:
    """Set the pubchem id for the inchikey and compute/print Tanimoto vs descriptor_data if match_mol provided."""
    print(f'Comparing {inchikey} {reason} with CID {pubchem_id}.')
    matched_fg = obtain_mordred_morgan_descriptors(match_mol)
    tc = tanimoto_similarity(reference_fg, matched_fg['morgan_fp_bit_array'])
    # For ions, set threshold to 0.0 as TC can vary a lot with small changes
    threshold = threshold if match_mol.GetNumAtoms() > 1 else 0.0
    if tc >= threshold:
        ligand_references[inchikey]['pubchem'] = str(pubchem_id)
        ligand_references[inchikey]['tc'] = tc
        result_str = f'    TC above threshold {threshold}. Match accepted.'
    else:
        result_str = f'    TC below threshold {threshold}. Ignoring match.'
    print(f'    Tanimoto coef. (morgan fingerprint): {tc:.3f}' + result_str)
    return tc >= threshold


def generate_ligand_references(
    structure: 'Structure',
    cache: 'Cache',
    input_ligands: Optional[list[dict]],
    pdb_ids: list[str],
    inchikeys: dict[str, 'InChIKeyData'],
    lipid_references: dict[str, dict],
    membrane_map: dict,
    mercy: list[str] = [],
    ) -> dict[str, dict]:
    """Generate a map of residues associated to ligands.

    This function identifies and maps ligands in the molecular structure through a
    multi-step matching process:

    1. **Direct InChIKey matching**: Extracts InChIKeys from structure fragments
       (excluding lipids and membrane components) and attempts direct matching with
       PubChem database.

    2. **Chemical similarity matching**: If direct matching fails, progressively
       modifies the molecular structure and calculates Tanimoto coefficient (TC)
       for similarity assessment:
       1. Neutralize charges
       2. Remove stereochemistry information
       3. Apply PubChem standardization (tautomer, protonation, etc.)
       4. Match against PDB-derived ligands (TC threshold ≥ 0.9)
       5. Perform similarity search in PubChem/ChEMBL (TC threshold ≥ 0.9)

    3. **Fallback handling**: Unmatched ligands are saved as-is with warnings.

    4. **User-forced selections**: Respects user-specified ligand selections from
       inputs.yaml, with warnings if TC compared to original fragment is insufficient.
    """
    # Check input ligands format validity
    input_ligands = input_ligands or []
    for i, ligand in enumerate(input_ligands):
        if type(ligand) is int or (type(ligand) is str and ligand.isnumeric()):
            print(f'A ligand number ID has been identified {ligand}, assuming that is a PubChem ID...')
            input_ligands[i] = {'pubchem': str(ligand)}
        elif type(ligand) is str:
            raise InputError(f'A name of ligand has been identified: {ligand}. Anyway, provide at least one of the following IDs: DrugBank, PubChem, ChEMBL.')
    # Get non-lipid inchikeys
    ligand_keys = {key: {} for key in inchikeys if key not in lipid_references.keys()}
    residx_2_inchikey = {resid: inchikey
                         for inchikey, data in inchikeys.items()
                         for resid in data.resindices}
    # Add non-membrane lipid inchikeys
    for resid in membrane_map['no_mem_lipid']:
        # Some residues may not have InChIKeys (e.g. CG)
        inchikey = residx_2_inchikey.get(resid, None)
        if inchikey and inchikey not in ligand_keys:
            ligand_keys[inchikey][resid] = {}
    # Create a dictionary to store the references generated without user input
    automatic_references = {key: {} for key in ligand_keys}
    pubchem_ids_from_pdb = []
    pdb_ids = pdb_ids or []
    pubchem_ids_from_pdb = pdbs_2_pubchems(pdb_ids, cache)
    for inchikey in ligand_keys:
        # 1. Direct InChIKey match
        pubchem_id = inchikey_2_pubchem(inchikey)
        inchi = inchikeys[inchikey].inchi
        mol = Chem.MolFromInchi(inchi, sanitize=False)
        # Original descriptors to calculate similarity later
        descriptor_data = obtain_mordred_morgan_descriptors(mol)
        ligand_fg = descriptor_data['morgan_fp_bit_array']
        automatic_references[inchikey].update(descriptor_data)
        if pubchem_id:
            automatic_references[inchikey]['pubchem'] = pubchem_id
            print(f'Found CID {pubchem_id} by direct InChIKey match.')
            continue
        # Try catch for some weird inchis (MCV1900695, no topology)
        try:
            Chem.SanitizeMol(mol, catchErrors=False)
            # print_molecule_terminal(mol)
        except Exception:
            warn(f'Could not sanitize molecule with InChIKey {inchikey}. Skipping further matching attempts.')
            # print_molecule_terminal(mol)
            continue
        # 2.1. Neutralize charges
        neutral_mol = rdMolStandardize.ChargeParent(mol)
        neutral_inchikey = Chem.MolToInchiKey(neutral_mol)
        if pubchem_id := inchikey_2_pubchem(neutral_inchikey):
            if record_pubchem_match(inchikey, pubchem_id, neutral_mol, ligand_fg,
                                    automatic_references, 'after neutralization'):
                continue
        # 2.2. Remove stereochemistry
        snon_inchikey = Chem.MolToInchiKey(mol, options='-SNon')
        if pubchem_id := inchikey_2_pubchem(snon_inchikey):
            Chem.RemoveStereochemistry(mol)
            if record_pubchem_match(inchikey, pubchem_id, mol, ligand_fg,
                                    automatic_references, 'after stereochemistry removal'):
                continue
        # Neutralize charges and remove stereochemistry
        snon_neutral_inchikey = Chem.MolToInchiKey(neutral_mol, options='-SNon')
        if pubchem_id := inchikey_2_pubchem(snon_neutral_inchikey):
            if record_pubchem_match(inchikey, pubchem_id, mol, ligand_fg,
                                    automatic_references, 'after neutralization + stereochemistry removal'):
                continue
        # 2.3. Standardize molecule
        standar_id = pubchem_standardization(inchi)
        if standar_id and len(standar_id) == 1:
            standar_pubchem = standar_id[0]['pubchem']
            standar_mol = Chem.MolFromInchi(standar_id[0]['inchi'])
            if record_pubchem_match(inchikey, standar_pubchem, standar_mol, ligand_fg,
                                    automatic_references, 'after standardization'):
                continue
        # 2.4. Match against PDB-derived ligands
        if pubchem_ids_from_pdb:
            # We remove Hs to be more permissive with the matching
            # This its okey since we can be biased towards matching to PDB ligands
            noH_mol = Chem.RemoveHs(mol)
            descriptor_data = obtain_mordred_morgan_descriptors(noH_mol)
            noH_ligand_fg = descriptor_data['morgan_fp_bit_array']
            for pdb_ligand in pubchem_ids_from_pdb:
                if not pdb_ligand.get('inchi', None):
                    ligand_data = obtain_pubchem_data_from_input(pdb_ligand)
                    pdb_ligand.update(ligand_data)

                pdb_mol = Chem.MolFromInchi(pdb_ligand['inchi'])
                success = record_pubchem_match(inchikey, pdb_ligand['pubchem'], pdb_mol, noH_ligand_fg,
                                        automatic_references, 'with PDB-derived ligand')
                if success:
                    break
        # RUBEN: disabled as it does not work very well and it is slow
        if False:
            # 2.5. Similarity search in PubChem using neutral molecule
            neutral_inchi = Chem.MolToInchi(neutral_mol)
            for threshold in [95, 90]:
                similarity_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/inchi/JSON?Threshold={threshold}&MaxRecords=5"
                r = _make_post_request(similarity_url, payload={'inchi': neutral_inchi})
                if r.status_code != 200:
                    continue
                data = r.json()
                tc = 0
                for compound in data['PC_Compounds']:
                    if compound['id']['id']['cid']:
                        cid = str(compound['id']['id']['cid'])
                    for prop in compound['props']:
                        if prop['urn']['label'] == 'InChI':
                            standard_inchi = prop['value']['sval']
                    if not compound or not standard_inchi:
                        assert False, "PubChem similarity search did not return valid compound data"
                    simil_mol = Chem.MolFromInchi(standard_inchi)
                    success = record_pubchem_match(inchikey, cid, simil_mol, ligand_fg,
                                        ligand_references, 'similarity search in PubChem')
                    if success: break
                if success: break

    for inchikey, ligand_data in automatic_references.items():
        if ligand_data.get('pubchem', None):
            automatic_references[inchikey].update(obtain_pubchem_data_from_input(ligand_data))

    # 4. User-forced selections
    input_ligands_data = []
    user_references = {}
    for input_ligand in input_ligands:
        ligand_data = obtain_pubchem_data_from_input(input_ligand)
        # If the user defined a ligand name, it will be respected and added to the metadata
        if forced_name := input_ligand.get('name', None):
            ligand_data['forced_name'] = forced_name
        forced_selection = input_ligand.get('selection', None)
        if forced_selection:
            # Could be a single residue or a list of residues
            selection_atoms = structure.select(forced_selection)
            residue_indices = structure.get_selection_residue_indices(selection_atoms)
            print(f'User forced selection of {len(residue_indices)} residue/s for ligand {input_ligand}.')
            if len(residue_indices) == 0 and LIGANDS_MATCH_FLAG not in mercy:
                raise InputError(f'Ligand with PubChem ID {pubchem_id} did not map with any residue')
            ligand_data['resindices'] = residue_indices
        # Search for this pubchem ID in automatic_references
        matched_inchikey = None
        for auto_inchikey, auto_data in automatic_references.items():
            if auto_data.get('pubchem') == ligand_data['pubchem']:
                matched_inchikey = auto_inchikey
                break
        # If found, move from automatic_references to user_references
        if matched_inchikey:
            user_references[matched_inchikey] = automatic_references.pop(matched_inchikey)
            if forced_selection:
                user_references[matched_inchikey]['resindices'] = residue_indices
        elif forced_selection:
            # If not found, create a new entry in user_references
            user_references[ligand_data['inchikey']] = ligand_data
        else:
            # Add to input ligands for similarity search
            input_fg = obtain_mordred_morgan_descriptors(Chem.MolFromInchi(ligand_data['inchi']))
            input_ligands_data.append({**input_ligand, **ligand_data, **input_fg})
    # Create similarity matrix (N x M)
    auto_inchikeys = list(automatic_references.keys())
    similarity_matrix = np.zeros((len(input_ligands_data), len(auto_inchikeys)))
    for i, input_fg in enumerate(input_ligands_data):
        for j, auto_inchikey in enumerate(auto_inchikeys):
            tc = tanimoto_similarity(
                input_fg['morgan_fp_bit_array'],
                automatic_references[auto_inchikey]['morgan_fp_bit_array'])
            similarity_matrix[i, j] = tc
    # Use Hungarian algorithm to find optimal assignment (maximize similarity)
    # Convert to minimization problem by negating similarities
    row_ind, col_ind = linear_sum_assignment(-similarity_matrix)
    # Process matched pairs
    for input_idx, ref_idx in zip(row_ind, col_ind):
        input_data = input_ligands_data[input_idx]
        matched_inchikey = auto_inchikeys[ref_idx]
        tc = similarity_matrix[input_idx, ref_idx]

        if tc < 0.6:
            warn(f'Ligand "{input_data.get("name", input_data["inchikey"])}" input selection has Tanimoto coefficient {tc:.3f} '
                 f'which is below the threshold of 0.6 when compared to the matched fragment (InChIKey: {matched_inchikey}).')

        # Update ligand reference with input data
        user_references[matched_inchikey] = input_data
        automatic_references.pop(matched_inchikey)
        print(f'Matched input ligand {input_data["pubchem"]} to reference {matched_inchikey} with TC={tc:.3f}')

    # Merge user and automatic references into ligand_references
    ligand_references = {**user_references, **automatic_references}
    not_matched_ligands = []
    for k, v in ligand_references.items():
        if 'pubchem' not in v:
            not_matched_ligands.append((inchikeys[k].molname, inchikeys[k].resindices))
            breakpoint()
    if not_matched_ligands:
        for k, v in not_matched_ligands:
            warn(f'Ligand {k} could not be matched to any PubChem ID. Residues: {v}')
        if LIGANDS_MATCH_FLAG not in mercy:
            raise InputError('Provide PubChem IDs for all ligands, a PDB code where it is '
                             f'present or use the flag "-m {LIGANDS_MATCH_FLAG}" to bypass this check.')

    if not ligand_references:
        print('No ligands were matched')

    return ligand_references


def pdb_ligand_2_prd(pdb_id: str) -> Optional[str]:
    """Given a PDB ID, get its PRD ligand code."""
    query = '''query structure($id: String!) {
        entry(entry_id: $id) {
            pdbx_molecule_features {
                prd_id
            }
        }
    }'''
    parsed_response = request_pdb_data(pdb_id, query)
    pdbx_id = parsed_response.get('pdbx_molecule_features', None)
    if pdbx_id is None:
        print(f'WARNING: Cannot find PRD ligand code for PDB ID {pdb_id}')
        return None
    prd_id = pdbx_id[0].get('prd_id', None)  # AGUS: podría haber casos donde haya más de uno?
    if prd_id is None:
        print(f'WARNING: Cannot find PRD ligand code for PDB ID {pdb_id}')
        return None
    # If the PRD ID is not empty then return it
    return prd_id


def get_pdb_ligand_codes(pdb_id: str) -> list[str]:
    """Given a PDB ID, get all its ligand codes.
    e.g. 2I3I -> 618, BTB, ZN, EDO, LI.
    """
    # Set the request query
    query = '''query structure($id: String!) {
        entry(entry_id: $id) {
            nonpolymer_entities { nonpolymer_comp { chem_comp { id } } }
        }
    }'''
    # Request PDB data
    parsed_response = request_pdb_data(pdb_id, query)
    # The response may be None
    # e.g. an obsolete entry with no replacement
    if parsed_response is None: return []
    # Mine data for nonpolymer entities
    nonpolymers = parsed_response['nonpolymer_entities']
    # If there are no nonpolymer entities, another type of entitie could be used
    # AGUS: esto es un caso muy concreto que me encontré con el PDB 1N3W
    # AGUS: el ligando en este caso es una 'Biologically Interesting Molecules' y se muestran como 'PRD_'
    prd_code = None
    if nonpolymers is None:
        # Get the prd ligand code
        prd_code = pdb_ligand_2_prd(pdb_id)
        if prd_code is not None:
            return [prd_code]

    if nonpolymers is None and prd_code is None: return []
    # Iterate nonpolymer entities to mine each PDB code
    ligand_codes = []
    for nonpolymer in nonpolymers:
        ligand_code = nonpolymer['nonpolymer_comp']['chem_comp']['id']
        ligand_codes.append(ligand_code)
    print(f' Ligand codes for PDB ID {pdb_id}: ' + ', '.join(ligand_codes))
    return ligand_codes


def pdb_ligand_2_pubchem(pdb_ligand_id: str) -> list:
    """Given a PDB ligand code, get its PubChem ID."""
    # Set the request query
    query = '''query molecule($id:String!){
        chem_comp(comp_id:$id) { rcsb_chem_comp_related{ resource_name resource_accession_code } }
    }'''
    # Request PDB data
    parsed_response = request_pdb_data(pdb_ligand_id, query)
    related_resources = parsed_response['rcsb_chem_comp_related']
    # It may happend that a ligand code has no related resources at all
    # e.g. ZN
    if not related_resources: return []
    pubchem_resource = [resource['resource_accession_code']
                        for resource in related_resources
                        if resource['resource_name'] == 'PubChem']
    return pubchem_resource


def pdb_ligand_2_pubchem_RAW(pdb_ligand_id: str) -> Optional[str]:
    """Given a PDB ligand code, get its PubChem ID.
    Use a web crawler to avoid having to use the PDB API.
    """
    # Set the request URL
    request_url = f'https://www.rcsb.org/ligand/{pdb_ligand_id}'
    # Run the query
    parsed_response = handle_http_request(request_url, "PDB ligand request")
    # Mine the PubChem ID out of the whole response
    pattern = re.compile('pubchem.ncbi.nlm.nih.gov\/compound\/([0-9]*)\"')
    match = re.search(pattern, parsed_response)
    # If there is no PubChem ID then return none
    # This is normal for some ligands such as counter ions (e.g. LI)
    if not match:
        return None
    pubchem_id = match[1]
    return pubchem_id


# DANI: No se ha provado a fondo
def pdb_ligand_2_pubchem_RAW_RAW(pdb_ligand_id: str) -> Optional[str]:
    """Given a PDB ligand code, get its PubChem ID.
    Ask to PubChem if the ligand exists and hope there is only one result.
    """
    # Set the request URL
    request_url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{pdb_ligand_id}/json'
    # Run the query
    response = handle_http_request(request_url, "PDB ligand request")
    # There may be no result for unkown atoms (e.g. UNK), ions (e.g. IOD) and some other ligands (e.g. DX9)
    if not response: return None
    parsed_response = json.loads(response)
    # Mine the PubChem ID
    compounds = parsed_response['PC_Compounds']
    if len(compounds) != 1:
        # There could be more than one result
        # For example, the case HEM: https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/HEM/json
        # In this case we picked the first result
        compound = compounds[0]
        # raise RuntimeError('We are not having 1 and only 1 result from PubChem: ' + request_url)
    compound = compounds[0]
    id1 = compound['id']
    id2 = id1['id']
    pubchem_id = str(id2['cid'])
    return pubchem_id


def pdb_2_pubchem(pdb_id: str) -> list[str]:
    """Given a PDB ID, get its PubChem IDs."""
    print(f'Searching PubChem IDs for PDB ID {pdb_id}')
    pubchem_ids = []
    # Iterate over pdb ligand codes
    ligand_codes = get_pdb_ligand_codes(pdb_id)
    for ligand_code in ligand_codes:
        # Ask the PDB API for the ligand
        pubchem_id = pdb_ligand_2_pubchem(ligand_code)
        # If this did not work then try mining the PDB client with a web crawler
        if not pubchem_id:
            pubchem_id = pdb_ligand_2_pubchem_RAW(ligand_code)
        # If this did not work then try it from PubChem
        if not pubchem_id:
            pubchem_id = pdb_ligand_2_pubchem_RAW_RAW(ligand_code)
        # Otherwise we surrender
        if not pubchem_id:
            print(f' {ligand_code} -> No PubChem ID')
            continue
        print(f' {ligand_code} -> {pubchem_id}')
        pubchem_ids.extend(pubchem_id)

    return pubchem_ids


def pdbs_2_pubchems(pdb_ids: list[str], cache: 'Cache') -> list[dict]:
    """Given a list of PDB IDs, get their PubChem IDs and add them to the ligands list."""
    pdb_ligands = []
    # Check we have cached pdb 2 PubChem values
    pdb_2_pubchem_cache = cache.retrieve(PDB_TO_PUBCHEM, {})
    new_data_to_cache = False
    # Get input ligands from the pdb IDs, if any
    for pdb_id in pdb_ids:
        # Check we have cached this specific pdb
        pubchem_ids_from_pdb = pdb_2_pubchem_cache.get(pdb_id, None)
        if pubchem_ids_from_pdb is not None:
            print(f' Retrieving from cache PubChem IDs for PDB ID {pdb_id}: ')
            if len(pubchem_ids_from_pdb) > 0:
                print('  PubChem IDs: ' + ', '.join(pubchem_ids_from_pdb))
            else:
                print('  This PDB ID has no PubChem IDs')

        # If we had no cached pdb 2 PubChem then ask for them
        if pubchem_ids_from_pdb is None:
            pubchem_ids_from_pdb = pdb_2_pubchem(pdb_id)
            # Save the result in the cache object so it is saved to cache later
            pdb_2_pubchem_cache[pdb_id] = pubchem_ids_from_pdb
            new_data_to_cache = True
        for pubchem_id in pubchem_ids_from_pdb:
            # Ligands in the structure (PDB) and the 'inputs.json' could be the same so it's not necessary to do it twice
            if not any('pubchem' in ligand and ligand['pubchem'] == pubchem_id for ligand in pdb_ligands):
                pdb_ligands.append({'pubchem': pubchem_id, 'pdb': True})
    # Save all pdb to PubChem results in cache, in case there is anything new
    if new_data_to_cache: cache.update(PDB_TO_PUBCHEM, pdb_2_pubchem_cache)

    return pdb_ligands


def get_pubchem_data(pubchem_id: str) -> Optional[dict]:
    """Given a PubChem ID, use the uniprot API to request its data and then mine what is needed for the database."""
    request_url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{pubchem_id}/JSON/'
    parsed_response = handle_http_request(request_url, "PubChem request")
    if parsed_response is None:
        warn('Cannot find PubChem entry for accession ' + pubchem_id)
        return None
    parsed_response = json.loads(parsed_response)
    # Set part of the error message, in case we find a problem
    error_message = f'Unexpected PubChem data structure in {request_url}\n  '
    # Mine target data: SMILES
    record = parsed_response.get('Record', None)
    if record is None:
        raise RuntimeError(error_message + 'no record')
    sections = record.get('Section', None)
    if sections is None:
        raise RuntimeError(error_message + 'no sections')
    names_and_ids_section = next((section for section in sections if section.get('TOCHeading', None) == 'Names and Identifiers'), None)
    if names_and_ids_section is None:
        raise RuntimeError(error_message + 'no name and IDs section')
    names_and_ids_subsections = names_and_ids_section.get('Section', None)
    if names_and_ids_subsections is None:
        raise RuntimeError(error_message + 'no name and IDs subsection')

    # Mine the name from the synonims section
    name_substance = None
    synonims = next((s for s in names_and_ids_subsections if s.get('TOCHeading', None) == 'Synonyms'), None)
    if synonims:
        synonims_subsections = synonims.get('Section', None)
        if synonims_subsections is None:
            raise RuntimeError(error_message + 'no synonims subsection')
        depositor_supplied_synonims = next((s for s in synonims_subsections if s.get('TOCHeading', None) == 'Depositor-Supplied Synonyms'), None)
        if depositor_supplied_synonims is None:
            removed_synonims = next((s for s in synonims_subsections if s.get('TOCHeading', None) == 'Removed Synonyms'), None)
            name_substance = removed_synonims.get('Information', None)[0].get('Value', {}).get('StringWithMarkup', None)[0].get('String', None)
        else:
            name_substance = depositor_supplied_synonims.get('Information', None)[0].get('Value', {}).get('StringWithMarkup', None)[0].get('String', None)
    # If we still don't have a name then try with the IUPAC name
    if name_substance is None:
        descriptors = next((s for s in names_and_ids_subsections if s.get('TOCHeading', None) == 'Computed Descriptors'), None)
        descriptors_subsections = descriptors.get('Section', None)
        if descriptors_subsections is None:
            raise RuntimeError(error_message + 'no descriptors subsection')
        iupac_name_section = next((s for s in descriptors_subsections if s.get('TOCHeading', None) == 'IUPAC Name'), None)
        if iupac_name_section:
            name_substance = iupac_name_section.get('Information', None)[0].get('Value', {}).get('StringWithMarkup', None)[0].get('String', None)
    if name_substance is None:
        # If we still do not have a name then we assume the compound has no name
        # This may happen (e.g. 57449604)
        name_substance = 'Unnamed'

    # Mine the SMILES
    computed_descriptors_subsection = next((s for s in names_and_ids_subsections if s.get('TOCHeading', None) == 'Computed Descriptors'), None)
    if computed_descriptors_subsection is None:
        raise RuntimeError(error_message + 'no computed descriptors')
    canonical_smiles_section = computed_descriptors_subsection.get('Section', None)
    if canonical_smiles_section is None:
        raise RuntimeError(error_message + 'no canonical SMILES section')
    canonical_smiles = next((s for s in canonical_smiles_section if s.get('TOCHeading', None) == 'Canonical SMILES'), None)
    if canonical_smiles is None:
        # In some cases there is no canonical SMILES but a non-canonical one could exists
        non_canonical_smiles_section = next((s for s in canonical_smiles_section if s.get('TOCHeading', None) == 'SMILES'), None)
        if non_canonical_smiles_section is None:
            raise RuntimeError(error_message + 'no canonical SMILES')

    if canonical_smiles:
        smiles = canonical_smiles.get('Information', None)[0].get('Value', {}).get('StringWithMarkup', None)[0].get('String', None)
    if non_canonical_smiles_section:
        smiles = non_canonical_smiles_section.get('Information', None)[0].get('Value', {}).get('StringWithMarkup', None)[0].get('String', None)

    if smiles is None:
        raise RuntimeError(error_message + 'no SMILES')

    # Mine target data: MOLECULAR FORMULA
    molecular_formula_subsection = next((s for s in names_and_ids_subsections if s.get('TOCHeading', None) == 'Molecular Formula'), None)
    if molecular_formula_subsection is None:
        raise RuntimeError(error_message + 'no molecular formula section')
    molecular_formula = molecular_formula_subsection.get('Information', None)[0].get('Value', {}).get('StringWithMarkup', None)[0].get('String', None)
    if molecular_formula is None:
        raise RuntimeError(error_message + 'no molecular formula')

    # Mine target data: PDB ID
    pdb_id = None
    pdb_id_subsection = next((s for s in sections if s.get('TOCHeading', None) == 'Interactions and Pathways'), None)
    # If this section is missing then it means this PubChem compound has no PDB ID
    if pdb_id_subsection:
        pdb_id_subsections = pdb_id_subsection.get('Section', None)
        if pdb_id_subsections is None:
            raise RuntimeError(error_message + 'no PDB IDs subsection')
        bond_structures = next((s for s in pdb_id_subsections if s.get('TOCHeading', None) == 'Protein Bound 3D Structures'), None)
        if bond_structures:
            bond_structures_section = bond_structures.get('Section', None)
            # If this section is missing then it means this PubChem compound has no PDB ID
            if bond_structures_section:
                ligands_structure = next((s for s in bond_structures_section if s.get('TOCHeading', None) == 'Ligands from Protein Bound 3D Structures'), None)
                if ligands_structure is None:
                    raise RuntimeError(error_message + 'no Ligands from Protein Bound 3D Structures section')
                ligands_structure_subsections = ligands_structure.get('Section', None)
                if ligands_structure_subsections is None:
                    raise RuntimeError(error_message + 'no Ligands from Protein Bound 3D Structures section')
                ligands_pdb = next((s for s in ligands_structure_subsections if s.get('TOCHeading', None) == 'PDBe Ligand Code'), None)
                if ligands_pdb is None:
                    raise RuntimeError(error_message + 'no PDBe ligand code')
                pdb_id = ligands_pdb.get('Information', None)[0].get('Value', {}).get('StringWithMarkup', None)[0].get('String', None)

    # Mine de INCHI and INCHIKEY
    inchi_section = next((s for s in canonical_smiles_section if s.get('TOCHeading', None) == 'InChI'), None)
    if inchi_section is None:
        raise RuntimeError(error_message + 'no InChI section')
    if inchi_section:
        inchi = inchi_section.get('Information', None)[0].get('Value', {}).get('StringWithMarkup', None)[0].get('String', None)

    inchikey_section = next((s for s in canonical_smiles_section if s.get('TOCHeading', None) == 'InChIKey'), None)
    if inchikey_section is None:
        raise RuntimeError(error_message + 'no InChI key')
    if inchikey_section:
        inchikey = inchikey_section.get('Information', None)[0].get('Value', {}).get('StringWithMarkup', None)[0].get('String', None)
    # Prepare the PubChem data to be returned
    return {'name': name_substance, 'smiles': smiles, 'formula': molecular_formula, 'pdbid': pdb_id, 'inchi': inchi, 'inchikey': inchikey}


def drugbank_2_pubchem(drugbank_id) -> Optional[str]:
    """Given a DrugBank ID, request its PubChem compound ID."""
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drugbank_id}/JSON'
    r = _make_get_request(url)
    if r.ok:
        data = r.json()
        return str(data['PC_Compounds'][0]['id']['id']['cid'])
    else:
        raise ValueError("Something went wrong with the DrugBank to PubChem request (error " + str(r.status_code) + ")")


def chembl_2_pubchem(chembl_id) -> Optional[str]:
    """Given a ChemBL ID, use the uniprot API to request its PubChem compound ID."""
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/RegistryID/{chembl_id}/cids/JSON'
    r = _make_get_request(url)
    if r.ok:
        data = r.json()
        return str(data['IdentifierList']['CID'][0])
    else:
        raise ValueError("Something went wrong with the ChEMBL to PubChem request (error " + str(r.status_code) + ")")


def obtain_mordred_morgan_descriptors(mol: Chem.Mol) -> dict:
    """Calculate the Morgan fingerprint and the Mordred descriptors from a molecule.

    References:
        - RDKit fingerprints: https://www.rdkit.org/docs/RDKit_Book.html#additional-information-about-the-fingerprints
        - RDKit tutorial: https://greglandrum.github.io/rdkit-blog/posts/2023-01-18-fingerprint-generator-tutorial.html
        - Fingerpring comparison: https://pmc.ncbi.nlm.nih.gov/articles/PMC10964529/

    """
    Chem.AllChem.Compute2DCoords(mol)
    mol_block = Chem.MolToMolBlock(mol)
    # We can select the different submodules of mordred descriptors, avaible in: 'https://mordred-descriptor.github.io/documentation/master/'
    calc = Calculator([
        descriptors.ABCIndex,  # Índice de ramificación
        descriptors.AcidBase.AcidicGroupCount,  # Grupos ácidos
        descriptors.AcidBase.BasicGroupCount,  # Grupos básicos
        descriptors.RingCount,  # Conteo de anillos
        descriptors.Constitutional,  # Propiedades generales como número de átomos, peso molecular
        descriptors.TopologicalCharge,  # Índices topológicos, Cargas parciales, polaridad
        descriptors.HydrogenBond,  # Donantes y aceptores de enlaces de hidrógeno
        descriptors.Lipinski,  # Reglas de Lipinski (drug-likeness)
        descriptors.FragmentComplexity,  # Identificación de subestructuras frecuentes
        descriptors.PathCount,  # Conteo de caminos moleculares
    ], ignore_3D=True)

    # Calculate Mordred results
    try:
        mordred_results = calc(mol).drop_missing().asdict()
    except Exception as e:
        print(f'WARNING: Mordred calculation failed with error: {e}. Retrying with ignore_3D=False.')
        mordred_results = {}

    # MORGAN FINGERPRINT
    morgan_fp_gen = Chem.rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    ao = Chem.rdFingerprintGenerator.AdditionalOutput()
    ao.AllocateAtomCounts()
    ao.AllocateAtomToBits()
    ao.AllocateBitInfoMap()
    fp = morgan_fp_gen.GetFingerprint(mol, additionalOutput=ao)
    morgan_fp_bit_array = list(fp)
    morgan_highlight_atoms = {}
    for bit, atoms in ao.GetBitInfoMap().items():
        morgan_highlight_atoms[bit] = list(set(atom for atom, radius in atoms))
    data = {'mordred': mordred_results,
            'morgan_fp_bit_array': morgan_fp_bit_array,
            'morgan_highlight_atoms': morgan_highlight_atoms,
            'mol_block': mol_block}
    return data


def tanimoto_similarity(fp1: list[int], fp2: list[int]) -> float:
    """Calculate the Tanimoto similarity between two fingerprints."""
    if len(fp1) != len(fp2):
        raise ValueError('Fingerprints must be of the same length')
    a = sum(1 for i in range(len(fp1)) if fp1[i] == 1 and fp2[i] == 1)
    b = sum(1 for i in range(len(fp1)) if fp1[i] == 1 and fp2[i] == 0)
    c = sum(1 for i in range(len(fp1)) if fp1[i] == 0 and fp2[i] == 1)
    if (a + b + c) == 0:
        return 0.0
    return a / (a + b + c)


def obtain_pubchem_data_from_input(ligand: dict) -> dict[str, Optional[str]]:
    """Given an input ligand, obtain all necessary data."""
    # Save in a dictionary all ligand data including its name and IDs
    # The ID can be of the databases: 'drugbank' , 'pubchem' , 'chembl'
    # Define the needed variables to check if the ligand has a database ID or it is None
    ligand_data = {field: None for field in LIGAND_DATA_FIELDS}
    # Set ligand data PubChem ID, even if the input ID is not from pubhcme (e.g. drugbank, chembl)
    if 'pubchem' in ligand:
        ligand_data['pubchem'] = str(ligand.get('pubchem'))
    elif 'drugbank' in ligand:
        ligand_data['drugbank'] = ligand.get('drugbank')
        ligand_data['pubchem'] = str(drugbank_2_pubchem(ligand_data['drugbank']))
    elif 'chembl' in ligand:
        ligand_data['chembl'] = ligand.get('chembl')
        ligand_data['pubchem'] = str(chembl_2_pubchem(ligand_data['chembl']))
    else:
        raise InputError('None of the ligand IDs are defined. Please provide at least one of the following IDs: DrugBank, PubChem, ChEMBL.')

    # Request ligand data from pubchem
    pubchem_data = get_pubchem_data(ligand_data['pubchem'])
    if not pubchem_data:
        raise RuntimeError('No PubChem data avilable')

    # Add PubChem data to ligand data
    ligand_data = {**ligand_data, **pubchem_data}
    return ligand_data


@lru_cache(maxsize=None)
def inchikey_2_pubchem(inchikey: str, conectivity_only=False) -> Optional[str]:
    """Given an InChIKey, get the PubChem CID."""
    inchikey = inchikey.split('-')[0] if conectivity_only else inchikey
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/cids/JSON"

    r = _make_get_request(url)
    if r.ok:
        data = r.json()
        if conectivity_only:
            return [str(cid) for cid in data['IdentifierList']['CID']]
        else:
            return str(data['IdentifierList']['CID'][0])
    return None


@lru_cache(maxsize=None)
def pubchem_standardization(inchi: str) -> Optional[list[dict]]:
    """Try PubChem standardization service to write things like the bonds from aromatic rings in the same order.

    You can see examples on the paper https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0293-8
    """
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/standardize/inchi/JSON"

    r = _make_post_request(url, payload={'inchi': inchi})
    if r.ok:
        data = r.json()
        results = []
        for compound in data['PC_Compounds']:
            if compound['id']['id']['cid']:
                cid = str(compound['id']['id']['cid'])
            else:
                continue
            standard_inchi = None
            for prop in compound['props']:
                if prop['urn']['label'] == 'InChI':
                    standard_inchi = prop['value']['sval']
            results.append({'pubchem': cid, 'inchi': standard_inchi})
        return results
    return None


def print_molecule_terminal(mol: Chem.Mol):
    """Print a molecule in the terminal. Useful for debugging."""
    from rdkit.Chem import Draw
    im = Draw.MolToImage(mol)
    im.show()
