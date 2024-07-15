import os
import sys
import json
from rdkit import Chem
from rdkit.Chem import AllChem
from mordred import Calculator, descriptors
from typing import List, Tuple, Optional, Callable

from model_workflow.utils.auxiliar import InputError, load_json, save_json
from urllib.request import Request, urlopen
from urllib.parse import urlencode
from urllib.error import HTTPError, URLError
import re

def get_drugbank_smiles (id_drugbank : str) -> Optional[str]:
    # Request Drugbank
    request_url = Request(
        url= f'https://go.drugbank.com/structures/small_molecule_drugs/{id_drugbank}.smiles',
        headers={'User-Agent': 'Mozilla/5.0'}
    )
    try:
        with urlopen(request_url) as response:
            smiles = response.read()
    # If the accession is not found in the database then we stop here
    except HTTPError as error:
        # If the drugbank ID is not yet in the Drugbank references then return None
        if error.code == 404:
            return None
        else:
            print('Error when requesting ' + request_url)
            raise ValueError('Something went wrong with the Drugbank request (error ' + str(error.code) + ')')
    # This error may occur if there is no internet connection
    except URLError as error:
        print('Error when requesting ' + request_url)
        raise ValueError('Something went wrong with the MDposit request')

    return smiles



# Given a ChemBL ID, use the uniprot API to request its data and then mine what is needed for the database
def get_chembl_smiles (id_chembl : str) -> Optional[str]:
    # Request ChemBL
    parsed_response = None
    request_url = Request(
        url= f'https://www.ebi.ac.uk/chembl/interface_api/es_proxy/es_data/get_es_document/chembl_molecule/{id_chembl}',
        headers={'User-Agent': 'Mozilla/5.0'}
    )
    try:
        with urlopen(request_url) as response:
            parsed_response = json.loads(response.read().decode("utf-8"))
            smiles = parsed_response['_source']['molecule_structures']['canonical_smiles']
            pubchem_id = parsed_response['_source']['_metadata']['unichem'][8]['id']
    # If the accession is not found in ChemBL then the id is not valid
    except HTTPError as error:
        if error.code == 404:
            print('WARNING: Cannot find ChemBL entry for accession ' + id_chembl)
            return None
        print('Error when requesting ' + request_url)
        raise ValueError('Something went wrong with the ChemBL request (error ' + str(error.code) + ')')
    # If we have not a response at this point then it may mean we are trying to access an obsolete entry (e.g. P01607)
    if parsed_response == None:
        print('WARNING: Cannot find ChemBL entry for accession ' + id_chembl)
        return None
    return smiles

# Given a PubChem ID, use the uniprot API to request its data and then mine what is needed for the database
def get_pubchem_data (id_pubchem : str) -> Optional[dict]:
    # Request PubChem
    parsed_response = None
    request_url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{id_pubchem}/JSON/'
    try:
        with urlopen(request_url) as response:
            #parsed_response = json.loads(response.read().decode("windows-1252"))
            parsed_response = json.loads(response.read().decode("utf-8", errors='ignore'))
    # If the accession is not found in PubChem then the id is not valid
    # This may happen with pubchem ids of non-discrete compounds (e.g. 483927498)
    except HTTPError as error:
        if error.code == 404:
            print('WARNING: Cannot find PubChem entry for accession ', id_pubchem)
            return None
        print('Error when requesting ', request_url)
        raise ValueError('Something went wrong with the PubChem request (error ', str(error.code), ')')
    # If we have not a response at this point then it may mean we are trying to access an obsolete entry (e.g. P01607)
    if parsed_response == None:
        print('WARNING: Cannot find PubChem entry for accession ' + id_pubchem)
        return None
    # Mine target data: SMILES
    record = parsed_response.get('Record', None)
    if record == None:
        raise RuntimeError('Wrong Pubchem data structure: no record')
    sections = record.get('Section', None)
    if sections == None:
        raise RuntimeError('Wrong Pubchem data structure: no sections')
    names_and_ids_section = next((section for section in sections if section.get('TOCHeading', None) == 'Names and Identifiers'), None)
    if names_and_ids_section == None:
        raise RuntimeError('Wrong Pubchem data structure: no name and ids section')
    names_and_ids_subsections = names_and_ids_section.get('Section', None)
    if names_and_ids_subsections == None:
        raise RuntimeError('Wrong Pubchem data structure: no name and ids subsections')

    # Mine the name
    synonims = next((s for s in names_and_ids_subsections if s.get('TOCHeading', None) == 'Synonyms'), None)
    if synonims == None:
        descriptors = next((s for s in names_and_ids_subsections if s.get('TOCHeading', None) == 'Computed Descriptors'), None)
        descriptors_subsections = descriptors.get('Section', None)
        if descriptors_subsections == None:
            raise RuntimeError('Wrong Pubchem data structure: no name and ids subsections')
        depositor_supplied_descriptors = next((s for s in descriptors_subsections if s.get('TOCHeading', None) == 'IUPAC Name'), None)
        name_substance = depositor_supplied_descriptors.get('Information', None)[0].get('Value', {}).get('StringWithMarkup', None)[0].get('String', None)
    else:
        synonims_subsections = synonims.get('Section', None)
        if synonims_subsections == None:
            raise RuntimeError('Wrong Pubchem data structure: no name and ids subsections')
        depositor_supplied_synonims = next((s for s in synonims_subsections if s.get('TOCHeading', None) == 'Depositor-Supplied Synonyms'), None)
        if depositor_supplied_synonims == None:
            removed_synonims = next((s for s in synonims_subsections if s.get('TOCHeading', None) == 'Removed Synonyms'), None)
            name_substance = removed_synonims.get('Information', None)[0].get('Value', {}).get('StringWithMarkup', None)[0].get('String', None)
        else:
            name_substance = depositor_supplied_synonims.get('Information', None)[0].get('Value', {}).get('StringWithMarkup', None)[0].get('String', None)
        
    # Mine the SMILES
    computed_descriptors_subsection = next((s for s in names_and_ids_subsections if s.get('TOCHeading', None) == 'Computed Descriptors'), None)
    if computed_descriptors_subsection == None:
        raise RuntimeError('Wrong Pubchem data structure: no computeed descriptors')
    canonical_smiles_section = computed_descriptors_subsection.get('Section', None)
    if canonical_smiles_section == None:
        raise RuntimeError('Wrong Pubchem data structure: no canonical SMILES section')
    canonical_smiles = next((s for s in canonical_smiles_section if s.get('TOCHeading', None) == 'Canonical SMILES'), None)
    if canonical_smiles == None:
        raise RuntimeError('Wrong Pubchem data structure: no canonical SMILES')
    smiles = canonical_smiles.get('Information', None)[0].get('Value', {}).get('StringWithMarkup', None)[0].get('String', None)
    if smiles == None:
        raise RuntimeError('Wrong Pubchem data structure: no SMILES')

    # Mine target data: MOLECULAR FORMULA
    molecular_formula_subsection = next((s for s in names_and_ids_subsections if s.get('TOCHeading', None) == 'Molecular Formula'), None)
    if molecular_formula_subsection == None:
        raise RuntimeError('Wrong Pubchem data structure: no molecular formula section')
    molecular_formula = molecular_formula_subsection.get('Information', None)[0].get('Value', {}).get('StringWithMarkup', None)[0].get('String', None)
    if molecular_formula == None:
        raise RuntimeError('Wrong Pubchem data structure: no molecular formula')

    # Prepare the pubchem data to be returned
    return { 'name': name_substance, 'smiles': smiles, 'formula': molecular_formula }


def find_drugbank_pubchem (drugbank_id):
    # Request Drugbank
    request_url = Request(
        url= f'https://go.drugbank.com/drugs/{drugbank_id}',
        headers={'User-Agent': 'Mozilla/5.0'}
    )
    pubchem_id = None
    try:
        with urlopen(request_url) as response:
            content = response.read().decode("utf-8")
            pattern = re.compile("http\:\/\/pubchem.ncbi.nlm.nih.gov\/summary\/summary.cgi\?cid\=([0-9]*)")
            match = re.search(pattern, str(content))
            if not match:
                raise ValueError("No se encontró información sobre el Pubchem compound en esta página.")
            pubchem_id = match[1]
            if not pubchem_id:
                raise ValueError("No se encontró información sobre el Pubchem compound en esta página.")
    # If the accession is not found in the database then we stop here
    except HTTPError as error:
        # If the drugbank ID is not yet in the Drugbank references then return None
        raise ValueError(f'Wrong request. Code: {error.code}')
    # This error may occur if there is no internet connection
    except URLError as error:
        print('Error when requesting ' + request_url)
        raise ValueError('Something went wrong with the DrugBank request')
    
    return pubchem_id


def find_chembl_pubchem (id_chembl):
    # Request ChemBL
    parsed_response = None
    request_url = Request(
         url= f'https://www.ebi.ac.uk/chembl/interface_api/es_proxy/es_data/get_es_document/chembl_molecule/{id_chembl}',
         headers={'User-Agent': 'Mozilla/5.0'}
    )
    try:
        with urlopen(request_url) as response:
            parsed_response = json.loads(response.read().decode("utf-8"))
            unichem = parsed_response['_source']['_metadata']['unichem']
            if unichem == None:
                raise RuntimeError('Wrong Pubchem data structure: no unichem section')
            pubchem_section = next((section for section in unichem if section.get('src_url', None) == "http://pubchem.ncbi.nlm.nih.gov"), None)
            if pubchem_section == None:
                raise RuntimeError('Wrong Pubchem data structure: no pubchem section')
            pubchem_id = pubchem_section.get('id', None)
            if pubchem_id == None:
                raise RuntimeError('Wrong Pubchem data structure: no pubchem ID')
    # If the accession is not found in ChemBL then the id is not valid
    except HTTPError as error:
        if error.code == 404:
            print('WARNING: Cannot find ChemBL entry for accession ' + id_chembl)
            return None
        print('Error when requesting ' + request_url)
        raise ValueError('Something went wrong with the ChemBL request (error ' + str(error.code) + ')')
    # If we have not a response at this point then it may mean we are trying to access an obsolete entry (e.g. P01607)
    if parsed_response == None:
        print('WARNING: Cannot find ChemBL entry for accession ' + id_chembl)
        return None
    
    return pubchem_id

# Calculate the Morgan fingerprint and the Mordred descriptors from a ligand SMILES
def obtain_mordred_morgan_descriptors (smiles : str) -> Tuple:
    smiles = Chem.MolFromSmiles(smiles)

    # We can select the different submodules of mordred descriptors, avaible in: 'https://mordred-descriptor.github.io/documentation/master/'
    #calc = Calculator(descriptors, ignore_3D=True)
    calc = Calculator([descriptors.ABCIndex, 
                       descriptors.AcidBase.AcidicGroupCount, 
                       descriptors.AcidBase.BasicGroupCount, 
                       descriptors.RingCount], ignore_3D=True)

    #print(f"Number of descriptors in calculator: {len(calc.descriptors)}")

    # Calculate Mordred results
    mordred_results = calc(smiles).drop_missing().asdict()

    #print(f"Number of calculated descriptors: {len(results_dict)}")
    #print(results_dict)

    ######## MORGAN FINGERPRINT ###########

    morgan_fp = list(AllChem.GetMorganFingerprintAsBitVect(smiles, radius=2, nBits=1024))

    #print("Morgan fingerprint:")
    #print(list(morgan_fp))

    return mordred_results, morgan_fp

# Generate a map of residues associated to ligands
def generate_ligand_mapping (
    input_ligands : Optional[List[dict]],
    input_pdb_ids : List[str],
    structure : 'Structure',
    output_ligands_filepath : str,
    ) -> dict:

    # Merge input ligands and pdb ligands
    ligands = []
    if input_ligands:
        ligands += input_ligands
    # Get input ligands from the pdb ids, if any
    if input_pdb_ids:
        for pdb_id in input_pdb_ids:
            pubchem_ids = pdb_to_pubchem(pdb_id)
            for pubchem_id in pubchem_ids:
                ligands.append({ 'pubchem': pubchem_id, 'pdb': True })
    # If no input ligands are passed then stop here
    if len(ligands) == 0:
        return [], {}
    
    # Save data from all ligands to be saved in a file
    json_ligands_data = []
    ligands_data = []

    # Load the ligands file if exists already
    if os.path.exists(output_ligands_filepath):
        json_ligands_data += import_ligands(output_ligands_filepath)

    # Save apart ligand names forced by the user
    ligand_names = {}

    # Visited pubchem ids
    visited_pubchem_ids = []
    # Visited formulas
    visited_formulas = []
    # Save the maps of every ligand
    ligand_maps = []
    # Iterate input ligands
    for ligand in ligands:
        # If input ligand is not a dict but a single int/string then handle it
        if type(ligand) == int:
            print(f'A ligand number ID has been identified {ligand}, assuming that is a PubChem ID...')
            ligand = { 'pubchem': str(ligand) }
        elif type(ligand) == str:
            raise InputError(f'A name of ligand has been identified: {ligand}. Anyway, provide at least one of the following IDs: DrugBank, PubChem, ChEMBL.')
        # Check if we already have this ligand data
        ligand_data = obtain_ligand_data_from_file(json_ligands_data, ligand)
        # If we do not have its data yet then get data from pubchem
        if not ligand_data:
            ligand_data = obtain_ligand_data_from_pubchem(ligand)
            # Get morgan and mordred descriptors, the SMILES is needed for this part
            smiles = ligand_data['smiles']
            mordred_results, morgan_fp = obtain_mordred_morgan_descriptors(smiles)
            ligand_data['mordred'] = mordred_results
            ligand_data['morgan'] = morgan_fp
        # Add current ligand data to the general list
        ligands_data.append(ligand_data)
        # Get pubchem id
        pubchem_id = ligand_data['pubchem']
        # If we already visited this then skip it
        # Ligands in the structure (PDB) and the 'inputs.json' could be the same so it's not necessary to do it twice
        if pubchem_id in visited_pubchem_ids:
            continue
        # Add it to the list of visited ids
        visited_pubchem_ids.append(pubchem_id)
        # If we already visited a different ligand but with identical formula then we skip this ligand as well
        # Note that the mapping will be identical thus overwritting the previous map
        # However, ligands forced by the user are processed before so we keep them as priority
        formula = ligand_data['formula']
        if formula in visited_formulas:
            print(f'WARNING: Ligand with PubChem Id {pubchem_id} has a formula which has been already matched')
            ligands_data.pop()
            continue
        # Add it to the list of visited formulas
        visited_formulas.append(formula)
        # If the user defined a ligand name, it will be respectedand added to the metadata
        # if there isn't a defined name, it will be mined from PubChem
        user_forced_ligand_name = ligand.get('name', None)
        if user_forced_ligand_name:
            ligand_names[pubchem_id] = user_forced_ligand_name        
        # Map structure residues with the current ligand
        ligand_map = map_ligand_residues(structure, ligand_data)
        # If the ligand did not map then discard it
        if len(ligand_map['residue_indices']) == 0:
            if not ligand.get('pdb', False):
                raise InputError(f'Ligand with PubChem Id {pubchem_id} did not map with any residue')
            ligands_data.pop()
            continue
        # Add current ligand map to the general list
        ligand_maps.append(ligand_map)

    # Export ligands to a file
    save_json(ligands_data, output_ligands_filepath)

    # print('Ligands: ', ligands)
    # print('A total of',len(atom_elements_per_residue), 'residues and', len(atom_elements_per_ligand), 'ligands were found.')
    # print('Matched:', len(ligands),'/', len(atom_elements_per_ligand),'of the ligands.')

    return ligand_maps, ligand_names
    
# Set the expected ligand data fields
LIGAND_DATA_FIELDS = set(['name', 'pubchem', 'drugbank', 'chembl', 'smiles', 'formula', 'morgan', 'mordred'])

# Import ligands json file so we do not have to rerun this process
def import_ligands (output_ligands_filepath : str) -> dict:
    # Read the file
    imported_ligands = load_json(output_ligands_filepath)
    # Format data as the process expects to find it
    for imported_ligand in imported_ligands:
        for expected_field in LIGAND_DATA_FIELDS:
            if expected_field not in imported_ligand:
                imported_ligand[expected_field] = None
    return imported_ligands

# Check if the current input ligand is already among the ligands we already have data for
def obtain_ligand_data_from_file ( ligands_data : List[dict], ligand: dict ) -> Optional[dict]:
    for ligand_data in ligands_data:
        if (ligand.get('pubchem') and ligand['pubchem'] == ligand_data.get('pubchem')) or \
        (ligand.get('drugbank') and ligand['drugbank'] == ligand_data.get('drugbank')) or \
        (ligand.get('chembl') and ligand['chembl'] == ligand_data.get('chembl')):
            return ligand_data
    return None

# Given an input ligand, obtain all necessary data
def obtain_ligand_data_from_pubchem (ligand : dict) -> dict:
    # Save in a dictionary all ligand data including its name and ids
    # The ID can be of the databases: 'drugbank' , 'pubchem' , 'chembl'
    # Define the needed variables to check if the ligand has a database ID or it is None
    ligand_data = {
        'name': None,
        'pubchem': None,
        'drugbank': None,
        'chembl': None,
        'smiles': None,
        'formula': None
    }

    # Set ligand data pubchem id, even if the input id is not from pubhcme (e.g. drugbank, chembl)
    if 'pubchem' in ligand:
        ligand_data['pubchem'] = str(ligand.get('pubchem'))
    elif 'drugbank' in ligand:
        ligand_data['drugbank'] = ligand.get('drugbank')
        ligand_data['pubchem'] = str(find_drugbank_pubchem(ligand_data['drugbank']))
    elif 'chembl' in ligand:
        ligand_data['chembl'] = ligand.get('chembl')
        ligand_data['pubchem'] = str(find_chembl_pubchem(ligand_data['chembl']))
    else:
        raise InputError('None of the ligand IDs are defined. Please provide at least one of the following IDs: DrugBank, PubChem, ChEMBL.')
    
    # Request ligand data from pubchem
    pubchem_data = get_pubchem_data(ligand_data['pubchem'])
    if not pubchem_data:
        raise RuntimeError('No PubChem data avilable')

    # Add pubchem data to ligand data
    ligand_data = { **ligand_data, **pubchem_data }
    
    return ligand_data


# For each residue, count the number of atom elements
def count_atom_elements_per_residue ( structure : 'Structure' ) -> dict:
    # Obtain a list of residues of ligand candidates
    residue_element_count = {}
    # Iterate over the residue(s)
    for residue in structure.residues:
        # Skip protein residues
        if residue.classification == 'protein':
            continue
        # Obtain the list of elements for each residue
        atom_elements = [ atom.element for atom in residue.atoms ]
        # Generate a dict with elements (keys) and their counting (values)
        atoms_count = {}
        for atom in atom_elements:
            if atom in atoms_count:
                atoms_count[atom] += 1
            else:
                atoms_count[atom] = 1
        # Define the key of the dictionary as the residue name and save it in a list with the all the residues
        residue_element_count[residue] = atoms_count

    return residue_element_count


def nest_brackets(tokens, i = 0):
    l = []
    while i < len(tokens):
        if tokens[i] == ')':
            return i,l
        elif tokens[i] == '(':
            i,subl = nest_brackets(tokens, i+1)
            l.append(subl)
        else:
            l.append(tokens[i])
        i += 1
    return i,l

def parse_compound(formula : str) -> List[str]:
    tokens = [''.join(t) for t in split_when(formula, lambda a,b: b.isupper() or b in '()' or (b.isdigit() and not a.isdigit()))]
    tokens = [(int(t) if t.isdigit() else t) for t in tokens]
    i, l = nest_brackets(tokens)
    assert(i == len(tokens)) # crash if unmatched ')'
    return l

# Split a string using a function
def split_when(string : str, func : Callable) -> List[str]:
    splits = []
    last_split = string[0]
    for s in string[1:]:
        if func(last_split, s):
            splits.append(last_split)
            last_split = s
        else:
            last_split += s
    splits.append(last_split)
    return splits

# Given a chemical formula, get the count of atoms per element
def count_atom_elements (molecular_formula : str) -> dict:
    parsed_molecular_formula = molecular_formula.replace("+", "").replace("-", "")
    l = parse_compound(parsed_molecular_formula)
    c = parse_splits(l)
    return c

# This function associates elements in a list
# If a string is followed by a number then they go together
# If a string has no number then the number is 1 for this specific string
def parse_splits (splits : List[str]) -> dict:
    parsed = {}
    previous_key = None
    for split in splits:
        if type(split) == str:
            if previous_key:
                parsed[previous_key] = 1
            previous_key = split
        elif type(split) == int:
            parsed[previous_key] = split
            previous_key = None
        else:
            raise ValueError('Not supported type ' + type(split))
    if previous_key:
        parsed[previous_key] = 1
    return parsed


def match_ligandsID_to_res (ligand_atom_element_count : dict, residue_atom_element_count : dict) -> bool:
        # Remove hydrogens since we only pay attention to heavy atoms
        # This is to avoid having mismatched_residues due to different protonation states
        if 'H' in ligand_atom_element_count:
            del ligand_atom_element_count['H']
        if 'H' in residue_atom_element_count:
            del residue_atom_element_count['H']
        #shared_items = {k: ligand_atom_element_count[k] for k in ligand_atom_element_count if k in residue_atom_element_count and ligand_atom_element_count[k] == residue_atom_element_count[k]}
        ligand_elements = set(ligand_atom_element_count.keys())
        residue_elements = set(residue_atom_element_count.keys())

        # Check if there are different elements
        # If so, we are done
        different_keys = ligand_elements.symmetric_difference(residue_elements)
        if len(different_keys) > 0:
        #    print("Elements that are not matching:")
        #     for atom in different_keys:
        #         if atom in ligand_atom_element_count:
        #             print(f"In Pubchem {list(ligand_atom_element_count.keys())[0]} : {atom}: {ligand_atom_element_count[atom]}")
        #         else:
        #             print(f"In PDB residue {list(residue_atom_element_count.keys())[0]}: {atom}: {residue_atom_element_count[atom]}")
            return False

        # Different values
        different_values = []
        for element in ligand_elements:
            if ligand_atom_element_count[element] != residue_atom_element_count[element]:
                different_values.append(element)

        # If the count of elements does not match then stop here
        if len(different_values) > 0:
        #     Print the differences found between the elements
        #     print("Number of atoms that are different:")
        #     for atom in different_values:
        #         print(f"In Pubchem {list(ligand_atom_element_count.keys())[0]} {atom}: {ligand_atom_element_count[atom]}, In PDB residue {list(pdb_dict.keys())[0]}: {atom}: {residue_atom_element_count[atom]}\n")
            return False
        
        return True

def map_ligand_residues (structure : 'Structure', ligand_data : dict) -> dict:
    # Load the structure where the atoms will be minresidue.indexed and obtain a residues dict
    atom_elements_per_residue = count_atom_elements_per_residue(structure)
    # Get the ligand formula
    ligand_formula = ligand_data['formula']
    if not ligand_formula:
        raise RuntimeError(f'Ligand with PubChem id {ligand_data["pubchem"]} is missing formula')
    # From pubchem mine the atoms of the ligands and save it in a dict
    ligand_atom_element_count = count_atom_elements(ligand_data['formula'])
    matched_residues = []
    # Get ligand pubchem id
    pubchem_id = ligand_data['pubchem']
    for residue, residue_atom_element_count in atom_elements_per_residue.items():
        if match_ligandsID_to_res(ligand_atom_element_count, residue_atom_element_count):
            matched_residues.append(residue.index)
    
    # Format the output as we expect it
    ligand_map = { 'name': pubchem_id, 'residue_indices': matched_residues, 'match': { 'ref': { 'pubchem': pubchem_id } } }
    return ligand_map

# Given a smiles, get the pubchem id
# e.g. CC1=C(C=NC=C1)NC(=O)CC2=CC(=CC(=C2)Cl)OC -> 154876006
# e.g. O=C4N3C(C(=O)Nc2cc(nn2c1ccccc1)C)C(SC3CC=CC4NC(=O)C(NC)C)(C)C -> None
def smiles_to_pubchem_id (smiles : str) -> Optional[str]:
    # Set the request URL
    request_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/JSON'
    # Set the POST data
    data = urlencode({ 'smiles': smiles }).encode()
    try:
        with urlopen(request_url, data=data) as response:
            parsed_response = json.loads(response.read().decode("utf-8"))
    # If the smiles is not found in pubchem then we can stop here
    except HTTPError as error:
        if error.code == 404:
            print(f' Smiles {smiles} not found')
            return None
        else:
            raise ValueError('Something went wrong with the PubChem request: ' + request_url)
    # Get the PubChem id
    compounds = parsed_response.get('PC_Compounds', None)
    if not compounds:
        raise RuntimeError(f'Something went wrong when mining pubchem data for SMILES {smiles}')
    # Make sure there is only 1 compound
    # DANI: Esto algún día podría no ser así, ya lo gestionaremos
    if len(compounds) != 1:
        raise RuntimeError(f'There should be one and only one compound for SMILES {smiles}')
    # Keep mining
    compound = compounds[0]
    first_id = compound.get('id', None)
    # If there is no first id then it means there is no direct match with pubchem
    if not first_id:
        return None
    second_id = first_id.get('id', None)
    if not second_id:
        raise RuntimeError(f'Missing second id when mining pubchem data for SMILES {smiles}')
    pubchem_id = second_id.get('cid', None)
    if not pubchem_id:
        raise RuntimeError(f'Missing pubchem id when mining pubchem data for SMILES {smiles}')
    return str(pubchem_id)

# Given a pdb Id, get its ligand smiles
# e.g. 4YDF -> O=S(=O)(N(C3C(N(Cc1ccccc1)S(=O)(=O)c2ccc([N+]([O-])=O)cc2)CNC3)Cc4ccccc4)c5cc([N+]([O-])=O)ccc5, [O-]S([O-])(=O)=O
def pdb_to_smiles (pdb_id : str) -> Optional[ List[str] ]:
    # Request the MMB service to retrieve pdb data
    request_url = f'http://mdb-login.bsc.es/api/pdb/{pdb_id}/entry'
    try:
        with urlopen(request_url) as response:
            parsed_response = json.loads(response.read().decode("utf-8"))
    # If the accession is not found in the PDB then we can stop here
    except HTTPError as error:
        if error.code == 404:
            print(f' PDB code {pdb_id} not found')
            return None
        else:
            raise ValueError('Something went wrong with the PDB request: ' + request_url)
    # Get the 'het' hits or ligands available
    smiles_list = []
    hets = parsed_response.get('hets', None)
    if not hets:
        raise RuntimeError('Failed to mine pdb data: missing hets section')
    # If there are more than one ligand, iterate over them
    for het in hets:
        descriptor = het.get('pdbx_chem_comp_descriptor', None)
        if descriptor == None:
            raise RuntimeError('Failed to mine pdb data: missing pdbx_chem_comp_descriptor section')
        depositor_descriptors = next((s for s in descriptor if s.get('type', None) == 'SMILES'), None)
        if depositor_descriptors == None:
            raise RuntimeError('Failed to mine pdb data: missing SMILES section')
        smiles = depositor_descriptors.get('descriptor', None)
        smiles_list.append(smiles)

    return smiles_list

# Given a PDB ligand code, get its pubchem
# DANI: Me rindo, he perdido demasiado tiempo con esto
# DANI: La api del pdb solo me devuelve errores 400 que no me dicen cual es el problema
def pdb_ligand_to_pubchem (pdb_ligand_id : str) -> str:
    # Set the request URL
    request_url = 'https://data.rcsb.org/graphql'
    request = Request(url = request_url, headers={
        'Authorization': 'Bearer 5959649b3b067a55a3c1ffad',
        'contentType': 'application/json'
    })
    # Set the POST data
    data = urlencode({
        "operationName": "molecule",
        #'query': 'query molecule($id:String!){chem_comp(comp_id:$id){rcsb_chem_comp_related{resource_name resource_accession_code}}}',
        "query": "query molecule($id:String!){\n chem_comp(comp_id:$id) {\n rcsb_chem_comp_related{\n resource_name\n resource_accession_code\n }\n }\n}\n",
        "variables": { "id": pdb_ligand_id }
    }).encode()
    # Run the query
    parsed_response = None
    try:
        with urlopen(request, data=data) as response:
            parsed_response = json.loads(response.read().decode("utf-8"))
    # If the accession is not found in the PDB then we can stop here
    except HTTPError as error:
        if error.code == 404:
            print(f' PDB ligand {pdb_ligand_id} not found')
            return None
        else:
            print(error.msg)
            raise ValueError('Something went wrong with the PDB ligand request: ' + request_url)
    print(parsed_response)

# Given a PDB ligand code, get its pubchem
# Use a web crawler to avoid having to use the PDB API
def pdb_ligand_to_pubchem_RAW (pdb_ligand_id : str) -> Optional[str]:
    # Set the request URL
    request_url = f'https://www.rcsb.org/ligand/{pdb_ligand_id}'
    # Run the query
    parsed_response = None
    try:
        with urlopen(request_url) as response:
            parsed_response = response.read().decode("utf-8")
    # If the accession is not found in the PDB then we can stop here
    except HTTPError as error:
        if error.code == 404:
            print(f' PDB ligand {pdb_ligand_id} not found')
            return None
        else:
            print(error.msg)
            raise ValueError('Something went wrong with the PDB ligand request: ' + request_url)
    # Mine the pubchem id out of the whole response
    pattern = re.compile('pubchem.ncbi.nlm.nih.gov\/compound\/([0-9]*)\"')
    match = re.search(pattern, parsed_response)
    # If there is no pubchem id then return none
    # This is normal for some ligands such as counter ions (e.g. LI)
    if not match:
        return None
    pubchem_id = match[1]
    return pubchem_id

# Given a PDB code, get all its ligand codes
def get_pdb_ligand_codes (pdb_id : str) -> List[str]:
    # Request the MMB service to retrieve pdb data
    request_url = f'http://mdb-login.bsc.es/api/pdb/{pdb_id}/entry'
    try:
        with urlopen(request_url) as response:
            parsed_response = json.loads(response.read().decode("utf-8"))
    # If the accession is not found in the PDB then we can stop here
    except HTTPError as error:
        if error.code == 404:
            print(f' PDB code {pdb_id} not found')
            return None
        else:
            raise ValueError('Something went wrong with the PDB request: ' + request_url)
    # Get the 'het' hits or ligands available
    ligand_codes = []
    hets = parsed_response.get('hets', None)
    # If the PDB has no ligands at all then there is no 'hets' section
    if not hets:
        return []
    # If there are more than one ligand, iterate over them
    for het in hets:
        ligand_code = het.get('_id', None)
        if ligand_code == None:
            raise RuntimeError('Failed to mine pdb data: missing id')
        ligand_codes.append(ligand_code)
    return ligand_codes

# Given a pdb Id, get its pubchem ids
# e.g. 4YDF -> 
def pdb_to_pubchem (pdb_id : str) -> List[str]:
    print(f'Searching PubChem Ids for PDB {pdb_id}')
    pubchem_ids = []
    # Iterate over pdb ligand codes
    ligand_codes = get_pdb_ligand_codes(pdb_id)
    for ligand_code in ligand_codes:
        pubchem_id = pdb_ligand_to_pubchem_RAW(ligand_code)
        if pubchem_id:
            pubchem_ids.append(pubchem_id)
    print(' Found ' + ', '.join(pubchem_ids))
    return pubchem_ids
    # DANI: De momento no usamos las SMILES que alguna vez me han dado problemas (e.g. 2I3I)
    # smiles_list = pdb_to_smiles(pdb_id)
    # # Iterate over the diferents SMILES to obtain the pubchem ID 
    # for smiles in smiles_list:
    #     pubchem_id = smiles_to_pubchem_id(smiles)
    #     if not pubchem_id:
    #         raise RuntimeError(f'SMILES {smiles} has no PubChem id')
    #     print(f' Found {pubchem_id}')
    #     pubchem_ids.append(pubchem_id)
    # return pubchem_ids