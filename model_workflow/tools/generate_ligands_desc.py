import os
import sys
import json
import urllib.request
from rdkit import Chem
from rdkit.Chem import AllChem
from mordred import Calculator, descriptors
from typing import List, Tuple, Optional, Callable

from model_workflow.utils.auxiliar import InputError, load_json, save_json
from urllib.request import Request, urlopen
import re

def get_drugbank_smiles (id_drugbank : str) -> Optional[str]:
    # Request Drugbank
    request_url = Request(
         url= f'https://go.drugbank.com/structures/small_molecule_drugs/{id_drugbank}.smiles', # f'https://go.drugbank.com/structures/small_molecule_drugs/{id_drugbank}.smiles', 'https://pubchem.ncbi.nlm.nih.gov/compound/1986#section=Canonical-SMILES&fullscreen=true'
         headers={'User-Agent': 'Mozilla/5.0'}
    )
    try:
        with urlopen(request_url) as response:
            smiles = response.read()
    # If the accession is not found in the database then we stop here
    except urllib.error.HTTPError as error:
        # If the drugbank ID is not yet in the Drugbank references then return None
        if error.code == 404:
            return None
        else:
            print('Error when requesting ' + request_url)
            raise ValueError('Something went wrong with the Drugbank request (error ' + str(error.code) + ')')
    # This error may occur if there is no internet connection
    except urllib.error.URLError as error:
        print('Error when requesting ' + request_url)
        raise ValueError('Something went wrong with the MDposit request')

    return smiles



# Given a ChemBL ID, use the uniprot API to request its data and then mine what is needed for the database
def get_chembl_smiles (id_chembl : str) -> Optional[str]:
    # Request ChemBL
    parsed_response = None
    request_url = Request(
         url= f'https://www.ebi.ac.uk/chembl/interface_api/es_proxy/es_data/get_es_document/chembl_molecule/{id_chembl}', # f'https://go.drugbank.com/structures/small_molecule_drugs/{id_drugbank}.smiles', 'https://pubchem.ncbi.nlm.nih.gov/compound/1986#section=Canonical-SMILES&fullscreen=true'
         headers={'User-Agent': 'Mozilla/5.0'}
    )
    try:
        with urlopen(request_url) as response:
            parsed_response = json.loads(response.read().decode("utf-8"))
            smiles = parsed_response['_source']['molecule_structures']['canonical_smiles']
            pubchem_id = parsed_response['_source']['_metadata']['unichem'][8]['id']
    # If the accession is not found in ChemBL then the id is not valid
    except urllib.error.HTTPError as error:
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

def purgeNonUtf8 (filename : str):
    with open(filename, mode="r+", encoding="utf-8", errors= 'ignore') as file:
        lines = file.readlines()
        file.seek(0)
        for line in lines:
            file.write(line)
        file.truncate()


# Given a PUBChem ID, use the uniprot API to request its data and then mine what is needed for the database
def get_pubchem_data (id_pubchem : str) -> Optional[str]:
    # Request PUBChem
    parsed_response = None
    request_url = Request(
        url= f'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{id_pubchem}/JSON/', #'https://pubchem.ncbi.nlm.nih.gov/compound/1986#section=Canonical-SMILES&fullscreen=true'
        headers={'User-Agent': 'Mozilla/5.0'}
    )
    try:
        with urlopen(request_url) as response:
            #parsed_response = json.loads(response.read().decode("windows-1252"))
            parsed_response = json.loads(response.read().decode("utf-8", errors='ignore'))
    # If the accession is not found in PUBChem then the id is not valid
    except urllib.error.HTTPError as error:
        if error.code == 404:
            print('WARNING: Cannot find PUBChem entry for accession ', id_pubchem)
            return None
        print('Error when requesting ', request_url)
        raise ValueError('Something went wrong with the PUBChem request (error ', str(error.code), ')')
    # If we have not a response at this point then it may mean we are trying to access an obsolete entry (e.g. P01607)
    if parsed_response == None:
        print('WARNING: Cannot find PUBChem entry for accession ' + id_pubchem)
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
    synonims_subsections = synonims.get('Section', None)
    if synonims_subsections == None:
        raise RuntimeError('Wrong Pubchem data structure: no name and ids subsections')
    depositor_supplied_synonims = next((s for s in synonims_subsections if s.get('TOCHeading', None) == 'Depositor-Supplied Synonyms'), None)
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

    return name_substance, smiles, molecular_formula


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
    except urllib.error.HTTPError as error:
        # If the drugbank ID is not yet in the Drugbank references then return None
        raise ValueError(f'Wrong request. Code: {error.code}')
    # This error may occur if there is no internet connection
    except urllib.error.URLError as error:
        print('Error when requesting ' + request_url)
        raise ValueError('Something went wrong with the DrugBank request')
    
    return pubchem_id


def find_chembl_pubchem (id_chembl):
    # Request ChemBL
    parsed_response = None
    request_url = Request(
         url= f'https://www.ebi.ac.uk/chembl/interface_api/es_proxy/es_data/get_es_document/chembl_molecule/{id_chembl}', # f'https://go.drugbank.com/structures/small_molecule_drugs/{id_drugbank}.smiles', 'https://pubchem.ncbi.nlm.nih.gov/compound/1986#section=Canonical-SMILES&fullscreen=true'
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
    except urllib.error.HTTPError as error:
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
    structure : 'Structure',
    output_ligands_filepath : str,
    ) -> dict:

    # If no input ligands are passed then stop here
    if not input_ligands:
        return []
    
    # Save data from all ligands to be saved in a file
    json_ligands_data = []
    ligands_data = []

    # Load the ligands file if exists already
    if os.path.exists(output_ligands_filepath):
        json_ligands_data += import_ligands(output_ligands_filepath)

    # Save apart ligand names forced by the user
    ligand_names = {}

    # Save the maps of every ligand
    ligand_maps = []
    # Iterate input ligands
    for input_ligand in input_ligands:
        # If input ligand is not a dict but a single int/string then handle it
        if type(input_ligand) == int:
            print(f'A ligand number ID has been identified {input_ligand}, assuming that is a PubChem ID...')
            input_ligand = { 'pubchem': str(input_ligand) }
        elif type(input_ligand) == str:
            raise InputError(f'A name of ligand has been identified: {input_ligand}. Anyway, provide at least one of the following IDs: DrugBank, PubChem, ChEMBL.')
        # Check if we already have this ligand data
        ligand_data = obtain_ligand_data_from_file(json_ligands_data, input_ligand)
        # If we do not have its data yet then get data from pubchem
        if not ligand_data:
            ligand_data = obtain_ligand_data_from_pubchem(input_ligand)
            # Get morgan and mordred descriptors, the SMILES is needed for this part
            smiles = ligand_data['smiles']
            mordred_results, morgan_fp = obtain_mordred_morgan_descriptors(smiles)
            ligand_data['mordred'] = mordred_results
            ligand_data['morgan'] = morgan_fp
        # Add current ligand data to the general list
        ligands_data.append(ligand_data)
        # Get pubchem id
        pubchem_id = ligand_data['pubchem']
        # If the user defined a ligand name, it will be respectedand added to the metadata
        # if there isn't a defined name, it will be mined from PubChem
        user_forced_ligand_name = input_ligand.get('name', None)
        if user_forced_ligand_name:
            ligand_names[pubchem_id] = user_forced_ligand_name        
        # Map structure residues with the current ligand
        ligand_map = map_ligand_residues(structure, ligand_data)
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
    print('imported ', imported_ligands)
    return imported_ligands

# Check if the current input ligand is already among the ligands we already have data for
def obtain_ligand_data_from_file ( ligands_data : List[dict], input_ligand: dict ) -> Optional[dict]:
    for ligand_data in ligands_data:
        if (input_ligand.get('pubchem') and input_ligand['pubchem'] == ligand_data.get('pubchem')) or \
        (input_ligand.get('drugbank') and input_ligand['drugbank'] == ligand_data.get('drugbank')) or \
        (input_ligand.get('chembl') and input_ligand['chembl'] == ligand_data.get('chembl')):
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

    # Set ligand data pubchem id, even if the input id is not from pubhcme (e.g. drugbakn, chembl)
    if 'pubchem' in ligand:
        ligand_data['pubchem'] = ligand.get('pubchem')
    elif 'drugbank' in ligand:
        ligand_data['drugbank'] = ligand.get('drugbank')
        ligand_data['pubchem'] = find_drugbank_pubchem(ligand_data['drugbank'])
    elif 'chembl' in ligand:
        ligand_data['chembl'] = ligand.get('chembl')
        ligand_data['pubchem'] = find_chembl_pubchem(ligand_data['chembl'])
    else:
        raise InputError('None of the ligand IDs are defined. Please provide at least one of the following IDs: DrugBank, PubChem, ChEMBL.')
 
    # Request ligand data from pubchem
    ligand_name, smiles, molecular_formula = get_pubchem_data(ligand_data['pubchem'])
    ligand_data['name'] = ligand_name
    ligand_data['smiles'] = smiles
    ligand_data['formula'] = molecular_formula
    #ligand_data[ligand_name] = ligand_data['pubchem'] # DANI: no se que hace
    
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
    # From pubchem mine the atoms of the ligands and save it in a dict
    ligand_atom_element_count = count_atom_elements(ligand_data['formula'])
    #print('Ligands: ', atom_elements_per_ligand)
    matched_residues = []
    # Get ligand pubchem id
    pubchem_id = ligand_data['pubchem']
    for residue, residue_atom_element_count in atom_elements_per_residue.items():
        if match_ligandsID_to_res(ligand_atom_element_count, residue_atom_element_count):
            matched_residues.append(residue.index)
    
    # Format the output as we expect it
    ligand_map = { 'name': pubchem_id, 'residue_indices': matched_residues, 'match': { 'ref': { 'pubchem': pubchem_id } } }
    return ligand_map