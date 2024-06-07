import os
import sys
import json
import urllib.request
from rdkit import Chem
from rdkit.Chem import AllChem
from mordred import Calculator, descriptors
from typing import List, Tuple, Optional, Union, Callable

from pathlib import Path
from model_workflow.tools.residues_library import residue_name_2_letter
from model_workflow.utils.auxiliar import load_json, save_json
from urllib.request import Request, urlopen
from itertools import chain
from collections import Counter
from model_workflow.utils.constants import STRUCTURE_FILENAME, OUTPUT_LIGANDS_FILENAME
import re

# Import ligands json file so we do not have to rerun this process
def import_ligands () -> dict:
    # Read the file
    file_content = load_json(OUTPUT_LIGANDS_FILENAME)
    # Format data as the process expects to find it
    ligands = {}
    drugbank_id = None
    chembl_id = None
    pubchem_id = None
    for ligand in file_content:

        if 'pubchem' in ligand:
            pubchem_id = ligand['pubchem']
        else:
            pubchem_id = None
        
        if 'name' in ligand:
            name = ligand['name']
        else:
            name = None

        if 'drugbank' in ligand:
            drugbank_id = ligand['drugbank']
        else:
            drugbank_id = None
        
        if 'chembl' in ligand:
            chembl_id = ligand['chembl']
        else:
            chembl_id = None
        
        ligands[name] = {'name': name,
                         'pubchem': pubchem_id,
                         'drugbank': drugbank_id,
                         'chembl': chembl_id
                        }
    return ligands, file_content


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
            parsed_response = json.loads(response.read().decode("utf-8"))
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

    return smiles, molecular_formula

def get_pubchem_name (id_pubchem : str) -> Optional[str]:
    # Request PUBChem
    parsed_response = None
    request_url = Request(
        url= f'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{id_pubchem}/JSON/', #'https://pubchem.ncbi.nlm.nih.gov/compound/1986#section=Canonical-SMILES&fullscreen=true'
        headers={'User-Agent': 'Mozilla/5.0'}
    )
    try:
        with urlopen(request_url) as response:
            #parsed_response = json.loads(response.read().decode("windows-1252"))
            parsed_response = json.loads(response.read().decode("utf-8"))
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
    # Mine target data: NAME
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
    computed_descriptors_subsection = next((s for s in names_and_ids_subsections if s.get('TOCHeading', None) == 'Synonyms'), None)
    names_and_ids_subsubsections = computed_descriptors_subsection.get('Section', None)
    if names_and_ids_subsubsections == None:
        raise RuntimeError('Wrong Pubchem data structure: no name and ids subsections')
    computed_descriptors_subsection = next((s for s in names_and_ids_subsubsections if s.get('TOCHeading', None) == 'Depositor-Supplied Synonyms'), None)
    name_substance = computed_descriptors_subsection.get('Information', None)[0].get('Value', {}).get('StringWithMarkup', None)[0].get('String', None)

    return name_substance


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
            print("Pubchem compound:", pubchem_id)
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

def obtain_ligands_id (input_ligands : list, ligands_to_include) -> dict:


    pubchem_id_list = []
    pubchem_name_dic = {}

    pubchem_dic_list = []
    ligands_to_add = []
    # for id in ligands_to_include:
    #     for ligand in input_ligands:
    #         if ligand.get('pubchem'):
    #             if ligand['pubchem'] == id:
    #                 ligands_to_add.append(ligand)
    #         else:
    #             continue
    # Save in a dictionary the name of the ligand as key and the ID as value
    # The ID can be of the databases: 'drugbank' , 'pubchem' , 'chembl'
    # Also, generate a list with all the pubchems IDs which will be used
    for ligand in ligands_to_include:
        drugbank_dic = {}
        pubchem_dic = {}
        chembl_dic = {}
        # Define the needed variables to check if the ligand has a database ID or it is None
        drugbank_id = None
        pubchem_id = None
        chembl_id = None
        ligand_name = None

        pubchem_dic['name'] = ligand_name
        pubchem_dic['drugbank'] = drugbank_id
        pubchem_dic['chembl'] = chembl_id
        pubchem_dic['pubchem'] = pubchem_id

        # If the user defined a ligand name, it will be respectedand added to the metadata 
        # if there isn't a defined name, it will be mined from PubChem
        if 'name' in ligand and 'pubchem' in ligand:
            ligand_name_metadata = ligand.get('name') # ligand.get('pubchem'): 
            ligand_name = get_pubchem_name(ligand.get('pubchem'))
            ligand_customized_name = ligand.get('name')
            pubchem_id = ligand.get('pubchem')
            pubchem_name_dic[pubchem_id] = ligand_name_metadata
            pubchem_dic[ligand_name] = pubchem_id
            pubchem_dic['name'] = ligand_name
            pubchem_dic['pubchem'] = pubchem_id
            pubchem_id_list.append(pubchem_id)
            pubchem_dic_list.append(pubchem_dic)
            print(f'A customized ligand name has been identified: {ligand_customized_name}. Adding to metadata.')
            print(f'Retrieved name of ligand ID {ligand}: {ligand_name}.')

        if 'pubchem' in ligand and not 'name' in ligand:
            ligand_name = get_pubchem_name(ligand.get('pubchem'))
            pubchem_id = ligand.get('pubchem')
            pubchem_dic[ligand_name] = pubchem_id
            pubchem_dic['name'] = ligand_name
            pubchem_dic['pubchem'] = pubchem_id
            pubchem_id_list.append(pubchem_id)
            pubchem_dic_list.append(pubchem_dic)
            print(f'A ligand number ID has been identified {ligand} without associated name, assuming that is a PubChem ID...')
            print(f'Retrieved name of ligand ID {ligand}: {ligand_name}.')

        if 'drugbank' in ligand and 'pubchem' in ligand:
            continue
        
        if 'drugbank' in ligand and 'name' in ligand:
            pubchem_id = find_drugbank_pubchem(ligand.get('drugbank'))
            ligand_name_metadata = ligand.get('name') # ligand.get('pubchem'): 
            ligand_name = get_pubchem_name(pubchem_id)
            pubchem_name_dic[pubchem_id] = ligand_name_metadata
            drugbank_id = ligand.get('drugbank')
            pubchem_dic[ligand_name] = pubchem_id
            pubchem_dic['drugbank'] = drugbank_id
            pubchem_dic['pubchem'] = pubchem_id
            pubchem_dic['name'] = ligand_name
            pubchem_dic_list.append(pubchem_dic)
            pubchem_id_list.append(pubchem_id)
        
        if 'drugbank' in ligand and not 'name' in ligand:
            pubchem_id = find_drugbank_pubchem(ligand.get('drugbank'))
            ligand_name = get_pubchem_name(pubchem_id)
            pubchem_dic[ligand_name] = pubchem_id
            pubchem_dic['pubchem'] = pubchem_id
            pubchem_dic['name'] = ligand_name
            pubchem_dic['drugbank'] = drugbank_id
            pubchem_id_list.append(pubchem_id)
            pubchem_dic_list.append(pubchem_dic)
            print(f'A ligand number ID has been identified {ligand} without associated name...')
            print(f'Retrieved name of ligand ID {ligand}: {ligand_name}.')

        if 'chembl' in ligand and 'name' in ligand:
            pubchem_id = find_chembl_pubchem(ligand.get('chembl'))
            ligand_name_metadata = ligand.get('name') # ligand.get('pubchem'): 
            pubchem_name_dic[pubchem_id] = ligand_name_metadata
            ligand_name = get_pubchem_name(pubchem_id)
            chembl_id = ligand.get('chembl')
            pubchem_dic[ligand_name] = pubchem_id
            pubchem_dic['chembl'] = chembl_id
            pubchem_dic['pubchem'] = pubchem_id
            pubchem_dic['name'] = ligand_name
            pubchem_dic_list.append(pubchem_dic)
            pubchem_id_list.append(pubchem_id)
        
        if 'chembl' in ligand and not 'name' in ligand:
            chembl_id = ligand.get('chembl')
            pubchem_id = find_chembl_pubchem(chembl_id)
            ligand_name = get_pubchem_name(pubchem_id)
            pubchem_dic[ligand_name] = pubchem_id
            pubchem_dic['chembl'] = chembl_id
            pubchem_dic['pubchem'] = pubchem_id
            pubchem_dic['name'] = ligand_name
            pubchem_dic_list.append(pubchem_dic)
            pubchem_id_list.append(pubchem_id)

        if (drugbank_id is None) and (pubchem_id is None) and (chembl_id is None):
            print('The ligand name or ID is confusing.')
            if type(ligand) == int:
                print(f'A ligand number ID has been identified {ligand}, assuming that is a PubChem ID...')
                pubchem_id_list.append(ligand)
                ligand_name = get_pubchem_name(ligand)
                pubchem_dic[ligand_name] = ligand
                print(f'Retrieved name of ligand ID {ligand}: {ligand_name}.')
            elif type(ligand) == str:
                raise ValueError(f'A name of ligand has been identified: {ligand}. Anyway, provide at least one of the following IDs: DrugBank, PubChem, ChEMBL.')
            else:
                raise ValueError('None of the ligand IDs are defined. Please provide at least one of the following IDs: DrugBank, PubChem, ChEMBL.')
    
    return drugbank_dic, pubchem_dic_list, chembl_dic, pubchem_id_list, pubchem_name_dic


def obtain_mordred_morgan_descriptors (smiles : str) -> dict:
    smiles = Chem.MolFromSmiles(smiles)

    # We can select the different submodules of mordred descriptors, avaible in: 'https://mordred-descriptor.github.io/documentation/master/'
    #calc = Calculator(descriptors, ignore_3D=True)
    calc = Calculator([descriptors.ABCIndex, 
                       descriptors.AcidBase.AcidicGroupCount, 
                       descriptors.AcidBase.BasicGroupCount, 
                       descriptors.RingCount], ignore_3D=True)

    #print(f"Number of descriptors in calculator: {len(calc.descriptors)}")

    results = calc(smiles)
    results_dict = results.drop_missing().asdict()

    #print(f"Number of calculated descriptors: {len(results_dict)}")
    #print(results_dict)

    ######## MORGAN FINGERPRINT ###########

    morgan_fp = AllChem.GetMorganFingerprintAsBitVect(smiles, radius=2, nBits=1024)

    #print("Morgan fingerprint:")
    #print(list(morgan_fp))

    return results_dict, morgan_fp

def generate_dict (name : str, results_dict : str, morgan_fp : list, pubchem_id : str, drugbank_id : str, chembl_id : str) -> dict:
    dic = {}
    dic['name'] = name
    dic['mordred'] = results_dict
    dic['morgan'] = list(morgan_fp)
    dic['pubchem'] = pubchem_id
    dic['drugbank'] = drugbank_id
    dic['chembl'] = chembl_id
    # Añadir aqui lo que se quiera
    return dic


def generate_ligand_mapping (
    input_ligands : Optional[List[dict]],
    structure : 'Structure',
    output_mordred_filepath : str,
    ) -> dict:
    ligands = {}
    imported_ligands = None
    if os.path.exists(OUTPUT_LIGANDS_FILENAME):
        imported_ligands, _ = import_ligands()
        # Append the imported references to the overall ligands pool
        for k,v in imported_ligands.items():
            ligands[k] = v

    # If no input ligands are passed then stop here
    if not input_ligands:
        return []
    matching_ligands, non_matching_ligands = check_ligands_done(ligands, input_ligands)

    # print(ligands_to_include)

    # if ligands_to_exclude:
    #     for ligand_exclude in ligands_to_exclude:
    #         for ligand_dict in pubchem_dict:
    #             if ligand_dict == pubchem_dict[ligand_dict]:
    #                 del pubchem_dict[ligand_dict]
    
    if len(non_matching_ligands) >= 1:
        drugbank_dict, pubchem_dict, chembl_dict, pubchem_id_list, pubchem_name_list = obtain_ligands_id(input_ligands, non_matching_ligands)
        ligand_list  = []
        # Obtain all the information here and append it to a list
        if len(pubchem_id_list) == 1:
            pubchem_dict = pubchem_dict[0]
            #for ligand in pubchem_dict:
                #print('LIGAND:', ligand)
                # Get information from PubChem database
            smiles, _ = get_pubchem_data(pubchem_dict['pubchem'])

            # Get morgan and mordred descriptors, the SMILES is needed for this part
            results_dict, morgan_fp = obtain_mordred_morgan_descriptors(smiles)

            # Create a diccionary with the info obtained and refered to a ligand
            ligand_dict = generate_dict(pubchem_dict['name'], 
                                        results_dict, 
                                        morgan_fp, 
                                        pubchem_dict['pubchem'], 
                                        pubchem_dict['drugbank'], 
                                        pubchem_dict['chembl'])

            # Append the diccionary to a list
            ligand_list.append(ligand_dict)
        
        if len(pubchem_id_list) > 1:
            for ligand_dic in pubchem_dict:
                for ligand_id in pubchem_id_list:
                    if ligand_id == ligand_dic['pubchem']:
                        # Get information from PubChem database
                        smiles, _ = get_pubchem_data(ligand_dic['pubchem'])

                        # Get morgan and mordred descriptors, the SMILES is needed for this part
                        results_dict, morgan_fp = obtain_mordred_morgan_descriptors(smiles)
                        print(ligand_dic['name'])
                        # Create a diccionary with the info obtained and refered to a ligand
                        ligand_dict = generate_dict(ligand_dic['name'], 
                                                    results_dict, 
                                                    morgan_fp, 
                                                    ligand_dic['pubchem'], 
                                                    ligand_dic['drugbank'], 
                                                    ligand_dic['chembl'])

                        # Append the diccionary to a list
                        ligand_list.append(ligand_dict)

        # Write the ligands reference file
        # If there was a previous ligands.json, we have to check if these ligands match with the ligands defined in the inputs file and respect it
        if os.path.exists(OUTPUT_LIGANDS_FILENAME):
            _, file_content = import_ligands()
            file_content = file_content + ligand_list
            # if ligands_to_exclude:
            #     for ligand in file_content:
            #         for id in ligands_to_exclude:
            #             if ligand['pubchem'] == id:
            #                 file_content.remove(ligand)
            save_json(file_content, output_mordred_filepath)
        else:
            save_json(ligand_list, output_mordred_filepath)

    
        # Match the residues of the protein (PDB) to the ligands to identify it
        # This function returns a list of ligands that matched some of the PDB residues (index list)
        ligands = ligands_residues(structure, pubchem_id_list)
        return ligands, pubchem_name_list
    
    drugbank_dict, pubchem_dict, chembl_dict, pubchem_id_list, pubchem_name_list = obtain_ligands_id(input_ligands, matching_ligands)
    ligands = ligands_residues(structure, pubchem_id_list)
    # Also, it's needed to return a list of custom names for the ligands, if there were defined
    return ligands, pubchem_name_list
    

# If a 'ligands.json' exists, we need to check if it has the ligands according to the inputs
def check_ligands_done( ligands: dict, input_ligands: list ) -> list:

    #ids_json = list(ligands.values())

    # Variables de entrada
    #input_ligands = [{'name': 'Heme group', 'ngl': ':C', 'pubchem': '162640511'}, 
    #                {'name': 'Abiraterone', 'ngl': ':C', 'drugbank': 'DB05812'}]

    #dic_ligands = {'154229-19-3': {'name': '154229-19-3', 'pubchem': '132971', 'drugbank': None, 'chembl': None}}

    # Inicializar listas para los que coinciden y los que no coinciden
    matching_ligands = []
    non_matching_ligands = []

    # Función para comprobar coincidencias
    def check_match(input_ligand, dic_ligands):
        for dic_ligand in dic_ligands.values():
            if (input_ligand.get('pubchem') and input_ligand['pubchem'] == dic_ligand.get('pubchem')) or \
            (input_ligand.get('drugbank') and input_ligand['drugbank'] == dic_ligand.get('drugbank')) or \
            (input_ligand.get('chembl') and input_ligand['chembl'] == dic_ligand.get('chembl')):
                return True
        return False

    # Comprobar cada ligando en la lista de entrada
    for ligand in input_ligands:
        if check_match(ligand, ligands):
            matching_ligands.append(ligand)
        else:
            non_matching_ligands.append(ligand)

    # Imprimir resultados
    if matching_ligands:
        print("Ligands already in 'ligands.json':", matching_ligands)
    if non_matching_ligands:
        print("Ligands that have to be added:", non_matching_ligands)

    
    # pubchem_ids_inputs = []
    # for ligand_id in input_ligands:
    #     if isinstance(ligand_id, int):
    #         # Assuming the number is pubchem 
    #         pubchem_ids_inputs.append(ligand_id)
    #     if isinstance(ligand_id, dict) and 'pubchem' in ligand_id:
    #         pubchem_ids_inputs.append(ligand_id['pubchem'])
    #     if isinstance(ligand_id, dict) and 'drugbank' in ligand_id:
    #         pubchem_ids_inputs.append(find_drugbank_pubchem(ligand_id['drugbank']))
    #         #ligand_id['pubchem'] = find_drugbank_pubchem(ligand_id['drugbank'])
    #     if isinstance(ligand_id, dict) and 'chembl' in ligand_id:
    #         pubchem_ids_inputs.append(find_chembl_pubchem(ligand_id['chembl']))
    #         #ligand_id['pubchem'] = find_chembl_pubchem(ligand_id['chembl'])
    
    # ligands_to_include = [id for id in pubchem_ids_inputs if id not in ids_json] 
    # ligands_to_exclude = [id for id in ids_json if id not in pubchem_ids_inputs]
    
    return matching_ligands, non_matching_ligands


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

def count_atom_elements_per_ligand (pubchem_id_list : List[str]) -> dict:
    atom_elements_per_ligand_dict = {}
    for id_pubchem in pubchem_id_list:
        _, molecular_formula = get_pubchem_data(id_pubchem)
        molecular_formula = molecular_formula.replace("+", "").replace("-", "")
        l = parse_compound(molecular_formula)
        c = parseSplits(l)
        atom_elements_per_ligand_dict[id_pubchem] = c
    return atom_elements_per_ligand_dict

# This function associates elements in a list
# If a string is followed by a number then they go together
# If a string has no number then the number is 1 for this specific string
def parseSplits (splits : List[str]) -> dict:
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
        # This is to avoid having mismatches due to different protonation states
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

def ligands_residues (structure : 'Structure', pubchem_id_list : List[str]) -> List[dict]:
    # Load the structure where the atoms will be minresidue.indexed and obtain a residues dict
    atom_elements_per_residue = count_atom_elements_per_residue(structure)
    # From pubchem mine the atoms of the ligands and save it in a dict
    atom_elements_per_ligand = count_atom_elements_per_ligand(pubchem_id_list)
    #print('Ligands: ', atom_elements_per_ligand)
    matches = {}

    for pubchem, ligand_atom_element_count in atom_elements_per_ligand.items():
        for residue, residue_atom_element_count in atom_elements_per_residue.items():
            if match_ligandsID_to_res(ligand_atom_element_count, residue_atom_element_count):
                previous_match = matches.get(pubchem, None)
                if previous_match:
                    previous_match.append(residue.index)
                else:
                    matches[pubchem] = [ residue.index ]
    
    # Format the output as we expect it
    ligands = [ { 'name': key, 'residue_indices': value, 'match': { 'ref': { 'pubchem': key } } } for key, value in matches.items() ]
    print('Ligands: ', ligands)
    print('A total of',len(atom_elements_per_residue), 'residues and', len(atom_elements_per_ligand), 'ligands were found.')
    print('Matched:', len(ligands),'/', len(atom_elements_per_ligand),'of the ligands.')

    return ligands