import os
import json
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem import inchi
from mordred import Calculator, descriptors
from model_workflow.utils.constants import LIGANDS_MATCH_FLAG, PDB_TO_PUBCHEM, NOT_MATCHED_LIGANDS
from model_workflow.utils.auxiliar import InputError, load_json, save_json, request_pdb_data
from model_workflow.utils.type_hints import *
from model_workflow.utils.structures import Structure
from urllib.request import Request, urlopen
from urllib.parse import urlencode
from urllib.error import HTTPError, URLError
import re
import requests

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
        raise RuntimeError('Wrong Pubchem data structure: no record: ' + request_url)
    sections = record.get('Section', None)
    if sections == None:
        raise RuntimeError('Wrong Pubchem data structure: no sections: ' + request_url)
    names_and_ids_section = next((section for section in sections if section.get('TOCHeading', None) == 'Names and Identifiers'), None)
    if names_and_ids_section == None:
        raise RuntimeError('Wrong Pubchem data structure: no name and ids section: ' + request_url)
    names_and_ids_subsections = names_and_ids_section.get('Section', None)
    if names_and_ids_subsections == None:
        raise RuntimeError('Wrong Pubchem data structure: no name and ids subsections: ' + request_url)

    # Mine the name
    synonims = next((s for s in names_and_ids_subsections if s.get('TOCHeading', None) == 'Synonyms'), None)
    if synonims == None:
        descriptors = next((s for s in names_and_ids_subsections if s.get('TOCHeading', None) == 'Computed Descriptors'), None)
        descriptors_subsections = descriptors.get('Section', None)
        if descriptors_subsections == None:
            raise RuntimeError('Wrong Pubchem data structure: no name and ids subsections: ' + request_url)
        depositor_supplied_descriptors = next((s for s in descriptors_subsections if s.get('TOCHeading', None) == 'IUPAC Name'), None)
        name_substance = depositor_supplied_descriptors.get('Information', None)[0].get('Value', {}).get('StringWithMarkup', None)[0].get('String', None)
    else:
        synonims_subsections = synonims.get('Section', None)
        if synonims_subsections == None:
            raise RuntimeError('Wrong Pubchem data structure: no name and ids subsections: ' + request_url)
        depositor_supplied_synonims = next((s for s in synonims_subsections if s.get('TOCHeading', None) == 'Depositor-Supplied Synonyms'), None)
        if depositor_supplied_synonims == None:
            removed_synonims = next((s for s in synonims_subsections if s.get('TOCHeading', None) == 'Removed Synonyms'), None)
            name_substance = removed_synonims.get('Information', None)[0].get('Value', {}).get('StringWithMarkup', None)[0].get('String', None)
        else:
            name_substance = depositor_supplied_synonims.get('Information', None)[0].get('Value', {}).get('StringWithMarkup', None)[0].get('String', None)
        
    # Mine the SMILES
    computed_descriptors_subsection = next((s for s in names_and_ids_subsections if s.get('TOCHeading', None) == 'Computed Descriptors'), None)
    if computed_descriptors_subsection == None:
        raise RuntimeError('Wrong Pubchem data structure: no computeed descriptors: ' + request_url)
    canonical_smiles_section = computed_descriptors_subsection.get('Section', None)
    if canonical_smiles_section == None:
        raise RuntimeError('Wrong Pubchem data structure: no canonical SMILES section: ' + request_url)
    canonical_smiles = next((s for s in canonical_smiles_section if s.get('TOCHeading', None) == 'Canonical SMILES'), None)
    if canonical_smiles == None:
        # In some cases there is no canonical SMILES but a non-canonical one could exists
        non_canonical_smiles_section = next((s for s in canonical_smiles_section if s.get('TOCHeading', None) == 'SMILES'), None)
        if non_canonical_smiles_section == None:
            raise RuntimeError('Wrong Pubchem data structure: no canonical SMILES: ' + request_url)
    
    if canonical_smiles:
        smiles = canonical_smiles.get('Information', None)[0].get('Value', {}).get('StringWithMarkup', None)[0].get('String', None)
    if non_canonical_smiles_section:
        smiles = non_canonical_smiles_section.get('Information', None)[0].get('Value', {}).get('StringWithMarkup', None)[0].get('String', None)
    
    if smiles == None:
        raise RuntimeError('Wrong Pubchem data structure: no SMILES: ' + request_url)

    # Mine target data: MOLECULAR FORMULA
    molecular_formula_subsection = next((s for s in names_and_ids_subsections if s.get('TOCHeading', None) == 'Molecular Formula'), None)
    if molecular_formula_subsection == None:
        raise RuntimeError('Wrong Pubchem data structure: no molecular formula section: ' + request_url)
    molecular_formula = molecular_formula_subsection.get('Information', None)[0].get('Value', {}).get('StringWithMarkup', None)[0].get('String', None)
    if molecular_formula == None:
        raise RuntimeError('Wrong Pubchem data structure: no molecular formula: ' + request_url)
    
    # Mine target data: PDB ID
    pdb_id = None
    pdb_id_subsection = next((s for s in sections if s.get('TOCHeading', None) == 'Interactions and Pathways'), None)
    # If this section is missing then it means this PubChem compound has no PDB id
    if pdb_id_subsection:
        pdb_id_subsections = pdb_id_subsection.get('Section', None)
        if pdb_id_subsections == None:
            raise RuntimeError('Wrong Pubchem data structure: no name and ids subsections: ' + request_url)
        bond_structures = next((s for s in pdb_id_subsections if s.get('TOCHeading', None) == 'Protein Bound 3D Structures'), None)
        if bond_structures:
            bond_structures_section = bond_structures.get('Section', None)
            # If this section is missing then it means this PubChem compound has no PDB id
            if bond_structures_section:
                ligands_structure = next((s for s in bond_structures_section if s.get('TOCHeading', None) == 'Ligands from Protein Bound 3D Structures'), None)
                if ligands_structure == None:
                    raise RuntimeError('Wrong Pubchem data structure: no Ligands from Protein Bound 3D Structures section: ' + request_url)
                ligands_structure_subsections = ligands_structure.get('Section', None)
                if ligands_structure_subsections == None:
                    raise RuntimeError('Wrong Pubchem data structure: no Ligands from Protein Bound 3D Structures subsections: ' + request_url)
                ligands_pdb = next((s for s in ligands_structure_subsections if s.get('TOCHeading', None) == 'PDBe Ligand Code'), None)
                if ligands_pdb == None:
                    raise RuntimeError('Wrong Pubchem data structure: no PDBe Ligand Code: ' + request_url)
                pdb_id = ligands_pdb.get('Information', None)[0].get('Value', {}).get('StringWithMarkup', None)[0].get('String', None)
    
    # Mine de INCHI and INCHIKEY
    inchi_section = next((s for s in canonical_smiles_section if s.get('TOCHeading', None) == 'InChI'), None)
    if inchi_section == None:
        raise RuntimeError('Wrong Pubchem data structure: no InChI: ' + request_url)
    if inchi_section:
        inchi = inchi_section.get('Information', None)[0].get('Value', {}).get('StringWithMarkup', None)[0].get('String', None)
    
    inchikey_section = next((s for s in canonical_smiles_section if s.get('TOCHeading', None) == 'InChIKey'), None)
    if inchikey_section == None:
        raise RuntimeError('Wrong Pubchem data structure: no InChIKey: ' + request_url)
    if inchikey_section:
        inchikey = inchikey_section.get('Information', None)[0].get('Value', {}).get('StringWithMarkup', None)[0].get('String', None)
    # Prepare the pubchem data to be returned
    return { 'name': name_substance, 'smiles': smiles, 'formula': molecular_formula , 'pdbid': pdb_id, 'inchi': inchi, 'inchikey': inchikey }


def find_drugbank_pubchem (drugbank_id):
    # Request Drugbank
    request_url = Request(
    url=f'https://go.drugbank.com/drugs/{drugbank_id}',
    headers={
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
        'Accept-Language': 'en-US,en;q=0.9',
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,*/*;q=0.8'
    }
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
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f'WARNING: Cannot generate a molecule from SMILES {smiles}')
        return None

    AllChem.Compute2DCoords(mol)
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
    mordred_results = calc(mol).drop_missing().asdict()

    ######## MORGAN FINGERPRINT ###########
    morgan_fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    ao = rdFingerprintGenerator.AdditionalOutput()
    ao.AllocateAtomCounts()
    ao.AllocateAtomToBits()
    ao.AllocateBitInfoMap()
    fp = morgan_fp_gen.GetFingerprint(mol,additionalOutput=ao)
    morgan_fp_bit_array = list(fp)
    morgan_highlight_atoms = {}
    for bit, atoms in ao.GetBitInfoMap().items():
        morgan_highlight_atoms[bit] = list(set(atom for atom, radius in atoms))

    return mordred_results, morgan_fp_bit_array, morgan_highlight_atoms, mol_block

# Check the cache to see if the ligand has already been matched or not, if not the ligand is skipped
def check_matched_ligand (ligand: dict, ligand_data: dict, cache: 'Cache') -> bool:
    # The pubchem id could be in the ligands data or the ligand input
    pubchem_id_1 = ligand.get('pubchem', None)
    if ligand_data != None:
        pubchem_id_2 = ligand_data.get('pubchem', None)
    else:
        pubchem_id_2 = None
    # If we have no pubchem id then we cannot check if the ligand matched
    if not pubchem_id_1 and not pubchem_id_2:
        return False
    # If we have a pubchem id then check if it is already in the cache
    cache_not_matched = cache.retrieve(NOT_MATCHED_LIGANDS, [])
    for pubchem in cache_not_matched:
        if pubchem == pubchem_id_1 or pubchem == pubchem_id_2:
            print(f" Ligand with PubChem id {pubchem} didn't match before, skipping it")
            return True
    
    return False
    
def generate_ligand_mapping (
    structure : 'Structure',
    cache : 'Cache',
    input_ligands : Optional[List[dict]],
    pdb_ids : List[str],
    output_filepath : str,
    mercy : List[str] = [],
    ) -> List[dict]:

    """Generate a map of residues associated to ligands."""
    # Merge input ligands and pdb ligands
    ligands = []
    if input_ligands:
        ligands += input_ligands
    # Check we have cached pdb 2 pubchem values
    pdb_to_pubchem_cache = cache.retrieve(PDB_TO_PUBCHEM, {})
    new_data_to_cache = False
    # Get input ligands from the pdb ids, if any
    for pdb_id in pdb_ids:
        # Check we have cached this specific pdb
        pubchem_ids_from_pdb = pdb_to_pubchem_cache.get(pdb_id, None)
        if pubchem_ids_from_pdb != None:
            print(f' Retriving from cache PubChem ids for PDB id {pdb_id}: ')
            if len(pubchem_ids_from_pdb) > 0:
                print('  PubChem ids: ' + ', '.join(pubchem_ids_from_pdb))
            else:
                print('  This PDB id has no PubChem ids')

        # If we had no cached pdb 2 pubchem then ask for them
        if pubchem_ids_from_pdb == None:
            pubchem_ids_from_pdb = pdb_to_pubchem(pdb_id)
            # Save the result in the cache object so it is saved to cache later
            pdb_to_pubchem_cache[pdb_id] = pubchem_ids_from_pdb
            new_data_to_cache = True
        for pubchem_id in pubchem_ids_from_pdb:
            # Ligands in the structure (PDB) and the 'inputs.json' could be the same so it's not necessary to do it twice
            if not any('pubchem' in ligand and ligand['pubchem'] == pubchem_id for ligand in ligands):
                ligands.append({ 'pubchem': pubchem_id, 'pdb': True })
    # Save all pdb to pubchem results in cache, in case there is anything new
    if new_data_to_cache: cache.update(PDB_TO_PUBCHEM, pdb_to_pubchem_cache)

    # If no input ligands are passed then stop here
    if len(ligands) == 0:
        return []
    
    # Save data from all ligands to be saved in a file
    json_ligands_data = []
    ligands_data = []

    # Load the ligands file if exists already
    if os.path.exists(output_filepath):
        json_ligands_data += import_ligands(output_filepath)
        # If json ligands exists and it is empty means that ligands analysis has been already done but no ligands were matched
        # so the file will contain an empty list []
        if len(json_ligands_data) == 0:
            print('No ligands have been matched yet.\nIf you want to force a ligand to be matched, please provide the field "residues" in the inputs.json file.')
            return []

    # Visited formulas
    visited_formulas = []
    # Save the maps of every ligand
    ligand_maps = []
    # Get cached ligand data
    # AGUS: esto lo creamos por alguna razón para que funcione sin internet (en el cluster) pero realmente
    # AGUS: llena el cache con datos que no son necesarios porque toda esa info está en el ligand_references.json
    #ligand_data_cache = cache.retrieve(LIGANDS_DATA, {})
    # Iterate input ligands
    for ligand in ligands:
        # Set the pubchem id which may be assigned in different steps
        pubchem_id = None
        # If input ligand is not a dict but a single int/string then handle it
        if type(ligand) == int:
            print(f'A ligand number ID has been identified {ligand}, assuming that is a PubChem ID...')
            ligand = { 'pubchem': str(ligand) }
        elif type(ligand) == str:
            raise InputError(f'A name of ligand has been identified: {ligand}. Anyway, provide at least one of the following IDs: DrugBank, PubChem, ChEMBL.')
        # Check if we already have this ligand data
        ligand_data = obtain_ligand_data_from_file(json_ligands_data, ligand)
        # Check if the ligand didn't match before
        if check_matched_ligand(ligand, ligand_data, cache): continue
        # If we do not have its data try to get from the cache
        # if not ligand_data:
        #     # If this is a ligand not in ligans.json but in cache it means it comes from PDB, never form user inputs
        #     # For this reason, this ligand will always have a PubChem id
        #     pubchem_id = ligand.get('pubchem', None)
        #     if pubchem_id:
        #         ligand_data = ligand_data_cache.get(pubchem_id, None)
        #     # If we still have no ligand data then request it to PubChem
        #     if not ligand_data:
        #         ligand_data = obtain_ligand_data_from_pubchem(ligand)
        #         # Save data mined from PubChem in the cache
        #         pubchem_id = ligand_data['pubchem']
        #         ligand_data_cache[pubchem_id] = { **ligand_data }
        #         cache.update(LIGANDS_DATA, ligand_data_cache)
        # Add current ligand data to the general list
        if not ligand_data:
            ligand_data = obtain_ligand_data_from_pubchem(ligand)
        ligands_data.append(ligand_data)
        # Get pubchem id
        pubchem_id = ligand_data.get('pubchem', None)
        # If we already visited a different ligand but with identical formula then we skip this ligand 
        # Note that the mapping will be identical thus overwritting the previous map
        # However, ligands forced by the user are processed before so we keep them as priority
        formula = ligand_data.get('formula', None)
        if formula in visited_formulas:
            print(f'WARNING: Ligand with PubChem Id {pubchem_id} has a formula which has been already matched')
            # AGUS: en este punto si el usuario ha definido el mismo ligando con diferente selección quiere decir que está repetido
            # AGUS: y que deberíamos mantener el ligando para diferentes residuos matcheados
            if forced_selection:
                visited_formulas.pop()
                ligands_data.pop()
                print(f'WARNING: Ligand with PubChem Id {pubchem_id} has been forced by the user and its repeated')
            else:
                ligands_data.pop()
                continue
        # Add it to the list of visited formulas
        visited_formulas.append(formula)
        # If the user defined a ligand name, it will be respectedand added to the metadata
        # if there isn't a defined name, it will be mined from PubChem
        user_forced_ligand_name = ligand.get('name', None)
        # Map structure residues with the current ligand
        # If the user forced the residues then use them
        forced_selection = ligand.get('selection', None)
        if forced_selection:
            # Could be a single residue or a list of residues
            selection_atoms = structure.select(forced_selection)
            residue_indices = structure.get_selection_residue_indices(selection_atoms)
            ligand_map = { 'name': pubchem_id, 'residue_indices': residue_indices, 'match': { 'ref': { 'pubchem': pubchem_id } } }
            if user_forced_ligand_name: ligand_map['forced_name'] = user_forced_ligand_name
        # If the user did not force the residues then map them
        else:
            ligand_map = map_ligand_residues(structure, ligand_data)
        # Add current ligand map to the general list
        ligand_maps.append(ligand_map)
        # If the ligand did not map then discard it
        did_not_match = len(ligand_map['residue_indices']) == 0
        if did_not_match:
            is_forced_by_user = not bool(ligand.get('pdb', False))
            if is_forced_by_user:
                # At this point we should have macthed all sequences
                # If not, kill the process unless mercy was given
                must_be_killed = LIGANDS_MATCH_FLAG not in mercy
                if must_be_killed:
                    raise InputError(f'Ligand with PubChem Id {pubchem_id} did not map with any residue')
            # Otherwise simply remove it from the final ligand references file and the forwarded maps 
            ligands_data.pop()
            ligand_maps.pop()
            # Update the cache with the pubchem id of the ligands that didn't match
            not_matched_pubchems = cache.retrieve(NOT_MATCHED_LIGANDS, [])
            not_matched_pubchems.append(ligand_map['name'])
            cache.update(NOT_MATCHED_LIGANDS, not_matched_pubchems)
        # If the ligand matched then calculate some additional data
        # Get morgan and mordred descriptors, the SMILES is needed for this part
        smiles = ligand_data['smiles']
        mordred_results, morgan_fp_bit_array, morgan_highlight_atoms, mol_block = obtain_mordred_morgan_descriptors(smiles)
        ligand_data['mordred'] = mordred_results
        ligand_data['morgan_fp_bit_array'] = morgan_fp_bit_array
        ligand_data['morgan_highlight_atoms'] = morgan_highlight_atoms
        ligand_data['mol_block'] = mol_block

    # At this point maybe some ligands in the structure may not match because they have not been included by either the author or the user. Alternatively, the authors may have added them
    # So those "ligands" are residues that are not matched to any ligand in the structure and the analyses will not run
    # LORE AGUS: en el dataset de MDBind (cineca) hay muchas simulaciones con ligandos modificados (incluso solo un átomo) y no matchean 
    # LORE AGUS: con los que están en el PDB, por lo que no se hace ningún análisis sobre estos (y se debería)
    # LORE AGUS: ahora se intenta buscar su pubchem ID para extraer la información a partir del inchikey (da menos problemas que el smiles) y se hacen los análisis pertinentes. 
    # AGUS: si no se puede extraer información de estos porque el inchikey no funciona, se muestra la información más básica posible y se genera una referencia vacía (será problemático)

    # Obtain a list of the possible residues in the structure that are not either protein/nucleic/water/lipid or matched ligands
    residue_element_count = count_atom_elements_per_residue(structure)
    # Iterate over the residues and check if they are not matched
    for residue, atoms_count in residue_element_count.items():
        # If the residue is not connected to other residues then it may be a ligand
        is_single_residue = len(residue.get_bonded_residue_indices()) == 0
        if not is_single_residue: continue
        # If the residue is not solvent, lipid or ion then it may be a ligand
        has_right_classification = residue.classification not in ['solvent', 'lipid', 'ion']
        if not has_right_classification: continue
        # Check if the residue is already matched to a ligand
        matched = False
        for ligand_map in ligand_maps:
            if residue.index in ligand_map['residue_indices']:
                matched = True
                break
        # If the residue is not matched then add it to the ligands data
        if matched: continue
        # Create a new ligand data with the residue name and its atoms count
        ligand_data = {
            'name': residue.name,
            'pubchem': None,
            'drugbank': None,
            'chembl': None,
            'smiles': None,
            'formula': None,
            'morgan': None,
            'mordred': None,
            'pdbid': None
        }
        # Add the residue as a ligand map
        ligand_map = { 
            'name': residue.name, 
            'residue_indices': [residue.index], 
            'match': { 'ref': { 'pubchem': None } } 
        }
        # Get a PDB only with the target residue
        residue_selection = structure.select_residue_indices([residue.index])
        filtered_structure = structure.filter(residue_selection)
        ligand_filename = "ligand.pdb"
        filtered_structure.generate_pdb_file(ligand_filename)

        # Obtain the InChIKey from the PDB file
        inchikey = obtain_inchikey_from_pdb(ligand_filename)
        # Remove the temporary PDB file
        os.remove(ligand_filename)
        # If the InChIKey is not None then try to obtain the PubChem CID
        if inchikey:
            cid = search_cid_by_inchikey(inchikey)
            if cid:
                ligand_data = obtain_ligand_data_from_pubchem( {'pubchem': cid} )
                smiles = ligand_data['smiles']
                mordred_results, morgan_fp_bit_array, morgan_highlight_atoms, mol_block = obtain_mordred_morgan_descriptors(smiles)
                ligand_data['mordred'] = mordred_results
                ligand_data['morgan_fp_bit_array'] = morgan_fp_bit_array
                ligand_data['morgan_highlight_atoms'] = morgan_highlight_atoms
                ligand_data['mol_block'] = mol_block
        # Add the ligand data to the list of ligands data
        ligands_data.append(ligand_data)
        # Add the ligand map to the list of ligand maps
        ligand_maps.append(ligand_map)
    
    # Export ligands to a file
    save_json(ligands_data, output_filepath)

    # Log matched ligands
    if not ligand_maps:
        print('No ligands were matched')
    else:
        print('Matched ligands:')
        for ligand_map in ligand_maps:
            residue_indices = ligand_map["residue_indices"]
            residue_count = len(residue_indices)
            plural_suffix = '' if residue_count == 1 else 's'
            print(f' - {ligand_map["name"]}: {residue_count} residue{plural_suffix}')

    return ligand_maps
    
# Set the expected ligand data fields
LIGAND_DATA_FIELDS = set(['name', 'pubchem', 'drugbank', 'chembl', 'smiles', 'formula', 'morgan', 'mordred', 'pdbid'])

# Import ligands json file so we do not have to rerun this process
def import_ligands (output_filepath : str) -> dict:
    # Read the file
    imported_ligands = load_json(output_filepath)
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
        'formula': None,
        'pdbid': None,
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
    # Clean the formula from charges
    # These charges include numbers which are not atom counts (e.g. 49867169 -> C18H16NO5PS-2)
    parsed_molecular_formula = re.sub('[+-][0-9]*', '', molecular_formula)
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

# Given a PDB ligand code, get its pubchem
def pdb_ligand_to_pubchem (pdb_ligand_id : str) -> Optional[str]:
    # Set the request query
    query = '''query molecule($id:String!){
        chem_comp(comp_id:$id) { rcsb_chem_comp_related{ resource_name resource_accession_code } }
    }'''
    # Request PDB data
    parsed_response = request_pdb_data(pdb_ligand_id, query)
    related_resources = parsed_response['rcsb_chem_comp_related']
    # It may happend that a ligand code has no related resources at all
    # e.g. ZN
    if not related_resources: return None
    pubchem_resource = next((resource for resource in related_resources if resource['resource_name'] == 'PubChem'), None)
    if not pubchem_resource: return None
    return pubchem_resource['resource_accession_code']

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

# Given a PDB ligand code, get its pubchem
# Ask to pubchem if the ligand exists and hope there is only one result
# DANI: No se ha provado a fondo
def pdb_ligand_to_pubchem_RAW_RAW (pdb_ligand_id : str) -> Optional[str]:
    # Set the request URL
    request_url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{pdb_ligand_id}/json'
    # Run the query
    parsed_response = None
    try:
        with urlopen(request_url) as response:
            parsed_response = json.loads(response.read().decode("utf-8"))
    # If the accession is not found in the PDB then we can stop here
    except HTTPError as error:
        # This may happen for weird things such as UNX (unknown atom or ion)
        if error.code == 404:
            return None
        print(error.msg)
        raise RuntimeError('Something went wrong with the PDB ligand request in PubChem: ' + request_url)
    # Mine the pubchem id
    compounds = parsed_response['PC_Compounds']
    if len(compounds) != 1:
        # There could be more than one result
        # For example, the case HEM: https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/HEM/json 
        # In this case we picked the first result
        compound = compounds[0]
        #raise RuntimeError('We are not having 1 and only 1 result from PubChem: ' + request_url)
    compound = compounds[0]
    id1 = compound['id']
    id2 = id1['id']
    pubchem_id = str(id2['cid'])
    return pubchem_id

# Given a PDB id, get all its ligand codes
# e.g. 2I3I -> 618, BTB, ZN, EDO, LI
def get_pdb_ligand_codes (pdb_id : str) -> List[str]:
    # Set the request query
    query = '''query structure($id: String!) {
        entry(entry_id: $id) {
            nonpolymer_entities { nonpolymer_comp { chem_comp { id } } }
        }
    }'''
    # Request PDB data
    parsed_response = request_pdb_data(pdb_id, query)
    # Mine data for nonpolymer entities
    nonpolymers = parsed_response['nonpolymer_entities']
    # If there are no nonpolymer entities, another type of entitie could be used
    # AGUS: esto es un caso muy concreto que me encontré con el PDB 1N3W
    # AGUS: el ligando en este caso es una 'Biologically Interesting Molecules' y se muestran como 'PRD_'
    prd_code = None
    if nonpolymers == None:
        # Get the prd ligand code
        prd_code = get_prd_ligand_code(pdb_id)
        if prd_code != None:
            return [prd_code]

    if nonpolymers == None and prd_code == None: return []
    # Iterate nonpolymer entities to mine each PDB code
    ligand_codes = []
    for nonpolymer in nonpolymers:
        ligand_code = nonpolymer['nonpolymer_comp']['chem_comp']['id']
        ligand_codes.append(ligand_code)
    print(f' Ligand codes for PDB id {pdb_id}: ' + ', '.join(ligand_codes))
    return ligand_codes

# Given a PDB id, get its PRD ligand code
def get_prd_ligand_code (pdb_id : str) -> Optional[str]:
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
        print(f'WARNING: Cannot find PRD ligand code for PDB id {pdb_id}')
        return None
    prd_id = pdbx_id[0].get('prd_id', None) # AGUS: podría haber casos donde haya más de uno? 
    if prd_id is None:
        print(f'WARNING: Cannot find PRD ligand code for PDB id {pdb_id}')
        return None
    # If the PRD id is not empty then return it
    return prd_id

# Given a pdb id, get its pubchem ids
# DANI: De momento no usamos las SMILES que alguna vez me han dado problemas (e.g. 2I3I)
# e.g. 4YDF -> 
def pdb_to_pubchem (pdb_id : str) -> List[str]:
    print(f'Searching PubChem ids for PDB {pdb_id}')
    pubchem_ids = []
    # Iterate over pdb ligand codes
    ligand_codes = get_pdb_ligand_codes(pdb_id)
    for ligand_code in ligand_codes:
        # Ask the PDB API for the ligand
        pubchem_id = pdb_ligand_to_pubchem(ligand_code)
        # If this did not work then try mining the PDB client with a web crawler
        if not pubchem_id:
            pubchem_id = pdb_ligand_to_pubchem_RAW(ligand_code)
        # If this did not work then try it from PubChem
        if not pubchem_id:
            pubchem_id = pdb_ligand_to_pubchem_RAW_RAW(ligand_code)
        # Otherwise we surrender
        if not pubchem_id:
            print(f' {ligand_code} -> No PubChem id')
            continue
        print(f' {ligand_code} -> {pubchem_id}')
        pubchem_ids.append(pubchem_id)
            
    return pubchem_ids

def obtain_inchikey_from_pdb(pdb_file : str) -> Optional[str]:
    mol = Chem.MolFromPDBFile(pdb_file, sanitize=True)
    if mol is None:
        print("❌ InchiKey couldn't be constructed from PDB file.")
        return None
    return inchi.MolToInchiKey(mol)

def search_cid_by_inchikey(inchikey : str) -> Optional[str]:
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/cids/JSON"
    r = requests.get(url)
    if r.ok:
        data = r.json()
        return data['IdentifierList']['CID'][0]
    else:
        return None
