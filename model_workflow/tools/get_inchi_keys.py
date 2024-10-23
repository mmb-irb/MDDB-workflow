from typing import Optional
import MDAnalysis
from rdkit import Chem
from model_workflow.utils.file import File
from model_workflow.utils.structures import Structure
from model_workflow.tools.residues_library import aminoacids
import requests

# This function uses MDAnalysis to create a RDkit Molecule
def get_inchi_keys (
    input_structure_file : 'File',
    input_topology_file : Optional['File'],
    structure : 'Structure'
) -> dict:
    # Set the MDAnalysis universe
    uni = MDAnalysis.Universe(input_topology_file, input_structure_file)
    # Get the residues that are not aminoacids and save the first found
    unq_res = {}
    for res in structure.residues:
        res_nm = res.name
        if res_nm not in aminoacids and res_nm not in unq_res:
            unq_res[res_nm] = res
        # TO-DO check the other residues are duplicates from the first one
        #if res_nm not in aminoacids and res_nm in uid_res:
        #    assert uid_res[res_nm] == res

    # Loop though the unique residues and find the InChi key
    inchi_keys = {}
    for res_nm,residue in unq_res.items():
        res_sele = residue.get_selection().to_mdanalysis()
        # Select residues atoms with MDAnalysis
        res_mda = uni.select_atoms(res_sele)
        # Convert to RDKIT and get InChi data
        res_RD = res_mda.atoms.convert_to("RDKIT")
        inchi = Chem.MolToInchi(res_RD)
        inchikey = Chem.InchiToInchiKey(inchi)
        inchi_keys[res_nm] = inchikey
    return inchi_keys

def is_in_LIPID_MAPS(inchikey) -> dict:
# search the InChi keys in LIPID MAPS
    headers = {'accept': 'json'}
    url = f"https://www.lipidmaps.org/rest/compound/inchi_key/{inchikey}/all"
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        js = response.text
        if js != '[]':
           return js        
    else:
        print(f"Error for {inchikey}: {response.status_code}")
