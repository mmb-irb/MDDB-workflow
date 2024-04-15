import os
import sys
import json
import urllib.request

from typing import List, Tuple, Optional, Union

from pathlib import Path

from model_workflow.tools.residues_library import residue_name_2_letter
from model_workflow.utils.auxiliar import load_json, save_json
from model_workflow.utils.constants import OUTPUT_MORDRED_FILENAME, REFERENCE_SEQUENCE_FLAG, DEFAULT_INPUTS_FILENAME 

import xmltodict

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Blast import NCBIWWW
from urllib.request import Request, urlopen

from model_workflow.utils.structures import Structure
from rdkit import Chem
from rdkit.Chem import AllChem
from mordred import Calculator, descriptors
import os
import pandas as pd
import numpy as np
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
import json
#import requests
import lxml
from bs4 import BeautifulSoup


def obtener_id_drugbank(archivo_json):
    with open(archivo_json, 'r') as f:
        data = json.load(f)
        for ligand in data.get('ligands', []):
            drugbank_id = ligand.get('drugbank')
            print(drugbank_id)
            if drugbank_id:
                return drugbank_id
        return None 

ruta_archivo = '/home/agarciad/workflow/model_workflow/analyses/inputs.json'
id_drugbank = obtener_id_drugbank(ruta_archivo)

if id_drugbank:
    print("ID de DrugBank obtenido:", id_drugbank)
else:
    print("No se encontró el ID de DrugBank en el archivo.")

# inputs_json = DEFAULT_INPUTS_FILENAME
# if not os.path.exists(inputs_json):
#     raise SystemExit('ERROR: The file does not exist')

req = Request(
    url= f'https://go.drugbank.com/structures/small_molecule_drugs/{id_drugbank}.smiles', # f'https://go.drugbank.com/structures/small_molecule_drugs/{id_drugbank}.smiles', 'https://pubchem.ncbi.nlm.nih.gov/compound/1986#section=Canonical-SMILES&fullscreen=true'
    headers={'User-Agent': 'Mozilla/5.0'}
)
webpage = urlopen(req).read()
print(webpage[1:])
smiles = Chem.MolFromSmiles(webpage)

calc = Calculator(descriptors, ignore_3D=True)

# calc = Calculator([descriptors.ABCIndex, 
#                    descriptors.AcidBase.AcidicGroupCount, 
#                    descriptors.AcidBase.BasicGroupCount, 
#                    descriptors.RingCount], ignore_3D=True)
# results = calc(mol)

print(f"Number of descriptors in calculator: {len(calc.descriptors)}")

#results = calc(smiles)
#results_dict = results.drop_missing().asdict()

#print(f"Number of calculated descriptors: {len(results_dict)}")
#print(results_dict)


########################################
######### MORGAN FINGERPRINT ###########
########################################

# morgan_fp = AllChem.GetMorganFingerprintAsBitVect(smiles, radius=2, nBits=1024)

# print("Morgan fingerprint:")
# print(list(morgan_fp))


########################################
# PRUEBA WEB-CRAWLER
# url = "https://pubchem.ncbi.nlm.nih.gov/compound/1986#section=Canonical-SMILES&fullscreen=true"
# headers = {
#   'User-Agent': 'Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/63.0.3239.132 Safari/537.36 QIHU 360SE'
# }
# f = requests.get(url, headers = headers)

#soup = BeautifulSoup(req,'lxml')
#print(soup)




