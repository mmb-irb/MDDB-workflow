# WARNING: This is very deprecated. Update the inputs if you want to reuse it

import sys
import json

# The mutation name
mutation = sys.argv[1]
# The ligand name
ligand = sys.argv[2]

# Set the default inputs
inputs = {
    "chainnames": {
        "A": "Protein",
        "B": "Ligand"
    },
    "ligands": [],
    "exceptions": [],
    "interactions": [
        {
            "name": "protein-ligand interaction",
            "agent_1": "protein",
            "selection_1": "chain A",
            "agent_2": "ligand",
            "selection_2": "chain B"
        }
    ],
    "membranes": [],
    "customs": [],
    "unit": "Other",
    "pdbIds": [],
    "name": "Unnamed",
    "description": None,
    "authors": "Adam Hospital",
    "groups": "Orozco lab, IRB Barcelona",
    "contact": None,
    "program": None,
    "version": None,
    "method": "Classical MD",
    "license": "This trajectory dataset is released under a Creative Commons Attribution 4.0 International Public License",
    "linkcense": "https://creativecommons.org/licenses/by/4.0/",
    "citation": None,
    "thanks": None,
    "length": 200,
    "temp": None,
    "ensemble": None,
    "timestep": None,
    "ff": None,
    "wat": None,
    "boxtype": None
}

# Set the available ligands
available_ligands = {
    'yy3': {
        "name": "YY3",
        "ngl": ":B",
        "drugbank": "DB09330",
        "chembl": "CHEMBL3353410"
    },
    'ire': {
        "name": "IRE",
        "ngl": ":B",
        "drugbank": "DB00317",
        "chembl": "CHEMBL939"
    },
    'aq4': {
        "name": "AQ4",
        "ngl": ":B",
        "drugbank": "DB00530",
        "chembl": "CHEMBL553"
    },
    'fmm': {
        "name": "FMM",
        "ngl": ":B",
        "drugbank": "DB01259",
        "chembl": "CHEMBL554"
    },
    'ico': {
        "name": "ICO",
        "ngl": ":B",
        "drugbank": "DB11737",
        "chembl": "CHEMBL2087361"
    }
}

# Find the corresponding ligand
ligand_key = ligand.split('_')[0]
input_ligand = available_ligands[ligand_key]

# Build the input name
input_name = 'EGFR - ' + mutation.replace('_',' ') + ' ' + ligand.replace('_',' ')

# Update default inputs
inputs['name'] = input_name
inputs['ligands'].append(input_ligand)

# Rewrite the inputs file in a pretty formatted way
inputs_filename = 'inputs.json'
with open(inputs_filename, 'w') as file:
    json.dump(inputs, file, indent=4)