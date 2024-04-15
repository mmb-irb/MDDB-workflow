import os
import sys
import json
import urllib.request
from rdkit import Chem
from rdkit.Chem import AllChem
from mordred import Calculator, descriptors
from typing import List, Tuple, Optional, Union

from pathlib import Path

from model_workflow.tools.residues_library import residue_name_2_letter
from model_workflow.utils.auxiliar import load_json, save_json
from urllib.request import Request, urlopen
from model_workflow.utils.structures import Structure


print(Structure)