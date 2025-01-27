import pytraj as pt

from os.path import exists
from json import load
from pathlib import Path
from subprocess import run, PIPE
from model_workflow.utils.type_hints import *

import re

from model_workflow.utils.constants import GROMACS_EXECUTABLE, STANDARD_TOPOLOGY_FILENAME, RAW_CHARGES_FILENAME

from MDAnalysis.topology.TPRParser import TPRParser
from MDAnalysis.topology.TOPParser import TOPParser

def get_charges (charges_source_file : 'File') -> List[float]:
    """
    Extract charges from a source file.

    Returns:
    List[float]: A list of atomic charges if extraction is successful, 
                 otherwise None if the file does not exist.

    """
    if not charges_source_file or not charges_source_file.exists:
        return None
    print('Charges in the "' + charges_source_file.filename + '" file will be used')
    charges = None
    # If we have the standard topology then get charges from it
    if charges_source_file.filename == STANDARD_TOPOLOGY_FILENAME:
        with open(charges_source_file.path, 'r') as file:
            standard_topology = load(file)
            charges = standard_topology['atom_charges']
    # In some ocasions, charges may come inside a raw charges file
    elif charges_source_file.filename == RAW_CHARGES_FILENAME:
        charges = get_raw_charges(charges_source_file.path)
    # In some ocasions, charges may come inside a topology which can be parsed through pytraj
    elif charges_source_file.is_pytraj_supported():
        charges = get_topology_charges(charges_source_file.path)
        # DANI: De momento ya no generaré más charges.txt ahora que las cargas estan en la topologia json
        #generate_raw_energies_file(charges)
    elif charges_source_file.format == 'tpr':
        charges = get_tpr_charges(charges_source_file.path)
    else:
        raise ValueError('Charges file (' + charges_source_file.filename + ') is in a non supported format')
    return charges

# Given a topology which includes charges
# Extract those charges and save them in a list to be returned
# Use different tools, since ones may fail where others success
# Note that this not a matter of the format only, but also of the format version
# There are different versions of .top files, for instance
def get_topology_charges (topology_filename : str) -> list:
    try:
        topology_charges = get_topology_charges_pytraj(topology_filename)
    except Exception as err:
        print(err)
        print('The canonical charges mining (pytraj) failed. Retrying with alternative mining (mdanalysis)')
        topology_charges = get_topology_charges_mdanalysis(topology_filename)
    return topology_charges

# Get topology charges using pytraj
# Supported formats (tested): prmtop, top, psf (standard psf, not from DESRES)
# Supported formats (not tested): mol2, cif, sdf
# Non supported formats: mae, tpr, pdb (has no charges)
def get_topology_charges_pytraj (topology_filename : str) -> list:
    topology = pt.load_topology(filename=topology_filename)
    # WARNING: We must convert this numpy ndarray to a normal list
    # Otherwise the search by index is extremly ineficient
    topology_charges = list(topology.charge)
    return topology_charges

# Get topology charges using mdanalysis
def get_topology_charges_mdanalysis (topology_filename : str) -> list:
    parser = TOPParser(topology_filename)
    topology = parser.parse()
    charges = list(topology.charges.values)
    return charges

# Write the raw charges file from a list of charges
def generate_raw_energies_file (charges : list, filename : str = RAW_CHARGES_FILENAME):
    with open(filename, 'w') as file:
        for charge in charges:
            file.write("{:.6f}".format(charge) + '\n')

# Given a raw file with listed charges
# Extract those charges and save them in a list to be returned
def get_raw_charges (topology_filename : str) -> list:
    charges = []
    with open(topology_filename, 'r') as file:
        lines = file.readlines()
        for line in lines:
            charges.append(float(line))
    return charges

# Given a tpr file, extract charges in a list
# Try 2 different methods and 1 of them should work
def get_tpr_charges (topology_filename : str) -> list:
    try:
        charges = get_tpr_charges_mdanalysis(topology_filename)
    except:
        print('WARNING: mdanalysis failed to extract charges. Using manual extraction...')
        charges = get_tpr_charges_manual(topology_filename)
    return charges

# This works for the new tpr format (tested in 122)
def get_tpr_charges_manual (topology_filename : str) -> list:
    charges = []
    # Read the tpr file making a 'dump'
    process = run([
        GROMACS_EXECUTABLE,
        "dump",
        "-s",
        topology_filename,
        "-quiet"
    ], stdout=PIPE, stderr=PIPE)
    readable_tpr = process.stdout.decode()
    # Mine the atomic charges
    for line in readable_tpr.split('\n'):
        # Skip everything which is not atomic charges data
        if line[0:16] != '            atom':
            continue
        # Parse the line to get only charges
        search = re.search(r"q=([0-9e+-. ]*),", line)
        if search:
            charges.append(float(search[1]))
    if len(charges) == 0:
        error_logs = process.stderr.decode()
        print(error_logs)
        raise SystemExit('Charges extraction from tpr file has failed')
    return charges

# This works for the old tpr format (tested in 112)
def get_tpr_charges_mdanalysis (topology_filename : str) -> list:
    parser = TPRParser(topology_filename)
    topology = parser.parse()
    charges = list(topology.charges.values)
    return charges