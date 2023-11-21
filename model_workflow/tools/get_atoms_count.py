import json
from pathlib import Path

# Load the reference for protein residue names
resources = str(Path(__file__).parent.parent / "resources") # / "utils"
residues_source = resources + '/residues.json'

# This script is used to count different type of atoms and residues in a pdb topology
# input_topology_filename - The name string of the input topology file (path)
# Tested supported formats are .pdb
def get_atoms_count (
    input_topology_filename : str,
    ) -> tuple:

    # Number of system atoms
    systats = 0

    # Number of protein atoms
    protats = 0

    # Number of protein residues
    prot = 0

    # Number of membrane phospholipid chains
    dppc = 0

    # Number of solven atoms
    sol = 0

    # Number of sodium atoms
    na = 0

    # Number of chlorine atoms
    cl = 0

    # List all possible aminoacids
    aminoacids = None
    with open(residues_source, 'r') as file:
        aminoacids = json.load(file)

    # Read the first frame pdb file and get the desired data
    # Residue count in the pdb file is reset after residue 9999
    # Atoms count in the pdb file is reset after atom 99999
    with open(input_topology_filename, 'r') as file:
        
        # Keep track of the current residue
        current_residue_number = None
        
        for line in file:
            if line.split()[0] == 'ATOM':
                systats += 1
                
                # Get the residue number in this line
                residue_number = int(line[22:26])
                
                # If the residue number in this line is different to the current residue...
                # Get its data and set it as the current
                if current_residue_number != residue_number:
                    current_residue_number = residue_number
                    
                    # Get the name of this residue
                    residue_name = line[17:20].strip().upper()
                    
                    # If this is an aminoacid
                    if (residue_name in aminoacids):
                        prot += 1
                        residue_name = 'PROT'
                    # DANI: Esto está fuertemente hardcodeado y fallará muchas veces
                    # DANI: Habría que montar una librería con nombres de residuos de membrana
                    elif(residue_name == 'DPP'):
                        dppc += 1
                        
                if residue_name == 'PROT':
                    protats += 1
                    
                if residue_name == 'SOL' or residue_name == 'WAT':
                    sol += 1
                        
                if residue_name == 'NA' or residue_name == 'NA+' or residue_name == 'K' or residue_name == 'K+':
                    na += 1
                    
                if residue_name == 'CL' or residue_name == 'CL-':
                    cl += 1

    # Display it
    print('Number of system atoms: ' + str(systats))
    print('Number of protein atoms: ' + str(protats))
    print('Number of protein residues: ' + str(prot))
    print('Number of membrane phospholipid chains: ' + str(dppc))
    print('Number of solvent atoms: ' + str(sol))
    print('Number of sodium atoms: ' + str(na))
    print('Number of chlorine atoms: ' + str(cl))

    # Return gromacs logs
    return (systats, protats, prot, dppc, sol, na, cl)