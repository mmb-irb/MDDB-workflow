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

    # Read the first frame pdb file and get the desired data
    # Elements count in the pdb file is reset after element 9999
    # Atoms count in the pdb file is reset after atom 99999
    with open(input_topology_filename, 'r') as file:
        
        # Keep track of the current element
        current_element = 0
        
        # List all possible aminoacids
        aminoacids = ['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE',
        'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
        
        for line in file:
            if(line.split()[0] == 'ATOM'):
                systats += 1
                
                # Get the element number in this line
                element = int(line[22:26])
                
                # If the element in this line is different to the current element...
                # Get its data and set it as the current
                if(current_element != element):
                    current_element = element
                    
                    # Get the name of this element
                    n = line[17:20]
                    
                    # If this is an aminoacid
                    if (n in aminoacids):
                        prot += 1
                        n = 'PROT'
                    elif(n == 'DPP'):
                        dppc += 1
                        
                if(n == 'PROT'):
                    protats += 1
                    
                if(n == 'SOL' or n == 'WAT'):
                    sol += 1
                        
                if(n == ' NA'):
                    na += 1
                    
                if(n == ' CL'):
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