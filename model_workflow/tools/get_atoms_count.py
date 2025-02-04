from model_workflow.utils.constants import RESIDUE_NAME_LETTERS

# Count different type of atoms and residues in the structure
def get_atoms_count (structure : 'Structure') -> tuple:

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

    # Iterate the structure residues
    for residue in structure.residues:
        systats += residue.atom_count
        # Add one to the corresponding residue counter
        residue_type = residue.classification
        if residue_type == 'protein':
            prot += 1
            protats += residue.atom_count
        if residue_type == 'fatty' or residue_type == 'steroid':
            dppc += 1
        if residue_type == 'solvent':
            sol += residue.atom_count
        if residue.name == 'NA' or residue.name == 'NA+' or residue.name == 'K' or residue.name == 'K+':
            na += 1
        if residue.name == 'CL' or residue.name == 'CL-':
            cl += 1

    # Display it
    print(f'Number of system atoms: {systats}')
    print(f'Number of protein atoms: {protats}')
    print(f'Number of protein residues: {prot}')
    print(f'Number of lipid residues: {dppc}')
    print(f'Number of solvent atoms: {sol}')
    print(f'Number of sodium atoms: {na}')
    print(f'Number of chlorine atoms: {cl}')

    # Return gromacs logs
    return systats, protats, prot, dppc, sol, na, cl