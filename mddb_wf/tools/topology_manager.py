# Topology/Structure management tools

# Import local utils
from mddb_wf.utils.structures import Structure

# Set a pair of independent functions to save and recover chains from a pdb file
# WARNING: These functions must be used only when the pdb has not changed in number of atoms

# Get a list with each atom chain from a pdb
def get_chains (pdb_filename : str) -> list:
    structure = Structure.from_pdb_file(pdb_filename)
    chains = structure.chains
    return chains

# Set each atom chain from a pdb from a list
def set_chains (pdb_filename : str, chains : list):
    structure = Structure.from_pdb_file(pdb_filename)
    structure.chains = chains
    structure.generate_pdb_file(pdb_filename)

# This is the new system to handle the topology
def setup_structure (pdb_filename : str) -> 'Structure':
    # Set the structure
    structure = Structure.from_pdb_file(pdb_filename)
    # Set a reference system to handle conversions to pytraj topology
    # First set the pytraj topology
    pytraj_topology = structure.get_pytraj_topology()
    pytraj_residues = list(pytraj_topology.residues)
    # Set functions to handle the conversions using the pytraj topology
    # Transform a structure residue to the pytraj residue numeration (1, 2, ... n)
    def residue_2_pytraj_residue_index (residue : 'Residue') -> int:
        residue_index = residue.index
        residue_number = residue.number
        residue_name = residue.name[0:3]
        # And check that this residue data matches the pytraj residues data
        pytraj_residue = pytraj_residues[residue_index]
        if (residue_number == pytraj_residue.original_resid and residue_name == pytraj_residue.name):
            return residue_index + 1
        # If not, we must iterate over all pytraj residues to find a match
        for index, pytraj_residue in enumerate(pytraj_residues):
            if (residue_number == pytraj_residue.original_resid and residue_name == pytraj_residue.name):
                return index + 1
        # Return None if there is no match
        return None
    structure.residue_2_pytraj_residue_index = residue_2_pytraj_residue_index
    # Transform a pytraj residue numeration (1, 2, ... n) to the topology residue number
    def pytraj_residue_index_2_residue (pytraj_residue_index : int) -> 'Residue':
        residue_index = pytraj_residue_index - 1
        pytraj_residue = pytraj_residues[residue_index]
        expected_number = pytraj_residue.original_resid
        expected_name = pytraj_residue.name
        # In the canonical way this index is equivalent to the structure resiude index
        if residue_index < len(structure.residues):
            residue = structure.residues[residue_index]
            if (residue.number == expected_number and residue.name[0:3] == expected_name):
                return residue
        # Pytraj index may not match the structure index in caotic topologies
        # (i.e. when heavy atoms and hydrogen are listed independently)
        # When this happens we can try to find the residue by comparing resnum and resname
        # WARNING: Note that this alternative method is nos sensitive to chains or icodes
        # WARNING: This is because pytraj does not deal with chains or icodes
        for residue in structure.residues:
            #print(str(residue.getResnum()) + ' -> ' + str(expectedNumber) + ' / ' + residue.getResname()[0:3] + ' -> ' + expectedName)
            if (residue.number == expected_number and residue.name[0:3] == expected_name):
                return residue
        # Return None if there are no results    
        return None
    structure.pytraj_residue_index_2_residue = pytraj_residue_index_2_residue
    # Set a function to find the absolute atom index of a residue atom provided the atom name
    def get_atom_index (residue : 'Residue', atom_name : str) -> int:
        return next(( atom.index for atom in residue.atoms if atom.name == atom_name ), None)
    structure.get_atom_index = get_atom_index
    # Return the modified structure instance
    return structure