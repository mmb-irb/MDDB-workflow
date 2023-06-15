from typing import List

# Given a structure and a VMD selection, return the list of residue indices
def get_pbc_residues (structure : 'Structure', input_pbc_selection : str) -> List[int]:
    selection = structure.select(input_pbc_selection, syntax='vmd')
    return structure.get_selection_residue_indices(selection)
