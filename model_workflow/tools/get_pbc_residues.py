from typing import List

# Given a structure and a VMD selection, return the list of residue indices
def get_pbc_residues (structure : 'Structure', input_pbc_selection : str) -> List[int]:
    if not input_pbc_selection:
        return []
    selection = structure.select(input_pbc_selection, syntax='vmd')
    pbc_residues = structure.get_selection_residue_indices(selection)
    print(' "' + input_pbc_selection + '" -> ' + str(len(pbc_residues)) + ' residues')
    return pbc_residues
