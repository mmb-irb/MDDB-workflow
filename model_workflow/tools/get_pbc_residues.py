from model_workflow.utils.type_hints import *

# Given a structure and a VMD selection, return the list of residue indices
def get_pbc_residues (structure : 'Structure', input_pbc_selection : str) -> List[int]:
    if not input_pbc_selection:
        return []
    selection = structure.select(input_pbc_selection, syntax='vmd')
    pbc_residues = structure.get_selection_residue_indices(selection)
    print(f'PBC residues "{input_pbc_selection}" -> {len(pbc_residues)} residues')
    return pbc_residues
