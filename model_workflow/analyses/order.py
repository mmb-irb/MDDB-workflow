from model_workflow.utils.topology_converter import to_MDAnalysis_topology
from model_workflow.utils.auxiliar import save_json
from model_workflow.utils.type_hints import *
from biobb_mem.gromacs.gmx_order import gmx_order
import MDAnalysis


def lipid_order (
    input_structure_filepath : str,
    input_trajectory_filepath : str,
    topology_file : 'File',
    output_analysis_filepath : str,
    membrane_map: dict,):
    print('-> Running lipid order analysis')

    if membrane_map['n_mems'] == 0:
        print(' No membranes found in the structure. Skipping analysis.')
        return

    mda_top = to_MDAnalysis_topology(topology_file.absolute_path)
    u = MDAnalysis.Universe(mda_top, input_structure_filepath)
    membrane_map['references'] = ['DPPC']  # Placeholder
    for resname in range(membrane_map['references']):
        # Select one residue, we have confirm before in the membrane mapping that all the lipids are the same
        res = u.select_atoms(f'resname {resname}').residues[0]
        carbon_groups = get_all_carbon_groups(res)

    # Save the data
    data = { 'data':{
        'order': 1,
        }
    }
    save_json(data, output_analysis_filepath)


def get_all_carbon_groups(residue: 'MDAnalysis.Residue') -> list:
    """
    Finds all groups of connected Carbon atoms within a residue, including cyclic structures.
    
    Parameters
    ----------
    residue : MDAnalysis.core.groups.Residue
        The residue to analyze
        
    Returns
    -------
    list
        A list of lists, where each inner list contains atom indices of 
        connected Carbon atoms forming a distinct group
    """
    def explore_carbon_group(start_atom, visited):
        """Helper function to explore a connected group of carbons"""
        to_visit = [start_atom]
        group = set()
        while to_visit:
            current_atom = to_visit.pop(0)
            if current_atom.index in visited:
                continue
            visited.add(current_atom.index)
            if current_atom.name.startswith('C'):
                group.add(current_atom.index)
                # Add all bonded carbon atoms to the visit queue
                for bond in current_atom.bonds:
                    for bonded_atom in bond.atoms:
                        if (bonded_atom.name.startswith('C') and 
                            bonded_atom.index not in visited):
                            to_visit.append(bonded_atom)
        return list(group)

    # Get all Carbon atoms in the residue
    carbon_atoms = residue.atoms.select_atoms('name C*')
    visited = set()
    carbon_groups = []
    # Find all distinct groups of connected carbons
    for carbon in carbon_atoms:
        if carbon.index not in visited:
            # Get all carbons connected to this one
            connected_carbons = explore_carbon_group(carbon, visited)
            # Only add groups > 6: Acyl tails and non-cyclic structures
            if len(connected_carbons)>6:  
                carbon_groups.append(connected_carbons)
    return carbon_groups