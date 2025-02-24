from model_workflow.utils.topology_converter import to_MDAnalysis_topology
from model_workflow.utils.auxiliar import save_json
from model_workflow.utils.type_hints import *
import numpy as np
import MDAnalysis


def lipid_order (
    input_trajectory_filepath : str,
    topology_file : 'File',
    output_analysis_filepath : str,
    membrane_map: dict,):
    print('-> Running lipid order analysis')

    if membrane_map['n_mems'] == 0:
        print(' No membranes found in the structure. Skipping analysis.')
        return
    
    mda_top = to_MDAnalysis_topology(topology_file.absolute_path)
    u = MDAnalysis.Universe(mda_top, input_trajectory_filepath)
    order_parameters_dict = {}
    for ref_data in membrane_map['references'].values():
        res = u.select_atoms(f'resname {ref_data["resname"]}').residues[0]
        carbon_groups = get_all_acyl_chains(res)
        order_parameters_dict[ref_data['resname']] = {}
        # For every 'tail'
        for chain_idx, group in enumerate(carbon_groups):
            atoms = res.universe.atoms[group]
            C_names = [atom.name for atom in atoms]
            # Find all C-H bonds indices
            ch_pairs = find_ch_bonds(u, ref_data["resname"], C_names)
            # Initialize the order parameters to sum over the trajectory
            order_parameters = []
            costheta_sums = {C_name: np.zeros(len(ch_pairs[C_name]['C'])) for C_name in C_names}
            n = 0
            for ts in u.trajectory[0:10001:100]:
                for C_name in C_names:
                    d = u.atoms[ch_pairs[C_name]['C']].positions - u.atoms[ch_pairs[C_name]['H']].positions
                    costheta_sums[C_name] += d[:,2]**2/np.linalg.norm(d, axis=1)**2
                n += 1
            order_parameters = 1.5 * np.array([costheta_sums[C_name].mean() for C_name in C_names])/n - 0.5
            serrors = 1.5 * np.array([costheta_sums[C_name].std() for C_name in C_names])/n
            order_parameters_dict[ref_data['resname']][str(chain_idx)] = {
                'atoms': C_names,
                'avg': order_parameters.tolist(),
                'std': serrors.tolist()}
    # Save the data
    data = { 'data': order_parameters_dict}
    save_json(data, output_analysis_filepath)


def get_all_acyl_chains(residue: 'MDAnalysis.Residue') -> list:
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
            if current_atom.element == 'C':
                group.add(current_atom.index)
                # Add all bonded carbon atoms to the visit queue
                for bond in current_atom.bonds:
                    for bonded_atom in bond.atoms:
                        if (bonded_atom.element == 'C' and 
                            bonded_atom.index not in visited):
                            to_visit.append(bonded_atom)
        return list(group)

    # Get all Carbon atoms in the residue
    carbon_atoms = residue.atoms.select_atoms('element C')
    visited = set()
    for at in carbon_atoms:
        if 'H' not in at.bonded_atoms.elements:
            visited.add(at.index)
    carbon_groups = []
    # Find all distinct groups of connected carbons
    for carbon in carbon_atoms:
        if carbon.index not in visited:
            # Get all carbons connected to this one
            connected_carbons = explore_carbon_group(carbon, visited)
            if len(connected_carbons)>6:  # Only add non-empty groups
                carbon_groups.append(sorted(connected_carbons))
    return carbon_groups


def find_ch_bonds(universe, lipid_resname, atom_names):
    """
    Find all carbon-hydrogen bonds in a lipid molecule using topology information.
    
    Args:
        universe: MDAnalysis Universe object with topology information
        lipid_resname: Residue name of the lipid (default: 'DPPC')
    
    Returns:
        dict: Dict for each carbon, of tuples containing (carbon_index, hydrogen_index) for each C-H bond
    """
    # Get all carbons in the lipid
    carbons = universe.select_atoms(f'resname {lipid_resname} and name ' + ' '.join(atom_names))
    # Get how many carbons with different name that we will average over later
    ch_pairs = {}
    for name in atom_names:
        ch_pairs[name] = { 'C':[], 'H':[] }

    # Iterate through each carbon
    for carbon in carbons:
        # Get all bonds for this carbon
        for bond in carbon.bonds:
            # Check if the bond is between C and H
            atoms = bond.atoms
            at1_nm = atoms[0].name
            at2_nm = atoms[1].name
            if ('C' in at1_nm and 'H' in at2_nm):
                # Store the indices making sure carbon is first
                ch_pairs[carbon.name]['C'].append(atoms[0].index)
                ch_pairs[carbon.name]['H'].append(atoms[1].index)
            elif 'H' in at1_nm and 'C' in at2_nm:
                ch_pairs[carbon.name]['C'].append(atoms[1].index)
                ch_pairs[carbon.name]['H'].append(atoms[0].index)
    # Convert lists to numpy arrays
    for name in ch_pairs:
        ch_pairs[name]['C'] = np.array(ch_pairs[name]['C'])
        ch_pairs[name]['H'] = np.array(ch_pairs[name]['H'])
    return ch_pairs

