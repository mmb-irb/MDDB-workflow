from mddb_workflow.tools.get_reduced_trajectory import calculate_frame_step
from mddb_workflow.utils.auxiliar import save_json, load_json, natural_sort_key, warn
from mddb_workflow.utils.constants import LIPIDS_RESIDUE_NAMES, OUTPUT_LIPID_ORDER_FILENAME
from mddb_workflow.utils.type_hints import *
import numpy as np
import MDAnalysis
import re
import os


def lipid_order(
    universe: 'MDAnalysis.Universe',
    topology_file: 'File',
    output_directory: str,
    membrane_map: dict,
    inchikey_map: list[dict],
    cg_residues: list[int],
    snapshots: int,
    frames_limit: int = 100
):
    """Calculate lipid order parameters for membranes.
    This function computes the order parameter (S) for lipid acyl chains, defined as:
    S = 0.5*(3*<cos²θ> - 1)
    where θ is the angle between the C-H bond and the membrane normal (z-axis).

    Notes:
    - For non-cholesterol lipids, analyzes all acyl chains separately.
    - For cholesterol, uses carbon atoms bonded to hydrogen.
    - The order parameter is calculated for each carbon atom in the chains.

    """
    if membrane_map is None or membrane_map['n_mems'] == 0:
        print('-> Skipping lipid order analysis')
        return
    frame_step, _ = calculate_frame_step(snapshots, frames_limit)
    if len(cg_residues) > 0:
        order_parameters_dict = cg_lipid_order(universe, topology_file, output_directory, frame_step)
    else:
        order_parameters_dict = aa_lipid_order(universe, inchikey_map, snapshots, frame_step)

    # Save the data
    data = {'data': order_parameters_dict}
    output_analysis_filepath = f'{output_directory}/{OUTPUT_LIPID_ORDER_FILENAME}'
    save_json(data, output_analysis_filepath)


def aa_lipid_order(universe, inchikey_map, snapshots, frame_step):
    """Calculate lipid order parameters for all-atom simulations."""
    order_parameters_dict = {}
    for ref in inchikey_map:
        if not ref['is_lipid']: continue
        inchikey = ref['generated_inchikey']
        # Take the first residue of the reference
        res = universe.residues[ref['residue_indices'][0]]
        # If not the cholesterol inchikey
        if inchikey != 'HVYWMOMLDIMFJA-DPAQBDIFSA-N':
            # Lipids can have multiple acyl chains
            carbon_groups = get_all_acyl_chains(res)
        else:
            carbon_groups = [res.atoms.select_atoms(
                'element C and bonded element H').indices]

        order_parameters_dict[inchikey] = {}
        # For every 'tail'
        for chain_idx, group in enumerate(carbon_groups):
            atoms = res.universe.atoms[group]
            C_names = sorted([atom.name for atom in atoms],
                             key=natural_sort_key)
            # Find all C-H bonds indices
            ch_pairs = find_CH_bonds(universe, ref["residue_indices"], C_names)
            # Initialize the order parameters to sum over the trajectory
            order_parameters = []
            costheta_sums = {C_name: np.zeros(
                len(ch_pairs[C_name]['C'])) for C_name in C_names}
            n = 0
            # Loop over the trajectory
            for ts in universe.trajectory[0:snapshots:frame_step]:
                for C_name in C_names:
                    d = universe.atoms[ch_pairs[C_name]['C']].positions - \
                        universe.atoms[ch_pairs[C_name]['H']].positions
                    costheta_sums[C_name] += d[:, 2]**2 / \
                        np.linalg.norm(d, axis=1)**2
                n += 1
            order_parameters = 1.5 * \
                np.array([costheta_sums[C_name].mean()
                         for C_name in C_names]) / n - 0.5
            serrors = 1.5 * np.array([costheta_sums[C_name].std()
                                     for C_name in C_names]) / n
            order_parameters_dict[inchikey][str(chain_idx)] = {
                'atoms': C_names,
                'avg': order_parameters.tolist(),
                'std': serrors.tolist()}
    return order_parameters_dict


def cg_lipid_order(universe, topology_file, output_directory, frame_step):
    """Calculate lipid order parameters for coarse-grain simulations using gorder.

    Runs gorder's CGOrder analysis on the universe trajectory, then parses the
    YAML output into the same format as the all-atom lipid_order output:
        { inchikey: { "0": { atoms: [...], avg: [...], std: [...] }, ... } }

    Each chain is identified by the trailing letter of the CG bead name
    (e.g. C1A, C2A → chain 'A'; C1B, C2B → chain 'B').
    """
    import yaml
    import gorder

    lipid_residues = [
        resname for resname in LIPIDS_RESIDUE_NAMES
        if universe.select_atoms(f'resname {resname}').n_atoms > 0
    ]

    if not lipid_residues:
        print('-> No lipid residues found for coarse-grain analysis')
        return {}

    trajectory_file = universe.trajectory.filename
    if topology_file.format not in ['tpr', 'gro']:
        # Supported formats: https://vachalab.github.io/gorder-manual/other_structure.html
        warn('Coarse-grain lipid order analysis requires a TPR o GRO topology file.')
        return {}

    yaml_out = os.path.join(output_directory, 'lorder.yaml')

    try:
        analysis = gorder.Analysis(
            structure=topology_file.absolute_path,
            trajectory=trajectory_file,
            analysis_type=gorder.analysis_types.CGOrder('@membrane'),
            estimate_error=gorder.estimate_error.EstimateError(),
            step=frame_step,
            output_yaml=yaml_out,
            handle_pbc=False,
        )
        results = analysis.run()
        results.write()

        with open(yaml_out, 'r') as f:
            gorder_data = yaml.safe_load(f)
    finally:
        if os.path.exists(yaml_out):
            os.remove(yaml_out)

    order_parameters_dict = {}

    # Pattern matching CG tail bead names: C<digit(s)><chain_letter>
    _tail_bond_re = re.compile(r'^[A-Z0-9]+\s+\w+\s+\(\d+\)\s+-\s+[A-Z0-9]+\s+(C\d+([A-Z]))\s+\(\d+\)')

    for resname, res_data in gorder_data.items():
        if resname in ('average order',):
            continue
        if not isinstance(res_data, dict) or 'order parameters' not in res_data:
            continue

        # Group bonds by chain letter, preserving bead order
        chains: dict[str, list[tuple[str, float]]] = {}
        for bond_label, bond_data in res_data['order parameters'].items():
            m = _tail_bond_re.match(bond_label)
            if m is None:
                continue
            bead_name = m.group(1)   # e.g. "C1A"
            chain_letter = m.group(2)  # e.g. "A"
            total = bond_data.get('total', {})
            value = total.get('mean', float('nan')) if isinstance(total, dict) else float(total)
            error = total.get('error', 0.0) if isinstance(total, dict) else 0.0
            chains.setdefault(chain_letter, []).append((bead_name, value, error))

        if not chains:
            continue

        order_parameters_dict[resname] = {}
        for chain_idx, (chain_letter, bonds) in enumerate(sorted(chains.items())):
            order_parameters_dict[resname][str(chain_idx)] = {
                'atoms': [b[0] for b in bonds],
                'avg': [b[1] for b in bonds],
                'std': [b[2] for b in bonds],
            }
    return order_parameters_dict


def get_all_acyl_chains(residue: 'MDAnalysis.Residue') -> list[list[int]]:
    """Find all groups of connected Carbon atoms within a residue, including cyclic structures.

    Args:
    residue (MDAnalysis.core.groups.Residue): The residue to analyze

    Returns:
        list
            A list of lists, where each inner list contains atom indices of
            connected Carbon atoms forming a distinct group

    """
    def explore_carbon_group(start_atom, visited):
        """Explore a connected group of carbons."""
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
    # Remove ester carbons (without hydrogens)
    for at in carbon_atoms:
        if 'H' not in at.bonded_atoms.elements:
            visited.add(at.index)
    carbon_groups = []
    # Find all distinct groups of connected carbons
    for carbon in carbon_atoms:
        if carbon.index not in visited:
            # Get all carbons connected to this one
            connected_carbons = explore_carbon_group(carbon, visited)
            if len(connected_carbons) > 6:  # Only add non-empty groups
                carbon_groups.append(sorted(connected_carbons))
    return carbon_groups


def find_CH_bonds(universe, lipid_resindices, atom_names):
    """Find all carbon-hydrogen bonds in a lipid molecule using topology information.

    Args:
        universe: MDAnalysis Universe object with topology information
        lipid_resindices: List of residue indices corresponding to the lipid
        atom_names: List of carbon atom names to search for C-H bonds

    Returns:
        dict: Dict for each carbon, of tuples containing (carbon_index, hydrogen_index) for each C-H bond

    """
    # Get all carbons in the lipid
    carbons = universe.residues[lipid_resindices].atoms.select_atoms(
        'name ' + ' '.join(atom_names))
    # Get how many carbons with different name that we will average over later
    ch_pairs = {}
    for name in atom_names:
        ch_pairs[name] = {'C': [], 'H': []}

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


def plot_lipid_order(data_filepath: str):
    """Visualize the lipid order analysis by plotting the order parameters.

    Args:
        data_filepath (str): Path to the JSON file containing the lipid order analysis results.

    """
    import matplotlib.pyplot as plt

    # Load the lipid order data
    data = load_json(data_filepath)
    order_parameters_dict = data.get('data', {})

    # Create output directory if it doesn't exist

    # Iterate over each lipid type
    for inchikey, chains in order_parameters_dict.items():
        plt.figure(figsize=(10, 6))
        for chain_idx, chain_data in chains.items():
            avg = chain_data['avg']
            std = chain_data['std']

            # Use the index of the atom to align both chains in the graph
            atom_indices = range(len(avg))

            # Plot the order parameters with error bars for each chain
            eb = plt.errorbar(atom_indices, avg, yerr=std, fmt='o-', label=f'Chain {chain_idx}')
            eb[-1][0].set_linestyle('--')

        plt.xlabel('Atom Index')
        plt.ylabel('Order Parameter (S)')
        plt.title(f'Lipid Order Parameters for {inchikey}')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.legend()
        plt.show()
