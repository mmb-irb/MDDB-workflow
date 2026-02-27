from mddb_workflow.tools.get_reduced_trajectory import calculate_frame_step
from mddb_workflow.utils.auxiliar import save_json, load_json, natural_sort_key, warn
from mddb_workflow.utils.constants import LIPIDS_RESIDUE_NAMES, OUTPUT_LIPID_ORDER_FILENAME
from mddb_workflow.utils.mda_spells import get_acyl_chain_atom_names, get_cg_acyl_chains
from mddb_workflow.utils.type_hints import *
import numpy as np
import MDAnalysis
import re
import os
import yaml
import gorder


def lipid_order(
    membrane_map: dict,
    universe: 'MDAnalysis.Universe',
    topology_file: 'File',
    output_directory: str,
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
        print('-> No membranes found, skipping lipid order analysis')
        return
    frame_step, _ = calculate_frame_step(snapshots, frames_limit)
    # Choose the method according to the type of lipids present
    if len(cg_residues) > 0:
        # Check supported formats: https://vachalab.github.io/gorder-manual/other_structure.html
        if topology_file.format not in ['tpr', 'gro']:
            warn('Lipid order analysis requires a TPR o GRO topology file when using cg_selection.')
            return {}
        # We will save one example residue used to extract the carbon chains
        lipids = []
        for resname in LIPIDS_RESIDUE_NAMES:
            lipid_ag = universe.select_atoms(f'resname {resname}')
            if lipid_ag.n_atoms > 0:
                lipids.append(lipid_ag.residues[0])
        # Heuristic: if the residues have more than 40 atoms, assume United Atom.
        is_united_atom = any(residue.atoms.n_atoms > 40 for residue in lipids)
        if is_united_atom:
            order_parameters_dict = ua_lipid_order(universe, topology_file, output_directory, frame_step, lipids)
        else:
            order_parameters_dict = cg_lipid_order(universe, topology_file, output_directory, frame_step, lipids)
    else:
        lipids = [ref for ref in inchikey_map if ref['is_lipid']]
        order_parameters_dict = aa_lipid_order(universe, snapshots, frame_step, lipids)
    # Output can be empty due to hardcoded lipid residue names for CG/UA
    if not order_parameters_dict:
        warn('Something went wrong with the lipid order analysis, no order parameters were calculated. '
             'This can be due to no lipids being found with the hardcoded residue names. '
             'Please check the logs for warnings.')
        return
    # Save the data
    data = {'data': order_parameters_dict, 'version': '0.1.0'}
    output_analysis_filepath = f'{output_directory}/{OUTPUT_LIPID_ORDER_FILENAME}'
    save_json(data, output_analysis_filepath)


def aa_lipid_order(universe, snapshots, frame_step, lipids):
    """Calculate lipid order parameters for all-atom simulations."""
    order_parameters_dict = {}
    for ref in lipids:
        inchikey = ref['generated_inchikey']
        # Take the first residue of the reference
        res = universe.residues[ref['residue_indices'][0]]
        # If not the cholesterol inchikey
        if inchikey != 'HVYWMOMLDIMFJA-DPAQBDIFSA-N':
            # Lipids can have multiple acyl chains
            tail_C_names = get_acyl_chain_atom_names(res, sort_key=natural_sort_key)
        else:
            tail_C_names = [sorted(
                [a.name for a in res.atoms.select_atoms('element C and bonded element H')],
                key=natural_sort_key)]

        order_parameters_dict[inchikey] = {}
        # For every 'tail'
        for tail_idx, C_names in enumerate(tail_C_names):
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
            order_parameters_dict[inchikey][str(tail_idx)] = {
                'atoms': C_names,
                'avg': order_parameters.tolist(),
                'std': serrors.tolist()}
    return order_parameters_dict


def cg_lipid_order(universe, topology_file, output_directory, frame_step, lipids):
    """Calculate lipid order parameters for coarse-grain simulations using gorder.

    Runs gorder's CGOrder analysis on the universe trajectory, then parses the
    YAML output into the same format as the all-atom lipid_order output:
        { inchikey: { "0": { atoms: [...], avg: [...], std: [...] }, ... } }

    Each chain is identified by the trailing letter of the CG bead name
    (e.g. C1A, C2A → chain 'A'; C1B, C2B → chain 'B').
    """
    yaml_out = os.path.join(output_directory, 'lorder.yaml')
    analysis = gorder.Analysis(
        structure=topology_file.absolute_path,
        trajectory=universe.trajectory.filename,
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
    if os.path.exists(yaml_out):
        os.remove(yaml_out)

    # Build tails from universe: resname -> list of ordered bead-name lists (one per chain)
    tails: dict[str, list[list[str]]] = {}
    for residue in lipids:
        tails[residue.resname] = get_cg_acyl_chains(residue)

    # Pattern for extracting bead name from gorder bond labels
    _tail_bond_re = re.compile(r'^[A-Z0-9]+\s+\w+\s+\(\d+\)\s+-\s+[A-Z0-9]+\s+(C\d+[A-Z])\s+\(\d+\)')
    order_parameters_dict = {}
    for resname, res_data in gorder_data.items():
        if resname in ('average order',):
            continue
        if not isinstance(res_data, dict) or 'order parameters' not in res_data:
            continue

        # Parse per-bead order parameters from bond labels
        parsed: dict[str, tuple[float, float]] = {}
        for bond_label, bond_data in res_data['order parameters'].items():
            m = _tail_bond_re.match(bond_label)
            if m is None:
                continue
            bead_name = m.group(1)   # e.g. "C1A"
            total = bond_data.get('total', {})
            value = total.get('mean', float('nan')) if isinstance(total, dict) else float(total)
            error = total.get('error', 0.0) if isinstance(total, dict) else 0.0
            parsed[bead_name] = (value, error)

        if not parsed:
            continue

        # Map beads to pre-built chains
        order_parameters_dict[resname] = {}
        for chain_idx, chain_beads in enumerate(tails.get(resname, [])):
            chain_atoms, chain_avgs, chain_stds = [], [], []
            for bead_name in chain_beads:
                if bead_name in parsed:
                    val, err = parsed[bead_name]
                    chain_atoms.append(bead_name)
                    chain_avgs.append(val)
                    chain_stds.append(err)
            if chain_atoms:
                order_parameters_dict[resname][str(chain_idx)] = {
                    'atoms': chain_atoms,
                    'avg': chain_avgs,
                    'std': chain_stds,
                }
    return order_parameters_dict


def ua_lipid_order(universe, topology_file, output_directory, frame_step, lipids):
    """Calculate lipid order parameters for united-atom simulations using gorder.

    Runs gorder's UAOrder analysis on the universe trajectory, then parses the
    YAML output into the same format as the other lipid_order outputs:
        { resname: { "0": { atoms: [...], avg: [...], std: [...] } } }

    All heavy atoms with C-H order contributions are collected per residue into
    a single chain (index "0"), preserving the order they appear in the YAML.
    """
    yaml_out = os.path.join(output_directory, 'lorder.yaml')
    analysis = gorder.Analysis(
        structure=topology_file.absolute_path,
        trajectory=universe.trajectory.filename,
        analysis_type=gorder.analysis_types.UAOrder('@membrane'),
        estimate_error=gorder.estimate_error.EstimateError(),
        step=frame_step,
        output_yaml=yaml_out,
        handle_pbc=False,
    )
    results = analysis.run()
    results.write()

    with open(yaml_out, 'r') as f:
        gorder_data = yaml.safe_load(f)
    if os.path.exists(yaml_out):
        os.remove(yaml_out)

    # Build tails: resname -> list of ordered atom-name arrays (one per acyl chain)
    tails: dict[str, list] = {}
    for residue in lipids:
        residue.atoms.select_atoms('name C*').elements = 'C'  # Override element to Carbon
        tails[residue.resname] = get_acyl_chain_atom_names(residue, removeHs=False)

    # Parse YAML: for each resname build a flat map atom_name -> (mean, error)
    order_parameters_dict = {}
    _atom_re = re.compile(r'^(\S+)\s+(\S+)\s+\((\d+)\)$')

    for resname, res_data in gorder_data.items():
        if resname == 'average order':
            continue
        if not isinstance(res_data, dict) or 'order parameters' not in res_data:
            continue

        parsed: dict[str, tuple[float, float]] = {}
        for atom_label, atom_data in res_data['order parameters'].items():
            m = _atom_re.match(atom_label)
            if m is None:
                continue
            atom_name = m.group(2)
            total = atom_data.get('total', {})
            value = total.get('mean', float('nan')) if isinstance(total, dict) else float(total)
            error = total.get('error', 0.0) if isinstance(total, dict) else 0.0
            parsed[atom_name] = (value, error)

        if not parsed:
            continue

        order_parameters_dict[resname] = {}
        for chain_idx, tail_atoms in enumerate(tails.get(resname, [])):
            chain_atoms, chain_avgs, chain_stds = [], [], []
            for atom_name in tail_atoms:
                if atom_name in parsed:
                    val, err = parsed[atom_name]
                    chain_atoms.append(atom_name)
                    chain_avgs.append(val)
                    chain_stds.append(err)
            if chain_atoms:
                order_parameters_dict[resname][str(chain_idx)] = {
                    'atoms': chain_atoms,
                    'avg': chain_avgs,
                    'std': chain_stds,
                }

    return order_parameters_dict


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
