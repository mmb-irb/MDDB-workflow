from mddb_workflow.tools.get_reduced_trajectory import calculate_frame_step
from mddb_workflow.utils.auxiliar import save_json, load_json, warn
from mddb_workflow.utils.constants import OUTPUT_LIPID_INTERACTIONS_FILENAME, LIPIDS_RESIDUE_NAMES
from mddb_workflow.utils.mda_spells import get_head_tail_split, get_cg_head_tail_split, get_all_acyl_chains
from mddb_workflow.utils.type_hints import *
import numpy as np
import multiprocessing
from functools import partial
import MDAnalysis


def lipid_interactions(
    membrane_map: dict,
    universe: 'MDAnalysis.Universe',
    output_directory: str,
    inchikey_map: list[dict],
    cg_residues: list[int],
    snapshots: int,
    frames_limit: int = 100,
):
    """Lipid-protein interactions analysis."""
    if membrane_map is None or membrane_map['n_mems'] == 0:
        print('-> No membranes found, skipping lipid interactions analysis')
        return
    if universe.select_atoms('protein').n_atoms == 0:
        print('-> No protein found, skipping lipid interactions analysis')
        return

    is_united_atom = False
    if len(cg_residues) > 0:
        lipids = universe.select_atoms(f'resname {" ".join(LIPIDS_RESIDUE_NAMES)}')
        # Heuristic: if the residues have more than 40 atoms, assume United Atom.
        is_united_atom = any(residue.atoms.n_atoms > 40 for residue in lipids.residues)
    else:
        lipids = [ref for ref in inchikey_map if ref['is_lipid']]

    # Fixed frame step iterator so we do not have to pass all the arguments down
    iterator = partial(
        _iterate_protein_residues,
        universe=universe,
        snapshots=snapshots,
        frames_limit=frames_limit,
        n_jobs=multiprocessing.cpu_count(),
    )

    # Choose the method according to the type of lipids present
    if is_united_atom:
        data = ua_lipid_interactions(universe, lipids, iterator)
    elif len(cg_residues) > 0:
        data = cg_lipid_interactions(universe, lipids, iterator)
    else:
        data = aa_lipid_interactions(universe, lipids, iterator)
    if not data:
        return

    # Wrap the data in a dictionary
    data = {'data': data, 'version': '0.1.0'}
    output_analysis_filepath = f'{output_directory}/{OUTPUT_LIPID_INTERACTIONS_FILENAME}'
    save_json(data, output_analysis_filepath)


def _per_frame_contacts(frame_index, universe, protein_residx, lipid_group):
    """Per-frame worker for parallel analysis.
    Seeks to *frame_index*, computes the lipid atoms near each protein residue,
    and returns a list of (i, residx, atom_indices) – one entry per protein residue.
    Atom indices are plain ints so the result is fully picklable.
    """
    universe.trajectory[frame_index]
    lipid_near_prot = universe.select_atoms(
        '(around 6 protein) and group lipid_group and not protein',
        lipid_group=lipid_group)
    frame_results = []
    for i, residx in enumerate(protein_residx):
        residue_atoms = universe.select_atoms(f'resindex {residx}')
        lipid_near_res = lipid_near_prot.select_atoms(
            'around 6 global group residuegroup',
            residuegroup=residue_atoms)
        frame_results.append((i, int(residx), lipid_near_res.indices.tolist()))
    return frame_results


def _iterate_protein_residues(universe, snapshots, frames_limit, lipid_group, accumulate_fn, n_jobs=1):
    """Shared scaffold: iterate frames and protein residues, calling accumulate_fn
    for each (frame, protein_residue_index, lipid_near_res) triplet.
    Returns (protein_residx, frame_count).

    When n_jobs != 1 frames are processed in parallel using multiprocessing.Pool.
    Set n_jobs=-1 to use all available CPU cores.
    """
    frame_step, frame_count = calculate_frame_step(snapshots, frames_limit)
    protein_residx = universe.select_atoms('protein').residues.resindices
    frame_indices = list(range(0, snapshots, frame_step))

    if n_jobs != 1:
        # Parallelization with multiprocessing (per-frame fashion)
        # https://userguide.mdanalysis.org/stable/examples/analysis/custom_parallel_analysis.html#Parallelization-in-a-simple-per-frame-fashion
        n_workers = multiprocessing.cpu_count() if n_jobs == -1 else n_jobs
        # Fixed version of run_per_frame where only the frame index changes
        run_per_frame = partial(
            _per_frame_contacts,
            universe=universe,
            protein_residx=protein_residx,
            lipid_group=lipid_group,
        )
        with multiprocessing.Pool(n_workers) as worker_pool:
            all_frame_results = worker_pool.map(run_per_frame, frame_indices)
        # Reconstruct AtomGroups from returned indices and call accumulate_fn
        for frame_results in all_frame_results:
            for i, residx, atom_indices in frame_results:
                lipid_near_res = universe.atoms[atom_indices]
                accumulate_fn(i, residx, lipid_near_res)
    else:
        # Serial fallback
        for frame_index in frame_indices:
            for i, residx, atom_indices in _per_frame_contacts(frame_index, universe, protein_residx, lipid_group):
                lipid_near_res = universe.atoms[atom_indices]
                accumulate_fn(i, residx, lipid_near_res)

    return protein_residx, frame_count


def aa_lipid_interactions(universe, lipid_map, iterator):
    """Atomistic lipid-protein interactions analysis."""
    lipids_residx = [residue_idx for lipid in lipid_map for residue_idx in lipid['residue_indices']]
    lipid_group = universe.select_atoms(f'resindex {" ".join(map(str, lipids_residx))}')

    # Pre-compute headgroup and tail AtomGroups per lipid residue (topology is static)
    head_atoms = {}  # lipid_residx -> AtomGroup of headgroup atoms
    tail_atoms = {}  # lipid_residx -> AtomGroup of acyl-chain atoms
    for lipid_residx in lipids_residx:
        residue = universe.select_atoms(f'resindex {lipid_residx}').residues[0]
        head_ag, tail_ag = get_head_tail_split(residue)
        head_atoms[lipid_residx] = head_ag
        tail_atoms[lipid_residx] = tail_ag
    # Map each lipid residue index to its lipid type (inchikey index) for binary-per-frame scoring
    lipid_type_map = {}
    for k, lipid in enumerate(lipid_map):
        for residx in lipid['residue_indices']:
            lipid_type_map[residx] = k

    # Contact probability arrays: fraction of frames where ANY atom of a lipid part
    # is within 6 Å of a given protein residue (binary per frame, per lipid type)
    n_prot = universe.select_atoms('protein').residues.n_residues
    n_lipid_types = len(lipid_map)
    contacts_head = np.zeros((n_prot, n_lipid_types))
    contacts_tail = np.zeros((n_prot, n_lipid_types))

    def accumulate(i, residx, lipid_near_res):
        # Use sets to score each lipid type at most once per frame
        seen_head = set()
        seen_tail = set()
        for lipid_residx in lipid_near_res.residues.resindices:
            k = lipid_type_map[lipid_residx]
            if len(head_atoms[lipid_residx] & lipid_near_res) > 0:
                seen_head.add(k)
            if len(tail_atoms[lipid_residx] & lipid_near_res) > 0:
                seen_tail.add(k)
        for k in seen_head:
            contacts_head[i, k] += 1
        for k in seen_tail:
            contacts_tail[i, k] += 1

    protein_residx, frame_count = iterator(
        lipid_group=lipid_group, accumulate_fn=accumulate)

    # Divide by frame count to get fraction-of-frames probability
    contacts_head /= frame_count
    contacts_tail /= frame_count

    # Save the data
    data = {'residue_indices': protein_residx.tolist()}
    for k, lipid in enumerate(lipid_map):
        data[lipid['inchikey']] = {
            'head': contacts_head[:, k].tolist(),
            'tail': contacts_tail[:, k].tolist(),
        }

    return data


def cg_lipid_interactions(universe, lipid_group, iterator):
    """Coarse-grain lipid-protein interactions analysis."""
    print('-> Performing coarse-grain lipid-protein interactions analysis')
    lipid_resnames = set(residue.resname for residue in lipid_group.residues)
    n_protein = universe.select_atoms('protein').residues.n_residues

    # Pre-compute headgroup and tail AtomGroups per CG lipid residue (topology is static)
    head_atoms = {}  # lipid_residx -> AtomGroup of headgroup atoms
    tail_atoms = {}  # lipid_residx -> AtomGroup of acyl-chain atoms
    for residue in lipid_group.residues:
        head_atoms[residue.resindex], tail_atoms[residue.resindex] = get_cg_head_tail_split(residue)
    # If no tail atoms found, skip the analysis (probably a false positive lipid)
    if all(len(tail) == 0 for tail in tail_atoms.values()):
        warn('Coarse-grain analysis not implemented for lipids without tails following the C1A, C2A, etc. naming convention.')
        return {}
    # Contact probability arrays: fraction of frames where ANY bead of a lipid part
    # is within 6 Å of a given protein residue (binary per frame, per lipid type)
    contacts_head = {resname: np.zeros(n_protein) for resname in lipid_resnames}
    contacts_tail = {resname: np.zeros(n_protein) for resname in lipid_resnames}

    def accumulate(i, residx, lipid_near_res):
        # Use sets to score each lipid type at most once per frame
        seen_head = set()
        seen_tail = set()
        for resname in lipid_resnames:
            near = lipid_near_res.select_atoms(f'resname {resname}')
            for lipid_res in near.residues:
                if len(head_atoms[lipid_res.resindex] & near) > 0:
                    seen_head.add(resname)
                if len(tail_atoms[lipid_res.resindex] & near) > 0:
                    seen_tail.add(resname)
        for resname in seen_head:
            contacts_head[resname][i] += 1
        for resname in seen_tail:
            contacts_tail[resname][i] += 1

    protein_residx, frame_count = iterator(
        lipid_group=lipid_group, accumulate_fn=accumulate)

    # Divide by frame count to get fraction-of-frames probability
    for resname in lipid_resnames:
        contacts_head[resname] /= frame_count
        contacts_tail[resname] /= frame_count

    # Save the data in the same format as atomistic analysis
    data = {'residue_indices': protein_residx.tolist()}
    for resname in lipid_resnames:
        data[resname] = {
            'head': contacts_head[resname].tolist(),
            'tail': contacts_tail[resname].tolist(),
        }

    return data


def ua_lipid_interactions(universe, lipid_group, iterator):
    """United-atom lipid-protein interactions analysis.

    Similar to CG interactions but uses bond-based acyl-chain detection
    with ``removeHs=False`` since united-atom topologies have no explicit
    hydrogen atoms on carbon chains.
    """
    print('-> Performing united-atom lipid-protein interactions analysis')
    lipid_resnames = set(residue.resname for residue in lipid_group.residues)
    n_protein = universe.select_atoms('protein').residues.n_residues

    # Pre-compute headgroup and tail AtomGroups per UA lipid residue
    head_atoms = {}  # lipid_resindex -> AtomGroup of headgroup atoms
    tail_atoms = {}  # lipid_resindex -> AtomGroup of acyl-chain atoms
    for residue in lipid_group.residues:
        # Override elements for C* atoms (UA topology may lack proper elements)
        residue.atoms.select_atoms('name C*').elements = 'C'
        # Use removeHs=False since UA has no explicit hydrogens
        chains = get_all_acyl_chains(residue, removeHs=False)
        tail_idx = [idx for chain in chains for idx in chain]
        if tail_idx:
            tail_ag = universe.atoms[tail_idx]
        else:
            tail_ag = residue.atoms[[]]  # empty AtomGroup
        head_ag = residue.atoms - tail_ag
        head_atoms[residue.resindex] = head_ag
        tail_atoms[residue.resindex] = tail_ag

    # If no tail atoms found, skip the analysis
    if all(len(tail) == 0 for tail in tail_atoms.values()):
        warn('United-atom analysis could not identify acyl chains. '
             'Please check if the lipid topology contains bond information.')
        return {}

    # Contact probability arrays: fraction of frames where ANY atom of a lipid part
    # is within 6 Å of a given protein residue (binary per frame, per lipid type)
    contacts_head = {resname: np.zeros(n_protein) for resname in lipid_resnames}
    contacts_tail = {resname: np.zeros(n_protein) for resname in lipid_resnames}

    def accumulate(i, residx, lipid_near_res):
        # Use sets to score each lipid type at most once per frame
        seen_head = set()
        seen_tail = set()
        for resname in lipid_resnames:
            near = lipid_near_res.select_atoms(f'resname {resname}')
            for lipid_res in near.residues:
                if len(head_atoms[lipid_res.resindex] & near) > 0:
                    seen_head.add(resname)
                if len(tail_atoms[lipid_res.resindex] & near) > 0:
                    seen_tail.add(resname)
        for resname in seen_head:
            contacts_head[resname][i] += 1
        for resname in seen_tail:
            contacts_tail[resname][i] += 1

    protein_residx, frame_count = iterator(
        lipid_group=lipid_group, accumulate_fn=accumulate)

    # Divide by frame count to get fraction-of-frames probability
    for resname in lipid_resnames:
        contacts_head[resname] /= frame_count
        contacts_tail[resname] /= frame_count

    # Save the data in the same format as other analyses
    data = {'residue_indices': protein_residx.tolist()}
    for resname in lipid_resnames:
        data[resname] = {
            'head': contacts_head[resname].tolist(),
            'tail': contacts_tail[resname].tolist(),
        }

    return data


def plot_lipid_interactions(output_analysis_filepath: str):
    """Plot lipid-protein interactions analysis results."""
    import matplotlib.pyplot as plt
    data = load_json(output_analysis_filepath)['data']

    protein_residues = data['residue_indices']
    lipid_keys = [key for key in data.keys() if key != 'residue_indices']

    # Create a figure for the plot
    plt.figure(figsize=(10, 6))

    for lipid_key in lipid_keys:
        interactions = data[lipid_key]
        if isinstance(interactions, dict):
            # New head/tail dict format
            for part, values in interactions.items():
                plt.plot(protein_residues, values, label=f'{lipid_key} ({part})')
        else:
            # Old flat format
            plt.plot(protein_residues, interactions, label=lipid_key)

    plt.xlabel('Protein Residue Indices')
    plt.ylabel('Interaction Count (Normalized)')
    plt.title('Lipid-Protein Interactions')
    plt.legend(); plt.grid(True); plt.tight_layout(); plt.show()
