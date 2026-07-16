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
        print('   No membranes found, skipping lipid interactions analysis')
        return
    if universe.select_atoms('protein').n_atoms == 0:
        print('   No protein found, skipping lipid interactions analysis')
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
    data = {'data': data, 'version': '0.2.0'}
    output_analysis_filepath = f'{output_directory}/{OUTPUT_LIPID_INTERACTIONS_FILENAME}'
    save_json(data, output_analysis_filepath)


def _build_frame_task(frame_index, universe, protein_residx, lipid_group):
    """Seek the frame in the main process and return a picklable snapshot.

    Uses MDAnalysis.Merge to create a minimal self-contained universe that holds
    only the protein atoms and the lipids already within 6 Å of the protein for
    this frame.  Workers receive this small object instead of the full universe,
    avoiding fork-COW memory copies proportional to trajectory size.
    """
    universe.trajectory[frame_index]
    protein_ag = universe.select_atoms('protein')
    lipid_near_prot = universe.select_atoms(
        '(around 6 protein) and group lipid_group and not protein',
        lipid_group=lipid_group
    )
    original_lipid_indices = lipid_near_prot.indices
    n_protein_atoms = protein_ag.n_atoms
    # Merge into a small single-frame universe with current positions
    if len(lipid_near_prot) > 0:
        merged = MDAnalysis.Merge(protein_ag, lipid_near_prot)
    else:
        merged = MDAnalysis.Merge(protein_ag)
    return (merged.atoms, protein_residx, original_lipid_indices, n_protein_atoms)


def _per_frame_contacts(task) -> list[tuple[int, int, list[int]]]:
    """Per-frame worker for parallel analysis.
    Receives a merged single-frame snapshot (no trajectory attached) and computes
    the lipid atoms near each protein residue.
    Returns a list of (i, residx, atom_indices) where atom_indices are indices
    in the *original* universe so the caller can reconstruct AtomGroups.
    """
    merged_atoms, protein_residx, original_lipid_indices, n_protein_atoms = task
    merged_protein = merged_atoms[:n_protein_atoms]
    merged_lipids = merged_atoms[n_protein_atoms:]
    frame_results = []
    for i, (residx, merged_res) in enumerate(zip(protein_residx, merged_protein.residues)):
        res_atoms = merged_res.atoms
        lipid_near_res = merged_lipids.select_atoms(
            'around 6 global group residuegroup',
            residuegroup=res_atoms)
        local_idx = lipid_near_res.indices - n_protein_atoms
        frame_results.append((i, int(residx), original_lipid_indices[local_idx].tolist()))
    return frame_results


def _accumulate_results(frame_iter, accumulate_fn):
        """Sequentially accumulate results from a frame iterator (either map or imap)."""
        for frame_results in frame_iter:
            for i, residx, atom_indices in frame_results:
                accumulate_fn(i, residx, frozenset(atom_indices))


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

    # Build tasks lazily: main process seeks each frame and merges relevant atoms.
    # Workers receive small picklable snapshots — no fork-COW on the full universe.
    print(f'   Iterating {len(frame_indices)} frames with step {frame_step} and {len(protein_residx)} protein residues')
    tasks = (_build_frame_task(fi, universe, protein_residx, lipid_group) for fi in frame_indices)

    if n_jobs != 1:
        n_workers = multiprocessing.cpu_count() if n_jobs == -1 else n_jobs
        with multiprocessing.Pool(n_workers) as worker_pool:
            _accumulate_results(worker_pool.imap(_per_frame_contacts, tasks, chunksize=1), accumulate_fn)
    else:
        _accumulate_results(map(_per_frame_contacts, tasks), accumulate_fn)
    print('   Finished iterating frames and accumulating results')
    return protein_residx, frame_count


def aa_lipid_interactions(universe: 'MDAnalysis.Universe', lipid_map: list, iterator: callable):
    """Atomistic lipid-protein interactions analysis."""
    lipids_residx = [residue_idx for lipid in lipid_map for residue_idx in lipid['residue_indices']]
    inchikeys = [lipid['inchikey'] for lipid in lipid_map for residue_idx in lipid['residue_indices']]
    lipid_group = universe.select_atoms(f'resindex {" ".join(map(str, lipids_residx))}')
    has_multiple_lipid_types = len(lipid_map) > 1

    # Pre-compute headgroup and tail index sets per lipid residue (topology is static)
    head_indices = {}  # lipid_residx -> frozenset of headgroup atom indices
    tail_indices = {}  # lipid_residx -> frozenset of acyl-chain atom indices
    for lipid_residx, inchikey in zip(lipids_residx, inchikeys):
        residue = universe.select_atoms(f'resindex {lipid_residx}').residues[0]
        if inchikey == 'HVYWMOMLDIMFJA-DPAQBDIFSA-N':
            # Special case for cholesterol: use the O-H atom as headgroup and the rest as tail
            head_ag = residue.atoms.select_atoms('element O and bonded element H')
            tail_ag = residue.atoms - head_ag
        else:
            head_ag, tail_ag = get_head_tail_split(residue)
        head_indices[lipid_residx] = frozenset(head_ag.indices.tolist())
        tail_indices[lipid_residx] = frozenset(tail_ag.indices.tolist())
    # Map each atom index to its lipid residue for fast per-frame contact lookup
    atom_to_lipid_resindex = {int(atom.index): atom.residue.resindex for atom in lipid_group.atoms}
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
    contacts_all_head = np.zeros(n_prot) if has_multiple_lipid_types else None
    contacts_all_tail = np.zeros(n_prot) if has_multiple_lipid_types else None

    def accumulate(i, residx, atom_idx_set):
        seen_head = set()
        seen_tail = set()
        seen_all_head = False
        seen_all_tail = False
        lipid_resindices = {atom_to_lipid_resindex[idx] for idx in atom_idx_set if idx in atom_to_lipid_resindex}
        for lipid_residx in lipid_resindices:
            k = lipid_type_map[lipid_residx]
            if head_indices[lipid_residx] & atom_idx_set:
                seen_head.add(k)
                seen_all_head = True
            if tail_indices[lipid_residx] & atom_idx_set:
                seen_tail.add(k)
                seen_all_tail = True
        for k in seen_head:
            contacts_head[i, k] += 1
        for k in seen_tail:
            contacts_tail[i, k] += 1
        if has_multiple_lipid_types:
            if seen_all_head:
                contacts_all_head[i] += 1
            if seen_all_tail:
                contacts_all_tail[i] += 1

    protein_residx, frame_count = iterator(
        lipid_group=lipid_group, accumulate_fn=accumulate)

    # Divide by frame count to get fraction-of-frames probability
    contacts_head /= frame_count
    contacts_tail /= frame_count
    if has_multiple_lipid_types:
        contacts_all_head /= frame_count
        contacts_all_tail /= frame_count

    # Save the data
    data = {'residue_indices': protein_residx.tolist()}
    if has_multiple_lipid_types:
        data['all_lipids'] = {
            'head': contacts_all_head.tolist(),
            'tail': contacts_all_tail.tolist(),
        }
    for k, lipid in enumerate(lipid_map):
        data[lipid['inchikey']] = {
            'head': contacts_head[:, k].tolist(),
            'tail': contacts_tail[:, k].tolist(),
        }

    return data


def cg_lipid_interactions(universe, lipid_group, iterator):
    """Coarse-grain lipid-protein interactions analysis."""
    print('   Performing coarse-grain lipid-protein interactions analysis')
    lipid_resnames = set(residue.resname for residue in lipid_group.residues)
    has_multiple_lipid_types = len(lipid_resnames) > 1
    n_protein = universe.select_atoms('protein').residues.n_residues

    # Pre-compute headgroup and tail index sets per CG lipid residue (topology is static)
    head_indices = {}  # lipid_residx -> frozenset of headgroup atom indices
    tail_indices = {}  # lipid_residx -> frozenset of acyl-chain atom indices
    lipid_resindex_to_resname = {}
    for residue in lipid_group.residues:
        head_ag, tail_ag = get_cg_head_tail_split(residue)
        head_indices[residue.resindex] = frozenset(head_ag.indices.tolist())
        tail_indices[residue.resindex] = frozenset(tail_ag.indices.tolist())
        lipid_resindex_to_resname[residue.resindex] = residue.resname
    if all(len(t) == 0 for t in tail_indices.values()):
        warn('Coarse-grain analysis not implemented for lipids without tails following the C1A, C2A, etc. naming convention.')
        return {}
    atom_to_lipid_resindex = {int(atom.index): atom.residue.resindex for atom in lipid_group.atoms}
    # Contact probability arrays: fraction of frames where ANY bead of a lipid part
    # is within 6 Å of a given protein residue (binary per frame, per lipid type)
    contacts_head = {resname: np.zeros(n_protein) for resname in lipid_resnames}
    contacts_tail = {resname: np.zeros(n_protein) for resname in lipid_resnames}
    contacts_all_head = np.zeros(n_protein) if has_multiple_lipid_types else None
    contacts_all_tail = np.zeros(n_protein) if has_multiple_lipid_types else None

    def accumulate(i, residx, atom_idx_set):
        seen_head = set()
        seen_tail = set()
        seen_all_head = False
        seen_all_tail = False
        lipid_resindices = {atom_to_lipid_resindex[idx] for idx in atom_idx_set if idx in atom_to_lipid_resindex}
        for lipid_residx in lipid_resindices:
            resname = lipid_resindex_to_resname[lipid_residx]
            if head_indices[lipid_residx] & atom_idx_set:
                seen_head.add(resname)
                seen_all_head = True
            if tail_indices[lipid_residx] & atom_idx_set:
                seen_tail.add(resname)
                seen_all_tail = True
        for resname in seen_head:
            contacts_head[resname][i] += 1
        for resname in seen_tail:
            contacts_tail[resname][i] += 1
        if has_multiple_lipid_types:
            if seen_all_head:
                contacts_all_head[i] += 1
            if seen_all_tail:
                contacts_all_tail[i] += 1

    protein_residx, frame_count = iterator(
        lipid_group=lipid_group, accumulate_fn=accumulate)

    # Divide by frame count to get fraction-of-frames probability
    for resname in lipid_resnames:
        contacts_head[resname] /= frame_count
        contacts_tail[resname] /= frame_count

    if has_multiple_lipid_types:
        contacts_all_head /= frame_count
        contacts_all_tail /= frame_count

    # Save the data in the same format as atomistic analysis
    data = {'residue_indices': protein_residx.tolist()}
    if has_multiple_lipid_types:
        data['all_lipids'] = {
            'head': contacts_all_head.tolist(),
            'tail': contacts_all_tail.tolist(),
        }
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
    print('   Performing united-atom lipid-protein interactions analysis')
    lipid_resnames = set(residue.resname for residue in lipid_group.residues)
    has_multiple_lipid_types = len(lipid_resnames) > 1
    n_protein = universe.select_atoms('protein').residues.n_residues

    # Pre-compute headgroup and tail index sets per UA lipid residue
    head_indices = {}  # lipid_resindex -> frozenset of headgroup atom indices
    tail_indices = {}  # lipid_resindex -> frozenset of acyl-chain atom indices
    lipid_resindex_to_resname = {}
    for residue in lipid_group.residues:
        residue.atoms.select_atoms('name C*').elements = 'C'
        chains = get_all_acyl_chains(residue, removeHs=False)
        tail_idx = [idx for chain in chains for idx in chain]
        if tail_idx:
            tail_ag = universe.atoms[tail_idx]
        else:
            tail_ag = residue.atoms[[]]
        head_ag = residue.atoms - tail_ag
        head_indices[residue.resindex] = frozenset(head_ag.indices.tolist())
        tail_indices[residue.resindex] = frozenset(tail_ag.indices.tolist())
        lipid_resindex_to_resname[residue.resindex] = residue.resname

    if all(len(t) == 0 for t in tail_indices.values()):
        warn('United-atom analysis could not identify acyl chains. '
             'Please check if the lipid topology contains bond information.')
        return {}
    atom_to_lipid_resindex = {int(atom.index): atom.residue.resindex for atom in lipid_group.atoms}

    # Contact probability arrays: fraction of frames where ANY atom of a lipid part
    # is within 6 Å of a given protein residue (binary per frame, per lipid type)
    contacts_head = {resname: np.zeros(n_protein) for resname in lipid_resnames}
    contacts_tail = {resname: np.zeros(n_protein) for resname in lipid_resnames}
    contacts_all_head = np.zeros(n_protein) if has_multiple_lipid_types else None
    contacts_all_tail = np.zeros(n_protein) if has_multiple_lipid_types else None

    def accumulate(i, residx, atom_idx_set):
        seen_head = set()
        seen_tail = set()
        seen_all_head = False
        seen_all_tail = False
        lipid_resindices = {atom_to_lipid_resindex[idx] for idx in atom_idx_set if idx in atom_to_lipid_resindex}
        for lipid_residx in lipid_resindices:
            resname = lipid_resindex_to_resname[lipid_residx]
            if head_indices[lipid_residx] & atom_idx_set:
                seen_head.add(resname)
                seen_all_head = True
            if tail_indices[lipid_residx] & atom_idx_set:
                seen_tail.add(resname)
                seen_all_tail = True
        for resname in seen_head:
            contacts_head[resname][i] += 1
        for resname in seen_tail:
            contacts_tail[resname][i] += 1
        if has_multiple_lipid_types:
            if seen_all_head:
                contacts_all_head[i] += 1
            if seen_all_tail:
                contacts_all_tail[i] += 1

    protein_residx, frame_count = iterator(
        lipid_group=lipid_group, accumulate_fn=accumulate)

    # Divide by frame count to get fraction-of-frames probability
    for resname in lipid_resnames:
        contacts_head[resname] /= frame_count
        contacts_tail[resname] /= frame_count

    if has_multiple_lipid_types:
        contacts_all_head /= frame_count
        contacts_all_tail /= frame_count

    # Save the data in the same format as other analyses
    data = {'residue_indices': protein_residx.tolist()}
    if has_multiple_lipid_types:
        data['all_lipids'] = {
            'head': contacts_all_head.tolist(),
            'tail': contacts_all_tail.tolist(),
        }
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
