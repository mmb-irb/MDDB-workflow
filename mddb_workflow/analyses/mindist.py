import itertools
import mdtraj as mdt
import numpy as np
from scipy.spatial import cKDTree
from tqdm import tqdm

from mddb_workflow.utils.auxiliar import TestFailure, warn
from mddb_workflow.utils.constants import CROSS_PBC_FLAG
from mddb_workflow.utils.type_hints import *

# All 26 neighbour cell offsets, in units of box vectors
NEIGHBOUR_CELL_SHIFTS = np.array(
    [shift for shift in itertools.product((-1, 0, 1), repeat=3) if shift != (0, 0, 0)],
    dtype=np.float64)


def find_frame_cross_contacts(
    positions : np.ndarray,
    box_vectors : np.ndarray,
    cutoff : float,
) -> tuple[np.ndarray, np.ndarray]:
    """Find contacts across periodic boundaries in a frame, for any (orthogonal or triclinic) box.

    Return the contact pairs as an (n, 2) array of position indices, sorted by
    distance, along with their distances through the periodic boundary.

    Positions are used as they are, without wrapping them into the unit cell,
    so the system does not need to be centred in the box and any imaging
    representation works (triclinic, compact, rectangular).
    Periodic copies of the atoms are generated for the 26 neighbour cells and
    only those falling within the cutoff of the atoms bounding box are kept.
    A pair is a cross-PBC contact when an atom is within the cutoff of the
    periodic copy of another atom but not within the cutoff of the atom itself.
    Images beyond the first neighbour cells are not checked.
    Molecules are expected to be whole and imaged together at this point.

    Made by Claude.
    """
    lower_bound = positions.min(axis=0) - cutoff
    upper_bound = positions.max(axis=0) + cutoff
    copy_positions = []
    copy_sources = []
    for shift in NEIGHBOUR_CELL_SHIFTS @ box_vectors:
        shifted = positions + shift
        near_system = np.all((shifted >= lower_bound) & (shifted <= upper_bound), axis=1)
        if not near_system.any(): continue
        copy_positions.append(shifted[near_system])
        copy_sources.append(np.where(near_system)[0])
    no_contacts = (np.empty((0, 2), dtype=np.int64), np.empty(0))
    if not copy_positions: return no_contacts
    copy_positions = np.concatenate(copy_positions)
    copy_sources = np.concatenate(copy_sources)
    # Find real atoms within the cutoff of every periodic copy
    neighbour_lists = cKDTree(positions).query_ball_point(copy_positions, cutoff, workers=-1)
    counts = [len(neighbours) for neighbours in neighbour_lists]
    total = sum(counts)
    if total == 0: return no_contacts
    pair_copies = np.repeat(np.arange(len(copy_positions)), counts)
    pair_sources = copy_sources[pair_copies]
    pair_neighbours = np.fromiter(itertools.chain.from_iterable(neighbour_lists), dtype=np.int64, count=total)
    # Discard pairs which are also in contact directly
    delta = positions[pair_sources] - positions[pair_neighbours]
    direct_distances = np.sqrt(np.einsum('ij,ij->i', delta, delta))
    cross = direct_distances >= cutoff
    if not cross.any(): return no_contacts
    pbc_delta = copy_positions[pair_copies[cross]] - positions[pair_neighbours[cross]]
    pbc_distances = np.sqrt(np.einsum('ij,ij->i', pbc_delta, pbc_delta))
    pairs = np.sort(np.stack([pair_sources[cross], pair_neighbours[cross]], axis=1), axis=1)
    # Every pair is found twice (from the copy of each atom), so keep only its shortest distance
    by_distance = np.argsort(pbc_distances, kind='stable')
    pairs, pbc_distances = pairs[by_distance], pbc_distances[by_distance]
    _, first_index = np.unique(pairs, axis=0, return_index=True)
    first_index = np.sort(first_index)
    return pairs[first_index], pbc_distances[first_index]

# This function has been validated agains different simulations
# Any project using -smp -> no box
# A01GX -> boxed but not centered
# A01KM -> boxed and ceneted, has no contacts
# A01GY -> boxed and centered, non orthogonal, has no contacts
# MCV1900439 -> boxed and centered, has contacts
def check_cross_periodic_contacts(
    input_structure_filename : str,
    input_trajectory_filename : str,
    structure : 'Structure',
    pbc_selection : 'Selection',
    mercy : list[str],
    trust: list[str],
    register : 'Register',
    check_selection : str,
    distance_cutoff : float,
    snapshots : int,
    simulation_box : Optional[tuple | str] = None,
) -> Optional[bool]:
    """Check if non-PBC atoms contact each other across periodic boundaries.

    Every frame is checked with find_frame_cross_contacts, which works for any
    box shape and does not need the system to be centred in the box.

    distance_cutoff is kept as a fixed parameter (default 5 Å in mwf.py).
    No dynamic scaling is needed: the comparison d_PBC < cutoff ≤ d_direct
    adapts naturally to the box geometry.
    """

    # If it is to be skipped
    if CROSS_PBC_FLAG in trust:
        return True

    # If the test was run already
    if register.tests.get(CROSS_PBC_FLAG, None):
        return True

    register.remove_warnings(CROSS_PBC_FLAG)

    # Skip if the system has no box at all for there is no way to run this test without it
    if simulation_box is None:
        print('Skipping cross-PBC contacts check: missing simulation box')
        register.update_test(CROSS_PBC_FLAG, 'na')
        return True

    full_selection = structure.select(check_selection, syntax='vmd')
    non_pbc_selection = full_selection - pbc_selection

    if len(non_pbc_selection) < 2:
        print('No non-PBC atoms found; skipping cross-PBC contacts check')
        register.update_test(CROSS_PBC_FLAG, 'na')
        return True

    cutoff_nm = distance_cutoff / 10.0   # Å → nm (MDtraj unit)
    non_pbc_indices = np.array(non_pbc_selection.atom_indices, dtype=np.int32)
    n_nonpbc = len(non_pbc_indices)

    print(f'Checking cross-PBC contacts ({n_nonpbc} non-PBC atoms, cutoff {distance_cutoff} Å)')

    violation_frames = []
    # Frame with the largest number of atoms in cross-PBC contact, to be reported as example
    example_frame = None
    example_atoms_count = 0
    example_pairs = example_distances = None

    trajectory = mdt.iterload(input_trajectory_filename, top=input_structure_filename, chunk=1)
    pbar = tqdm(trajectory, total=snapshots, desc=' Frame', unit='frame')

    # Iterate frames
    for frame_idx, frame in enumerate(pbar, 1):
        # Ignore frames with no unit vectors
        # DANI: Si ya hemos comprobado que tiene caja esto puede pasar?
        if frame.unitcell_vectors is None: raise RuntimeError('Missing unitcell vectors')
        # Check if there are atoms too close to its period images in this frame
        box_vectors = frame.unitcell_vectors[0].astype(np.float64)   # (3, 3) nm, rows are vectors
        positions = frame.xyz[0][non_pbc_indices].astype(np.float64)  # (n_nonpbc, 3) nm
        pairs, distances = find_frame_cross_contacts(positions, box_vectors, cutoff_nm)
        if len(pairs) == 0: continue
        violation_frames.append(frame_idx)
        atoms_count = len(np.unique(pairs))
        if atoms_count > example_atoms_count:
            example_frame = frame_idx
            example_atoms_count = atoms_count
            example_pairs, example_distances = pairs, distances

    # Report any problem
    if violation_frames:
        n_violations = len(violation_frames)
        message = (
            f'Cross-PBC contacts detected in {n_violations}/{snapshots} frames. '
            'Non-PBC molecules appear to contact each other across periodic boundaries.\n'
            f' Frame {example_frame} has the most atoms in contact ({example_atoms_count}).'
            f' Closest pairs in this frame:'
        )
        MAX_REPORTED_PAIRS = 5
        for (atom_a, atom_b), distance in zip(example_pairs[:MAX_REPORTED_PAIRS], example_distances):
            label_a = structure.atoms[int(non_pbc_indices[atom_a])].label
            label_b = structure.atoms[int(non_pbc_indices[atom_b])].label
            message += f'\n  {label_a} — {label_b} ({distance * 10:.2f} Å)'
        if len(example_pairs) > MAX_REPORTED_PAIRS:
            message += f'\n  etc. ({len(example_pairs)} pairs)'
        if CROSS_PBC_FLAG in mercy:
            register.add_warning(CROSS_PBC_FLAG, message)
            register.update_test(CROSS_PBC_FLAG, False)
            return False
        raise TestFailure(message)

    print(' Test passed: no cross-PBC contacts detected')
    register.update_test(CROSS_PBC_FLAG, True)
    return True
