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


def frame_has_triclinic_cross_contacts(
    positions : np.ndarray,
    box_vectors : np.ndarray,
    cutoff : float,
) -> bool:
    """Check a frame for contacts across periodic boundaries in any (e.g. triclinic) box.

    Positions are used as they are, without wrapping them into the unit cell,
    so any imaging representation works (triclinic, compact, rectangular).
    Periodic copies of the atoms are generated for the 26 neighbour cells and
    only those falling within the cutoff of the atoms bounding box are kept.
    A pair is a cross-PBC contact when an atom is within the cutoff of the
    periodic copy of another atom but not within the cutoff of the atom itself.
    Images beyond the first neighbour cells are not checked.

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
    if not copy_positions: return False
    copy_positions = np.concatenate(copy_positions)
    copy_sources = np.concatenate(copy_sources)
    # Find real atoms within the cutoff of every periodic copy
    neighbour_lists = cKDTree(positions).query_ball_point(copy_positions, cutoff, workers=-1)
    counts = [len(neighbours) for neighbours in neighbour_lists]
    total = sum(counts)
    if total == 0: return False
    pair_sources = np.repeat(copy_sources, counts)
    pair_neighbours = np.fromiter(itertools.chain.from_iterable(neighbour_lists), dtype=np.int64, count=total)
    # Discard pairs which are also in contact directly
    delta = positions[pair_sources] - positions[pair_neighbours]
    direct_distances = np.sqrt(np.einsum('ij,ij->i', delta, delta))
    return bool(np.any(direct_distances >= cutoff))

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
    is_system_centered : Optional[bool] = None,
    is_simulation_box_orthogonal : Optional[bool] = None,
) -> Optional[bool]:
    """Check if non-PBC atoms contact each other across periodic boundaries.

    Algorithm (per frame):
      1. Build a periodic cKDTree of all non-PBC atom positions.
      2. Query only "boundary atoms" (those within distance_cutoff of any box
         face) against the full tree using the minimum-image metric.  Any
         atom found by the tree has d_PBC < cutoff by construction.
      3. For the found pairs compute the direct Euclidean distance.  If
         d_direct >= cutoff the contact only exists through the periodic
         image → cross-PBC contact.

    The cKDTree approach is O(N log N) per frame, replacing the previous
    O(K × N) batched numpy loop that caused memory and speed issues.

    Non-orthogonal boxes use frame_has_triclinic_cross_contacts instead.

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
    if is_system_centered is None:
        print('Skipping cross-PBC contacts check: missing simulation box')
        register.update_test(CROSS_PBC_FLAG, 'na')
        return True

    # Non-orthogonal boxes are checked with periodic copies of the atoms, which needs no centering
    is_triclinic = is_simulation_box_orthogonal == False

    # Skip if the system is not centred in the box
    # The orthogonal check wraps atoms into the box, so this would produce a lot of false positives
    if not is_triclinic and is_system_centered == False:
        print('Skipping cross-PBC contacts check: system not centred in the simulation box')
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
    MAX_REPORTED = 5
    reported = 0

    trajectory = mdt.iterload(input_trajectory_filename, top=input_structure_filename, chunk=1)
    pbar = tqdm(trajectory, total=snapshots, desc=' Frame', unit='frame')

    for frame_idx, frame in enumerate(pbar, 1):

        box = frame.unitcell_lengths
        if box is None:
            continue
        box = box[0].astype(np.float64)   # (3,) nm

        pos = frame.xyz[0]                                       # (n_atoms, 3) float32
        pos_nonpbc = pos[non_pbc_indices].astype(np.float64)    # (n_nonpbc, 3)

        if is_triclinic:
            box_vectors = frame.unitcell_vectors[0].astype(np.float64)   # (3, 3) nm, rows are vectors
            if frame_has_triclinic_cross_contacts(pos_nonpbc, box_vectors, cutoff_nm):
                violation_frames.append(frame_idx)
            continue

        # Positions of non-PBC atoms, wrapped to [0, box) as required by cKDTree
        pos_wrapped = pos_nonpbc % box                           # [0, box)

        # Boundary atoms: within cutoff_nm of any box face
        near_boundary = np.zeros(n_nonpbc, dtype=bool)
        for dim in range(3):
            near_boundary |= pos_wrapped[:, dim] < cutoff_nm
            near_boundary |= pos_wrapped[:, dim] > box[dim] - cutoff_nm

        boundary_local = np.where(near_boundary)[0]   # indices into non_pbc_indices
        if len(boundary_local) == 0:
            continue

        # Build periodic KD-tree: O(N log N)
        tree = cKDTree(pos_wrapped, boxsize=box)

        # Query: for each boundary atom find all non-PBC atoms within cutoff via PBC
        # Returns list[array]: one array of LOCAL indices per boundary atom
        # O(K log N + M) where K = boundary count, M = total found pairs
        pbc_neighbor_lists = tree.query_ball_point(
            pos_wrapped[boundary_local], cutoff_nm, workers=-1
        )

        # Flatten to pair arrays without Python loops over individual atoms
        counts = [len(ns) for ns in pbc_neighbor_lists]
        total  = sum(counts)
        if total == 0:
            continue

        pair_bi = np.repeat(np.arange(len(boundary_local), dtype=np.int32), counts)
        pair_nj = np.fromiter(
            itertools.chain.from_iterable(pbc_neighbor_lists),
            dtype=np.int32, count=total
        )

        # Remove self-pairs
        not_self = boundary_local[pair_bi] != pair_nj
        pair_bi  = pair_bi[not_self]
        pair_nj  = pair_nj[not_self]
        if len(pair_bi) == 0:
            continue

        # Direct (Euclidean) distance between each pair using wrapped positions
        delta    = pos_wrapped[boundary_local[pair_bi]] - pos_wrapped[pair_nj]
        d_direct = np.sqrt(np.einsum('ij,ij->i', delta, delta))

        # The cKDTree guarantees d_PBC < cutoff; flag pairs where d_direct >= cutoff
        cross = d_direct >= cutoff_nm
        if cross.any():
            violation_frames.append(frame_idx)
            # DANI: Esto va bien para debugear pero rompe el tqdm
            # if reported < MAX_REPORTED:
            #     idx    = int(np.argmax(cross))
            #     ai_loc = boundary_local[pair_bi[idx]]
            #     aj_loc = pair_nj[idx]
            #     # Compute PBC distance only for the reported pair
            #     dv     = pos_wrapped[ai_loc] - pos_wrapped[aj_loc]
            #     dv_mic = dv - box * np.round(dv / box)
            #     d_pbc_A  = float(np.sqrt(dv_mic @ dv_mic) * 10)
            #     d_dir_A  = float(d_direct[idx] * 10)
            #     rep_ai = int(non_pbc_indices[ai_loc])
            #     rep_aj = int(non_pbc_indices[aj_loc])
            #     name_a = structure.atoms[rep_ai].label
            #     name_b = structure.atoms[rep_aj].label
            #     print(
            #         f' FAIL: Frame {frame_idx}: cross-PBC contact '
            #         f'"{name_a}" — "{name_b}" '
            #         f'(PBC {d_pbc_A:.2f} Å, direct {d_dir_A:.2f} Å)'
            #     )
            #     reported += 1
            # elif reported == MAX_REPORTED:
            #     print(' etc...')
            #     reported += 1

    if violation_frames:
        n_violations = len(violation_frames)
        message = (
            f'Cross-PBC contacts detected in {n_violations}/{snapshots} frames. '
            'Non-PBC molecules appear to contact each other across periodic boundaries.'
        )
        if CROSS_PBC_FLAG in mercy:
            register.add_warning(CROSS_PBC_FLAG, message)
            register.update_test(CROSS_PBC_FLAG, False)
            return False
        raise TestFailure(message)

    print(' Test passed: no cross-PBC contacts detected')
    register.update_test(CROSS_PBC_FLAG, True)
    return True
