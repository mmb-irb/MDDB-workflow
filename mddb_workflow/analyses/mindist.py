import itertools
from os.path import abspath
import numpy as np
from scipy.spatial import cKDTree
from tqdm import tqdm

from mddb_workflow.utils.auxiliar import TestFailure, warn
from mddb_workflow.utils.constants import CROSS_PBC_FLAG
from mddb_workflow.utils.file import File
from mddb_workflow.utils.mdt_spells import iterate_xtc_frames
from mddb_workflow.utils.selections import Selection
from mddb_workflow.tools.generate_topology import format_selection
from mddb_workflow.utils.type_hints import *

# All 26 neighbour cell offsets, in units of box vectors
NEIGHBOUR_CELL_SHIFTS = np.array(
    [shift for shift in itertools.product((-1, 0, 1), repeat=3) if shift != (0, 0, 0)],
    dtype=np.float64)

# Distance (nm) below which atoms in cross-PBC contact are considered to overlap
OVERLAP_DISTANCE = 0.1

# Minimum size in Angstroms of the focused region, to avoid zooming too much on a few atoms
MIN_FOCUS_SIZE = 10

# Distance (nm) a fragment center must be closer to a periodic image of the reference than to the reference itself
# to be considered escaped from the rest of the system
ESCAPE_MARGIN = 1.0


def find_frame_escaped_fragments(
    reference_center : np.ndarray,
    fragment_centers : np.ndarray,
    box_vectors : np.ndarray,
    margin : float,
) -> list[tuple[int, float, float]]:
    """Find fragments whose center is considerably closer to a periodic image of the reference than to the reference itself.

    The reference is the geometric center of all non-PBC atoms, so it does not depend on where the system is in the box.
    Return a list of (fragment index, distance to the reference, distance to the closest reference image).
    A fragment is escaped when its distance to a reference image is shorter than the distance to the reference minus the margin.
    Only images in the first neighbour cells are checked, but a fragment further away is still closer to one of them.
    """
    shifts = NEIGHBOUR_CELL_SHIFTS @ box_vectors
    deltas = fragment_centers - reference_center
    own_distances = np.linalg.norm(deltas, axis=1)
    # The distance to a reference image shifted by s is the norm of delta - s
    image_distances = np.linalg.norm(deltas[:, np.newaxis, :] - shifts[np.newaxis, :, :], axis=2).min(axis=1)
    escaped_indices = np.where(own_distances - image_distances > margin)[0]
    return [ (int(index), own_distances[index], image_distances[index]) for index in escaped_indices ]


def find_frame_cross_contacts(
    positions : np.ndarray,
    box_vectors : np.ndarray,
    cutoff : float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Find contacts across periodic boundaries in a frame, for any (orthogonal or triclinic) box.

    Return the contact pairs as an (n, 2) array of position indices, sorted by
    distance, along with their distances through the periodic boundary and the
    cell shifts as an (n, 3) array in units of box vectors.
    The second atom of a pair shifted by its cell shift is in contact with the first atom.

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
    copy_cell_shifts = []
    for cell_shift, shift in zip(NEIGHBOUR_CELL_SHIFTS, NEIGHBOUR_CELL_SHIFTS @ box_vectors):
        shifted = positions + shift
        near_system = np.all((shifted >= lower_bound) & (shifted <= upper_bound), axis=1)
        if not near_system.any(): continue
        copy_positions.append(shifted[near_system])
        copy_sources.append(np.where(near_system)[0])
        copy_cell_shifts.append(np.tile(cell_shift, (near_system.sum(), 1)))
    no_contacts = (np.empty((0, 2), dtype=np.int64), np.empty(0), np.empty((0, 3), dtype=np.int64))
    if not copy_positions: return no_contacts
    copy_positions = np.concatenate(copy_positions)
    copy_sources = np.concatenate(copy_sources)
    copy_cell_shifts = np.concatenate(copy_cell_shifts).astype(np.int64)
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
    # The copied atom goes second, so the cell shift moves the second atom towards the first one
    pairs = np.stack([pair_neighbours[cross], pair_sources[cross]], axis=1)
    cell_shifts = copy_cell_shifts[pair_copies[cross]]
    # Sort atoms in every pair, so the cell shift is reversed when atoms are swapped
    swapped = pairs[:, 0] > pairs[:, 1]
    pairs[swapped] = pairs[swapped][:, ::-1]
    cell_shifts[swapped] *= -1
    # Every pair is found twice (from the copy of each atom), so keep only its shortest distance
    by_distance = np.argsort(pbc_distances, kind='stable')
    pairs, pbc_distances, cell_shifts = pairs[by_distance], pbc_distances[by_distance], cell_shifts[by_distance]
    _, first_index = np.unique(pairs, axis=0, return_index=True)
    first_index = np.sort(first_index)
    return pairs[first_index], pbc_distances[first_index], cell_shifts[first_index]

# This function has been validated agains different simulations
# Any project using -smp -> no box
# A01GX -> boxed but not centered
# A01KM -> boxed and ceneted, has no contacts
# A01GY -> boxed and centered, non orthogonal, has no contacts
# MCV1900439 -> boxed and centered, has contacts
def check_cross_periodic_contacts(
    input_structure_filename : str,
    input_trajectory_filename : str,
    view_script_filename : str,
    structure : 'Structure',
    pbc_selection : 'Selection',
    mercy : list[str],
    trust: list[str],
    register : 'Register',
    check_selection : str,
    distance_cutoff : float,
    snapshots : int,
    simulation_box : Optional[tuple | str] = None,
    ignore_box : bool = False,
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

    # Skip if the box is to be ignored
    if ignore_box:
        print('Skipping cross-PBC contacts check: simulation box is ignored')
        register.remove_warnings(CROSS_PBC_FLAG)
        remove_view_script(view_script_filename)
        register.update_test(CROSS_PBC_FLAG, 'na')
        return True

    # If the test was passed already
    # Note that a not applicable result is checked again, since it may be applicable now (e.g. box is not ignored anymore)
    if register.tests.get(CROSS_PBC_FLAG, None) is True:
        return True

    register.remove_warnings(CROSS_PBC_FLAG)
    remove_view_script(view_script_filename)

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

    # # Split non-PBC atoms in fragments to check none of them goes away from the rest of the system
    # # The geometric center of every fragment is compared against the geometric center of all non-PBC atoms
    # # If bonds are missing then we can not find fragments so we use chains instead
    # if structure.is_missing_any_bonds():
    #     fragments = [ chain.get_selection() & non_pbc_selection for chain in structure.chains ]
    #     fragments = [ fragment for fragment in fragments if fragment ]
    # else:
    #     fragments = list(structure.find_fragments(non_pbc_selection))
    # fragments_atom_indices = [ np.array(fragment.atom_indices) for fragment in fragments ]
    # # Frames where every fragment escaped, and the distances in the first escaped frame
    # escape_frames = [ [] for _ in fragments ]
    # first_escape_distances = [ None for _ in fragments ]

    violation_frames = []
    # Frame with the largest number of atoms in cross-PBC contact, to be reported as example
    example_frame = None
    example_atoms_count = 0
    example_pairs = example_distances = example_cell_shifts = example_box_vectors = None
    # Shortest cross-PBC contact distance along the whole trajectory
    closest_distance = np.inf

    if not input_trajectory_filename.endswith('.xtc'):
        raise ValueError(f'Cross-PBC contacts check is only supported for xtc trajectories: {input_trajectory_filename}')
    frames = iterate_xtc_frames(input_trajectory_filename)
    pbar = tqdm(frames, total=snapshots, desc=' Frame', unit='frame')

    # Iterate frames
    for frame_idx, (frame_positions, frame_box) in enumerate(pbar):
        # Ignore frames with no unit vectors
        # DANI: Si ya hemos comprobado que tiene caja esto puede pasar?
        if not frame_box.any(): raise RuntimeError('Missing unitcell vectors')
        # Check if there are atoms too close to its period images in this frame
        box_vectors = frame_box.astype(np.float64)   # (3, 3) nm, rows are vectors
        all_positions = frame_positions.astype(np.float64)
        positions = all_positions[non_pbc_indices]  # (n_nonpbc, 3) nm
        # # Check if any fragment is closer to a periodic image of the non-PBC center than to the center itself
        # fragment_centers = np.array([ all_positions[atom_indices].mean(axis=0) for atom_indices in fragments_atom_indices ])
        # escaped = find_frame_escaped_fragments(positions.mean(axis=0), fragment_centers, box_vectors, ESCAPE_MARGIN)
        # for fragment_index, own_distance, image_distance in escaped:
        #     escape_frames[fragment_index].append(frame_idx)
        #     if first_escape_distances[fragment_index] is None:
        #         first_escape_distances[fragment_index] = (own_distance, image_distance)
        # Check if there are atoms too close to its period images in this frame
        pairs, distances, cell_shifts = find_frame_cross_contacts(positions, box_vectors, cutoff_nm)
        if len(pairs) == 0: continue
        violation_frames.append(frame_idx)
        closest_distance = min(closest_distance, distances[0])
        atoms_count = len(np.unique(pairs))
        if atoms_count > example_atoms_count:
            example_frame = frame_idx
            example_atoms_count = atoms_count
            example_pairs, example_distances, example_cell_shifts = pairs, distances, cell_shifts
            example_box_vectors = box_vectors

    # Report any problem
    # escaped_fragment_indices = [ index for index, frames in enumerate(escape_frames) if frames ]
    # if not violation_frames and not escaped_fragment_indices:
    if not violation_frames:
        print(' Test passed: no cross-PBC contacts detected')
        register.update_test(CROSS_PBC_FLAG, True)
        return True
    # Each problem has a short title and a description with the details
    titles, descriptions = [], []
    # Frames and atoms to be shown in the client
    data = {}
    # # Fragments which went away from the rest of the system
    # if escaped_fragment_indices:
    #     titles.append('Some non-PBC fragments go away from the rest of the system')
    #     description = ''
    #     for fragment_index in escaped_fragment_indices:
    #         fragment_name = name_fragment(structure, fragments[fragment_index])
    #         frames = escape_frames[fragment_index]
    #         own_distance, image_distance = first_escape_distances[fragment_index]
    #         description += (f'{fragment_name} in {len(frames)}/{snapshots} frames.'
    #             f' First in frame {frames[0]}\n')
    #     description += 'If these fragments diffuse freely (e.g. an unbound ligand) then they must be part of the PBC selection'
    #     descriptions.append(description)
    # Contacts across periodic boundaries
    if violation_frames:
        n_violations = len(violation_frames)
        titles.append(f'Cross-PBC contacts detected in {n_violations}/{snapshots} frames')
        description = (
            'Non-PBC molecules appear to contact each other across periodic boundaries.\n'
            f'Frame {example_frame} has the most atoms in contact ({example_atoms_count}).'
            f' Closest pairs in this frame:'
        )
        MAX_REPORTED_PAIRS = 5
        for (atom_a, atom_b), distance in zip(example_pairs[:MAX_REPORTED_PAIRS], example_distances):
            label_a = structure.atoms[int(non_pbc_indices[atom_a])].label
            label_b = structure.atoms[int(non_pbc_indices[atom_b])].label
            description += f'\n  {label_a} — {label_b} ({distance * 10:.2f} Å)'
        if len(example_pairs) > MAX_REPORTED_PAIRS:
            description += f'\n  etc. ({len(example_pairs)} pairs)'
        # Atoms overlapping through the boundary are not physical, so the box is likely wrong
        if closest_distance < OVERLAP_DISTANCE:
            description += (f'\nSome atoms overlap their periodic neighbours (down to {closest_distance * 10:.2f} Å).'
                ' This is not physical: the simulation box in the trajectory is probably wrong.')
        descriptions.append(description)
        # Atoms in contact in the example frame
        structure_pairs = non_pbc_indices[example_pairs]
        contact_selection = Selection([ int(index) for index in np.unique(structure_pairs) ])
        # Atoms whose periodic image is in contact with other atoms, grouped by the cell shift of the image
        images = []
        for cell_shift in np.unique(example_cell_shifts, axis=0):
            shift_pairs = structure_pairs[np.all(example_cell_shifts == cell_shift, axis=1)]
            image_selection = Selection([ int(index) for index in np.unique(shift_pairs[:, 1]) ])
            images.append((cell_shift, image_selection))
        data['frame'] = example_frame
        data['selection'] = format_selection(structure, contact_selection)
        data['images'] = [ { 'shift': cell_shift.tolist(), 'selection': format_selection(structure, image_selection) }
            for cell_shift, image_selection in images ]
        # Write a VMD script for the user to see the contacts
        # This is not part of the warning since it is only useful for the workflow user
        write_view_script(view_script_filename, input_structure_filename, input_trajectory_filename,
            example_frame, example_box_vectors, contact_selection, images)
    message = '. '.join(titles)
    description = '\n'.join(descriptions)
    view_hint = f'To see the atoms in contact run:\n  vmd -e {view_script_filename}' if violation_frames else ''
    if CROSS_PBC_FLAG in mercy:
        register.add_warning(CROSS_PBC_FLAG, message, description, data)
        if view_hint: print(f' {view_hint}')
        register.update_test(CROSS_PBC_FLAG, False)
        return False
    raise TestFailure(f'{message}\n{description}\n{view_hint}')


def write_view_script(
    script_filename : str,
    structure_filename : str,
    trajectory_filename : str,
    frame : int,
    box_vectors : np.ndarray,
    contact_selection : 'Selection',
    images : list[tuple[np.ndarray, 'Selection']],
):
    """Write a VMD script to see the atoms in cross-PBC contact in a frame.

    Atoms in contact are shown in CPK and their periodic images are shown in transparent yellow.
    The system is loaded again for every image and moved by the cell shift, and the box is drawn with lines.
    The whole image is shown in yellow lines.
    Note that the VMD periodic display is not used since VMD would misplace images and the box if the box is rotated.
    """
    load_commands = [
        f'mol new {{{abspath(structure_filename)}}} type pdb waitfor all',
        # Frames from the structure file are removed so frame numbers match the trajectory
        'animate delete all',
        f'mol addfile {{{abspath(trajectory_filename)}}} type xtc first {frame} last {frame} waitfor all',
    ]
    commands = [ f'# Cross-PBC contacts in frame {frame}', *load_commands, 'set main [molinfo top]' ]
    # Box vectors in Angstroms
    box_vectors = box_vectors * 10
    commands.append('draw color blue')
    # Every corner of the box is connected to the three next corners
    for shift in itertools.product((0, 1), repeat=3):
        corner = np.array(shift) @ box_vectors
        for axis in range(3):
            if shift[axis] == 1: continue
            next_corner = corner + box_vectors[axis]
            commands.append(f'draw line {to_tcl_vector(corner)} {to_tcl_vector(next_corner)}')
    # Atoms in contact
    commands += [
        'mol representation CPK',
        'mol color Name',
        f'mol selection {{{contact_selection.to_vmd()}}}',
        'mol material Opaque',
        'mol addrep $main',
    ]
    # Periodic images of the atoms in contact, with one representation per cell shift
    for cell_shift, image_selection in images:
        commands += [
            *load_commands,
            'mol modcolor 0 top ColorID 4',
            'mol representation CPK',
            'mol color ColorID 4',
            f'mol selection {{{image_selection.to_vmd()}}}',
            'mol material Transparent',
            'mol addrep top',
            f'[atomselect top all] moveby {to_tcl_vector(cell_shift @ box_vectors)}',
        ]
    # Center the view in the atoms in contact and zoom in
    # Note that every molecule has its own center, so all of them are centered together
    commands += [
        'mol top $main',
        f'set contact [atomselect $main {{{contact_selection.to_vmd()}}}]',
        'proc get_size {atoms} {',
        '    lassign [measure minmax $atoms] lowest highest',
        '    return [tcl::mathfunc::max {*}[vecsub $highest $lowest]]',
        '}',
        'display resetview',
        'set center [list [transoffset [vecscale -1 [measure center $contact]]]]',
        'foreach molecule [molinfo list] { molinfo $molecule set center_matrix $center }',
        f'set zoom [expr {{[get_size [atomselect $main all]] / max([get_size $contact], {MIN_FOCUS_SIZE})}}]',
        'if {$zoom > 1} { scale by $zoom }',
        'puts "Atoms in contact are shown in CPK and their periodic images in transparent yellow"',
    ]
    with open(script_filename, 'w') as script_file:
        script_file.write('\n'.join(commands) + '\n')


def to_tcl_vector(vector : np.ndarray) -> str:
    """Format a vector as a Tcl list, rounded to thousandths."""
    return '{' + ' '.join(f'{value:.3f}' for value in vector) + '}'


def remove_view_script(script_filename : str):
    """Remove the VMD script from a previous run, if any."""
    script_file = File(script_filename)
    if script_file.exists: script_file.remove()


def name_fragment(structure : 'Structure', fragment : 'Selection') -> str:
    """Name a fragment by its residue if it has only one, or by its chains otherwise."""
    residue_indices = structure.get_selection_residue_indices(fragment)
    if len(residue_indices) == 1: return structure.residues[residue_indices[0]].label
    return structure.name_selection(fragment)
