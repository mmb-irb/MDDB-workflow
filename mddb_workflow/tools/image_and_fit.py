# This script is used to fit a trajectory (i.e. eliminate translation and rotation)
# This process is carried by Gromacs

from collections import deque

import mdtraj as mdt
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from mddb_workflow.utils.structures import Structure
from mddb_workflow.utils.auxiliar import MISSING_TOPOLOGY, MISSING_BONDS, InputError, warn
from mddb_workflow.utils.constants import AUTOMATIC_FLAG
from mddb_workflow.utils.gmx_spells import run_gromacs
from mddb_workflow.utils.mdt_spells import rotate_boxes, iterate_xtc_frames
from mddb_workflow.utils.mdt_spells import read_xtc_sample_frames, is_box_rotated
from mddb_workflow.utils.mdt_spells import guess_translation_fitted, guess_rotation_fitted
from mddb_workflow.analyses.mindist import find_frame_cross_contacts
from mddb_workflow.tools.get_bonds import mine_topology_bonds
from mddb_workflow.utils.type_hints import *

# Set a pair of independent auxiliar functions to save and recover chains from a pdb file
# WARNING: These functions must be used only when the pdb has not changed in number of atoms


def get_chains(pdb_filename: str) -> list:
    """Get a list with each atom chain from a PDB."""
    structure = Structure.from_pdb_file(pdb_filename)
    chains = structure.chains
    return chains


def set_chains(pdb_filename: str, chains: list):
    """Set the chains for a PDB file."""
    structure = Structure.from_pdb_file(pdb_filename)
    structure.chains = chains
    structure.generate_pdb_file(pdb_filename)

# Set some constant values to handle NDX atom index selections
SYSTEM_SELECTION_NAME = 'system'
CENTER_SELECTION_NAME = 'center'
INDEX_FILEPATH = '.index.ndx'

# Distance (nm) below which non-PBC atoms are considered part of the same object when finding the automatic translation
AUTO_TRANSLATION_CUTOFF = 0.5

# Number of frames read at once in the no-jump step
NO_JUMP_CHUNK_SIZE = 100


def image_and_fit(
    input_structure_file: 'File',
    input_trajectory_file: 'File',
    input_topology_file: Union['File', Exception],
    output_structure_file: 'File',
    output_trajectory_file: 'File',
    image: bool,
    fit: bool,
    translation: list[float] | str,
    structure: 'Structure',
    pbc_selection: 'Selection',
    ) -> str:
    """Image and fit a trajectory.

    Args:
        input_structure_file:
            The input structure file. Tested supported formats are .pdb and .tpr.
        input_trajectory_file:
            The input trajectory file. Tested supported formats are .xtc and .trr.
        input_topology_file:
            The input topology file. Tested supported format is .tpr.
            This is optional if there are no PBC residues.
        output_structure_file:
            The output structure file.
        output_trajectory_file:
            The output trajectory file.
        image:
            Set if it must me imaged
        fit:
            Set if it must be fitted
        translation:
            Set how it must be translated during imaging, in nm.
            Use 'auto' to find a translation which keeps non-PBC atoms away from the box boundaries.
        structure:
            The structure of the system. Note that this is an early provisional structure.
        pbc_selection:
            Set a selection to exclude PBC residues from both the centering and the fitting focuses.
            Selection of atoms which are to stay in PBC conditions after imaging.
            Note that this is an early provisional atom selection.

    """
    print(f' Image: {image} | Fit: {fit}')

    if not image and not fit:
        return

    # Using a TPR topology may provide some advantage in the imaging process
    # This includes connectivity and atom mass data
    # This is important to keep atoms together when using -pbc res, mol or whole
    is_tpr_available = input_topology_file != MISSING_TOPOLOGY and input_topology_file.format == 'tpr'
    if not is_tpr_available:
        warn('Since we are missing a TPR file the automatic imaging protocol is more prone to fail.')
    # Check if we have PBC atoms
    # LORE: We were running the '-pbc nojump' step only when no pbc atoms were present in the system
    # LORE: However now we always run a no-jump followed by a '-pbc res' step thus recovering PBC atoms
    has_pbc_atoms = bool(pbc_selection)

    # First of all save chains
    # Gromacs will delete chains so we need to recover them after
    chains_backup = get_chains(input_structure_file.path)

    # Set a custom index file (.ndx) to select reference atoms for centering and fitting steps
    # We will use as center all atoms which are not under PBC
    system_selection = structure.select('all', syntax='vmd')
    center_selection = structure.invert_selection(pbc_selection)
    if not center_selection:
        warn(f'The default selection to center and fit (non PBC) is empty. These steps will be skipped.')
    # Warn the user if the trajectory was already fitted, since it may be a problem for the imaging
    if center_selection:
        warn_broken_pbc(input_trajectory_file, structure, center_selection)
    # Convert both selections to a single ndx file which gromacs can read
    with open(INDEX_FILEPATH, 'w') as file:
        system_selection_ndx = system_selection.to_ndx(SYSTEM_SELECTION_NAME)
        file.write(system_selection_ndx)
        if center_selection:
            center_selection_ndx = center_selection.to_ndx(CENTER_SELECTION_NAME)
            file.write(center_selection_ndx)

    # Use the first frame of the trajectory as structure
    # The input structure may be a different frame, but gromacs uses it as reference for the no-jump and the fit
    reset_structure(input_structure_file.path, input_trajectory_file.path, output_structure_file.path)

    # Set the latest structure and trajectory files
    # They are the ones to be used as inputs for the following steps
    latest_structure = output_structure_file
    latest_trajectory = input_trajectory_file

    # Imaging --------------------------------------------------------------------------------------

    if image:
        print(' Running imaging steps')

        # Find the translation automatically if requested
        if translation == AUTOMATIC_FLAG:
            translation = find_automatic_translation(latest_trajectory, center_selection)

        # Check if coordinates are to be translated
        must_translate = translation != [0, 0, 0]

        # Translate the system if requested
        # Then put residues back in the box, moving each residue as a whole
        # Gromacs applies the translation before putting residues in the box
        if must_translate:
            print(' Running translation')
            run_gromacs(f'trjconv -s {latest_structure.path} \
                -f {latest_trajectory.path} -o {output_trajectory_file.path} \
                -trans {translation[0]} {translation[1]} {translation[2]} -pbc res',
                user_input=SYSTEM_SELECTION_NAME, show_error_logs=True)
            latest_trajectory = output_trajectory_file

            # Select the first frame of the recently imaged trayectory as the new structure
            reset_structure(
                latest_structure.path,
                latest_trajectory.path,
                output_structure_file.path)
            latest_structure = output_structure_file

        # Run an intial -pbc res step just to keep residue atoms together before the no-jump step
        # Note that this can be done only if we have the TPR file
        if is_tpr_available:
            print(' Running initial -pbc res')
            run_gromacs(f'trjconv -s {input_topology_file.path} \
                -f {latest_trajectory.path} -o {output_trajectory_file.path} \
                -pbc res -n {INDEX_FILEPATH}',
                user_input=SYSTEM_SELECTION_NAME, show_error_logs=True)
            latest_trajectory = output_trajectory_file
            # Select the first frame of the recently imaged trayectory as the new structure
            reset_structure(
                latest_structure.path,
                latest_trajectory.path,
                output_structure_file.path)
            latest_structure = output_structure_file

        # !!! We assume that at this point (after a possible translation) the first frame is "correct"
        # This means that the different molecules are placed as desired
        # e.g. two protein monomers are together, and not interacting across boundaries
        # To keep things like this we will run a no-jump step
        # DANI: The nojump step may cause artifacts, specially in already pre-imaged simulations
        if center_selection:
            print(' Running no-jump')
            # LORE: We were using the Gromacs no-jump, but it decides every atom jump independently
            # LORE: Atoms of the same molecule could end in different periodic images when they moved around half box
            # run_gromacs(f'trjconv -s {latest_structure.path} \
            #     -f {latest_trajectory.path} -o {output_trajectory_file.path} \
            #     -pbc nojump', user_input=SYSTEM_SELECTION_NAME, show_error_logs=True)
            # Input and output trajectories may be the same file so we write an intermediate file
            atom_groups = find_no_jump_groups(input_topology_file, structure)
            no_jump_trajectory_file = output_trajectory_file.get_prefixed_file('nojump.')
            no_jump(latest_trajectory, no_jump_trajectory_file, atom_groups)
            no_jump_trajectory_file.rename_to(output_trajectory_file)
            latest_trajectory = output_trajectory_file

            # Place all residues in the box again in case we have PBC residues
            # This is critical to recover all PBC residues which were diluted because of the no-jump step
            # Also center the non-PBC region if there is a center selection
            # We try to do both processes in a single step just to be more efficient
            # WARNING: If there is no center selection then there is no need with the initial -pbc res we should be good
            if has_pbc_atoms:
                print(' Running center -pbc res')
                run_gromacs(f'trjconv -s {latest_structure.path} \
                    -f {latest_trajectory.path} -o {output_trajectory_file.path} \
                    -pbc res -center -n {INDEX_FILEPATH}',
                    user_input=f'{CENTER_SELECTION_NAME} {SYSTEM_SELECTION_NAME}',
                    show_error_logs=True)
                latest_trajectory = output_trajectory_file
            # If there are no PBC residues then just center the system
            else:
                print(' Running center')
                run_gromacs(f'trjconv -s {latest_structure.path} \
                    -f {latest_trajectory.path} -o {output_trajectory_file.path} \
                    -center -n {INDEX_FILEPATH}',
                    user_input=f'{CENTER_SELECTION_NAME} {SYSTEM_SELECTION_NAME}',
                    show_error_logs=True)
                latest_trajectory = output_trajectory_file
        # If there is no center and no PBC then do nothing, but this should never happen

        # Select the first frame of the recently imaged trayectory as the new structure
        reset_structure(latest_structure.path, latest_trajectory.path, output_structure_file.path)
        latest_structure = output_structure_file

    # Fitting --------------------------------------------------------------------------------------

    # Remove translation and rotation in the trajectory using the custom selection as reference
    # WARNING: Sometimes, simulations with membranes do not require this step and they may get worse when fitted
    if fit:
        print(' Running fitting step')

        # Here we use the structure, not the topology, even if it is a TPR
        # Note that using the topology would make useless the last 'nojump' process

        # Gromacs rotates the coordinates but it keeps the original box
        # Then periodic images are no longer valid after the fit
        # To fix it we write the fit in an intermediate file and then we rotate boxes in the output file
        fitted_trajectory_file = output_trajectory_file.get_prefixed_file('fitted.')
        rotated_trajectory_file = output_trajectory_file.get_prefixed_file('rotated.')

        # Run Gromacs
        run_gromacs(f'trjconv -s {latest_structure.path} \
            -f {latest_trajectory.path} -o {fitted_trajectory_file.path} \
            -fit rot+trans -n {INDEX_FILEPATH}',
            user_input=f'{CENTER_SELECTION_NAME} {SYSTEM_SELECTION_NAME}',
            show_error_logs=True)

        # Rotate boxes as coordinates were rotated by the fit
        print(' Rotating boxes')
        rotate_boxes(latest_trajectory, fitted_trajectory_file, rotated_trajectory_file)
        fitted_trajectory_file.remove()
        rotated_trajectory_file.rename_to(output_trajectory_file)

        # If there is no output structure at this time (i.e. there was no imaging) then create it now
        # If there was imaging, we reset again to it its not far from the origin
        # Note that the input structure is not necessarily the output structure in this scenario
        # The fit may have been done using the topology and there may be an offset, so better dump it
        reset_structure(
            input_structure_file.path,
            output_trajectory_file.path,
            output_structure_file.path)

    # Recover chains
    set_chains(output_structure_file.path, chains_backup)

    # Clean up the index file
    # if exists(INDEX_FILEPATH):
    #     remove(INDEX_FILEPATH)

# Auxiliar functions to aid in the automatic imaging process
# Mostly made by Claude

def find_no_jump_groups(topology_file: Union['File', Exception], structure: 'Structure') -> np.ndarray:
    """Set the group of every atom for the no-jump step, as atoms in the same group move together.
    Molecules are used when there are bonds in the topology, and residues otherwise.
    Bonds are not guessed here since guessing them requires the trajectory to be imaged already.
    """
    atom_bonds = mine_topology_bonds(topology_file)
    if atom_bonds and any(bonds == MISSING_BONDS for bonds in atom_bonds):
        print('  Some atoms are missing bonds in the topology')
        atom_bonds = None
    if not atom_bonds:
        print('  No-jump groups: residues')
        atom_residue_indices = np.empty(structure.atom_count, dtype=int)
        for residue in structure.residues:
            atom_residue_indices[residue.atom_indices] = residue.index
        return atom_residue_indices
    if len(atom_bonds) != structure.atom_count:
        raise InputError(f'Topology bonds have {len(atom_bonds)} atoms, but the structure has {structure.atom_count} atoms')
    print('  No-jump groups: molecules from topology bonds')
    pairs = np.array([(atom_index, bonded) for atom_index, bonds in enumerate(atom_bonds) for bonded in bonds], dtype=int)
    pairs = pairs.reshape(-1, 2)
    atom_count = structure.atom_count
    graph = coo_matrix((np.ones(len(pairs), dtype=bool), (pairs[:, 0], pairs[:, 1])), shape=(atom_count, atom_count))
    _, atom_molecules = connected_components(graph, directed=False)
    return atom_molecules


def periodic_mean(fractional: np.ndarray, groups: np.ndarray, group_count: int) -> np.ndarray:
    """Get the mean of fractional displacements in every group, where displacements are periodic.
    Every displacement is converted to an angle, so displacements which differ in whole boxes are the same.
    The output is in fractional coordinates in the range [-0.5, 0.5].
    """
    angles = 2 * np.pi * fractional
    means = np.empty((group_count, 3))
    for axis in range(3):
        sines = np.bincount(groups, weights=np.sin(angles[:, axis]), minlength=group_count)
        cosines = np.bincount(groups, weights=np.cos(angles[:, axis]), minlength=group_count)
        means[:, axis] = np.arctan2(sines, cosines) / (2 * np.pi)
    return means


def closest_images(fractional: np.ndarray, box_vectors: np.ndarray) -> np.ndarray:
    """Get the whole boxes to remove from every fractional displacement to make it as short as possible.
    Rounding fractional coordinates is not enough in triclinic boxes, so neighbour images are also checked.
    This is only done for long displacements, since short ones are already the shortest after rounding.
    """
    images = np.round(fractional)
    remainders = fractional - images
    neighbours = np.array([(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)], dtype=np.float64)
    neighbour_vectors = neighbours @ box_vectors
    shortest_vector = np.linalg.norm(neighbour_vectors, axis=1)[np.any(neighbours != 0, axis=1)].min()
    # A displacement shorter than half the shortest box vector is always the shortest among its images
    long_ones = np.linalg.norm(remainders @ box_vectors, axis=1) >= shortest_vector / 2
    if long_ones.any():
        candidates = remainders[long_ones, None, :] @ box_vectors - neighbour_vectors[None, :, :]
        best_neighbours = np.argmin(np.linalg.norm(candidates, axis=2), axis=1)
        images[long_ones] += neighbours[best_neighbours]
    return images


def no_jump(
    input_trajectory_file: 'File',
    output_trajectory_file: 'File',
    atom_groups: np.ndarray,
):
    """Remove jumps across the box boundaries along the trajectory, keeping the first frame as it is.
    In every frame, each atom is placed in the periodic image closest to where it would be if it moved like its group.
    The movement of every group is found from the previous frame, and it is chosen as the closest to the movement of the whole system.
    This way atoms of a group stay together even if they move around half box, and molecules broken in the input are fixed.
    Any box shape is supported since everything is done in box vector units.
    """
    for trajectory_file in [input_trajectory_file, output_trajectory_file]:
        if trajectory_file.format != 'xtc':
            raise ValueError(f'No-jump is only supported for xtc trajectories: {trajectory_file.path}')
    group_count = atom_groups.max() + 1
    all_atoms = np.zeros(len(atom_groups), dtype=int)
    previous_positions = None
    with mdt.formats.XTCTrajectoryFile(input_trajectory_file.path) as reader, \
        mdt.formats.XTCTrajectoryFile(output_trajectory_file.path, 'w') as writer:
        while True:
            xyz, time, step, boxes = reader.read(n_frames=NO_JUMP_CHUNK_SIZE)
            if len(xyz) == 0: break
            output_xyz = np.empty_like(xyz)
            for i, (positions, box_vectors) in enumerate(zip(xyz, boxes)):
                positions = positions.astype(np.float64)
                if previous_positions is None:
                    previous_positions = positions
                    output_xyz[i] = positions
                    continue
                if not box_vectors.any():
                    raise InputError(f'Frame {step[i]} has no box, so jumps can not be removed')
                box_vectors = box_vectors.astype(np.float64)
                displacements = (positions - previous_positions) @ np.linalg.inv(box_vectors)
                # Movement of the whole system, which may be large if the input was centered in every frame
                system_displacement = periodic_mean(displacements, all_atoms, 1)
                system_displacement = (system_displacement - closest_images(system_displacement, box_vectors))[0]
                # Movement of every group, as close as possible to the movement of the whole system
                group_displacements = periodic_mean(displacements - system_displacement, atom_groups, group_count)
                group_displacements = system_displacement + group_displacements - closest_images(group_displacements, box_vectors)
                # Remove the whole boxes which separate every atom from where its group says it should be
                jumps = closest_images(displacements - group_displacements[atom_groups], box_vectors)
                positions = positions - jumps @ box_vectors
                previous_positions = positions
                output_xyz[i] = positions
            writer.write(output_xyz, time=time, step=step, box=boxes)


def warn_broken_pbc(trajectory_file: 'File', structure: 'Structure', selection: 'Selection'):
    """Warn the user if the trajectory seems to be fitted already in rotation and periodic boundary conditions are thus broken.
    Periodic boundary conditions are kept only if boxes were rotated along with the coordinates.
    Note that some tools (e.g. Gromacs) do not rotate boxes when fitting.
    Note that a translation (e.g. centering) does not break periodic boundary conditions.
    Only a few frames are checked so this is just a guess.
    """
    frames = read_xtc_sample_frames(trajectory_file)
    selected_atom_indices = set(selection.atom_indices)
    residue_atom_indices = [ [ index for index in residue.atom_indices if index in selected_atom_indices ] for residue in structure.residues ]
    residue_atom_indices = [ atom_indices for atom_indices in residue_atom_indices if atom_indices ]
    translation_fitted = guess_translation_fitted(frames, selection.atom_indices)
    rotation_fitted = guess_rotation_fitted(frames, residue_atom_indices)
    print(f' The trajectory seems to be fitted in translation: {translation_fitted} | in rotation: {rotation_fitted}')
    if not rotation_fitted: return
    boxes = [ box for _, box in frames ]
    if not all(box is not None and box.any() for box in boxes):
        warn('The trajectory seems to be fitted already, but it has no box.\n'
            ' Periodic boundary conditions are probably broken, so the imaging will probably fail.')
    elif not any(is_box_rotated(box) for box in boxes):
        warn('The trajectory seems to be fitted already, but its box is not rotated accordingly.\n'
            ' Periodic boundary conditions are probably broken, so the imaging will probably fail.')


def reset_structure(
    input_structure_filepath: str,
    input_trajectory_filepath: str,
    output_structure_filepath: str,
):
    """Get the first frame of a trajectory and use it to reset the structure."""
    print(' Reseting structure after imaging step')
    run_gromacs(f'trjconv -s {input_structure_filepath} \
        -f {input_trajectory_filepath} -o {output_structure_filepath} \
        -dump 0', user_input=SYSTEM_SELECTION_NAME)


def find_automatic_translation(trajectory_file: 'File', selection: 'Selection') -> list[float]:
    """Find a translation which places the box boundaries where no non-PBC object is split.

    Only the first frame of the trajectory is used.
    Atoms in the selection which are close to each other, directly or across the box boundaries, are grouped in objects.
    Then every object is made whole and the boundary along every box vector is placed in the middle of the largest gap between objects.
    """
    if not selection: return [0, 0, 0]
    print(' Finding the automatic translation')
    positions, box_vectors = next(iterate_xtc_frames(trajectory_file.path))
    positions = positions[selection.atom_indices].astype(np.float64)
    box_vectors = box_vectors.astype(np.float64)
    inverse_box = np.linalg.inv(box_vectors)
    # Group atoms in fragments using direct contacts only
    atom_count = len(positions)
    pairs = cKDTree(positions).query_pairs(AUTO_TRANSLATION_CUTOFF, output_type='ndarray')
    graph = coo_matrix((np.ones(len(pairs), dtype=bool), (pairs[:, 0], pairs[:, 1])), shape=(atom_count, atom_count))
    fragment_count, atom_fragments = connected_components(graph, directed=False)
    # Find contacts across the box boundaries and the periodic shift (in box vectors) which brings the second atom next to the first one
    cross_pairs, _, _ = find_frame_cross_contacts(positions, box_vectors, AUTO_TRANSLATION_CUTOFF)
    # If there are no contacts across the box boundaries then nothing is split and there is no need to translate
    if len(cross_pairs) == 0:
        print('  No atoms are split by the box boundaries, no translation is needed')
        return [0, 0, 0]
    deltas = positions[cross_pairs[:, 1]] - positions[cross_pairs[:, 0]]
    shifts = np.rint(deltas @ inverse_box).astype(int)
    links = np.unique(np.column_stack([atom_fragments[cross_pairs[:, 0]], atom_fragments[cross_pairs[:, 1]], shifts]), axis=0)
    neighbours = [[] for _ in range(fragment_count)]
    for fragment, other_fragment, *shift in links:
        neighbours[fragment].append((other_fragment, np.array(shift)))
        neighbours[other_fragment].append((fragment, -np.array(shift)))
    # Join fragments in objects, giving every fragment the periodic shift which makes the object whole
    self_contact_error = InputError('Non-PBC atoms are in contact with their own periodic image in the first frame.\n'
        '  The automatic translation does not work in these cases.')
    fragment_shifts = [None] * fragment_count
    fragment_objects = np.empty(fragment_count, dtype=int)
    object_count = 0
    for seed in range(fragment_count):
        if fragment_shifts[seed] is not None: continue
        fragment_shifts[seed] = np.zeros(3, dtype=int)
        fragment_objects[seed] = object_count
        queue = deque([seed])
        while queue:
            fragment = queue.popleft()
            for other_fragment, shift in neighbours[fragment]:
                expected_shift = fragment_shifts[fragment] + shift
                if fragment_shifts[other_fragment] is None:
                    fragment_shifts[other_fragment] = expected_shift
                    fragment_objects[other_fragment] = object_count
                    queue.append(other_fragment)
                elif not np.array_equal(fragment_shifts[other_fragment], expected_shift):
                    raise self_contact_error
        object_count += 1
    # Get the fractional coordinates of every atom once its object is whole
    fractional = positions @ inverse_box - np.array(fragment_shifts)[atom_fragments]
    atom_objects = fragment_objects[atom_fragments]
    translation = np.zeros(3)
    for axis in range(3):
        values = fractional[:, axis]
        lows = np.full(object_count, np.inf)
        highs = np.full(object_count, -np.inf)
        np.minimum.at(lows, atom_objects, values)
        np.maximum.at(highs, atom_objects, values)
        if np.any(highs - lows >= 1): raise self_contact_error
        # Bring every object interval into [0, 1) and split the ones going beyond 1
        offsets = np.floor(lows)
        lows, highs = lows - offsets, highs - offsets
        crossing = highs > 1
        lows = np.concatenate([lows, np.zeros(crossing.sum())])
        highs = np.concatenate([np.minimum(highs, 1), highs[crossing] - 1])
        # Merge the intervals and find the largest gap between them, going around the box
        order = np.argsort(lows)
        merged = [[lows[order[0]], highs[order[0]]]]
        for low, high in zip(lows[order[1:]], highs[order[1:]]):
            if low <= merged[-1][1]: merged[-1][1] = max(merged[-1][1], high)
            else: merged.append([low, high])
        gaps = [(merged[i][1], merged[i + 1][0]) for i in range(len(merged) - 1)]
        gaps.append((merged[-1][1], merged[0][0] + 1))
        gap_low, gap_high = max(gaps, key=lambda gap: gap[1] - gap[0])
        if gap_high <= gap_low:
            raise InputError(f'Non-PBC atoms fill the whole box along the box vector {axis + 1} in the first frame.\n'
                '  The automatic translation does not work in these cases.')
        # Move the middle of the gap to the boundary, choosing the shortest way
        middle = ((gap_low + gap_high) / 2) % 1
        translation[axis] = -middle if middle < 0.5 else 1 - middle
    translation = translation @ box_vectors
    print(f'  Automatic translation: {translation[0]:.3f} {translation[1]:.3f} {translation[2]:.3f} nm')
    return [round(float(value), 3) for value in translation]
