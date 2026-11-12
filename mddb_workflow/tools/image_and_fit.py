# This script is used to fit a trajectory (i.e. eliminate translation and rotation)
# This process is carried by Gromacs

from collections import deque

import mdtraj as mdt
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components, breadth_first_order
from scipy.spatial import cKDTree
from tqdm import tqdm

from mddb_workflow.utils.structures import Structure
from mddb_workflow.utils.auxiliar import MISSING_TOPOLOGY, MISSING_BONDS, InputError, warn
from mddb_workflow.utils.constants import AUTOMATIC_FLAG
from mddb_workflow.utils.gmx_spells import run_gromacs
from mddb_workflow.utils.mdt_spells import rotate_boxes, iterate_xtc_frames
from mddb_workflow.utils.mdt_spells import read_xtc_sample_frames, is_box_rotated, has_box
from mddb_workflow.utils.mdt_spells import guess_translation_fitted, guess_rotation_fitted
from mddb_workflow.analyses.mindist import find_frame_cross_contacts
from mddb_workflow.tools.get_bonds import mine_topology_bonds
from mddb_workflow.utils.type_hints import *

# Set some constant values to handle NDX atom index selections
SYSTEM_SELECTION_NAME = 'system'
CENTER_SELECTION_NAME = 'center'
INDEX_FILEPATH = '.index.ndx'

# Distance (nm) below which non-PBC atoms are considered part of the same object when finding the automatic translation
AUTO_TRANSLATION_CUTOFF = 0.5

# Number of frames read at once in the no-jump step
NO_JUMP_CHUNK_SIZE = 100

# Number of frames along the trajectory where split molecules are searched when the image is automatic
SPLIT_CHECK_FRAMES = 10


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
    find_bonds: Callable[['File', int], list[list[int]]],
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
            Set if it must be fitted.
            Note that the automatic flag must be solved before calling this function.
        translation:
            Set how it must be translated during imaging, in nm.
            Use 'auto' to find a translation which keeps non-PBC atoms away from the box boundaries.
            If all atoms are PBC then lipids are used instead, so membranes are not split along their normal.
        structure:
            The structure of the system. Note that this is an early provisional structure.
        pbc_selection:
            Set a selection to exclude PBC residues from both the centering and the fitting focuses.
            Selection of atoms which are to stay in PBC conditions after imaging.
            Note that this is an early provisional atom selection.
        find_bonds:
            Function to find the bonds of the system, given a trajectory and its number of frames.
            Bonds may be guessed from atom distances, so they are found once the trajectory has no jumps.

    """
    print(f' Image: {image} | Fit: {fit}')

    if not image and not fit:
        return

    # Without a simulation box there are no periodic images, so imaging makes no sense
    # Fitting is still possible, but there is no box to rotate
    trajectory_has_box = has_box(input_trajectory_file.path)
    if image and not trajectory_has_box:
        raise InputError('The trajectory has no simulation box, so it can not be imaged.\n'
            ' Please do not ask for the imaging (i.e. remove the "image" input or the "-img" argument).')

    # Using a TPR topology may provide some advantage in the imaging process
    # Bonds are taken from the topology instead of being guessed from atom distances
    is_tpr_available = input_topology_file != MISSING_TOPOLOGY and input_topology_file.format == 'tpr'
    if not is_tpr_available:
        warn('Since we are missing a TPR file the automatic imaging protocol is more prone to fail.')

    # Gromacs can not read mmCIF so we work with PDB files, where chains do not matter
    # At the end, the output structure is written with the input structure data and the new coordinates
    input_structure = Structure.from_file(input_structure_file.path)
    gromacs_input_structure_file = output_structure_file.get_prefixed_file('input.').get_reformated_file('pdb')
    input_structure.generate_pdb_file(gromacs_input_structure_file.path, show_warnings=False, accept_limitations=True)
    gromacs_structure_file = output_structure_file.get_reformated_file('pdb')

    # Set a custom index file (.ndx) to select reference atoms for centering and fitting steps
    # We will use as center all atoms which are not under PBC
    system_selection = structure.select('all', syntax='vmd')
    center_selection = structure.invert_selection(pbc_selection)
    if not center_selection:
        warn(f'The default selection to center and fit (non PBC) is empty. These steps will be skipped.')
    # Convert both selections to a single ndx file which gromacs can read
    with open(INDEX_FILEPATH, 'w') as file:
        system_selection_ndx = system_selection.to_ndx(SYSTEM_SELECTION_NAME)
        file.write(system_selection_ndx)
        if center_selection:
            center_selection_ndx = center_selection.to_ndx(CENTER_SELECTION_NAME)
            file.write(center_selection_ndx)

    # Use the first frame of the trajectory as structure
    # The input structure may be a different frame, but gromacs uses it as reference for the no-jump and the fit
    reset_structure(gromacs_input_structure_file.path, input_trajectory_file.path, gromacs_structure_file.path)

    # Set the latest structure and trajectory files
    # They are the ones to be used as inputs for the following steps
    latest_structure = gromacs_structure_file
    latest_trajectory = input_trajectory_file

    # Imaging --------------------------------------------------------------------------------------

    if image:
        print(' Running imaging steps')

        # Find the translation automatically if requested
        # If all atoms are PBC then use lipids, so membranes are not split by the box boundaries (e.g. lipid-only systems)
        # Note that lipids are still PBC atoms in the rest of the imaging, so they may be periodic (e.g. membranes)
        # Periodic non-PBC atoms (e.g. a nanoparticle touching its own image) still fail, since they can not be kept whole
        if translation == AUTOMATIC_FLAG:
            if center_selection:
                translation = find_automatic_translation(latest_trajectory, center_selection)
            else:
                lipid_selection = structure.select_lipids()
                if lipid_selection: print(' There are only PBC atoms, so lipids are used to find the automatic translation')
                translation = find_automatic_translation(latest_trajectory, lipid_selection, ignore_periodic=True)

        # !!! We assume that at this point (after a possible translation) the first frame is "correct"
        # This means that the different molecules are placed as desired
        # e.g. two protein monomers are together, and not interacting across boundaries
        # To keep things like this we will run a no-jump step
        # The translation is applied in the no-jump step as well
        # In the first frame residues are made whole and, if translated, non-PBC residues are put in the box
        # This way non-PBC objects are together, since the automatic translation places the box boundaries between them
        # PBC residues are not put in the box here since molecules of many residues would be split
        # LORE: We were translating with '-trans -pbc res' and running an initial '-pbc res' step with the TPR
        # LORE: However '-pbc res' split PBC molecules of many residues, and the no-jump kept them split
        # DANI: The nojump step may cause artifacts, specially in already pre-imaged simulations
        print(' Running no-jump')
        # LORE: We were using the Gromacs no-jump, but it decides every atom jump independently
        # LORE: Atoms of the same molecule could end in different periodic images when they moved around half box
        # run_gromacs(f'trjconv -s {latest_structure.path} \
        #     -f {latest_trajectory.path} -o {output_trajectory_file.path} \
        #     -pbc nojump', user_input=SYSTEM_SELECTION_NAME, show_error_logs=True)
        atom_groups = find_no_jump_groups(input_topology_file, structure)
        must_translate = translation != [0, 0, 0]
        if must_translate: print(f'  Translation: {translation[0]} {translation[1]} {translation[2]} nm')
        atom_residues = np.array([ atom.residue_index for atom in structure.atoms ])
        prepare_first_frame = lambda positions, box_vectors: place_first_frame(
            positions, box_vectors, atom_residues, center_selection.atom_indices, must_translate)
        # Input and output trajectories may be the same file so we write an intermediate file
        no_jump_trajectory_file = output_trajectory_file.get_prefixed_file('nojump.')
        no_jump(latest_trajectory, no_jump_trajectory_file, atom_groups,
            translation=np.array(translation, dtype=np.float64), prepare_first_frame=prepare_first_frame)
        no_jump_trajectory_file.rename_to(output_trajectory_file)
        latest_trajectory = output_trajectory_file

        # Now that there are no jumps the bonds may be guessed from atom distances
        # Then every bonded fragment is made whole, the non-PBC atoms are centered and PBC fragments are put in the box
        # LORE: We were using '-pbc res -center', but it split PBC molecules of many residues (e.g. RNA, PEG lipids)
        # LORE: Gromacs '-pbc mol' is not used either since its molecules may be much bigger (e.g. a whole membrane)
        print(' Putting PBC fragments in the box')
        with mdt.formats.XTCTrajectoryFile(latest_trajectory.path) as reader:
            frames_count = len(reader)
        atom_bonds = find_bonds(latest_trajectory, frames_count)
        boxed_trajectory_file = output_trajectory_file.get_prefixed_file('boxed.')
        put_fragments_in_box(latest_trajectory, boxed_trajectory_file, atom_bonds, structure, center_selection)
        boxed_trajectory_file.rename_to(output_trajectory_file)
        latest_trajectory = output_trajectory_file

        # Select the first frame of the recently imaged trayectory as the new structure
        reset_structure(latest_structure.path, latest_trajectory.path, gromacs_structure_file.path)
        latest_structure = gromacs_structure_file

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
        if trajectory_has_box:
            print(' Rotating boxes')
            rotate_boxes(latest_trajectory, fitted_trajectory_file, rotated_trajectory_file)
            fitted_trajectory_file.remove()
            rotated_trajectory_file.rename_to(output_trajectory_file)
        else:
            print(' There is no box to rotate')
            fitted_trajectory_file.rename_to(output_trajectory_file)

        # If there is no output structure at this time (i.e. there was no imaging) then create it now
        # If there was imaging, we reset again to it its not far from the origin
        # Note that the input structure is not necessarily the output structure in this scenario
        # The fit may have been done using the topology and there may be an offset, so better dump it
        reset_structure(
            gromacs_input_structure_file.path,
            output_trajectory_file.path,
            gromacs_structure_file.path)

    # Write the output structure with the new coordinates
    output_coordinates = [ atom.coords for atom in Structure.from_pdb_file(gromacs_structure_file.path).atoms ]
    input_structure.set_new_coordinates(output_coordinates)
    input_structure.generate_file(output_structure_file.path)
    gromacs_input_structure_file.remove()
    gromacs_structure_file.remove()

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
    translation: np.ndarray = np.zeros(3),
    prepare_first_frame: Optional[Callable[[np.ndarray, np.ndarray], np.ndarray]] = None,
):
    """Remove jumps across the box boundaries along the trajectory, keeping the first frame as it is.
    All frames may be translated, and the first frame may be changed before (e.g. to make residues whole).
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
    with mdt.formats.XTCTrajectoryFile(input_trajectory_file.path) as reader:
        frames_count = len(reader)
    with mdt.formats.XTCTrajectoryFile(input_trajectory_file.path) as reader, \
        mdt.formats.XTCTrajectoryFile(output_trajectory_file.path, 'w') as writer, \
        tqdm(total=frames_count, desc='  Frame', unit='frame') as pbar:
        while True:
            xyz, time, step, boxes = reader.read(n_frames=NO_JUMP_CHUNK_SIZE)
            if len(xyz) == 0: break
            output_xyz = np.empty_like(xyz)
            for i, (positions, box_vectors) in enumerate(zip(xyz, boxes)):
                pbar.update(1)
                positions = positions.astype(np.float64) + translation
                if previous_positions is None:
                    if prepare_first_frame: positions = prepare_first_frame(positions, box_vectors.astype(np.float64))
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


def place_first_frame(
    positions: np.ndarray,
    box_vectors: np.ndarray,
    atom_residues: np.ndarray,
    selected_atom_indices: list[int],
    put_in_box: bool,
) -> np.ndarray:
    """Make every residue whole, by placing every atom in the image closest to the first atom of its residue.
    If requested, also put the selected atoms in the box, moving every residue as a whole.
    """
    inverse_box = np.linalg.inv(box_vectors)
    residues, residue_first_atoms = np.unique(atom_residues, return_index=True)
    first_atoms = residue_first_atoms[np.searchsorted(residues, atom_residues)]
    fractional = (positions - positions[first_atoms]) @ inverse_box
    positions = positions - closest_images(fractional, box_vectors) @ box_vectors
    if not put_in_box or len(selected_atom_indices) == 0: return positions
    selected_residues = np.unique(atom_residues[selected_atom_indices], return_inverse=True)[1]
    fractional = positions[selected_atom_indices] @ inverse_box
    sizes = np.bincount(selected_residues)
    centers = np.stack([ np.bincount(selected_residues, weights=fractional[:, axis]) / sizes for axis in range(3) ], axis=1)
    positions[selected_atom_indices] -= np.floor(centers)[selected_residues] @ box_vectors
    return positions


def find_bond_tree(atom_bonds: list[list[int]], structure: 'Structure') -> tuple[np.ndarray, np.ndarray, list[np.ndarray]]:
    """Find the bonded fragments and a path along bonds from the first atom of every fragment to the rest of atoms.
    Atoms with missing bonds (e.g. coarse grain) are joined to the first atom of their residue instead.
    Return the fragment of every atom, the previous atom in the path of every atom and the atoms at every step of the paths.
    """
    atom_count = structure.atom_count
    if len(atom_bonds) != atom_count:
        raise ValueError(f'There are bonds for {len(atom_bonds)} atoms, but the structure has {structure.atom_count} atoms')
    pairs = []
    for atom_index, bonds in enumerate(atom_bonds):
        if bonds == MISSING_BONDS:
            pairs.append((structure.atoms[atom_index].residue.atom_indices[0], atom_index))
            continue
        pairs += [ (atom_index, bonded) for bonded in bonds if atom_index < bonded ]
    pairs = np.array(pairs, dtype=int).reshape(-1, 2)
    graph = coo_matrix((np.ones(len(pairs), dtype=bool), (pairs[:, 0], pairs[:, 1])), shape=(atom_count, atom_count))
    _, atom_fragments = connected_components(graph, directed=False)
    # An extra atom is bonded to the first atom of every fragment, so a single search covers all fragments
    roots = np.unique(atom_fragments, return_index=True)[1]
    extra_atom = atom_count
    sources = np.concatenate([pairs[:, 0], np.full(len(roots), extra_atom)])
    targets = np.concatenate([pairs[:, 1], roots])
    graph = coo_matrix((np.ones(len(sources), dtype=bool), (sources, targets)), shape=(atom_count + 1, atom_count + 1))
    order, parents = breadth_first_order(graph, extra_atom, directed=False, return_predecessors=True)
    depths = np.zeros(atom_count + 1, dtype=int)
    for atom_index in order[1:]:
        depths[atom_index] = depths[parents[atom_index]] + 1
    # Roots are at depth 1, so their atoms are not moved
    steps = [ order[depths[order] == depth] for depth in range(2, depths.max() + 1) ]
    return atom_fragments, parents[:atom_count], steps


def put_fragments_in_box(
    input_trajectory_file: 'File',
    output_trajectory_file: 'File',
    atom_bonds: list[list[int]],
    structure: 'Structure',
    center_selection: 'Selection',
):
    """In every frame, make every bonded fragment whole, center the selection and put fragments out of the selection in the box.
    Fragments are put in the box as a whole, according to their center.
    Fragments with any atom in the selection are never put in the box, since they are expected to have no jumps already.
    """
    for trajectory_file in [input_trajectory_file, output_trajectory_file]:
        if trajectory_file.format != 'xtc':
            raise ValueError(f'Putting fragments in the box is only supported for xtc trajectories: {trajectory_file.path}')
    atom_fragments, parents, steps = find_bond_tree(atom_bonds, structure)
    center_atom_indices = center_selection.atom_indices
    centered_fragments = np.unique(atom_fragments[center_atom_indices])
    boxed_atom_indices = np.flatnonzero(~np.isin(atom_fragments, centered_fragments))
    boxed_atom_fragments = np.unique(atom_fragments[boxed_atom_indices], return_inverse=True)[1]
    boxed_fragment_sizes = np.bincount(boxed_atom_fragments)
    print(f'  {len(centered_fragments)} fragments are centered and {len(boxed_fragment_sizes)} fragments are put in the box')
    with mdt.formats.XTCTrajectoryFile(input_trajectory_file.path) as reader:
        frames_count = len(reader)
    with mdt.formats.XTCTrajectoryFile(input_trajectory_file.path) as reader, \
        mdt.formats.XTCTrajectoryFile(output_trajectory_file.path, 'w') as writer, \
        tqdm(total=frames_count, desc='  Frame', unit='frame') as pbar:
        while True:
            xyz, time, step, boxes = reader.read(n_frames=NO_JUMP_CHUNK_SIZE)
            if len(xyz) == 0: break
            output_xyz = np.empty_like(xyz)
            for i, (positions, box_vectors) in enumerate(zip(xyz, boxes)):
                pbar.update(1)
                if not box_vectors.any():
                    raise InputError(f'Frame {step[i]} has no box, so fragments can not be put in the box')
                positions, box_vectors = positions.astype(np.float64), box_vectors.astype(np.float64)
                inverse_box = np.linalg.inv(box_vectors)
                for atom_indices in steps:
                    fractional = (positions[atom_indices] - positions[parents[atom_indices]]) @ inverse_box
                    positions[atom_indices] -= closest_images(fractional, box_vectors) @ box_vectors
                if len(center_atom_indices) > 0:
                    positions += box_vectors.sum(axis=0) / 2 - positions[center_atom_indices].mean(axis=0)
                if len(boxed_atom_indices) > 0:
                    fractional = positions[boxed_atom_indices] @ inverse_box
                    centers = np.stack([ np.bincount(boxed_atom_fragments, weights=fractional[:, axis]) / boxed_fragment_sizes
                        for axis in range(3) ], axis=1)
                    positions[boxed_atom_indices] -= np.floor(centers)[boxed_atom_fragments] @ box_vectors
                output_xyz[i] = positions
            writer.write(output_xyz, time=time, step=step, box=boxes)


def check_fitted_trajectory(trajectory_file: 'File', structure: 'Structure', selection: 'Selection') -> tuple[bool, bool]:
    """Check if the trajectory seems to be fitted already, both in translation and rotation.
    Check also if it is fitted in rotation and periodic boundary conditions are thus broken.
    Return both values: fitted and broken PBC.
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
    if translation_fitted and rotation_fitted:
        print(' The trajectory seems to be fitted in translation and rotation')
    elif translation_fitted:
        print(' The trajectory seems to be fitted in translation, but NOT in rotation')
    elif rotation_fitted:
        print(' The trajectory seems to be fitted in rotation, but NOT in translation')
    fitted = translation_fitted and rotation_fitted
    if not rotation_fitted: return fitted, False
    boxes = [ box for _, box in frames ]
    if not all(box is not None and box.any() for box in boxes):
        warn('The trajectory seems to be fitted already, but it has no box.\n'
            ' Periodic boundary conditions are probably broken, so the imaging will probably fail.')
        return fitted, False
    if any(is_box_rotated(box) for box in boxes): return fitted, False
    warn('The trajectory seems to be fitted already, but its box is not rotated accordingly.\n'
        ' Periodic boundary conditions are probably broken, so the simulation box will be ignored and the imaging will probably fail.')
    return fitted, True


def find_split_molecules(
    trajectory_file: 'File',
    topology_file: Union['File', Exception],
    structure: 'Structure',
    selection: 'Selection',
) -> bool:
    """Check if molecules in the selection are split by the box boundaries in some frames.
    Two bonded atoms are split when they are closer through the box boundaries than directly.
    Bonds are taken from the topology. If they are missing then every atom is paired with the first atom of its residue.
    Only a few frames along the trajectory are checked, so short problems may be missed.
    """
    if trajectory_file.format != 'xtc':
        raise ValueError(f'Finding split molecules is only supported for xtc trajectories: {trajectory_file.path}')
    selected_atom_indices = set(selection.atom_indices)
    atom_bonds = mine_topology_bonds(topology_file)
    if atom_bonds and any(bonds == MISSING_BONDS for bonds in atom_bonds):
        atom_bonds = None
    if atom_bonds:
        if len(atom_bonds) != structure.atom_count:
            raise InputError(f'Topology bonds have {len(atom_bonds)} atoms, but the structure has {structure.atom_count} atoms')
        pairs = [ (atom_index, bonded) for atom_index, bonds in enumerate(atom_bonds) for bonded in bonds
            if atom_index < bonded and atom_index in selected_atom_indices and bonded in selected_atom_indices ]
    else:
        pairs = [ (residue.atom_indices[0], atom_index) for residue in structure.residues for atom_index in residue.atom_indices[1:]
            if residue.atom_indices[0] in selected_atom_indices and atom_index in selected_atom_indices ]
    if not pairs: return False
    pairs = np.array(pairs, dtype=int)
    with mdt.formats.XTCTrajectoryFile(trajectory_file.path) as reader:
        frames_count = len(reader)
        for frame in np.unique(np.linspace(0, frames_count - 1, SPLIT_CHECK_FRAMES).astype(int)):
            reader.seek(frame)
            xyz, _, _, boxes = reader.read(n_frames=1)
            if boxes is None or not boxes[0].any(): return False
            positions, box_vectors = xyz[0].astype(np.float64), boxes[0].astype(np.float64)
            fractional = (positions[pairs[:, 1]] - positions[pairs[:, 0]]) @ np.linalg.inv(box_vectors)
            split_count = int(closest_images(fractional, box_vectors).any(axis=1).sum())
            if split_count:
                print(f' {split_count} bonds are split by the box boundaries in frame {frame}')
                return True
    return False


def reset_structure(
    input_structure_filepath: str,
    input_trajectory_filepath: str,
    output_structure_filepath: str,
):
    """Get the first frame of a trajectory and use it to reset the structure."""
    print(' Reseting structure')
    run_gromacs(f'trjconv -s {input_structure_filepath} \
        -f {input_trajectory_filepath} -o {output_structure_filepath} \
        -dump 0', user_input=SYSTEM_SELECTION_NAME)


def find_automatic_translation(trajectory_file: 'File', selection: 'Selection', ignore_periodic: bool = False) -> list[float]:
    """Find a translation which places the box boundaries where no object of the selection is split.

    Only the first frame of the trajectory is used.
    Atoms in the selection which are close to each other, directly or across the box boundaries, are grouped in objects.
    Then every object is made whole and the boundary along every box vector is placed in the middle of the largest gap between objects.
    Objects in contact with their own periodic image or filling the whole box can not be kept whole, so it fails.
    If periodic objects are to be ignored (e.g. membranes of PBC lipids) then every box vector is handled on its own instead.
    Objects are split by any boundary along the box vectors they are periodic along, so they are ignored along them only.
    e.g. a flat membrane is ignored along the box vectors in its plane, but not along its normal.
    If there is no gap along a box vector then its boundary is not moved.
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
    # A fragment reached again with a different shift means the object is in contact with its own periodic image
    # Then the object is periodic (e.g. a membrane) along the box vectors where both shifts differ
    # Along the other box vectors shifts are always the same, so the object can still be made whole along them
    self_contact_error = InputError('Non-PBC atoms are in contact with their own periodic image in the first frame.\n'
        '  The automatic translation does not work in these cases.')
    fragment_shifts = [None] * fragment_count
    fragment_objects = np.empty(fragment_count, dtype=int)
    periodic_objects = []
    object_count = 0
    for seed in range(fragment_count):
        if fragment_shifts[seed] is not None: continue
        fragment_shifts[seed] = np.zeros(3, dtype=int)
        fragment_objects[seed] = object_count
        periodic_axes = np.zeros(3, dtype=bool)
        queue = deque([seed])
        while queue:
            fragment = queue.popleft()
            for other_fragment, shift in neighbours[fragment]:
                expected_shift = fragment_shifts[fragment] + shift
                if fragment_shifts[other_fragment] is None:
                    fragment_shifts[other_fragment] = expected_shift
                    fragment_objects[other_fragment] = object_count
                    queue.append(other_fragment)
                else:
                    periodic_axes |= fragment_shifts[other_fragment] != expected_shift
        periodic_objects.append(periodic_axes)
        object_count += 1
    periodic_objects = np.array(periodic_objects)
    if not ignore_periodic and periodic_objects.any(): raise self_contact_error
    # Get the fractional coordinates of every atom once its object is whole
    # Coordinates along the box vectors where the object is periodic are not used
    fractional = positions @ inverse_box - np.array(fragment_shifts)[atom_fragments]
    atom_objects = fragment_objects[atom_fragments]
    translation = np.zeros(3)
    for axis in range(3):
        # Objects which are periodic along this box vector are split by any boundary, so they are ignored
        periodic = periodic_objects[:, axis]
        if periodic.all():
            print(f'  Objects are periodic along the box vector {axis + 1} (e.g. membranes), so there is no translation along it')
            continue
        if periodic.any():
            print(f'  Some objects are periodic along the box vector {axis + 1} (e.g. membranes), so they are ignored along it')
        values = fractional[:, axis]
        lows = np.full(object_count, np.inf)
        highs = np.full(object_count, -np.inf)
        np.minimum.at(lows, atom_objects, values)
        np.maximum.at(highs, atom_objects, values)
        lows, highs = lows[~periodic], highs[~periodic]
        if not ignore_periodic and np.any(highs - lows >= 1): raise self_contact_error
        # Bring every object interval into [0, 1) and split the ones going beyond 1
        # Objects which are longer than the box along this box vector fill it all
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
            if not ignore_periodic:
                raise InputError(f'Non-PBC atoms fill the whole box along the box vector {axis + 1} in the first frame.\n'
                    '  The automatic translation does not work in these cases.')
            warn(f'Atoms to find the automatic translation fill the whole box along the box vector {axis + 1} in the first frame.\n'
                ' Some of them will be split by the box boundary anyway, so there is no translation along it.')
            continue
        # Move the middle of the gap to the boundary, choosing the shortest way
        middle = ((gap_low + gap_high) / 2) % 1
        translation[axis] = -middle if middle < 0.5 else 1 - middle
    translation = translation @ box_vectors
    print(f'  Automatic translation: {translation[0]:.3f} {translation[1]:.3f} {translation[2]:.3f} nm')
    return [round(float(value), 3) for value in translation]
