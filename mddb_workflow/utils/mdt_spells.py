import os
from os.path import exists, getsize
from subprocess import run, PIPE

import mdtraj as mdt
import numpy as np

from mddb_workflow.utils.auxiliar import ToolError, InputError, reprint
from mddb_workflow.utils.file import File
from mddb_workflow.utils.gmx_spells import merge_xtc_files
from mddb_workflow.utils.constants import GREY_HEADER, COLOR_END
from mddb_workflow.utils.type_hints import *

CURSOR_UP_ONE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'
ERASE_PREVIOUS_LINE = CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE

# Set mdtraj supported formats
mdtraj_supported_structure_formats = {
    'pdb', 'pdb.gz' 'h5', 'lh5', 'prmtop', 'parm7', 'prm7', 'psf', 'mol2', 'hoomdxml', 'gro', 'arc', 'hdf5', 'gsd'
}
mdtraj_supported_trajectory_formats = {'dcd', 'xtc', 'trr', 'nc', 'h5', 'binpos', 'mdcrd', 'xyz', 'pdb'}

# Use mdtraj 'mdconvert' command-line script (there is no python version for this tool apparently)
# Multiple files may be selected with bash syntax
def merge_and_convert_trajectories (
    input_trajectory_filenames : list[str],
    output_trajectory_filename : str
    ):

    # Run MDtraj
    print(GREY_HEADER, end='\r')
    process = run([
        "mdconvert",
        "-o",
        output_trajectory_filename,
        *input_trajectory_filenames,
    ], stderr=PIPE)
    error_logs = process.stderr.decode()
    print(COLOR_END, end='\r')

    # If output has not been generated then warn the user
    if not exists(output_trajectory_filename) or getsize(output_trajectory_filename) == 0:
        print(error_logs)
        raise ToolError('Something went wrong with MDTraj')

# NEVER FORGET: mdconvert does not handle mdcrd format
mdconvert_supported_formats = {'dcd', 'xtc', 'trr', 'nc', 'h5', 'binpos'}
merge_and_convert_trajectories.format_sets = [
    {
        'inputs': {
            'input_trajectory_filenames': mdconvert_supported_formats
        },
        'outputs': {
            'output_trajectory_filename': mdconvert_supported_formats
        }
    },
]

# Merge and convert trajectories without the mdconvert command
# WARNING: This process is slow since we must iterate and merge each frame
# WARNING: This process is restricted to trr/xtc output files given the merger
def merge_and_convert_trajectories_alternative (
    input_structure_filename : str,
    input_trajectory_filenames : list[str],
    output_trajectory_filename : str
    ):

    # Assert we have input values
    if not input_structure_filename:
        raise InputError('Missing input structure filenames')
    if not input_trajectory_filenames:
        raise InputError('Missing input trajectory filenames')
    if not output_trajectory_filename:
        raise InputError('Missing output trajectory filename')

    # Load the topology, which is used further
    topology = mdt.load_topology(input_structure_filename)

    # If the output trajectory filename matches any input file then we must rename the input filename to not overwrite it
    if output_trajectory_filename in input_trajectory_filenames:
        backup_filename = 'backup.' + output_trajectory_filename
        os.rename(output_trajectory_filename, backup_filename)
        repeated_input_filename_index = input_trajectory_filenames.index(output_trajectory_filename)
        input_trajectory_filenames[repeated_input_filename_index] = backup_filename

    # If the output trajectory file already exists at this point then we must stop here
    # The raw trjcat implementation will not join things to the end of it
    if exists(output_trajectory_filename):
        raise InputError('The output file already exists and overwrite is not supported for this functionality')
    # Print an empty line for the first 'ERASE_PREVIOUS_LINE' to not delete a previous log
    print()
    # Iterate over the different input trajectory filenames
    frame_filename = '.frame.xtc'
    for input_trajectory_filename in input_trajectory_filenames:
        # Load the trajectory frame by frame
        trajectory = mdt.iterload(input_trajectory_filename, top=input_structure_filename, chunk=1)
        for f, frame in enumerate(trajectory):
            # Update the current frame log
            print(ERASE_PREVIOUS_LINE)
            print('Frame ' + str(f))
            frame.save(frame_filename)
            # Join current frame to the output trajectory
            merge_xtc_files(output_trajectory_filename, frame_filename)
    # Remove the residual file
    # WARNING: It may not exist if the trajectory has 1 frame
    if exists(frame_filename):
        os.remove(frame_filename)
merge_and_convert_trajectories_alternative.format_sets = [
    {
        'inputs': {
            'input_structure_filename': mdtraj_supported_structure_formats,
            'input_trajectory_filenames': mdtraj_supported_trajectory_formats
        },
        'outputs': {
            'output_trajectory_filename': {'xtc', 'trr'}
        }
    },
]

# Merge and convert trajectories without the mdconvert command
# WARNING: Note that this process is not memory efficient so beware the size of trajectories to be converted
# DEPRECATED: This function was meant to convert trajectories to mdcrd, which is not supported by mdconvert
# DEPRECATED: However the export to mdcrd in mdtraj does not allow to remove the periodic box, which may be a problem
# DEPRECTAED: Use VMD instead
def merge_and_convert_trajectories_unefficient (
    input_structure_filename : str,
    input_trajectory_filenames : list[str],
    output_trajectory_filename : str,
):
    print('WARNING: You are using a not memory efficient tool. If the trajectory is too big your system may not hold it.')
    # Load all trajectories together
    sample_trajectory = input_trajectory_filenames[0]
    if input_structure_filename:
        trajectory = mdt.load(sample_trajectory, top=input_structure_filename)
    else:
        trajectory = mdt.load(sample_trajectory)
    for extra_trajectory in input_trajectory_filenames[1:]:
        if input_structure_filename:
            extra = mdt.load(sample_trajectory, top=input_structure_filename)
        else:
            extra = mdt.load(sample_trajectory)
        trajectory = mdt.join([trajectory, extra], check_topology=False)

    # Write the whole trajectory
    trajectory.save(output_trajectory_filename)

merge_and_convert_trajectories_unefficient.format_sets = [
    {
        'inputs': {
            'input_structure_filename': mdtraj_supported_structure_formats,
            'input_trajectory_filenames': mdtraj_supported_trajectory_formats
        },
        'outputs': {
            'output_trajectory_filename': mdtraj_supported_trajectory_formats
        }
    },
]

# Get specific frames from a trajectory
# WARNING: This function is time unefficient
def get_trajectory_subset (
    input_structure_filename : str,
    input_trajectory_filename : str,
    output_trajectory_filename : str,
    start : int = 0,
    end : int = None,
    step : int = 1,
    frames : list[int] = [],
    skip : list[int] = [],
):

    # Load the trajectory frame by frame and get only the desired frames
    if input_structure_filename:
        trajectory = mdt.iterload(input_trajectory_filename, top=input_structure_filename, chunk=1)
    else:
        trajectory = mdt.iterload(input_trajectory_filename, chunk=1)
    # Set the reduced trajectory to be returned
    reduced_trajectory = None
    # Specific frames scenario
    # DANI: No se ha provado
    if frames and len(frames) > 0:
        # Get the first and the last frames
        sorted_frames = sorted(frames)
        first_frame = sorted_frames[0]
        last_frame = sorted_frames[-1]
        # Get the list of frames as a set to speed up the searches
        selected_frames = set(frames)
        # Get the first chunk/frame
        for frame, chunk in enumerate(trajectory):
            if frame == first_frame:
                reduced_trajectory = chunk
                break
        # If we have nothing at this point then it means our start is out of the frames range
        if not reduced_trajectory:
            raise InputError('None of the selected frames are in range of the trajectory frames number')
        # Get further chunks/frames
        for frame, chunk in enumerate(trajectory, first_frame + 1):
            if frame in selected_frames:
                reduced_trajectory = mdt.join([reduced_trajectory, chunk], check_topology=False)
            # If this is the last frame then stop here to avoid reading the rest of the trajectory
            if frame >= last_frame:
                break
    # Start / End / Step + Skip scenario
    else:
        # Get the first chunk
        frame_count = 0 # This count works only in case the start frame is out of range, for the logs
        for i, chunk in enumerate(trajectory):
            frame_count = i
            if i == start:
                reduced_trajectory = chunk
                break
        # If we have nothing at this point then it means our start is out of the frames range
        if not reduced_trajectory:
            frame_count += 1
            raise InputError(f'The trajectory has {frame_count} frames so we can not start at frame {start}')
        # Get further chunks
        for i, chunk in enumerate(trajectory, 1): # Start the count at 1
            frame = start + i
            if frame in skip:
                continue
            if frame == end:
                break
            if i % step == 0:
                reduced_trajectory = mdt.join([reduced_trajectory, chunk], check_topology=False)

    # WARNING: This print here is not just a log. DO NOT REMOVE IT
    # WARNING: It fixes an error (ValueError: unitcell angle < 0) which happens sometimes
    # WARNING: It fixes an error (ValueError: Only rectilinear boxes can be saved to mdcrd files) which happens sometimes
    # DANI: Esto es azaroso: a veces funciona y a veces no
    #print(reduced_trajectory)

    # This is necessary sometimes to avoid the following error:
    #     ValueError: Only rectilinear boxes can be saved to mdcrd files
    # DANI: Esto suele funcionar, aunque a veces mata el proceso
    # e.g. ValueError: unitcell_lengths must be shape (5, 3). You supplied  (1, 3)
    # Este ejemplo fue observado con las trayectorias raw de 'NAFlex_1d11'
    # mdtb subset -is structure.stripped.pdb -it structure.stripped.mdcrd -ot test.mdcrd -start 2 -end 17 -step 3
    # DANI: Va bien que esto esté desactivado, así cuando falle podré compartir la trayectoria que da problemas a los de mdtraj:
    # https://github.com/mdtraj/mdtraj/issues/1730
    #reduced_trajectory.unitcell_lengths = [0,0,0]
    #reduced_trajectory.unitcell_angles = [90,90,90]

    # Write reduced trajectory to output file
    reduced_trajectory.save(output_trajectory_filename)
get_trajectory_subset.format_sets = [
    {
        'inputs': {
            'input_structure_filename': mdtraj_supported_structure_formats,
            'input_trajectory_filename': mdtraj_supported_trajectory_formats
        },
        'outputs': {
            'output_trajectory_filename': mdtraj_supported_trajectory_formats
        }
    },
    {
        'inputs': {
            'input_structure_filename': None,
            'input_trajectory_filename': mdtraj_supported_structure_formats
        },
        'outputs': {
            'output_trajectory_filename': mdtraj_supported_structure_formats
        }
    }
]

# Get first frame from a trajectory
def get_first_frame (
    input_structure_filename : str,
    input_trajectory_filename : str,
    output_frame_filename : str,
):
    get_trajectory_subset(input_structure_filename, input_trajectory_filename, output_frame_filename, 0)
# Set function supported formats
single_frame_getter_format_sets = [
    {
        'inputs': {
            'input_structure_filename': None,
            'input_trajectory_filename': mdtraj_supported_trajectory_formats
        },
        'outputs': {
            'output_frame_filename': mdtraj_supported_trajectory_formats
        }
    },
    {
        'inputs': {
            'input_structure_filename': mdtraj_supported_structure_formats,
            'input_trajectory_filename': mdtraj_supported_trajectory_formats
        },
        'outputs': {
            'output_frame_filename': { *mdtraj_supported_structure_formats, *mdtraj_supported_trajectory_formats }
        }
    },
]
get_first_frame.format_sets = single_frame_getter_format_sets

# Split a trajectory which is actually a merge of independent trajectories back to the original pieces
# Run an RMSD analysis to guess where the pieces are
# The cutoff sets the minimum RMSD change to consider it is a different trajectory
def split_merged_trajectories (
    input_structure_filename : str,
    input_trajectory_filename : str,
    sudden_jump_cutoff : float = 0.2,
    output_trajectory_prefix : str = 'split'):
    # Get the input trajectory format
    input_trajectory_file = File(input_trajectory_filename)
    input_trajectory_format = input_trajectory_file.format
    # The cutoff must be a negative number since independent trajectories RMSD sudden jumps will be negative
    cutoff = -abs(sudden_jump_cutoff)
    # Load the trajectory
    trajectory = mdt.load(input_trajectory_filename, top=input_structure_filename)
    # Run a RMSD analysis
    rmsd_data = mdt.rmsd(trajectory, trajectory, 0)
    # Find sudden jumps in RMSD values
    sudden_jumps = []
    previous_value = 0
    for i, value in enumerate(rmsd_data):
        diff = value - previous_value
        # In case this is a new trajectory the jump will be negative
        if diff < cutoff:
            print('New trajectory at frame ' + str(i))
            sudden_jumps.append(i)
        previous_value = value
    # In there was no jumps then stop here
    if len(sudden_jumps) == 0:
        print('Apparently there is a single trajectory')
        return
    # Generate a trajectory subset for each cut
    cut_points = [ 0, *sudden_jumps, len(rmsd_data) ]
    for i in range(len(cut_points) -1):
        start = cut_points[i]
        end = cut_points[i+1]
        trajectory_split = trajectory[start:end]
        split_filename = output_trajectory_prefix + '_' + str(i+1) + '.' + input_trajectory_format
        print('Writting from frame ' + str(start) + ' to frame ' + str(end) + ' to "' + split_filename + '"')
        trajectory_split.save(split_filename)

# Sort atoms in a trajectory file by sorting its coordinates
# A new atom indices list is expected in order to do the sorting
def sort_trajectory_atoms (
    input_structure_file : 'File',
    input_trajectory_file : 'File',
    output_trajectory_file : 'File',
    new_atom_indices : list[int]
):

    # Load the topology, which is used further
    topology = mdt.load_topology(input_structure_file.path)

    # Check the new atom indices to match the number of atoms in the topology
    if len(new_atom_indices) != topology.n_atoms:
        raise ValueError('The number of atom indices (' + str(len(new_atom_indices)) + ') does not match the number of atoms in topology(' + str(topology.n_atoms) + ')' )

    # If the input and the output trajectory filenames are equal then we must rename the input filename to not overwrite it
    if input_trajectory_file.path == output_trajectory_file.path:
        backup_file = input_trajectory_file.get_prefixed_file('backup.')
        os.rename(input_trajectory_file.path, backup_file.path)
        input_trajectory_file = backup_file

    # If the output trajectory file already exists at this point then we must stop here
    # The raw trjcat implementation will not join things to the end of it
    if output_trajectory_file.exists:
        raise InputError('The output file already exists and overwrite is not supported for this functionality')

    # Read and write the xtc files directly instead of using MDtraj trajectories
    # MDtraj trajectories rebuild box vectors with the first vector along the x axis
    # This would lose the rotation of boxes in fitted trajectories (see image_and_fit)
    if input_trajectory_file.format != 'xtc' or output_trajectory_file.format != 'xtc':
        raise ValueError('Sorting trajectory atoms is only supported for xtc trajectories')
    new_atom_indices = np.array(new_atom_indices)

    # Print an empty line for the first 'ERASE_PREVIOUS_LINE' to not delete a previous log
    print()

    with mdt.formats.XTCTrajectoryFile(input_trajectory_file.path) as reader, \
        mdt.formats.XTCTrajectoryFile(output_trajectory_file.path, 'w') as writer:
        frame_count = 0
        while True:
            xyz, time, step, box = reader.read(n_frames=1000)
            if len(xyz) == 0: break
            frame_count += len(xyz)
            # Update the current frame log
            print(ERASE_PREVIOUS_LINE)
            print('Frame ' + str(frame_count))
            # Sort coordinates according to new atom indices
            writer.write(xyz[:, new_atom_indices], time=time, step=step, box=box)

# Get a frame in MDtraj
# Most specifically, get a trajectry with a single frame
def get_frame (structure_file : 'File', trajectory_file : 'File', frame : int) -> 'mdt.Trajectory':
    # Read the trajectory until we get the desired frame
    trajectory = mdt.iterload(trajectory_file.path, top=structure_file.path, chunk=1)
    # Iterate frames until we get the one we want
    for f, mdt_frame in enumerate(trajectory):
        if f != frame: continue
        return mdt_frame

# Get the simulation box size from any trajectory format
# WARNING: There is a function for XTC format only which is way faster
def get_simulation_box(
    input_trajectory_filename : str,
    input_structure_filename : str,
    verbose: bool = True,
) -> tuple[Optional[tuple[float, float, float]], bool]:
    """Read the simulation box across all trajectory frames.
    Return also if it is constant or dynamic.
    If the simulation box is dynamic then the maximum size for every dimension is returned.
    """
    print()
    rows = []
    trajectory = mdt.iterload(input_trajectory_filename, top=input_structure_filename, chunk=1)
    # Iterate trajectory frames while mining every box size
    for frame_idx, frame in enumerate(trajectory, 1):
        reprint(f' Reading simulation box size: frame {frame_idx}')
        if frame.unitcell_lengths is not None:
            rows.append(frame.unitcell_lengths[0] * 10.0)

    # Handle when there is no box at all
    if not rows:
        return None, None

    # Check if the box is constant or dynamic
    boxes_arr = np.array(rows)
    if boxes_arr.shape[0] > 1:
        variation = np.std(boxes_arr, axis=0) / np.mean(boxes_arr, axis=0)
        is_constant = bool(np.all(variation < 0.001))
    else:
        is_constant = True

    # Set the final box and make it a tuple
    if is_constant:
        b = boxes_arr[0]
        box = (float(b[0]), float(b[1]), float(b[2]))
    else:
        m = boxes_arr.max(axis=0)
        box = (float(m[0]), float(m[1]), float(m[2]))

    if box is None:
        return None, None

    # Make a summary of the final box
    if verbose:
        if is_constant:
            print(f' Simulation box is constant ({box[0]:.2f} × {box[1]:.2f} × {box[2]:.2f} Å)')
        else:
            print(f' Simulation box varies. Max per dimension: '
                f'x={box[0]:.2f} Å, y={box[1]:.2f} Å, z={box[2]:.2f} Å')

    return box, is_constant


def iterate_xtc_frames(xtc_filename : str):
    """Iterate the frames in a xtc file as (positions, box vectors), both in nm.
    Boxes are read as they are in the file, since MDtraj rebuilds box vectors with the first vector along the x axis.
    This would be wrong for fitted trajectories, where boxes are rotated (see image_and_fit).
    Trajectories with no box give a box full of zeros, which means there are no periodic boundaries.
    """
    with mdt.formats.XTCTrajectoryFile(xtc_filename) as reader:
        while True:
            xyz, _, _, box = reader.read(n_frames=1)
            if len(xyz) == 0: return
            yield xyz[0], box[0] if box is not None else np.zeros((3, 3), dtype=np.float32)


# Maximum value (nm) of box vector components which must be zero in a not rotated box
# Note that the xtc precision is 0.001 nm
ROTATED_BOX_TOLERANCE = 0.001


def is_box_rotated(box_vectors: np.ndarray) -> bool:
    """Check if a box is rotated, i.e. the first vector is not along the x axis or the second vector is not in the xy plane.
    Most tools (e.g. Gromacs, MDtraj, VMD) expect boxes not rotated.
    Note that boxes are rotated when the trajectory is fitted (see rotate_boxes).
    """
    return bool(np.abs(box_vectors[[0, 0, 1], [1, 2, 2]]).max() > ROTATED_BOX_TOLERANCE)


def unrotate_box(positions: np.ndarray, box_vectors: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotate the coordinates and the box of a frame so the box is not rotated anymore (see is_box_rotated).
    This way tools which expect boxes not rotated (e.g. MDtraj) can be used with fitted trajectories.
    """
    first_vector, second_vector = box_vectors[0], box_vectors[1]
    x_axis = first_vector / np.linalg.norm(first_vector)
    y_axis = second_vector - (second_vector @ x_axis) * x_axis
    y_axis /= np.linalg.norm(y_axis)
    basis = np.stack([x_axis, y_axis, np.cross(x_axis, y_axis)], axis=1)
    return positions @ basis, box_vectors @ basis


# Thresholds to guess if a trajectory is fitted
# These thresholds were set by looking at 40 fitted trajectories in the database, a few unfitted trajectories
# and fitted trajectories which were rotated and translated on purpose
# Maximum translation (nm) of the center of the selected atoms
FITTED_MAX_TRANSLATION = 0.5
# Maximum fraction of the RMSD between centered coordinates which is removed by rotating them
# Note that a fitted trajectory still has some gain since its conformation changes
FITTED_MAX_ROTATION_GAIN = 0.1
# Fraction of residues which move the least, without any alignment, to be considered the region the trajectory was fitted to
# Some trajectories are fitted to a part of the system only, so the rest of the system may rotate a bit
FITTED_REGION_FRACTION = 0.1
# Minimum width (nm) of the region, since the rotation of a small region is not reliable
FITTED_REGION_MIN_WIDTH = 1
# Maximum rotation (degrees) of the region
FITTED_REGION_MAX_ROTATION = 2


def read_xtc_sample_frames(trajectory_file: 'File') -> list[tuple[np.ndarray, Optional[np.ndarray]]]:
    """Read the coordinates and the box of the first, middle and last frames of a xtc trajectory.
    The box is None if the trajectory has no box.
    """
    if trajectory_file.format != 'xtc':
        raise ValueError(f'Reading sample frames is only supported for xtc trajectories: {trajectory_file.path}')
    frames = []
    with mdt.formats.XTCTrajectoryFile(trajectory_file.path) as reader:
        frames_count = len(reader)
        for frame in sorted(set([0, frames_count // 2, frames_count - 1])):
            reader.seek(frame)
            xyz, _, _, box = reader.read(n_frames=1)
            frames.append((xyz[0].astype(np.float64), None if box is None else box[0].astype(np.float64)))
    return frames


def guess_translation_fitted(frames: list[tuple[np.ndarray, Optional[np.ndarray]]], atom_indices: list[int]) -> bool:
    """Guess if a trajectory is fitted (or centered) in translation, i.e. the center of the selected atoms barely moves.
    The first frame is compared with the rest of frames.
    """
    if len(frames) < 3: return False
    reference = frames[0][0][atom_indices]
    for coordinates, _ in frames[1:]:
        translation = np.linalg.norm(coordinates[atom_indices].mean(axis=0) - reference.mean(axis=0))
        if translation > FITTED_MAX_TRANSLATION: return False
    return True


def guess_rotation_fitted(frames: list[tuple[np.ndarray, Optional[np.ndarray]]], residue_atom_indices: list[list[int]]) -> bool:
    """Guess if a trajectory is fitted in rotation.
    The first frame is compared with the rest of frames, ignoring the translation.

    A trajectory is fitted in rotation if rotating the selected atoms barely improves their RMSD.
    Note that the angle of the rotation is not reliable when the conformation changes a lot (e.g. flexible molecules).

    A trajectory is also fitted in rotation if the region which moves the least is wide and it barely rotates.
    This is for trajectories which are fitted to a part of the system only (e.g. the base of a spike).
    """
    if len(frames) < 3: return False
    atom_indices = [ index for residue_indices in residue_atom_indices for index in residue_indices ]
    reference = frames[0][0]
    if all(get_rotation_gain(reference[atom_indices], coordinates[atom_indices]) < FITTED_MAX_ROTATION_GAIN
        for coordinates, _ in frames[1:]): return True
    # Find the region which moves the least
    if len(residue_atom_indices) < 3: return False
    def get_residue_centers(coordinates: np.ndarray) -> np.ndarray:
        return np.array([ coordinates[residue_indices].mean(axis=0) for residue_indices in residue_atom_indices ])
    reference_centers = get_residue_centers(reference)
    displacements = np.max([ np.linalg.norm(get_residue_centers(coordinates) - reference_centers, axis=1)
        for coordinates, _ in frames[1:] ], axis=0)
    region_residues_count = max(3, int(len(residue_atom_indices) * FITTED_REGION_FRACTION))
    region_residues = np.argsort(displacements)[:region_residues_count]
    region = np.concatenate([ residue_atom_indices[residue] for residue in region_residues ])
    # The width of the region is its spread along its second principal axis, so a line has no width
    centered_region = reference[region] - reference[region].mean(axis=0)
    spreads = np.sqrt(np.clip(np.linalg.eigvalsh(centered_region.T @ centered_region / len(region)), 0, None))
    if spreads[1] < FITTED_REGION_MIN_WIDTH: return False
    for coordinates, _ in frames[1:]:
        rotation, _ = find_rotation(reference[region], coordinates[region])
        angle = np.degrees(np.arccos(np.clip((np.trace(rotation) - 1) / 2, -1, 1)))
        if angle > FITTED_REGION_MAX_ROTATION: return False
    return True


def get_rotation_gain(reference: np.ndarray, coordinates: np.ndarray) -> float:
    """Get the fraction of the RMSD between centered coordinates which is removed by rotating them."""
    centered_rmsd = np.sqrt((((coordinates - coordinates.mean(axis=0)) - (reference - reference.mean(axis=0))) ** 2).sum(axis=1).mean())
    # Identical coordinates are not rotated
    if centered_rmsd == 0: return 0
    _, rotated_rmsd = find_rotation(coordinates, reference)
    return 1 - rotated_rmsd / centered_rmsd


# Maximum RMSD (nm) between fitted coordinates and original coordinates once rotated
# Since the fit is a rigid transformation this should be only the xtc precision noise
MAX_ROTATION_RESIDUAL = 0.02
# Number of frames to be read and written at once when rotating boxes
ROTATION_CHUNK_SIZE = 1000


def find_rotation(original_coordinates: np.ndarray, fitted_coordinates: np.ndarray) -> tuple[np.ndarray, float]:
    """Find the rotation matrix R which turns original coordinates into fitted coordinates (Kabsch).
    Coordinates are row vectors, so centered fitted coordinates ~ centered original coordinates @ R.
    Translation is ignored since both sets of coordinates are centered.
    Return the rotation matrix and the RMSD between fitted and rotated original coordinates.
    """
    original = original_coordinates - original_coordinates.mean(axis=0)
    fitted = fitted_coordinates - fitted_coordinates.mean(axis=0)
    u, _, vt = np.linalg.svd(original.T @ fitted)
    # Make sure we get a proper rotation and not a reflection
    sign = np.sign(np.linalg.det(u @ vt))
    rotation = u @ np.diag([1, 1, sign]) @ vt
    residual = np.sqrt(((original @ rotation - fitted) ** 2).sum(axis=1).mean())
    return rotation, residual


def rotate_boxes(
    original_trajectory_file: 'File',
    fitted_trajectory_file: 'File',
    output_trajectory_file: 'File',
):
    """Write the fitted trajectory with every box rotated as its frame coordinates were rotated by the fit.
    The rotation in every frame is found by comparing the original and the fitted coordinates.
    Note that most tools (e.g. Gromacs, MDtraj) expect the first box vector along the x axis.
    Thus they do not support rotated boxes and periodic images may be wrong when using them.
    """
    for trajectory_file in [original_trajectory_file, fitted_trajectory_file, output_trajectory_file]:
        if trajectory_file.format != 'xtc':
            raise ValueError(f'Rotating boxes is only supported for xtc trajectories: {trajectory_file.path}')
    with mdt.formats.XTCTrajectoryFile(original_trajectory_file.path) as original_reader, \
        mdt.formats.XTCTrajectoryFile(fitted_trajectory_file.path) as fitted_reader, \
        mdt.formats.XTCTrajectoryFile(output_trajectory_file.path, 'w') as writer:
        while True:
            original_xyz, _, _, _ = original_reader.read(n_frames=ROTATION_CHUNK_SIZE)
            fitted_xyz, time, step, boxes = fitted_reader.read(n_frames=ROTATION_CHUNK_SIZE)
            if len(original_xyz) != len(fitted_xyz):
                raise ValueError('Original and fitted trajectories have a different number of frames')
            if len(fitted_xyz) == 0: break
            rotated_boxes = np.empty_like(boxes)
            for i, (original_coordinates, fitted_coordinates) in enumerate(zip(original_xyz, fitted_xyz)):
                rotation, residual = find_rotation(original_coordinates.astype(np.float64), fitted_coordinates.astype(np.float64))
                if residual > MAX_ROTATION_RESIDUAL:
                    raise ValueError(f'Fitted coordinates are not a rigid rotation of the original coordinates (RMSD {residual:.4f} nm)')
                # Box vectors are rows, just like coordinates
                rotated_boxes[i] = boxes[i].astype(np.float64) @ rotation
            writer.write(fitted_xyz, time=time, step=step, box=rotated_boxes)
