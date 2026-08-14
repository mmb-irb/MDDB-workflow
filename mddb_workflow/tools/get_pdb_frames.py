from mddb_workflow.utils.pyt_spells import get_pytraj_trajectory, get_reduced_pytraj_trajectory
from mddb_workflow.utils.auxiliar import reprint, get_auxiliar_filepath
from mddb_workflow.utils.type_hints import *
from tqdm import tqdm
import os

import pytraj as pt

# Build a generator which returns frames along the trajectory in pdb format
# The frames limit is the maximum number of frames to be iterated
# Note that the number of frames iterated may be less than the specified number
def get_pdb_frames (
    structure : 'Structure',
    trajectory_filepath : str,
    snapshots : int,
    frames_limit : Optional[int] = None,
    output_frames_prefix : str = 'frame',
    pbar_bool : bool = False,
) -> tuple[ Generator[str, None, None], int, int ]:
    # WARNING: Do not set a pytraj iterload trajectory and read its 'n_frames' to get the snapshots
    # WARNING: Trying to read the number o frames of a xtc trajectory will read the whole trajectory
    # WARNING: This may be a lot of time for a huge trajectory. Use the snapshots input instead

    # In case we are missing a frames limit set the limit as the number os snapshots
    if frames_limit == None:
        frames_limit = snapshots

    # Make a copy of the structure to be used as a template, so we do not modify the original coordinates
    reference_structure = structure.copy()

    # Generate an auxiliar PDB to be used to reduce the trajectory
    reference_structure_sample_pdb = get_auxiliar_filepath('.reference.pdb')
    reference_structure.generate_pdb_file(reference_structure_sample_pdb)

    # Set a maximum number of frames
    # If trajectory has more frames than the limit create a reduced trajectory
    reduced_trajectory, frames_step, frames_count = get_reduced_pytraj_trajectory(
        reference_structure_sample_pdb,
        trajectory_filepath,
        snapshots,
        frames_limit
    )

    # Delete the auxiliar PDB sine it is not needed anymore
    os.remove(reference_structure_sample_pdb)

    def frames_generator():
        # Create a progress bar
        if pbar_bool: pbar = tqdm(initial=0, desc=' Frames', total=frames_count, unit='frame')
        # Or print an empty line for the reprint to not delete a previous log
        else: print()
        # Extract each frame in pdb format
        for frame_number, frame in enumerate(reduced_trajectory.iterframe()):
            # Update the current frame log
            if pbar_bool: pbar.update(1); pbar.refresh()
            else: reprint(f'Frame {frame_number+1} ({frame_number} / {frames_count})')
            current_frame_filepath = get_auxiliar_filepath(f'{output_frames_prefix}{frame_number+1}.pdb')
            reference_structure.set_new_coordinates(frame.xyz)
            reference_structure.generate_pdb_file(current_frame_filepath)
            yield current_frame_filepath
            # Delete current frame file before going for the next frame
            os.remove(current_frame_filepath)

    return frames_generator(), frames_step, frames_count

# Build a generator which returns frames at the start of the trajectory in pdb format
# The frames limit is the maximum number of frames to be iterated
# Note that the number of frames iterated may be less than the specified number
def get_starting_pdb_frames (
    structure : 'Structure',
    trajectory_filepath : str,
    snapshots : int,
    frames_limit : Optional[int] = None,
    output_frames_prefix : str = 'frame',
    pbar_bool : bool = False,
) -> Generator[str, None, None]:

    # In case we are missing a frames limit set the limit as the number os snapshots
    if frames_limit == None:
        frames_limit = snapshots

    # Make a copy of the structure to be used as a template, so we do not modify the original coordinates
    reference_structure = structure.copy()

    # Generate an auxiliar PDB to be used to reduce the trajectory
    reference_structure_sample_pdb = get_auxiliar_filepath('.reference.pdb')
    reference_structure.generate_pdb_file(reference_structure_sample_pdb)

    # Load the trajectory using pytraj
    trajectory = get_pytraj_trajectory(reference_structure_sample_pdb, trajectory_filepath)

    # Delete the auxiliar PDB sine it is not needed anymore
    os.remove(reference_structure_sample_pdb)

    def frames_generator():
        # Create a progress bar
        if pbar_bool: pbar = tqdm(initial=0, desc=' Frames', total=frames_limit, unit='frame')
        # Or print an empty line for the reprint to not delete a previous log
        else: print()
        # Extract each frame in pdb format
        for frame_number, frame in enumerate(trajectory.iterframe(), 1):
            if frame_number >= frames_limit: break
            # Update the current frame log
            if pbar_bool: pbar.update(1); pbar.refresh()
            else: reprint(f'Frame {frame_number+1} ({frame_number} / {frames_limit})')
            current_frame_filepath = get_auxiliar_filepath(f'{output_frames_prefix}{frame_number+1}.pdb')
            reference_structure.set_new_coordinates(frame.xyz)
            reference_structure.generate_pdb_file(current_frame_filepath)
            yield current_frame_filepath
            # Delete current frame file before going for the next frame
            os.remove(current_frame_filepath)

    return frames_generator()

# Get a specific trajectory frame in pdb format
# Return the name of the generated file
def get_pdb_frame (
    topology_filename : str,
    trajectory_filepath : str,
    frame : int
) -> str:

    # Load the trajectory using pytraj
    trajectory = get_pytraj_trajectory(topology_filename, trajectory_filepath)
    trajectory_frame = trajectory[frame:frame+1]
    trajectory_frame_filepath = get_auxiliar_filepath(f'frame{frame}.pdb')
    pt.write_traj(trajectory_frame_filepath, trajectory_frame, overwrite=True)
    return trajectory_frame_filepath
