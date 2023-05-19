# Parse the trajectory file to a pytraj trajectory which is used further in pytraj analyses
# Also create a reduced trajectory for heavy analyses

import math
import pytraj as pt
from distutils.version import StrictVersion

# Get the whole trajectory
def get_pytraj_trajectory (
    input_topology_filename : str,
    input_trajectory_filename : str):
    
    # Set the pytraj trayectory and get the number of frames
    pt_trajectory = pt.iterload(input_trajectory_filename, input_topology_filename)
    
    # WARNING: This extra line prevents the error "Segment violation (core dumped)" in some pdbs
    # This happens with some random pdbs which pytraj considers to have 0 Mols
    # More info: https://github.com/Amber-MD/cpptraj/pull/820
    # DANI: Esto es útil en pytraj <= 2.0.5 pero hace fallar el código a partir de pytraj 2.0.6
    if StrictVersion(pt.__version__) <= StrictVersion('2.0.5'):
        pt_trajectory.top.start_new_mol()

    return pt_trajectory

# Get the reduced trajectory
# WARNING: The final number of frames may be the specifided or less
def get_reduced_pytraj_trajectory (
    input_topology_filename : str,
    input_trajectory_filename : str,
    snapshots : int,
    reduced_trajectory_frames_limit : int):
    
    # Set the pytraj trayectory and get the number of frames
    pt_trajectory = get_pytraj_trajectory(input_topology_filename, input_trajectory_filename)
    # WARNING: Do not read pt_trajectory.n_frames to get the number of snapshots or you will read the whole trajectory
    # WARNING: This may be a lot of time for a huge trajectory. Use the snapshots input instead

    # Set a reduced trajectory used for heavy analyses
    # If the current trajectory has already less or the same frames than the limit
    # Then do nothing and use it also as reduced
    if snapshots <= reduced_trajectory_frames_limit:
        reduced_pt_trajectory = pt_trajectory
        # Add a step value which will be required later
        reduced_pt_trajectory.step = 1
    # Otherwise, create a reduced trajectory with as much frames as specified above
    # These frames are picked along the trajectory
    else:
        # Calculate the step between frames in the reduced trajectory to match the final number of frames
        # WARNING: Since the step must be an integer the thorical step must be rounded
        # This means the specified final number of frames may not be accomplished, but it is okey
        # WARNING: Since the step is rounded with the math.ceil function it will always be rounded up
        # This means the final number of frames will be the specified or less
        # CRITICAL WARNING:
        # This formula is exactly the same that the client uses to request stepped frames to the API
        # This means that the client and the workflow are coordinated and these formulas must not change
        # If you decide to change this formula (in both workflow and client)...
        # You will have to run again all the database analyses with reduced trajectories
        step = math.ceil(snapshots / reduced_trajectory_frames_limit)
        reduced_pt_trajectory = pt_trajectory[0:snapshots:step]
        # Add the step value to the reduced trajectory explicitly. It will be required later
        reduced_pt_trajectory.step = step

    return reduced_pt_trajectory