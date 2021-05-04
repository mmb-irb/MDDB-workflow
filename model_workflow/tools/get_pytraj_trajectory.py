# Parse the trajectory file to a pytraj trajectory which is used further in pytraj analyses
# Also create a reduced trajectory for heavy analyses

import math
import pytraj as pt

# Get the whole trajectory
def get_pytraj_trajectory (
    input_topology_filename : str,
    input_trajectory_filename : str) -> int:
    
    # Set the pytraj trayectory and get the number of frames
    pt_trajectory = pt.iterload(input_trajectory_filename, input_topology_filename)
    
    # WARNING: This extra line prevents the error "Segment violation (core dumped)" in some pdbs
    # This happens with some random pdbs which pytraj considers to have 0 Mols
    # More info: https://github.com/Amber-MD/cpptraj/pull/820
    # DANI: Cuando pytraj sea actualizado con la versión actual de cpptraj este error debería desaparecer
    pt_trajectory.top.start_new_mol()

    return pt_trajectory

# Set the maximum number of frames for the reduced trajectory
# WARNING: The final number of frames in some analyses may be +1
reduced_trajectory_frames = 200

# Get the reduced trajectory
def get_reduced_pytraj_trajectory (
    input_topology_filename : str,
    input_trajectory_filename : str) -> int:
    
    # Set the pytraj trayectory and get the number of frames
    pt_trajectory = get_pytraj_trajectory(input_topology_filename, input_trajectory_filename)

    # Get the number of frames
    trajectory_frames = pt_trajectory.n_frames

    # Set a reduced trajectory used for heavy analyses
    reduced_pt_trajectory = None
    # If the current trajectory has already less frames than the maximum then use it also as reduced
    if trajectory_frames < reduced_trajectory_frames:
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
        step = math.floor(trajectory_frames / reduced_trajectory_frames)
        reduced_pt_trajectory = pt_trajectory[0:trajectory_frames:step]
        # Add the step value to the reduced trajectory explicitly. It will be required later
        reduced_pt_trajectory.step = step

    return reduced_pt_trajectory