import pytraj as pt
 
# Get the average structure
# This process is carried by pytraj, since the Gromacs average may be displaced
# Output may be both pdb or xtc
def get_average (
    pytraj_trajectory,
    output_average_filename : str):

    # Create a new frame with the average positions
    # WARNING: Do not pass the argument 'autoimage=True'
    # WARNING: Autoimage makes some trajectories get displaced the same as in Gromacs
    average_frame = pt.mean_structure(pytraj_trajectory())

    # In order to export it, first create an empty trajectory only with the topology
    # Then add the average frame and write it to 'xtc' format
    average = pt.Trajectory(top=pytraj_trajectory.top)
    average.append(average_frame)
    pt.write_traj(output_average_filename, average, overwrite=True)