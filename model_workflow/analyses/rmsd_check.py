import mdtraj as mdt
from numpy import mean, std

CURSOR_UP_ONE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'
ERASE_PREVIOUS_LINE = CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE

# LORE
# This test was originaly intended to use a RMSD jump cutoff based on number of atoms and timestep
# However, after a deep study, it was observed that simulations with similar features may show very different RMSD jumps
# For this reason now we comptue RMSD jumps along the whole trajectory and check that the biggest jump is not an outlier
# The outlier is defined according to how many times the standard deviation far from the mean is a value
# DANI: He visto saltos 'correctos' pasar de 6
# DANI: He visto saltos 'incorrectos' no bajar de 10
standard_deviations_cutoff = 9

# Look for sudden raises of RMSd values from one frame to another
# To do so, we check the RMSD of every frame using its previous frame as reference
def check_sudden_jumps (
    input_structure_filename : str,
    input_trajectory_filename : str,
    structure : str,
    time_length : float,
    snapshots : int,
    check_selection : str = 'protein or nucleic',
    ) -> bool:

    # If the trajectory has only 1 value then there is no test to do
    if snapshots <= 1:
        return False

    # Parse the selection in VMD selection syntax
    parsed_selection = structure.select(check_selection, syntax='vmd')

    # If there is nothing to check then warn the user and stop here
    if not parsed_selection or len(parsed_selection.atom_indices) == 0:
        print('WARNING: There are not atoms to be analyzed for the RMSD analysis')
        return

    print('Checking trajectory integrity')

    # Load the trajectory frame by frame
    trajectory = mdt.iterload(input_trajectory_filename, top=input_structure_filename, chunk=1)

    # Save the previous frame any time
    previous_frame = next(trajectory)

    # Save all RMSD jumps
    rmsd_jumps = []    

    # Print an empty line for the first 'ERASE_PREVIOUS_LINE' to not delete a previous log
    print()

    for f, frame in enumerate(trajectory, 1):
        # Update the current frame log
        print(ERASE_PREVIOUS_LINE)
        print(' Frame ' + str(f))

        # Calculate RMSD value between previous and current frame
        rmsd_value = mdt.rmsd(frame, previous_frame, atom_indices=parsed_selection.atom_indices)[0]
        rmsd_jumps.append(rmsd_value)

        # Update the previous frame as the current one
        previous_frame = frame

    # Get the maximum RMSD value and check it is a reasonable deviation from the average values
    # Otherwise, if it is an outlier, the test fails
    mean_rmsd_jump = mean(rmsd_jumps)
    stdv_rmsd_jump = std(rmsd_jumps)

    # First frames may not be perfectly equilibrated and thus have stronger RMSD jumps
    # For this reason we allow the first frames to bypass the check
    # As soon as one frame is below the cutoff the bypass is finished for the following frames
    # Count the number of bypassed frames and warn the user in case there are any
    bypassed_frames = 0

    # Capture outliers
    # If we capture more than 5 we stop searching
    outliers_count = 0
    for i, rmsd_jump in enumerate(rmsd_jumps):
        z_score = (rmsd_jump - mean_rmsd_jump) / stdv_rmsd_jump
        if abs(z_score) > standard_deviations_cutoff:
            # If there are as many bypassed frames as the index then it means no frame has passed the cutoff yet
            if i == bypassed_frames:
                bypassed_frames += 1
                continue
            if outliers_count >= 4:
                print(' etc...')
                break
            print(' FAIL: Sudden RMSD jump between frames ' + str(i) + ' and ' + str(i+1))
            outliers_count += 1

    # If there were any outlier then the check has failed
    if outliers_count > 0:
        return True

    # Warn the user if we had bypassed frames
    if bypassed_frames > 0:
        print(' WARNING: First ' + str(bypassed_frames) + ' frames may be not equilibrated')

    return False