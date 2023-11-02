import mdtraj as mdt
from numpy import mean, std

from typing import List

test_name = 'intrajrity'

# LORE
# This test was originaly intended to use a RMSD jump cutoff based on number of atoms and timestep
# However, after a deep study, it was observed that simulations with similar features may show very different RMSD jumps
# For this reason now we comptue RMSD jumps along the whole trajectory and check that the biggest jump is not an outlier
# The outlier is defined according to how many times the standard deviation far from the mean is a value

# Look for sudden raises of RMSd values from one frame to another
# To do so, we check the RMSD of every frame using its previous frame as reference
def check_trajectory_integrity (
    input_structure_filename : str,
    input_trajectory_filename : str,
    structure : 'Structure',
    pbc_residues : List[int],
    mercy : List[str],
    trust: List[str],
    register : dict,
    #time_length : float,
    check_selection : str,
    # DANI: He visto saltos 'correctos' pasar de 6
    # DANI: He visto saltos 'incorrectos' no bajar de 10
    standard_deviations_cutoff : float) -> bool:

    # Skip the test if we trust
    if test_name in trust:
        return True

    # Parse the selection in VMD selection syntax
    parsed_selection = structure.select(check_selection, syntax='vmd')

    # If there is nothing to check then warn the user and stop here
    if not parsed_selection:
        raise Exception('WARNING: There are not atoms to be analyzed for the RMSD analysis')

    # Discard PBC residues from the selection to be checked
    pbc_selection = structure.select_residue_indices(pbc_residues)
    parsed_selection -= pbc_selection

    # If there is nothing to check then warn the user and stop here
    if not parsed_selection:
        raise Exception('WARNING: There are not atoms to be analyzed after PBC substraction for the RMSD analysis')

    print('Checking trajectory integrity')

    # Load the trajectory frame by frame
    trajectory = mdt.iterload(input_trajectory_filename, top=input_structure_filename, chunk=1)

    # Save the previous frame any time
    previous_frame = next(trajectory)

    # Save all RMSD jumps
    rmsd_jumps = []    

    for f, frame in enumerate(trajectory, 1):
        # Update the current frame log
        print(' Frame ' + str(f), end='\r')

        # Calculate RMSD value between previous and current frame
        rmsd_value = mdt.rmsd(frame, previous_frame, atom_indices=parsed_selection.atom_indices)[0]
        rmsd_jumps.append(rmsd_value)

        # Update the previous frame as the current one
        previous_frame = frame

    # If the trajectory has only 1 or 2 frames then there is no test to do
    if len(rmsd_jumps) <= 1:
        return True

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
    max_z_score = 0
    max_z_score_frame = 0
    for i, rmsd_jump in enumerate(rmsd_jumps):
        z_score = abs( (rmsd_jump - mean_rmsd_jump) / stdv_rmsd_jump )
        # Keep track of the maixmum z score
        if z_score > max_z_score:
            max_z_score = z_score
            max_z_score_frame = i
        # If z score bypassed the limit then report it
        if z_score > standard_deviations_cutoff:
            # If there are as many bypassed frames as the index then it means no frame has passed the cutoff yet
            if i == bypassed_frames:
                bypassed_frames += 1
                continue
            if outliers_count >= 4:
                print(' etc...')
                break
            print(' FAIL: Sudden RMSD jump between frames ' + str(i) + ' and ' + str(i+1))
            outliers_count += 1

    # Always print the maximum z score and its frames
    print(' Maximum z score (' + str(max_z_score) + ') reported between frames ' + str(max_z_score_frame) + ' and ' + str(max_z_score_frame + 1))

    # If there were any outlier then the check has failed
    if outliers_count > 0:
        # Add a warning an return True since the test failed in case we have mercy
        message = 'RMSD check has failed: there may be sudden jumps along the trajectory'
        if test_name in mercy:
            register['warnings'].append(message)
            return False
        # Otherwise kill the process right away
        raise Exception(message)

    # Warn the user if we had bypassed frames
    if bypassed_frames > 0:
        print(' WARNING: First ' + str(bypassed_frames) + ' frames may be not equilibrated')
        register['warnings'].append('First ' + str(bypassed_frames) + ' frames may be not equilibrated')

    print(' Test has passed successfully')
    return True