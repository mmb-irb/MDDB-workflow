import mdtraj as mdt
import pytraj as pt
import numpy as np
import time

from model_workflow.utils.auxiliar import delete_previous_log, reprint, TestFailure, warn
from model_workflow.utils.constants import TRAJECTORY_INTEGRITY_FLAG
from model_workflow.utils.pyt_spells import get_pytraj_trajectory
from model_workflow.utils.type_hints import *
from tqdm import tqdm

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
    pbc_selection : 'Selection',
    mercy : List[str],
    trust: List[str],
    register : 'Register',
    #time_length : float,
    check_selection : str,
    # DANI: He visto saltos 'correctos' pasar de 6
    # DANI: He visto saltos 'incorrectos' no bajar de 10
    standard_deviations_cutoff : float,
    snapshots : int) -> bool:

    # Skip the test if we trust
    if TRAJECTORY_INTEGRITY_FLAG in trust:
        return True

    # Skip the test if it is already passed according to the register
    if register.tests.get(TRAJECTORY_INTEGRITY_FLAG, None):
        return True

    # Remove old warnings
    register.remove_warnings(TRAJECTORY_INTEGRITY_FLAG)

    # Parse the selection in VMD selection syntax
    parsed_selection = structure.select(check_selection, syntax='vmd')

    # If there is nothing to check then warn the user and stop here
    if not parsed_selection:
        raise Exception('WARNING: There are not atoms to be analyzed for the RMSD analysis')

    # Discard PBC residues from the selection to be checked
    parsed_selection -= pbc_selection

    # If there is nothing to check then warn the user and stop here
    if not parsed_selection:
        warn('There are no atoms to be analyzed for the RMSD checking after PBC substraction')
        register.update_test(TRAJECTORY_INTEGRITY_FLAG, 'na')
        return True

    print('Checking trajectory integrity')

    # Load the trajectory frame by frame
    trajectory = mdt.iterload(input_trajectory_filename, top=input_structure_filename, chunk=1)

    # Save the previous frame any time
    previous_frame = next(trajectory)

    # Save all RMSD jumps
    rmsd_jumps = []

    # Initialize progress bar first
    pbar = tqdm(enumerate(trajectory, 1), total=snapshots, desc=' Frame', unit='frame', initial=1)
    
    for f, frame in pbar:
        # Calculate RMSD value between previous and current frame
        rmsd_value = mdt.rmsd(frame, previous_frame, atom_indices=parsed_selection.atom_indices)[0]
        rmsd_jumps.append(rmsd_value)

        # Update the previous frame as the current one
        previous_frame = frame
    time.sleep(0.1) # Needed for progress bar to be print at the correct place

    # If the trajectory has only 1 or 2 frames then there is no test to do
    if len(rmsd_jumps) <= 1:
        register.update_test(TRAJECTORY_INTEGRITY_FLAG, True)
        return True

    # Get the maximum RMSD value and check it is a reasonable deviation from the average values
    # Otherwise, if it is an outlier, the test fails
    mean_rmsd_jump = np.mean(rmsd_jumps)
    stdv_rmsd_jump = np.std(rmsd_jumps)

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
            print(f' FAIL: Sudden RMSD jump between frames {i} and {i+1}. RMSD jump: {rmsd_jump:4f}')
            outliers_count += 1

    # Always print the maximum z score and its frames
    print(f' Maximum z score {max_z_score:4f} reported between frames {max_z_score_frame} and {max_z_score_frame + 1}.\n'
          f' Mean: {mean_rmsd_jump:4f}. Stdv: {stdv_rmsd_jump:4f}')

    # If there were any outlier then the check has failed
    if outliers_count > 0:
        # Add a warning an return True since the test failed in case we have mercy
        message = 'RMSD check has failed: there may be sudden jumps along the trajectory'
        if TRAJECTORY_INTEGRITY_FLAG in mercy:
            register.add_warning(TRAJECTORY_INTEGRITY_FLAG, message)
            register.update_test(TRAJECTORY_INTEGRITY_FLAG, False)
            return False
        # Otherwise kill the process right away
        raise TestFailure(message)

    # Warn the user if we had bypassed frames
    if bypassed_frames > 0:
        register.add_warning(TRAJECTORY_INTEGRITY_FLAG, f'First {bypassed_frames} frames may be not equilibrated')

    print(' Test has passed successfully')
    register.update_test(TRAJECTORY_INTEGRITY_FLAG, True)
    return True

# Look for sudden raises of RMSd values from one frame to another
# To do so, we check the RMSD of every frame using its previous frame as reference
# DANI: Dependemos de que cambien el comportamiento de MDtraj para que esto funcione
# DANI: https://github.com/mdtraj/mdtraj/issues/1966
def check_trajectory_integrity_per_fragment (
    input_structure_filename : str,
    input_trajectory_filename : str,
    structure : 'Structure',
    pbc_selection : 'Selection',
    mercy : List[str],
    trust: List[str],
    register : 'Register',
    #time_length : float,
    check_selection : str,
    # DANI: He visto saltos 'correctos' pasar de 6
    # DANI: He visto saltos 'incorrectos' no bajar de 10
    standard_deviations_cutoff : float) -> bool:

    # Skip the test if we trust
    if TRAJECTORY_INTEGRITY_FLAG in trust:
        return True

    # Skip the test if it is already passed according to the register
    if register.tests.get(TRAJECTORY_INTEGRITY_FLAG, None):
        return True

    # Remove old warnings
    register.remove_warnings(TRAJECTORY_INTEGRITY_FLAG)

    # Parse the selection in VMD selection syntax
    parsed_selection = structure.select(check_selection, syntax='vmd')

    # If there is nothing to check then warn the user and stop here
    if not parsed_selection:
        raise Exception('WARNING: There are not atoms to be analyzed for the RMSD analysis')

    # Discard PBC residues from the selection to be checked
    parsed_selection -= pbc_selection

    # If there is nothing to check then warn the user and stop here
    if not parsed_selection:
        warn('There are no atoms to be analyzed for the RMSD checking after PBC substraction')
        register.update_test(TRAJECTORY_INTEGRITY_FLAG, 'na')
        return True
    
    # Get fragments out of the parsed selection
    # Fragments will be analyzed independently
    # A small fragment may cause a small RMSD perturbation when it is part of a large structure
    # Note that jumps of partial fragments along boundaries are rare imaging problems
    # Usually are whole fragments the ones which jump
    fragments = list(structure.find_fragments(parsed_selection))

    print(f'Checking trajectory integrity ({len(fragments)} fragments)')

    # Load the trajectory frame by frame
    trajectory = mdt.iterload(input_trajectory_filename, top=input_structure_filename, chunk=1)

    # Save the previous frame any time
    previous_frame = next(trajectory)

    # Save all RMSD jumps
    fragment_rmsd_jumps = { fragment: [] for fragment in fragments }

    # Add an extra breakline before the first log
    print()

    # Iterate trajectory frames
    for f, frame in enumerate(trajectory, 1):
        # Update the current frame log
        reprint(f' Frame {f}')

        # Iterate over the different fragments
        for fragment in fragment_rmsd_jumps:
            # Calculate RMSD value between previous and current frame
            # DANI: El centrado de MDtraj elimina el salto a través de las boundaries
            # DANI: El precentered=True debería evitarlo, pero es ignorado si hay atom_indices
            rmsd_value = mdt.rmsd(frame, previous_frame, atom_indices=fragment.atom_indices)[0]
            fragment_rmsd_jumps[fragment].append(rmsd_value)

        # Update the previous frame as the current one
        previous_frame = frame

    # First frames may not be perfectly equilibrated and thus have stronger RMSD jumps
    # For this reason we allow the first frames to bypass the check
    # As soon as one frame is below the cutoff the bypass is finished for the following frames
    # Count the maximum number of bypassed frames and warn the user in case there are any
    max_bypassed_frames = 0

    # Capture outliers
    # If we capture more than 5 we stop searching
    outliers_count = 0
    max_z_score = 0
    max_z_score_frame = 0

    # Iterate over the different fragments
    for fragment, rmsd_jumps in fragment_rmsd_jumps.items():

        # If the trajectory has only 1 or 2 frames then there is no test to do
        if len(rmsd_jumps) <= 1:
            register.update_test(TRAJECTORY_INTEGRITY_FLAG, True)
            return True

        # Get the maximum RMSD value and check it is a reasonable deviation from the average values
        # Otherwise, if it is an outlier, the test fails
        mean_rmsd_jump = np.mean(rmsd_jumps)
        stdv_rmsd_jump = np.std(rmsd_jumps)

        # Count the number of bypassed frames for this fragment
        bypassed_frames = 0

        fragment_max_z_score = 0
        fragment_max_z_score_frame = 0

        for i, rmsd_jump in enumerate(rmsd_jumps):
            z_score = abs( (rmsd_jump - mean_rmsd_jump) / stdv_rmsd_jump )
            # Keep track of the maixmum z score for this fragment
            if z_score > fragment_max_z_score:
                fragment_max_z_score = z_score
                fragment_max_z_score_frame = i
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
                print(f' FAIL: Sudden RMSD jump between frames {i} and {i+1}')
                outliers_count += 1

        print(f'Fragment {structure.name_selection(fragment)} -> {fragment_max_z_score} in frame {fragment_max_z_score_frame}')

        # Update the maximum number of bypassed frames
        if bypassed_frames > max_bypassed_frames:
            max_bypassed_frames = bypassed_frames

        # If there were any outlier then the check has failed
        if outliers_count > 0:
            # Add a warning an return True since the test failed in case we have mercy
            message = 'RMSD check has failed: there may be sudden jumps along the trajectory'
            if TRAJECTORY_INTEGRITY_FLAG in mercy:
                register.add_warning(TRAJECTORY_INTEGRITY_FLAG, message)
                register.update_test(TRAJECTORY_INTEGRITY_FLAG, False)
                return False
            # Otherwise kill the process right away
            raise TestFailure(message)
        
    # Always print the maximum z score and its frames
    print(f' Maximum z score {max_z_score} reported between frames {max_z_score_frame} and {max_z_score_frame + 1}')

    # Warn the user if we had bypassed frames
    if max_bypassed_frames > 0:
        register.add_warning(TRAJECTORY_INTEGRITY_FLAG, f'First {max_bypassed_frames} frames may be not equilibrated')

    print(' Test has passed successfully')
    register.update_test(TRAJECTORY_INTEGRITY_FLAG, True)
    return True

# Look for sudden raises of RMSd values from one frame to another
# To do so, we check the RMSD of every frame using its previous frame as reference
# DANI: Dependemos de que cambien el comportamiento de MDtraj para que esto funcione
# DANI: https://github.com/mdtraj/mdtraj/issues/1966
def check_trajectory_integrity_per_fragment_2 (
    input_structure_filename : str,
    input_trajectory_filename : str,
    structure : 'Structure',
    pbc_selection : 'Selection',
    mercy : List[str],
    trust: List[str],
    register : 'Register',
    #time_length : float,
    check_selection : str,
    # DANI: He visto saltos 'correctos' pasar de 6
    # DANI: He visto saltos 'incorrectos' no bajar de 10
    standard_deviations_cutoff : float) -> bool:

    # Skip the test if we trust
    if TRAJECTORY_INTEGRITY_FLAG in trust:
        return True

    # Skip the test if it is already passed according to the register
    if register.tests.get(TRAJECTORY_INTEGRITY_FLAG, None):
        return True

    # Remove old warnings
    register.remove_warnings(TRAJECTORY_INTEGRITY_FLAG)

    # Parse the selection in VMD selection syntax
    parsed_selection = structure.select(check_selection, syntax='vmd')

    # If there is nothing to check then warn the user and stop here
    if not parsed_selection:
        raise Exception('WARNING: There are not atoms to be analyzed for the RMSD analysis')

    # Discard PBC residues from the selection to be checked
    parsed_selection -= pbc_selection

    # If there is nothing to check then warn the user and stop here
    if not parsed_selection:
        warn('There are no atoms to be analyzed for the RMSD checking after PBC substraction')
        register.update_test(TRAJECTORY_INTEGRITY_FLAG, 'na')
        return True
    
    # First frames may not be perfectly equilibrated and thus have stronger RMSD jumps
    # For this reason we allow the first frames to bypass the check
    # As soon as one frame is below the cutoff the bypass is finished for the following frames
    # Count the maximum number of bypassed frames and warn the user in case there are any
    max_bypassed_frames = 0
    
    # Get fragments out of the parsed selection
    # Fragments will be analyzed independently
    # A small fragment may cause a small RMSD perturbation when it is part of a large structure
    # Note that jumps of partial fragments along boundaries are rare imaging problems
    # Usually are whole fragments the ones which jump
    fragments = list(structure.find_fragments(parsed_selection))

    print(f'Checking trajectory integrity')

    # Iterate over the different fragments
    for f, fragment in enumerate(fragments, 1):
        fragment_name = structure.name_selection(fragment)
        print(f'Fragment {f}/{len(fragments)} -> {fragment_name}')

        # Load the trajectory frame by frame
        trajectory = mdt.iterload(
            input_trajectory_filename,
            top=input_structure_filename,
            atom_indices=fragment.atom_indices,
            chunk=1)
        
        # Save all RMSD jumps
        rmsd_jumps = []

        # Save the previous frame any time
        previous_frame = next(trajectory)

        # Add an extra breakline before the first log
        print()

        # Iterate trajectory frames
        for f, frame in enumerate(trajectory, 1):
            # Update the current frame log
            reprint(f' Frame {f}')
        
            # Calculate RMSD value between previous and current frame
            # DANI: El centrado de MDtraj elimina el salto a través de las boundaries
            # DANI: El precentered=True debería evitarlo, pero es ignorado si hay atom_indices
            rmsd_value = mdt.rmsd(frame, previous_frame, precentered=True)[0]
            rmsd_jumps.append(rmsd_value)

            # Update the previous frame as the current one
            previous_frame = frame

        # If the trajectory has only 1 or 2 frames then there is no test to do
        if len(rmsd_jumps) <= 1:
            register.update_test(TRAJECTORY_INTEGRITY_FLAG, True)
            return True

        # Get the maximum RMSD value and check it is a reasonable deviation from the average values
        # Otherwise, if it is an outlier, the test fails
        mean_rmsd_jump = np.mean(rmsd_jumps)
        stdv_rmsd_jump = np.std(rmsd_jumps)

        # Count the number of bypassed frames for this fragment
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
                print(f' FAIL: Sudden RMSD jump between frames {i} and {i+1}')
                outliers_count += 1

        # Always print the maximum z score and its frames
        print(f' Maximum z score {max_z_score} reported between frames {max_z_score_frame} and {max_z_score_frame + 1}')

        # Update the maximum number of bypassed frames
        if bypassed_frames > max_bypassed_frames:
            max_bypassed_frames = bypassed_frames

        # If there were any outlier then the check has failed
        if outliers_count > 0:
            # Add a warning an return True since the test failed in case we have mercy
            message = 'RMSD check has failed: there may be sudden jumps along the trajectory'
            if TRAJECTORY_INTEGRITY_FLAG in mercy:
                register.add_warning(TRAJECTORY_INTEGRITY_FLAG, message)
                register.update_test(TRAJECTORY_INTEGRITY_FLAG, False)
                return False
            # Otherwise kill the process right away
            raise TestFailure(message)

    # Warn the user if we had bypassed frames
    if max_bypassed_frames > 0:
        register.add_warning(TRAJECTORY_INTEGRITY_FLAG, f'First {max_bypassed_frames} frames may be not equilibrated')

    print(' Test has passed successfully')
    register.update_test(TRAJECTORY_INTEGRITY_FLAG, True)
    return True

# Compute every residue RMSD to check if there are sudden jumps along the trajectory
# HARDCODE: This function is not fully implemented but enabled manually for specific cases
def check_trajectory_integrity_per_residue (
    input_structure_filename : str,
    input_trajectory_filename : str,
    structure : 'Structure',
    pbc_selection : 'Selection',
    mercy : List[str],
    trust: List[str],
    register : 'Register',
    #time_length : float,
    check_selection : str,
    # DANI: He visto saltos 'correctos' pasar de 11
    # DANI: He visto saltos 'incorrectos' no bajar de 14
    standard_deviations_cutoff : float):

    # HARDCODE: The default value does not work for a single residue
    standard_deviations_cutoff = 12

    # Skip the test if we trust
    if TRAJECTORY_INTEGRITY_FLAG in trust:
        return True

    # Skip the test if it is already passed according to the register
    if register.tests.get(TRAJECTORY_INTEGRITY_FLAG, None):
        return True

    # Remove old warnings
    register.remove_warnings(TRAJECTORY_INTEGRITY_FLAG)

    # Parse the selection in VMD selection syntax
    parsed_selection = structure.select(check_selection, syntax='vmd')

    # If there is nothing to check then warn the user and stop here
    if not parsed_selection:
        raise Exception('WARNING: There are not atoms to be analyzed for the RMSD analysis')

    # Discard PBC residues from the selection to be checked
    parsed_selection -= pbc_selection

    # If there is nothing to check then warn the user and stop here
    if not parsed_selection:
        warn('There are no atoms to be analyzed for the RMSD checking after PBC substraction')
        register.update_test(TRAJECTORY_INTEGRITY_FLAG, 'na')
        return True

    # We must filter out residues which only have 1 atom (e.g. ions)
    # This is because sometimes pytraj does not return results for them and then the number of results and residues does not match
    # More info: https://github.com/Amber-MD/pytraj/issues/1580
    ion_atom_indices = []
    for residue in structure.residues:
        if len(residue.atom_indices) == 1:
            ion_atom_indices += residue.atom_indices
    ions_selection = structure.select_atom_indices(ion_atom_indices)
    parsed_selection -= ions_selection

    # Parse the selection to a pytraj mask
    pytraj_selection = parsed_selection.to_pytraj()

    # Calculate the residue indices of the overall structure remaining in the filtered trajectory
    residue_indices = structure.get_selection_residue_indices(parsed_selection)
    n_residues = len(residue_indices)

    print('Checking trajectory integrity per residue')

    # Parse the trajectory into pytraj and apply the mask
    # NEVER FORGET: The pytraj iterload does not accept a mask, but we apply the mask later in the analysis
    pt_trajectory = get_pytraj_trajectory(input_structure_filename, input_trajectory_filename, atom_selection = parsed_selection)

    # Make sure the expected output number of residues to match with the pytraj results
    # These numbers may not match when ions are included so we better check
    # NEVER FORGET: The pytraj TrajectoryIterator is not an iterator
    first_frame = pt_trajectory[0:1]
    # DANI: When the 'resname' argument is missing it prints "Error: Range::SetRange(None): Range is -1 for None"
    # DANI: However there is no problem and the analysis runs flawlessly
    # DANI: For this reason we call this function with no resname and then we remove the log
    data_sample = pt.rmsd_perres(first_frame, ref=first_frame, perres_mask=pytraj_selection)
    # We remove the previous error log
    delete_previous_log()
    # We remove the first result, which is meant to be the whole rmsd and whose key is 'RMSD_00001'
    del data_sample[0]
    if n_residues != len(data_sample):
        raise ValueError(f'Number of target residues ({n_residues}) does not match number of residues in data ({len(data_sample)})')

    # Saving all RMSD jumps may take a lot of memory
    # Instead we will store the sum of values and the maximum
    # This way we can caluclate the average value at the end and check if the maximum is too far from it
    rmsd_per_residue_per_frame = []

    # Add an extra breakline before the first log
    print()

    # Iterate trajectory frames
    previous_frame_trajectory = first_frame
    frame_number = 1
    for frame in pt_trajectory:
        # Update the current frame log
        reprint(f' Frame {frame_number}')
        # Set a pytraj trajectory out of a pytraj frame
        frame_trajectory = pt.Trajectory(top=pt_trajectory.top)
        frame_trajectory.append(frame)
        # Run the analysis in pytraj
        # The result data is a custom pytraj class: pytraj.datasets.datasetlist.DatasetList
        # This class has keys but its attributes can not be accessed through the key
        # They must be accessed thorugh the index
        # DANI: When the 'resname' argument is missing it prints "Error: Range::SetRange(None): Range is -1 for None"
        # DANI: However there is no problem and the analysis runs flawlessly
        # DANI: Adding resrage as a list/range was tried and did not work, only string works
        # DANI: Adding a string resrange however strongly impacts the speed when this function is called repeatedly
        # DANI: For this reason we call this function with no resname and then we remove the log
        rmsd_per_residue = pt.rmsd_perres(frame_trajectory, ref=previous_frame_trajectory, mask=pytraj_selection)
        # We remove the previous error log
        delete_previous_log()
        # We remove the first result, which is meant to be the whole rmsd and whose key is 'RMSD_00001'
        del rmsd_per_residue[0]
        # Check we have no NaNs
        if np.isnan(rmsd_per_residue[0][0]):
            raise ValueError(f'We are having NaNs at frame {frame_number}')
        # Add last values to the list
        rmsd_per_residue_per_frame.append(rmsd_per_residue)
        # rmsd_per_residue_per_frame[:, frame]
        # Now update data for every residue
        # for index, residue_rmsd in enumerate(rmsd_per_residue):
        #     current_rmsd = residue_rmsd[0]
        #     # Get the current residue rmsd data
        #     total_rmsd = rmsd_per_residue_per_frame[index]
        #     total_rmsd['accumulated'] += current_rmsd
        #     total_rmsd['maximum'] = max(total_rmsd['maximum'], current_rmsd)
        # Update previous coordinates
        previous_frame_trajectory = frame_trajectory
        # Update the frame_number
        frame_number += 1

    # If the trajectory has only 1 or 2 frames then there is no test to do
    n_jumps = len(rmsd_per_residue_per_frame)
    if n_jumps <= 1:
        register.update_test(TRAJECTORY_INTEGRITY_FLAG, True)
        return True

    # Keep the overall maximum z score, and its residue and frame for the logs
    overall_max_z_score = 0
    overall_max_z_score_frame = None
    overall_max_z_score_residue = None

    # Keep the overall maximum bypassed frames number
    overall_bypassed_frames = 0

    # Keep the overall count of residues with outliers
    overall_outliered_residues = 0

    # Add an extra breakline before the next log
    print()

    # Now check there are not sudden jumps for each residue separattely
    for residue_number in range(n_residues):
        reprint(f' Residue {residue_number+1}')
        # Get the rmsd jumps for each frame for this specific residue
        rmsd_jumps = [ frame[residue_number] for frame in rmsd_per_residue_per_frame ]

        # Get the maximum RMSD value and check it is a reasonable deviation from the average values
        # Otherwise, if it is an outlier, the test fails
        mean_rmsd_jump = np.mean(rmsd_jumps)
        stdv_rmsd_jump = np.std(rmsd_jumps)

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
                # Otherwise we consider this as an outlier and thus the test has failed
                # However we keep checking just to find and report the highest outlier
                outliers_count += 1

        # Update the overall bypassed frames if we overcomed it
        if overall_bypassed_frames < bypassed_frames:
            overall_bypassed_frames = bypassed_frames

        # Update the overall max z score if we overcome it
        if max_z_score > overall_max_z_score:
            overall_max_z_score = max_z_score
            overall_max_z_score_frame = max_z_score_frame
            overall_max_z_score_residue = residue_number

        # If there were any outlier then add one to the overall count
        if outliers_count > 0:
            overall_outliered_residues += 1

    # Always print the overall maximum z score and its frames and residue
    overall_max_z_score_residue_label = pt_trajectory.top.residue(overall_max_z_score_residue)
    print(f' Maximum z score {overall_max_z_score} reported for residue {overall_max_z_score_residue_label} between frames {overall_max_z_score_frame} and {overall_max_z_score_frame + 1}')

    # If there were any outlier then the check has failed
    if overall_outliered_residues > 0:
        # Add a warning an return True since the test failed in case we have mercy
        message = 'RMSD check has failed: there may be sudden jumps along the trajectory'
        if TRAJECTORY_INTEGRITY_FLAG in mercy:
            register.add_warning(TRAJECTORY_INTEGRITY_FLAG, message)
            register.update_test(TRAJECTORY_INTEGRITY_FLAG, False)
            return False
        # Otherwise kill the process right away
        raise TestFailure(message)

    # Warn the user if we had bypassed frames
    if overall_bypassed_frames > 0:
        register.add_warning(TRAJECTORY_INTEGRITY_FLAG, f'First {overall_bypassed_frames} frames may be not equilibrated')

    print(' Test has passed successfully')
    register.update_test(TRAJECTORY_INTEGRITY_FLAG, True)
    return True