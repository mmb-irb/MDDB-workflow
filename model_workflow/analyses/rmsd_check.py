import mdtraj as mdt

CURSOR_UP_ONE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'
ERASE_PREVIOUS_LINE = CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE

# Set how big can be the RMSD jump between 2 frames
# If the jump is bigger then we return True
rmsd_cutoff = 1

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

    # Parse the selection in VMD selection syntax
    parsed_selection = structure.select(check_selection, syntax='vmd')

    # If there is nothing to check then warn the user and stop here
    if not parsed_selection or len(parsed_selection.atom_indices) == 0:
        print('WARNING: There are not atoms to be analyzed for the RMSD analysis')
        return

    print('Checking trajectory integrity')

    # Set the RMSD cutoff based in the number of atoms: the greater the structure the more flexible the RMSD cutoff
    natoms = len(parsed_selection.atom_indices)
    print(' Number of atoms evaluated: ' + str(natoms))
    # Set the RMSD cutoff based in the time step between frames: the longest time the more flexible the RMSD cutoff
    timestep = round((time_length / snapshots) * 100) / 100
    print(' Frame timestep: ' + str(timestep) + ' ns')
    # Set the RMSD cutoff
    rmsd_cutoff = round((timestep * natoms * 0.01) * 1000) / 1000 
    print(' RMSD jump cutoff -> ' + str(rmsd_cutoff) + ' Å')

    # Load the trajectory frame by frame
    trajectory = mdt.iterload(input_trajectory_filename, top=input_structure_filename, chunk=1)

    # Print an empty line for the first 'ERASE_PREVIOUS_LINE' to not delete a previous log
    print()

    # Save the previous frame any time
    previous_frame = next(trajectory)

    for f, frame in enumerate(trajectory, 1):
        # Update the current frame log
        print(ERASE_PREVIOUS_LINE)
        print(' Frame ' + str(f))

        # Calculate RMSD value between previous and current frame
        rmsd_value = mdt.rmsd(frame, previous_frame, atom_indices=parsed_selection.atom_indices)[0]
        #print(rmsd_value)

        # If the RMSD value is bigger than the cutoff then stop here
        if rmsd_value > rmsd_cutoff:
            print('FAIL: High RMSd (' + str(rmsd_value)  + ') values between frames ' + str(f) + ' and ' + str(f+1))
            return True

        # Update the previous frame as the current one
        previous_frame = frame

    return False