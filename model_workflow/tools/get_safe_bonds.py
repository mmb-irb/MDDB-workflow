from model_workflow.tools.get_pdb_frames import get_pdb_frames
from mdtoolbelt.vmd_spells import get_covalent_bonds

from typing import List, Optional


# A few frames should be enough
frames_limit = 10

# Check if two sets of bonds match perfectly
def do_bonds_match (bonds_1 : List[ List[int] ], bonds_2 : List[ List[int] ], verbose : bool = False) -> bool:
    # If the number of atoms in both lists is not matching then there is something very wrong
    if len(bonds_1) != len(bonds_2):
        raise ValueError('The number of atoms is not matching in both bond lists')
    # For each atom, check bonds to match perfectly
    # Order is not important
    for atom_index, (atom_bonds_1, atom_bonds_2) in enumerate(zip(bonds_1, bonds_2)):
        if len(atom_bonds_1) != len(atom_bonds_2):
            if verbose:
                print('Missmatch in atom ' + str(atom_index) + ': ' + str(atom_bonds_1) + ' -> ' + str(atom_bonds_2))
                print() # Extra line for the frame logs to not erase the previous log
            return False
        if any(bond not in atom_bonds_2 for bond in atom_bonds_1):
            if verbose:
                print('Missmatch in atom ' + str(atom_index) + ': ' + str(atom_bonds_1) + ' -> ' + str(atom_bonds_2))
                print() # Extra line for the frame logs to not erase the previous log
            return False
    return True

# Get covalent bonds using VMD along different frames
# This way we avoid having false positives because 2 atoms are very close in one frame by accident
# This way we avoid having false negatives because 2 atoms are very far in one frame by accident
def get_safe_bonds (structure_filename : str, trajectory_filename : str) -> List[ List[int] ]:

    # Get each frame in pdb format to run VMD
    print('Finding safe bonds')
    frames, step, count = get_pdb_frames(structure_filename, trajectory_filename, frames_limit)

    # Track bonds along frames
    frame_bonds = []

    # Iterate over the different frames
    for current_frame_pdb in frames:

        # Find the covalent bonds for the current frame
        current_frame_bonds = get_covalent_bonds(current_frame_pdb)
        frame_bonds.append(current_frame_bonds)

    # Then keep those bonds which are respected in the majority of frames
    # Usually wrongs bonds (both false positives and negatives) are formed only one frame
    # It should not happend that a bond is formed around half of times
    majority_cut = count / 2
    atom_count = len(frame_bonds[0])
    safe_bonds = []
    for atom in range(atom_count):
        total_bonds = []
        # Accumulate all bonds
        for frame in frame_bonds:
            atom_bonds = frame[atom]
            total_bonds += atom_bonds
        # Keep only those bonds with more occurrences than half the number of frames
        unique_bonds = set(total_bonds)
        atom_safe_bonds = []
        for bond in unique_bonds:
            occurrences = total_bonds.count(bond)
            if occurrences > majority_cut:
                atom_safe_bonds.append(bond)
        # Add current atom safe bonds to the total
        safe_bonds.append(atom_safe_bonds)

    return safe_bonds

# Return a canonical frame number where all bonds are exactly as they should
def get_safe_bonds_canonical_frame (
    structure_filename : str,
    trajectory_filename : str,
    safe_bonds : List[ List[int] ],
    patience : int = 100, # Limit of frames to check before we surrender
) -> Optional[int]:

    # Now that we have the safe bonds, we must find a frame where bonds are exactly the canonical ones
    print('Searching safe bonds canonical frame')
    frames, step, count = get_pdb_frames(structure_filename, trajectory_filename)

    # We check all frames but we stop as soon as we find a match
    safe_bonds_frame = None
    for frame_number, frame_pdb in enumerate(frames):
        bonds = get_covalent_bonds(frame_pdb)
        if do_bonds_match(bonds, safe_bonds):
            safe_bonds_frame = frame_number
            break
        # If we didn't find a canonical frame at this point we probablty won't
        if frame_number > patience:
            return None
    # If no frame has the canonical bonds then we return None
    if safe_bonds_frame == None:
        return None

    print(' Got it -> Frame ' + str(safe_bonds_frame + 1))

    return safe_bonds_frame