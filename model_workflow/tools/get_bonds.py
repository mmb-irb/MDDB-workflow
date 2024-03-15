from model_workflow.tools.get_pdb_frames import get_pdb_frames
from model_workflow.utils.auxiliar import load_json
from model_workflow.utils.constants import TOPOLOGY_FILENAME, SUPPORTED_ION_ELEMENTS
from model_workflow.utils.vmd_spells import get_covalent_bonds

from typing import List, Optional

import pytraj as pt


# Check if two sets of bonds match perfectly
def do_bonds_match (
    bonds_1 : List[ List[int] ],
    bonds_2 : List[ List[int] ],
    # Atom elements must be passed since we do not evaluate bonds with ions
    atom_elements : List[str],
    # Set verbose as true to show which are the atoms preventing the match
    verbose : bool = False
) -> bool:
    # If the number of atoms in both lists is not matching then there is something very wrong
    if len(bonds_1) != len(bonds_2):
        raise ValueError('The number of atoms is not matching in both bond lists')
    # Find ion atom indices
    ion_atom_indices = set()
    for atom_index, atom_element in enumerate(atom_elements):
        if atom_element in SUPPORTED_ION_ELEMENTS:
            ion_atom_indices.add(atom_index)
    # For each atom, check bonds to match perfectly
    # Order is not important
    for atom_index, (atom_bonds_1, atom_bonds_2) in enumerate(zip(bonds_1, bonds_2)):
        # Skip ion bonds
        if atom_index in ion_atom_indices:
            continue
        atom_bonds_set_1 = set(atom_bonds_1) - ion_atom_indices
        atom_bonds_set_2 = set(atom_bonds_2) - ion_atom_indices
        # Check atom bonds to match
        if len(atom_bonds_set_1) != len(atom_bonds_set_2) or any(bond not in atom_bonds_set_2 for bond in atom_bonds_set_1):
            if verbose:
                print('Missmatch in atom ' + str(atom_index) + ': it is ' + str(atom_bonds_set_1) + ' -> it must be ' + str(atom_bonds_set_2))
                print() # Extra line for the frame logs to not erase the previous log
            return False
    return True

# Get covalent bonds using VMD along different frames
# This way we avoid having false positives because 2 atoms are very close in one frame by accident
# This way we avoid having false negatives because 2 atoms are very far in one frame by accident
def get_most_stable_bonds (
    structure_filename : str,
    trajectory_filename : str,
    snapshots : int,
    frames_limit : int = 10
) -> List[ List[int] ]:

    # Get each frame in pdb format to run VMD
    print('Finding most stable bonds')
    frames, step, count = get_pdb_frames(structure_filename, trajectory_filename, snapshots, frames_limit)

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
    most_stable_bonds = []
    for atom in range(atom_count):
        total_bonds = []
        # Accumulate all bonds
        for frame in frame_bonds:
            atom_bonds = frame[atom]
            total_bonds += atom_bonds
        # Keep only those bonds with more occurrences than half the number of frames
        unique_bonds = set(total_bonds)
        atom_most_stable_bonds = []
        for bond in unique_bonds:
            occurrences = total_bonds.count(bond)
            if occurrences > majority_cut:
                atom_most_stable_bonds.append(bond)
        # Add current atom safe bonds to the total
        most_stable_bonds.append(atom_most_stable_bonds)

    return most_stable_bonds

# Return a canonical frame number where all bonds are exactly as they should
def get_bonds_canonical_frame (
    structure_filename : str,
    trajectory_filename : str,
    snapshots : int,
    reference_bonds : List[ List[int] ],
    atom_elements : List[str],
    patience : int = 100, # Limit of frames to check before we surrender
) -> Optional[int]:

    # Now that we have the reference bonds, we must find a frame where bonds are exactly the canonical ones
    print('Searching reference bonds canonical frame. Only first ' + str(patience) + ' frames will be checked.')
    frames, step, count = get_pdb_frames(structure_filename, trajectory_filename, snapshots)

    # We check all frames but we stop as soon as we find a match
    reference_bonds_frame = None
    for frame_number, frame_pdb in enumerate(frames):
        bonds = get_covalent_bonds(frame_pdb)
        if do_bonds_match(bonds, reference_bonds, atom_elements):
            reference_bonds_frame = frame_number
            break
        # If we didn't find a canonical frame at this point we probablty won't
        if frame_number > patience:
            return None
    # If no frame has the canonical bonds then we return None
    if reference_bonds_frame == None:
        return None

    print(' Got it -> Frame ' + str(reference_bonds_frame + 1))

    return reference_bonds_frame

# Extract bonds from a source file
def mine_topology_bonds (bonds_source_file : 'File') -> list:
    if not bonds_source_file or not bonds_source_file.exists:
        return None
    # If we have the standard topology then get bonds from it
    if bonds_source_file.filename == TOPOLOGY_FILENAME:
        print('Bonds in the "' + bonds_source_file.filename + '" file will be used')
        standard_topology = load_json(bonds_source_file.path)
        bonds = standard_topology.get('atom_bonds', None)
        if bonds:
            return bonds
        print('  There were no bonds in the topology file. Is this an old file?')
    # In some ocasions, bonds may come inside a topology which can be parsed through pytraj
    if bonds_source_file.is_pytraj_supported:
        print('Bonds will be mined from "' + bonds_source_file.path + '"')
        pt_topology = pt.load_topology(filename=bonds_source_file.path)
        atom_bonds = [ [] for i in range(pt_topology.n_atoms) ]
        for bond in pt_topology.bonds:
            a,b = bond.indices
            # Make sure atom indices are regular integers so they are JSON serializables
            atom_bonds[a].append(int(b))
            atom_bonds[b].append(int(a))
        # If there is any bonding data then return bonds
        if any(len(bonds) > 0 for bonds in atom_bonds):
            return atom_bonds
        # If all bonds are empty then it means the parsing failed or the pytraj topology has no bonds
        # We must guess them
        print(' Bonds could not be mined')
    # If we can not mine bonds then return None and they will be guessed further
    return None

# Get safe bonds
# First try to mine bonds from a topology files
# If the mining fails then search for the most stable bonds
# If we turst in stable bonds then simply return the structure bonds
def get_safe_bonds (
    input_topology_file : 'File',
    input_structure_file : 'File',
    input_trajectory_file : 'File',
    must_check_stable_bonds : bool,
    snapshots : int,
    structure : 'Structure'
) -> List[List[int]]:
    # Try to get bonds from the topology before guessing
    safe_bonds = mine_topology_bonds(input_topology_file)
    if safe_bonds:
        return safe_bonds
    # If failed to mine topology bonds then guess stable bonds
    print('Bonds will be guessed by atom distances and radius')
    # Find stable bonds if necessary
    if must_check_stable_bonds:
        # Using the trajectory, find the most stable bonds
        safe_bonds = get_most_stable_bonds(input_structure_file.path, input_trajectory_file.path, snapshots)
        return safe_bonds
    # If we trust stable bonds then simply use structure bonds
    safe_bonds = structure.bonds
    return safe_bonds