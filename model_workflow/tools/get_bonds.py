from model_workflow.tools.get_pdb_frames import get_pdb_frames
from model_workflow.utils.auxiliar import load_json, warn, MISSING_TOPOLOGY
from model_workflow.utils.auxiliar import MISSING_BONDS, JSON_SERIALIZABLE_MISSING_BONDS
from model_workflow.utils.constants import STANDARD_TOPOLOGY_FILENAME
from model_workflow.utils.vmd_spells import get_covalent_bonds
from model_workflow.utils.gmx_spells import get_tpr_bonds as get_tpr_bonds_gromacs
from model_workflow.utils.gmx_spells import get_tpr_atom_count
from model_workflow.utils.mda_spells import get_tpr_bonds_mdanalysis
from model_workflow.utils.type_hints import *
import pytraj as pt
from collections import Counter

FAILED_BOND_MINING_EXCEPTION = Exception('Failed to mine bonds')

# Set some atoms which are to be skipped from bonding tests given their "fake" nature
def get_excluded_atoms_selection (structure : 'Structure', pbc_selection : 'Selection') -> 'Selection':
    # Get a selection of ion atoms which are not in PBC
    # These ions are usually "tweaked" to be bonded to another atom although there is no real covalent bond
    # They are not taken in count when testing coherent bonds or looking for the reference frame
    non_pbc_ions_selection = structure.select_ions() - pbc_selection
    # We also exclude coarse grain atoms since their bonds will never be found by a distance/radius guess
    excluded_atoms_selection = non_pbc_ions_selection + structure.select_cg()
    return excluded_atoms_selection

# Check if two sets of bonds match perfectly
def do_bonds_match (
    bonds_1 : list[ list[int] ],
    bonds_2 : list[ list[int] ],
    # A selection of atoms whose bonds are not evaluated
    excluded_atoms_selection : 'Selection',
    # Set verbose as true to show which are the atoms preventing the match
    verbose : bool = False,
    # The rest of inputs are just for logs and debug
    atoms : Optional[ list['Atom'] ] = None,
    counter_list : Optional[ list[int] ] = None
) -> bool:
    # If the number of atoms in both lists is not matching then there is something very wrong
    if len(bonds_1) != len(bonds_2):
        raise ValueError(f'The number of atoms is not matching in both bond lists ({len(bonds_1)} and {len(bonds_2)})')
    # Find ion atom indices
    excluded_atom_indices = set(excluded_atoms_selection.atom_indices)
    # For each atom, check bonds to match perfectly
    # Order is not important
    for atom_index, (atom_bonds_1, atom_bonds_2) in enumerate(zip(bonds_1, bonds_2)):
        # Skip ion bonds
        if atom_index in excluded_atom_indices:
            continue
        atom_bonds_set_1 = set(atom_bonds_1) - excluded_atom_indices
        atom_bonds_set_2 = set(atom_bonds_2) - excluded_atom_indices
        # Check atom bonds to match
        if len(atom_bonds_set_1) != len(atom_bonds_set_2) or any(bond not in atom_bonds_set_2 for bond in atom_bonds_set_1):
            if verbose:
                if atoms:
                    mismatch_atom_label = atoms[atom_index].label
                    print(f' Mismatch in atom {mismatch_atom_label}:')
                    it_is_atom_labels = ', '.join([ atoms[index].label for index in atom_bonds_set_1 ])
                    print(f' It is bonded to atoms {it_is_atom_labels}')
                    it_should_be_atom_labels = ', '.join([ atoms[index].label for index in atom_bonds_set_2 ])
                    print(f' It should be bonded to atoms {it_should_be_atom_labels}')
                else:
                    print(f' Mismatch in atom with index {atom_index}:')
                    it_is_atom_indices = ','.join([ str(index) for index in atom_bonds_set_1 ])
                    print(f' It is bonded to atoms with indices {it_is_atom_indices}')
                    it_should_be_atom_indices = ','.join([ str(index) for index in atom_bonds_set_2 ])
                    print(f' It should be bonded to atoms with indices {it_should_be_atom_indices}')
            # Save for failure analysis
            if counter_list is not None:
                counter_list.append((atom_index,tuple(atom_bonds_set_1),tuple(atom_bonds_set_2)))
            return False
    return True

# Get covalent bonds using VMD along different frames
# This way we avoid having false positives because 2 atoms are very close in one frame by accident
# This way we avoid having false negatives because 2 atoms are very far in one frame by accident
def get_most_stable_bonds (
    structure_filepath : str,
    trajectory_filepath : str,
    snapshots : int,
    frames_limit : int = 10
) -> list[ list[int] ]:

    # Get each frame in pdb format to run VMD
    print('Finding most stable bonds')
    frames, step, count = get_pdb_frames(structure_filepath, trajectory_filepath, 
                                         snapshots, frames_limit, pbar_bool=True)

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


def get_bonds_reference_frame (
    structure_file : 'File',
    trajectory_file : 'File',
    snapshots : int,
    reference_bonds : list[ list[int] ],
    structure : 'Structure',
    pbc_selection : 'Selection',
    patience : int = 100,  # Limit of frames to check before we surrender
    verbose : bool = False,
) -> Optional[int]:
    """Return a reference frame number where all bonds are exactly as they should (by VMD standards).
    This is the frame used when representing the MD."""
    # Set some atoms which are to be skipped from these test given their "fake" nature
    excluded_atoms_selection = get_excluded_atoms_selection(structure, pbc_selection)

    # If all atoms are to be excluded then set the first frame as the reference frame and stop here
    if len(excluded_atoms_selection) == len(structure.atoms): return 0

    # Now that we have the reference bonds, we must find a frame where bonds are exactly the reference ones
    # IMPORTANT: Note that we do not set a frames limit here, so all frames will be read and the step will be 1
    frames, step, count = get_pdb_frames(structure_file.path, trajectory_file.path, snapshots,patience=patience)
    if step != 1: raise ValueError('If we are skipping frames then the code below will silently return a wrong reference frame')
    print(f'Searching the reference frame for the bonds. Only first {min(patience, count)} frames will be checked.')
    # We check all frames but we stop as soon as we find a match
    bonds_reference_frame = None
    counter_list = []
    for frame_number, frame_pdb in enumerate(frames):
        # Get the actual frame number
        bonds = get_covalent_bonds(frame_pdb)
        if do_bonds_match(bonds, reference_bonds, excluded_atoms_selection, counter_list=counter_list, verbose=verbose):
            bonds_reference_frame = frame_number
            break
    frames.close()
    # If no frame has the reference bonds then we return None
    if bonds_reference_frame == None:
        # Print the first clashes table
        print(' First clash stats:')
        headers = ['Count', 'Atom', 'Is bonding with', 'Should bond with']
        count = Counter(counter_list).most_common(10)
        # Calculate column widths
        table_data = []
        for (at, bond, should), n in count:
            table_data.append([n, at, bond, should])
        col_widths = [max(len(str(item)) for item in col) for col in zip(*table_data, headers)]
        # Format rows
        def format_row(row):
            return " | ".join(f"{str(item):>{col_widths[i]}}" for i, item in enumerate(row))
        # Print table
        print(format_row(headers))
        print("-+-".join('-' * width for width in col_widths))
        for row in table_data:
            print(format_row(row))
        return None
    print(f' Got it -> Frame {bonds_reference_frame + 1}')

    return bonds_reference_frame

# Extract bonds from a source file and format them per atom
def mine_topology_bonds (bonds_source_file : 'File' | Exception) -> list[ list[int] ]:
    # If there is no topology then return no bonds at all
    if bonds_source_file == MISSING_TOPOLOGY or not bonds_source_file.exists:
        return None
    print('Mining atom bonds from topology file')
    # If we have the standard topology then get bonds from it
    if bonds_source_file.filename == STANDARD_TOPOLOGY_FILENAME:
        print(f' Bonds in the "{bonds_source_file.filename}" file will be used')
        standard_topology = load_json(bonds_source_file.path)
        standard_atom_bonds = standard_topology.get('atom_bonds', None)
        # Convert missing bonds flags
        # These come from coarse grain (CG) simulations with no topology
        atom_bonds = []
        for bonds in standard_atom_bonds:
            if bonds == JSON_SERIALIZABLE_MISSING_BONDS:
                atom_bonds.append(MISSING_BONDS)
                continue
            atom_bonds.append(bonds)
        if atom_bonds: return atom_bonds
        print('  There were no bonds in the topology file. Is this an old file?')
    # In some ocasions, bonds may come inside a topology which can be parsed through pytraj
    elif bonds_source_file.is_pytraj_supported():
        print(f' Bonds will be mined from "{bonds_source_file.path}"')
        pt_topology = pt.load_topology(filename=bonds_source_file.path)
        raw_bonds = [ bonds.indices for bonds in pt_topology.bonds ]
        # Sort bonds
        atom_count = pt_topology.n_atoms
        atom_bonds = sort_bonds(raw_bonds, atom_count)
        # If there is any bonding data then return atom bonds
        if any(len(bonds) > 0 for bonds in atom_bonds): return atom_bonds
    # If we have a TPR then use our own tool
    elif bonds_source_file.format == 'tpr':
        print(f' Bonds will be mined from TPR file "{bonds_source_file.path}"')
        raw_bonds = get_tpr_bonds(bonds_source_file.path)
        if raw_bonds != FAILED_BOND_MINING_EXCEPTION:
            # Sort bonds
            atom_count = get_tpr_atom_count(bonds_source_file.path)
            atom_bonds = sort_bonds(raw_bonds, atom_count)
            return atom_bonds
    # If we failed to mine bonds then return None and they will be guessed further
    print (' Failed to mine bonds -> They will be guessed from atom distances and radius')
    return None

# Get TPR bonds
# Try 2 different methods and hope 1 of them works
def get_tpr_bonds (tpr_filepath : str) -> list[ tuple[int, int] ]:
    try:
        bonds = get_tpr_bonds_gromacs(tpr_filepath)
    except:
        print(' Our tool failed to extract bonds. Using MDAnalysis extraction...')
        try:
            bonds = get_tpr_bonds_mdanalysis(tpr_filepath)
            # This happens when the filter selection is what we want to exclude, and not to keep
            if len(bonds) == 0:
                warn('Bonds were mined successfully but it looks like there are no bonds. Is it all ions or CG solvent?')
        except Exception as err:
            print(f' MDAnalysis failed to extract bonds: {err}. Relying on guess...')
            bonds = FAILED_BOND_MINING_EXCEPTION
    return bonds

# Sort bonds according to our format: a list with the bonded atom indices for each atom
# Source data is the usual format to store bonds: a list of tuples with every pair of bonded atoms
def sort_bonds (source_bonds : list[ tuple[int, int] ], atom_count : int) -> list[ list[int] ]:
    # Set a list of lists with an empty list for every atom
    atom_bonds = [ [] for i in range(atom_count) ]
    for bond in source_bonds:
        a,b = bond
        # Make sure atom indices are regular integers so they are JSON serializables
        atom_bonds[a].append(int(b))
        atom_bonds[b].append(int(a))
    return atom_bonds

# Get safe bonds
# First try to mine bonds from a topology files
# If the mining fails then search for the most stable bonds
# If we trust in stable bonds then simply return the structure bonds
def find_safe_bonds (
    topology_file : 'File' | Exception,
    structure_file : 'File',
    trajectory_file : 'File',
    must_check_stable_bonds : bool,
    snapshots : int,
    structure : 'Structure',
    guess_bonds : bool = False,
    # Optional file with bonds sorted according a new atom order
    resorted_bonds_file : Optional['File'] = None
) -> list[list[int]]:
    """Find reference safe bonds in the system."""
    # If we have a resorted file then use it
    # Note that this is very excepcional
    if resorted_bonds_file != None and resorted_bonds_file.exists:
        warn('Using resorted safe bonds')
        return load_json(resorted_bonds_file.path)
    # Get a selection including coarse grain atoms in the structure
    cg_selection = structure.select_cg()
    # Try to get bonds from the topology before guessing
    if guess_bonds:
        print('Bonds were forced to be guessed, so bonds in the topology file will be ignored')
    else:
        safe_bonds = mine_topology_bonds(topology_file)
        if safe_bonds:
            return safe_bonds
    # If all bonds are in coarse grain then set all bonds "wrong" already
    if len(cg_selection) == structure.atom_count:
        safe_bonds = [ MISSING_BONDS for atom in range(structure.atom_count) ]
        return safe_bonds
    # If failed to mine topology bonds then guess stable bonds
    print('Bonds will be guessed by atom distances and radii along different frames in the trajectory')
    # Find stable bonds if necessary
    if must_check_stable_bonds:
        # Using the trajectory, find the most stable bonds
        print('Checking bonds along trajectory to determine which are stable')
        safe_bonds = get_most_stable_bonds(structure_file.path, trajectory_file.path, snapshots)
        discard_coarse_grain_bonds(safe_bonds, cg_selection)
        return safe_bonds
    # If we trust stable bonds then simply use structure bonds
    print('Default structure bonds will be used since they have been marked as trusted')
    safe_bonds = structure.bonds
    discard_coarse_grain_bonds(safe_bonds, cg_selection)
    return safe_bonds

# Given a list of bonds, discard the ones in the coarse grain selection
# Note that the input list will be mutated
def discard_coarse_grain_bonds (bonds : list, cg_selection : 'Selection'):
    # For every atom in CG, replace its bonds with a class which will raise and error when read
    # Thus we make sure using these wrong bonds anywhere further will result in failure
    for atom_index in cg_selection.atom_indices:
        bonds[atom_index] = MISSING_BONDS