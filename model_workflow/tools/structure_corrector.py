from os import remove
from os.path import exists
from json import load

# Import local tools
from model_workflow.tools.get_bonds import find_safe_bonds, do_bonds_match, get_bonds_canonical_frame
from model_workflow.tools.get_pdb_frames import get_pdb_frame
from model_workflow.utils.auxiliar import TestFailure, get_new_letter, save_json, warn
from model_workflow.utils.constants import CORRECT_ELEMENTS, STABLE_BONDS_FLAG, COHERENT_BONDS_FLAG
from model_workflow.utils.structures import Structure
from model_workflow.utils.type_hints import *

# Analyze the structure looking for irregularities and then modify the structure to standarize the format
#
# Supported cases:
#
# * Missing/Non-standard elements -> Atom elements are guessed when missing and standarized (e.g. ZN -> Zn)
# * Unstable atom bonds: A bond is formed / broken because its 2 atoms are too close / far in the structure
#   -> Coordinates in the pdb file are replaced by those of the first frame in the trajectory where bonds are stable
# * Incoherent atom bonds -> If an atom has unexpected bonds then the simulation may be wrong
#   -> Stop the workflow and warn the user
# * Missing chains -> Chains are added through VMD
# * Splitted chain -> Kill the process and warn the user. This should never happen by just reading a pdb, but after correcting
# * Repeated chains -> Chains are renamed (e.g. A, G, B, G, C, G -> A, G, B, H, C, I)
# * Splitted residues -> Atoms are sorted together by residue. Trajectory coordinates are also sorted
# * Repeated residues -> Residues are renumerated (e.g. 1, 2, 3, 1, 2, 1, 2 -> 1, 2, 3, 4, 5, 6, 7)
# * Repeated atoms -> Atoms are renamed with their numeration (e.g. C, C, C, O, O -> C1, C2, C3, O1, O2)

# Note that the 'mercy' flag may be passed for critical checkings to not kill the process on fail
# Note that the 'trust' flag may be passed for critical checkings to skip them

# This function also sets some values in the passed MD

def structure_corrector (
    input_structure_file : 'File',
    input_trajectory_file : Optional['File'],
    input_topology_file : Optional['File'],
    output_structure_file : 'File',
    output_trajectory_file : Optional['File'],
    MD : 'MD'
) -> dict:

    # Extract some MD features
    snapshots = MD._snapshots
    register = MD.register
    mercy = MD.project.mercy
    trust = MD.project.trust

    # Import the pdb file and parse it to a structure object
    structure = Structure.from_pdb_file(input_structure_file.path)

    # Write the inital output structure file which will be overwritten several times further
    structure.generate_pdb_file(output_structure_file.path)

    # ------------------------------------------------------------------------------------------
    # Missing/Non-standard elements ------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # It is important to fix elements before trying to fix bonds, since elements have an impact on bonds
    # VMD logic to find bonds relies in the atom element to set the covalent bond distance cutoff
    if structure.fix_atom_elements(trust=(CORRECT_ELEMENTS in trust)):
        # Update the structure file using the corrected structure
        print(' The structure file has been modified -> ' + output_structure_file.filename)
        structure.generate_pdb_file(output_structure_file.path)

    # ------------------------------------------------------------------------------------------
    # Unstable atom bonds ----------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # Set if stable bonds have to be checked
    # Note that we must not skip this even if the test already passed
    # It may be a corrected structure the one which passed the structure, while this structure comes from the raw input
    must_check_stable_bonds = STABLE_BONDS_FLAG not in trust

    # Get safe bonds
    # Use topology bonds if possible
    # Otherwise guess bonds by guessing bonds according to coordinates and atom radius for 10 frames along the trajectory
    safe_bonds = find_safe_bonds(
        input_topology_file,
        output_structure_file,
        input_trajectory_file,
        must_check_stable_bonds,
        snapshots,
        structure
    )
    # If safe bonds do not match structure bonds then we have to fix it
    safe_bonds_frame = None
    atom_elements = [ atom.element for atom in structure.atoms ]
    def check_stable_bonds ():
        # If we have been requested to skip this test then we are done
        if not must_check_stable_bonds:
            # Set the safe bonds as the structure bonds, just in case
            structure.bonds = safe_bonds
            return
        # Reset warnings related to this analysis
        register.remove_warnings(STABLE_BONDS_FLAG)
        # If bonds match from the begining we are done as well
        print('Checking default structure bonds')
        if do_bonds_match(structure.bonds, safe_bonds, atom_elements, verbose=True, atoms=structure.atoms):
            register.update_test(STABLE_BONDS_FLAG, True)
            print(' They are good')
            return
        print(' They are wrong')
        # Set the safe bonds as the structure bonds
        structure.bonds = safe_bonds
        # Find the first frame in the whole trajectory where safe bonds are respected
        safe_bonds_frame = get_bonds_canonical_frame(
            structure_filepath = output_structure_file.path,
            trajectory_filepath = input_trajectory_file.path,
            snapshots = snapshots,
            reference_bonds = safe_bonds,
            atom_elements = atom_elements
        )
        # If there is no canonical frame then stop here since there must be a problem
        if safe_bonds_frame == None:
            print('There is no canonical frame for safe bonds. Is the trajectory not imaged?')
            must_be_killed = STABLE_BONDS_FLAG not in mercy
            if must_be_killed:
                raise TestFailure('Failed to find stable bonds')
            register.update_test(STABLE_BONDS_FLAG, False)
            register.add_warning(STABLE_BONDS_FLAG, ('Could not find a frame in the trajectory respecting all bonds if bonds were guessed according to atom coordinates and radius.\n'
                'The main PDB structure is a default structure and it would be considered to have wrong bonds if they were predicted as previously stated.'))
            return
        # If we found a canonical frame then we are good
        # Save this frame as the reference frame for the current MD
        MD._reference_frame = safe_bonds_frame
        # Set also the safe bonds frame structure to mine its coordinates
        safe_bonds_frame_filename = get_pdb_frame(output_structure_file.path, input_trajectory_file.path, safe_bonds_frame)
        safe_bonds_frame_structure = Structure.from_pdb_file(safe_bonds_frame_filename)
        # Set all coordinates in the main structure by copying the safe bonds frame coordinates
        for atom_1, atom_2 in zip(structure.atoms, safe_bonds_frame_structure.atoms):
            atom_1.coords = atom_2.coords
        # Remove the safe bonds frame since it is not required anymore
        remove(safe_bonds_frame_filename)
        # Set the modified variable as true since we have changes the structure
        # Update the structure file using the corrected structure
        print(' The structure file has been modified -> ' + output_structure_file.filename)
        structure.generate_pdb_file(output_structure_file.path)
        # Set the test as passed
        register.update_test(STABLE_BONDS_FLAG, True)
    check_stable_bonds()

    # Write safe bonds back to the MD
    MD.project._safe_bonds = safe_bonds
    MD.project._safe_bonds_frame = safe_bonds_frame

    # ------------------------------------------------------------------------------------------
    # Incoherent atom bonds ---------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # Set if coherent bonds have to be checked
    must_check_coherent_bonds = COHERENT_BONDS_FLAG not in trust

    # If this analysis has been already passed then we skip the process
    if register.tests.get(COHERENT_BONDS_FLAG, None) == True:
        must_check_stable_bonds = False

    # Run the coherent bonds analysis if necessary
    if must_check_coherent_bonds:
        # Reset warnings related to this analysis
        register.remove_warnings(COHERENT_BONDS_FLAG)
        # If the test is not passed then report it
        if structure.check_incoherent_bonds():
            print('FAIL: Uncoherent bonds were found')
            must_be_killed = COHERENT_BONDS_FLAG not in mercy
            if must_be_killed:
                raise TestFailure('Failed to find coherent bonds')
            register.update_test(COHERENT_BONDS_FLAG, False)
            register.add_warning(COHERENT_BONDS_FLAG, 'Some atoms may have a higher or lower number of bonds than they should according to their element.')
        # Set the test as succeed if all was good
        else:
            register.update_test(COHERENT_BONDS_FLAG, True)

    # ------------------------------------------------------------------------------------------
    # Missing chains ---------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # Check if chains are missing. If so, create a new chainned structure and set it as the reference
    chains = structure.chains
    
    # In case there are not chains at all
    if len(chains) == 1 and ( chains[0].name == ' ' or chains[0].name == 'X' ):
        print('WARNING: chains are missing and they will be added')

        # Run the chainer
        #structure.chainer()
        structure.auto_chainer()
        # Update the structure file using the corrected structure
        print(' The structure file has been modified -> ' + output_structure_file.filename)
        structure.generate_pdb_file(output_structure_file.path)

    else:
        # In case there are some missing chains
        # Note that atoms with no chain are not a problem for the workflow but they are for the web client
        unlettered_chain = next((chain for chain in chains if chain.name == ' '), None)
        if unlettered_chain:
            current_letters = set([ chain.name for chain in chains ])
            new_letter = get_new_letter(current_letters)
            # If we run out of letters there may be some problematic chain configuration
            # In this cases we cannot respect the original chains
            if new_letter == None:
                warn('No more letters in the alphabel to fill missing chains -> All chains will be assigned from scratch')
                structure.auto_chainer()
            else:
                warn(f'Some chains are missing -> Unchained regions will be chained as {new_letter}')
                unlettered_chain.name = new_letter
            # Update the structure file using the corrected structure
            print(f' The structure file has been modified -> {output_structure_file.filename}')
            structure.generate_pdb_file(output_structure_file.path)

    # ------------------------------------------------------------------------------------------
    # Splitted chains --------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # Check if chains are splitted. If so, stop here and warn the user
    checked_chains = []
    last_chain = None
    for atom in structure.atoms:
        chain_name = atom.chain.name
        if chain_name == last_chain:
            continue
        last_chain = chain_name
        if chain_name in checked_chains:
            # Note that this will never happen because of the pdb itself, duplicated chains are handled further
            # This will only happen if chains were missing and guessed, and there is something very wrong in the structure
            # Check fragments with the VMD and searh for wrong bonds
            raise TestFailure(f'We are having splitted chains (e.g. {chain_name})')
        checked_chains.append(chain_name)

    # ------------------------------------------------------------------------------------------
    # Repeated chains --------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    if structure.check_repeated_chains(fix_chains=True, display_summary=True):
        # Update the structure file using the corrected structure
        print(' The structure file has been modified -> ' + output_structure_file.filename)
        structure.generate_pdb_file(output_structure_file.path)

    # ------------------------------------------------------------------------------------------
    # Splitted residues ------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    if structure.check_repeated_residues(fix_residues=True, display_summary=True):
        # Update the structure file using the corrected structure
        print(' The structure file has been modified -> ' + output_structure_file.filename)
        structure.generate_pdb_file(output_structure_file.path)

        # Sort trajectory coordinates in case atoms were sorted
        if input_trajectory_file.path and structure.trajectory_atom_sorter:
            # Save a warning in the register
            print('WARNING: Atoms have been sorted to solve splitted residues')
            register.update_cache('resorted_atoms', True)
            print('Creating resorted files for atom bonds and charges')
            # Bonds are already resorted
            save_json(safe_bonds, MD.project.resorted_bonds_file.path, indent=4)
            # Charges are to be resorted
            resorted_charges = [ MD.project._charges[index] for index in structure.new_atom_order ]
            MD.project._charges = resorted_charges
            save_json(resorted_charges, MD.project.resorted_charges_file.path, indent=4)
            print('Sorting trajectory coordinates to fit the new structure atom sort...')
            structure.trajectory_atom_sorter(output_structure_file, input_trajectory_file, output_trajectory_file)

    # ------------------------------------------------------------------------------------------
    # Repeated atoms ---------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    if structure.check_repeated_atoms(fix_atoms=True, display_summary=True):
        # Update the structure file using the corrected structure
        print(' The structure file has been modified -> ' + output_structure_file.filename)
        structure.generate_pdb_file(output_structure_file.path)