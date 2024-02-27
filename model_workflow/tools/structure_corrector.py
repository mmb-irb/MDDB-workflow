from typing import Optional, List
from os import remove
from os.path import exists
from json import load

# Import local tools
from model_workflow.utils.structures import Structure
from model_workflow.utils.auxiliar import TestFailure
from model_workflow.tools.get_bonds import get_safe_bonds, do_bonds_match, get_bonds_canonical_frame
from model_workflow.tools.get_pdb_frames import get_pdb_frame
from model_workflow.utils.auxiliar import get_new_letter, save_json
from model_workflow.utils.constants import CORRECT_ELEMENTS, STABLE_BONDS_FLAG, COHERENT_BONDS_FLAG

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

    # Track if there has been any modification and then the structure must be rewritten
    modified = False

    # Import the pdb file and parse it to a structure object
    structure = Structure.from_pdb_file(input_structure_file.path)

    # ------------------------------------------------------------------------------------------
    # Missing/Non-standard elements ------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # It is important to fix elements before trying to fix bonds, since elements have an impact on bonds
    # VMD logic to find bonds relies in the atom element to set the covalent bond distance cutoff
    if structure.fix_atom_elements(trust=(CORRECT_ELEMENTS in trust)):
        modified = True

    # ------------------------------------------------------------------------------------------
    # Unstable atom bonds ----------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # Set if stable bonds have to be checked
    must_check_stable_bonds = STABLE_BONDS_FLAG not in trust

    # If this analysis has been already passed then we skip the process
    if register.tests.get(STABLE_BONDS_FLAG, None) == True:
        must_check_stable_bonds = False

    # Get safe bonds
    # Use topology bonds if possible
    # Otherwise guess bonds by guessing bonds according to coordinates and atom radius for 10 frames along the trajectory
    safe_bonds = get_safe_bonds(
        input_topology_file,
        input_structure_file,
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
            return
        # Reset warnings related to this analysis
        register.remove_warnings(STABLE_BONDS_FLAG)
        # If bonds match from the begining we are done as well
        if do_bonds_match(structure.bonds, safe_bonds, atom_elements):
            register.tests[STABLE_BONDS_FLAG] = True
            return
        print('WARNING: Default structure has wrong bonds')
        # Set the safe bonds as the structure bonds
        structure.bonds = safe_bonds
        # Find the first frame in the whole trajectory where safe bonds are respected
        safe_bonds_frame = get_bonds_canonical_frame(
            structure_filename = input_structure_file.path,
            trajectory_filename = input_trajectory_file.path,
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
            register.tests[STABLE_BONDS_FLAG] = False
            register.add_warning(STABLE_BONDS_FLAG, ('Could not find a frame in the trajectory respecting all bonds if bonds were guessed according to atom coordinates and radius.\n'
                'The main PDB structure is a default structure and it would be considered to have wrong bonds if they were predicted as previously stated.'))
            return
        # If we found a canonical frame then we are good
        # Set also the safe bonds frame structure to mine its coordinates
        safe_bonds_frame_filename = get_pdb_frame(input_structure_file.path, input_trajectory_file.path, safe_bonds_frame)
        safe_bonds_frame_structure = Structure.from_pdb_file(safe_bonds_frame_filename)
        # Set all coordinates in the main structure by copying the safe bonds frame coordinates
        for atom_1, atom_2 in zip(structure.atoms, safe_bonds_frame_structure.atoms):
            atom_1.coords = atom_2.coords
        # Remove the safe bonds frame since it is not required anymore
        remove(safe_bonds_frame_filename)
        # Tag the test as succeed if we did not skip it
        register.tests[STABLE_BONDS_FLAG] = True
        # Set the modified variable as true since we have changes the structure
        modified = True
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
            register.tests[COHERENT_BONDS_FLAG] = False
            register.add_warning(COHERENT_BONDS_FLAG, 'Some atoms may have a higher or lower number of bonds than they should according to their element.')
        # Tag the test as succeed if all was good
        else:
            register.tests[COHERENT_BONDS_FLAG] = True
        

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
        modified = True

    else:
        # In case there are some missing chains
        # Note that atoms with no chain are not a problem for the workflow but they are for the web client
        unlettered_chain = next((chain for chain in chains if chain.name == ' '), None)
        if unlettered_chain:
            current_letters = [ chain.name for chain in chains ]
            new_letter = get_new_letter(current_letters)
            unlettered_chain.name = new_letter
            modified = True
            print('WARNING: Some missing chains -> Now chain ' + new_letter)

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
            raise SystemExit('We are having splitted chains (e.g. ' + chain_name + ')')
        checked_chains.append(chain_name)

    # ------------------------------------------------------------------------------------------
    # Repeated chains --------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    if structure.check_repeated_chains(fix_chains=True, display_summary=True):
        modified = True

    # ------------------------------------------------------------------------------------------
    # Splitted residues ------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    if structure.check_repeated_residues(fix_residues=True, display_summary=True):
        modified = True

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
            structure.trajectory_atom_sorter(input_structure_file, input_trajectory_file, output_trajectory_file)

    # ------------------------------------------------------------------------------------------
    # Repeated atoms ---------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    if structure.check_repeated_atoms(fix_atoms=True, display_summary=True):
        modified = True

    # ------------------------------------------------------------------------------------------
    # Final output -----------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # Write a new structure if any modification was done
    if modified:
        print(' The structure file has been modified -> ' + output_structure_file.filename)
    else:
        print(' Everything is fine')
    # Generate the file anyway so this new structure is used and not reclaulcated
    structure.generate_pdb_file(output_structure_file.path)