from os import remove
from mddb_workflow.tools.get_bonds import find_safe_bonds, do_bonds_match, get_bonds_reference_frame
from mddb_workflow.tools.get_bonds import get_excluded_atoms_selection
from mddb_workflow.tools.get_pdb_frames import get_pdb_frame
from mddb_workflow.tools.get_charges import get_charges
from mddb_workflow.utils.auxiliar import InputError, TestFailure, MISSING_BONDS
from mddb_workflow.utils.auxiliar import get_new_letter, save_json, warn
from mddb_workflow.utils.constants import CORRECT_ELEMENTS, STABLE_BONDS_FLAG, COHERENT_BONDS_FLAG
from mddb_workflow.utils.structures import Structure
from mddb_workflow.utils.type_hints import *


def structure_corrector(
    # Note that this is an early provisional structure
    structure: 'Structure',
    input_trajectory_file: Optional['File'],
    input_topology_file: Union['File', Exception],
    output_structure_file: 'File',
    output_trajectory_file: Optional['File'],
    MD: 'MD',
    # Note that this is an early provisional atom selection
    pbc_selection: 'Selection',
    snapshots: int,
    register: 'Register',
    mercy: list[str],
    trust: list[str],
    guess_bonds: bool,
    ignore_bonds: bool,
) -> dict:
    """Analyze the structure looking for irregularities and then modify the structure to standarize the format.

    Supported cases:

    * Missing/Non-standard elements -> Atom elements are guessed when missing and standarized (e.g. ZN -> Zn)
    * Unstable atom bonds: A bond is formed / broken because its 2 atoms are too close / far in the structure
      -> Coordinates in the pdb file are replaced by those of the first frame in the trajectory where bonds are stable
    * Incoherent atom bonds -> If an atom has unexpected bonds then the simulation may be wrong
      -> Stop the workflow and warn the user
    * Missing chains -> Chains are added through VMD
    * Splitted chain -> Kill the process and warn the user. This should never happen by just reading a pdb, but after correcting
    * Repeated chains -> Chains are renamed (e.g. A, G, B, G, C, G -> A, G, B, H, C, I)
    * Splitted residues -> Atoms are sorted together by residue. Trajectory coordinates are also sorted
    * Repeated residues -> Residues are renumerated (e.g. 1, 2, 3, 1, 2, 1, 2 -> 1, 2, 3, 4, 5, 6, 7)
    * Repeated atoms -> Atoms are renamed with their numeration (e.g. C, C, C, O, O -> C1, C2, C3, O1, O2)

    Note that the 'mercy' flag may be passed for critical checkings to not kill the process on fail.
    Note that the 'trust' flag may be passed for critical checkings to skip them.

    This function also sets some values in the passed MD.

    """
    # Write the inital output structure file which will be overwritten several times further
    print(' The structure file has been copied -> ' + output_structure_file.filename)
    structure.generate_pdb_file(output_structure_file.path)

    # ------------------------------------------------------------------------------------------
    # Missing/Non-standard elements ------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # It is important to fix elements before trying to fix bonds, since elements have an impact on bonds
    # VMD logic to find bonds relies in the atom element to set the covalent bond distance cutoff
    if structure.fix_atom_elements(trust=(CORRECT_ELEMENTS in trust)):
        # Update the structure file using the corrected structure
        print(' The structure file has been modified (fix atom elements) -> ' + output_structure_file.filename)
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
        structure,
        register,
        guess_bonds,
        ignore_bonds
    )

    # Check if the structure is missing bonds
    missing_any_bonds = any(bond == MISSING_BONDS for bond in safe_bonds)
    if missing_any_bonds: must_check_stable_bonds = False

    # If safe bonds do not match structure bonds then we have to fix it
    def check_stable_bonds():
        # Save the current structure bonds to further compare with the safe bonds
        current_bonds = structure.bonds
        # Now force safe bonds in the structure even if they "match" already
        # Note that some bonds are excluded from the check and they may be wrong in the structure
        # e.g. coarse grain atom bonds
        structure.bonds = safe_bonds
        # If we have been requested to skip this test then we are done
        if not must_check_stable_bonds: return
        # Reset warnings related to this analysis
        register.remove_warnings(STABLE_BONDS_FLAG)
        # Set some atoms which are to be skipped from these test given their "fake" nature
        excluded_atoms_selection = get_excluded_atoms_selection(structure, pbc_selection)
        # If bonds match from the begining we are done as well
        print(f'Checking default structure bonds ({STABLE_BONDS_FLAG})')
        if do_bonds_match(current_bonds, safe_bonds, excluded_atoms_selection, verbose=True, atoms=structure.atoms):
            register.update_test(STABLE_BONDS_FLAG, True)
            print(' They are good')
            return
        print(' They are wrong')
        # Find the first frame in the whole trajectory where safe bonds are respected
        bonds_reference_frame = get_bonds_reference_frame(
            structure_file=output_structure_file,
            trajectory_file=input_trajectory_file,
            snapshots=snapshots,
            reference_bonds=safe_bonds,
            structure=structure,
            pbc_selection=pbc_selection
        )
        # Update the task output so it does not have to be repeated further
        # IMPORTANT: Note that this is not always run, but only when default structure bonds are wrong
        # IMPORTANT: Thus it is NORMAL when the reference frame is calculated further in the workflow
        MD.get_reference_frame.prefill(MD, bonds_reference_frame, {
            'structure_file': output_structure_file,
            'trajectory_file': input_trajectory_file,
            'snapshots': snapshots,
            'reference_bonds': safe_bonds,
            'structure': structure,
            'pbc_selection': pbc_selection
        })
        # If there is no reference frame then stop here since there must be a problem
        if bonds_reference_frame == None:
            print('There is no reference frame for the safe bonds. Is the trajectory not imaged?')
            must_be_killed = STABLE_BONDS_FLAG not in mercy
            if must_be_killed:
                raise TestFailure('Failed to find stable bonds')
            register.update_test(STABLE_BONDS_FLAG, False)
            register.add_warning(STABLE_BONDS_FLAG, ('Could not find a frame in the trajectory respecting all bonds if bonds were guessed according to atom coordinates and radius.\n'
                'The main PDB structure is a default structure and it would be considered to have wrong bonds if they were predicted as previously stated.'))
            return
        # Set also the safe bonds frame structure to mine its coordinates
        safe_bonds_frame_filename = get_pdb_frame(output_structure_file.path, input_trajectory_file.path, bonds_reference_frame)
        safe_bonds_frame_structure = Structure.from_pdb_file(safe_bonds_frame_filename)
        # Set all coordinates in the main structure by copying the safe bonds frame coordinates
        for atom_1, atom_2 in zip(structure.atoms, safe_bonds_frame_structure.atoms):
            atom_1.coords = atom_2.coords
        # Remove the safe bonds frame since it is not required anymore
        remove(safe_bonds_frame_filename)
        # Set the modified variable as true since we have changes the structure
        # Update the structure file using the corrected structure
        print(' The structure file has been modified (stable bonds) -> ' + output_structure_file.filename)
        structure.generate_pdb_file(output_structure_file.path)
        # Set the test as passed
        register.update_test(STABLE_BONDS_FLAG, True)
    check_stable_bonds()

    # Write safe bonds back to the MD
    MD.project.get_reference_bonds.prefill(MD.project, safe_bonds, {
        'topology_file': input_topology_file,
        'structure_file': output_structure_file,
        'trajectory_file': input_trajectory_file,
        'must_check_stable_bonds': must_check_stable_bonds,
        'snapshots': snapshots,
        'structure': structure,
        'register': register,
        'guess_bonds': guess_bonds,
        'ignore_bonds': ignore_bonds,
    })

    # ------------------------------------------------------------------------------------------
    # Incoherent residue bonds ---------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # Make sure there are no disconnected groups of atoms in every residue
    for residue in structure.residues:
        # If the residue is missing bonds then we can not check if it is coherent
        # We assume it is and continue
        if residue.is_missing_any_bonds(): continue
        # If the residue is coherent then continue
        if residue.is_coherent(): continue
        # Otherwise we report the problem
        residue_selection = residue.get_selection()
        fragments = list(structure.find_fragments(residue_selection))
        if len(fragments) == 0: raise RuntimeError('Do we have an empty residue?')
        if len(fragments) == 1: raise RuntimeError('Test failed but should not')
        warn(f'Multiple fragments in residue {residue.index}: {residue}')
        for f, fragment in enumerate(fragments, 1):
            fragment_atoms = [structure.atoms[index] for index in fragment.atom_indices]
            atom_names = [atom.label for atom in fragment_atoms]
            print(f' Fragment {f}: {", ".join(atom_names)}')
        raise TestFailure(f'Residue {residue.index}: {residue} is not coherent: some atoms are disconnected')

    # ------------------------------------------------------------------------------------------
    # Incoherent atom bonds ---------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # Set if coherent bonds have to be checked
    must_check_coherent_bonds = COHERENT_BONDS_FLAG not in trust
    if missing_any_bonds: must_check_coherent_bonds = False

    # If this analysis has been already passed then we skip the process
    if register.tests.get(COHERENT_BONDS_FLAG, None) is True:
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
        warn('Chains are missing and they will be added')
        # Stop here if we have bonds guessed from coarse grain (i.e. we have no topology)
        # Note that we rely in fragments (and thus in bonds) to guess chains
        if structure.is_missing_any_bonds():
            raise InputError('We cannot guess chains with bonds guessed from coarse grain.\n'
                ' Please either provide a topology including bonds or set chains in the structure PDB file.')
        # Run the chainer
        structure.auto_chainer()
        # Update the structure file using the corrected structure
        print(' The structure file has been modified (no chains) -> ' + output_structure_file.filename)
        structure.generate_pdb_file(output_structure_file.path)

    else:
        # In case there are some missing chains
        # Note that atoms with no chain are not a problem for the workflow but they are for the web client
        unlettered_chain = next((chain for chain in chains if chain.name == ' '), None)
        if unlettered_chain:
            current_letters = set([chain.name for chain in chains])
            new_letter = get_new_letter(current_letters)
            # If we run out of letters there may be some problematic chain configuration
            # In this cases we cannot respect the original chains
            if new_letter is None:
                warn('No more letters in the alphabet to fill missing chains -> All chains will be assigned from scratch')
                # Stop here if we have bonds guessed from coarse grain (i.e. we have no topology)
                # Note that we rely in fragments (and thus in bonds) to guess chains
                if structure.is_missing_any_bonds():
                    raise InputError('We cannot guess chains with bonds guessed from coarse grain.\n'
                        ' Please either provide a topology including bonds or set chains in the structure PDB file.')
                structure.auto_chainer()
            else:
                warn(f'Some chains are missing -> Unchained regions will be chained as {new_letter}')
                unlettered_chain.name = new_letter
            # Update the structure file using the corrected structure
            print(f' The structure file has been modified (missing chains) -> {output_structure_file.filename}')
            structure.generate_pdb_file(output_structure_file.path)

    # ------------------------------------------------------------------------------------------
    # Splitted chains --------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # Note that this will never happen because of the pdb itself, duplicated chains are handled further
    # This will only happen if chains were missing and guessed
    # This may mean there is something wrong in the structure
    # Check fragments with the VMD and searh for wrong bonds
    if structure.check_splitted_chains(fix_chains=True, display_summary=True):
        # Update the structure file using the corrected structure
        print(' The structure file has been modified (splitted chains) -> ' + output_structure_file.filename)
        structure.generate_pdb_file(output_structure_file.path)

    # ------------------------------------------------------------------------------------------
    # Repeated chains --------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    if structure.check_repeated_chains(fix_chains=True, display_summary=True):
        # Update the structure file using the corrected structure
        print(' The structure file has been modified (repeated chains) -> ' + output_structure_file.filename)
        structure.generate_pdb_file(output_structure_file.path)

    # ------------------------------------------------------------------------------------------
    # Coherent chains --------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # Makes sure polymers with sequence are all bonded
    # In other words, make sure there are not multiple proteins/nucleics labeled with the same chain
    # Note that their sequences otherwise will be merged silently, thus not reflecting reality
    # This would make protein mapping impossible for instance
    sequence_polymers_selection = structure.select_protein() + structure.select_nucleic()

    # We also want to support having a solution of free aminoacids/nucleotides
    # So if residues in the chain are not bonded or bonded in very small fragments we accept it
    N_RESIDUES_CUTOFF = 3

    # Save if we splitted chains
    had_to_split_chains = False

    # Iterate sequence polymer chains
    for chain in structure.get_selection_chains(sequence_polymers_selection):
        # If the chain is missing bonds then we can not check if it is coherent
        # We assume it is and continue
        if chain.is_missing_any_bonds(): continue
        # If there is only one fragment we are good
        chain_selection = chain.get_selection()
        fragments = list(structure.find_fragments(chain_selection, atom_bonds=safe_bonds))
        if len(fragments) == 0: raise ValueError(f'No fragments found in chain {chain.name}')
        if len(fragments) == 1: continue
        # We also want to support having a solution of free aminoacids/nucleotides
        # So if residues in the chain are not bonded or bonded in very small fragments we accept it
        # Make a single fragments out of all small fragments
        independent_fragments = []
        merged_fragment = None
        for fragment in fragments:
            # If the fragments reaches the size cutoff then send it to the independent fragments list
            residues = structure.get_selection_residues(fragment)
            if len(residues) > N_RESIDUES_CUTOFF:
                independent_fragments.append(fragment)
                continue
            # Otherwise merge it with other small fragments
            if merged_fragment: merged_fragment += fragment
            else: merged_fragment = fragment
        # Add the merged fragment to the list as well
        if merged_fragment:
            independent_fragments.append(merged_fragment)
        # If we have only one fragment after the merging then we are good
        if len(independent_fragments) == 1: continue
        # Otherwise rename fragments
        warn(f'Chain {chain.name} has sequence polymer(s) but not all atoms are bonded.')
        # We can skip the first fragment, so it can keep the original chain name
        for independent_fragment in independent_fragments[1:]:
            print(f' A fragment of this chain including {len(independent_fragment)} atoms will be splitted appart.')
            had_to_split_chains = True
            # Make a new chain for the current fragment
            new_chain_name = structure.get_next_available_chain_name(chain.name)
            print(f' This fragment will be renamed as chain "{new_chain_name}"')
            structure.set_selection_chain_name(independent_fragment, new_chain_name)

    # If we did any change then save the structure
    if had_to_split_chains:
        print(' The structure file has been modified (coherent chains) -> ' + output_structure_file.filename)
        structure.generate_pdb_file(output_structure_file.path)

    # ------------------------------------------------------------------------------------------
    # Merged residues ------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # There may be residues which contain unconnected (unbonded) atoms. They are not allowed.
    # They may come from a wrong parsing and be indeed duplicated residues.
    # NEVER FORGET: Merged residues may be generated when calling the structure.auto_chainer
    # Note that we need bonds to check for merged residues
    # If bonds are missing then we skip this check
    if not missing_any_bonds and structure.check_merged_residues(fix_residues=True, display_summary=True):
        # Update the structure file using the corrected structure
        print(' The structure file has been modified (merged residues) -> ' + output_structure_file.filename)
        structure.generate_pdb_file(output_structure_file.path)

    # ------------------------------------------------------------------------------------------
    # Splitted residues ------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # Note that we need bonds to check for splitted residues
    # If bonds are missing then we skip this check
    if not missing_any_bonds and structure.check_repeated_residues(fix_residues=True, display_summary=True):
        # Update the structure file using the corrected structure
        print(' The structure file has been modified (repeated residues) -> ' + output_structure_file.filename)
        structure.generate_pdb_file(output_structure_file.path)

        # Sort trajectory coordinates in case atoms were sorted
        if input_trajectory_file.path and structure.trajectory_atom_sorter:
            # Save a warning in the register
            warn('Atoms have been sorted to solve splitted residues')
            print('Creating resorted files for atom bonds and charges')
            # Bonds are already resorted
            save_json(safe_bonds, MD.project.resorted_bonds_file.path, indent=4)
            # Charges are to be resorted
            # DANI: Esto no se prueba desde hace tiempo y ha sufrido cambios
            # DANI: Debería funcionar pero puede tener algun bug
            charges = get_charges(input_topology_file)
            resorted_charges = [charges[index] for index in structure.new_atom_order]
            MD.project.get_charges.prefill(MD.project, resorted_charges, {
                'topology_file': input_topology_file
            })
            save_json(resorted_charges, MD.project.resorted_charges_file.path, indent=4)
            print('Sorting trajectory coordinates to fit the new structure atom sort...')
            structure.trajectory_atom_sorter(output_structure_file, input_trajectory_file, output_trajectory_file)

    # ------------------------------------------------------------------------------------------
    # Repeated atoms ---------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    if structure.check_repeated_atoms(fix_atoms=True, display_summary=True):
        # Update the structure file using the corrected structure
        print(' The structure file has been modified (repeated atoms) -> ' + output_structure_file.filename)
        structure.generate_pdb_file(output_structure_file.path)
