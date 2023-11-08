from typing import Optional, List
from os import remove
from os.path import exists
from json import load

# Import mdtoolbelt tools
from mdtoolbelt.structures import Structure

# Import local tools
from model_workflow.constants import TOPOLOGY_FILENAME
from model_workflow.tools.get_safe_bonds import do_bonds_match, get_safe_bonds, get_safe_bonds_canonical_frame
from model_workflow.tools.get_pdb_frames import get_pdb_frame

# Import other software
import pytraj as pt

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

# Note that the 'mercy' flag may be passed for crticial checkings to not kill the process on fail

def structure_corrector (
    input_structure_file : 'File',
    input_trajectory_file : Optional['File'],
    input_topology_file : Optional['File'],
    output_structure_file : 'File',
    output_trajectory_file : Optional['File'],
    snapshots : int,
    register : 'Register',
    mercy : List[str],
    trust : List[str]
):

    # Track if there has been any modification and then the structure must be rewritten
    modified = False

    # Import the pdb file and parse it to a structure object
    structure = Structure.from_pdb_file(input_structure_file.path)

    # ------------------------------------------------------------------------------------------
    # Missing/Non-standard elements ------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # It is important to fix elements before trying to fix bonds, since elements have an impact on bonds
    # VMD logic to find bonds relies in the atom element to set the covalent bond distance cutoff
    if structure.fix_atom_elements(trust=False):
        modified = True

    # ------------------------------------------------------------------------------------------
    # Unstable atom bonds ----------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # Set if stable bonds have to be checked
    stable_bonds_flag = 'stabonds'
    must_check_stable_bonds = stable_bonds_flag not in trust

    # If this analysis has been already passed then we skip the process
    if register.tests.get(stable_bonds_flag, None) == True:
        must_check_stable_bonds = False

    # Try to get bonds from the topology before guessing
    safe_bonds = get_bonds(input_topology_file)
    # If failed to mine topology bonds then guess stable bonds
    if not safe_bonds:
        print('Bonds will be guessed by atom distances and radius')
        # Find stable bonds if necessary
        if must_check_stable_bonds:
            # Using the trajectory, find the safe bonds (i.e. bonds stable along several frames)
            safe_bonds = get_safe_bonds(input_structure_file.path, input_trajectory_file.path, snapshots)
        # If we trust stable bonds then simply use structure bonds
        else:
            safe_bonds = structure.bonds
    # If the safe bonds do not match the structure bonds then we have to fix it
    if must_check_stable_bonds and not do_bonds_match(structure.bonds, safe_bonds):
        modified = True
        print('WARNING: Default structure has wrong bonds')
        # Set the safe bonds as the structure bonds
        structure.bonds = safe_bonds
        # Find the first frame in the whole trajectory where safe bonds are respected
        safe_bonds_frame = get_safe_bonds_canonical_frame(input_structure_file.path, input_trajectory_file.path, snapshots, safe_bonds)
        # If there is no canonical frame then stop here since there must be a problem
        if safe_bonds_frame == None:
            print('There is no canonical frame for safe bonds. Is the trajectory not imaged?')
            must_be_killed = stable_bonds_flag not in mercy
            if must_be_killed:
                raise SystemExit('Failed to find stable bonds')
            register.tests[stable_bonds_flag] = False
            register.warnings.append(('Could not find a frame in the trajectory respecting all bonds if bonds were predicted according to atom coordinates.\n'
            'The main PDB structure is the default structure and it would be considered to have wrong bonds if they were predicted as previously stated.'))
        else:
            # Set also the safe bonds frame structure to mine its coordinates
            safe_bonds_frame_filename = get_pdb_frame(input_structure_file.path, input_trajectory_file.path, safe_bonds_frame)
            safe_bonds_frame_structure = Structure.from_pdb_file(safe_bonds_frame_filename)
            # Set all coordinates in the main structure by copying the safe bonds frame coordinates
            for atom_1, atom_2 in zip(structure.atoms, safe_bonds_frame_structure.atoms):
                atom_1.coords = atom_2.coords
            # Remove the safe bonds frame since it is not required anymore
            remove(safe_bonds_frame_filename)

    # Tag the test as succeed if we did not skip it
    if must_check_stable_bonds and register.tests.get(stable_bonds_flag, None) != False:
        register.tests[stable_bonds_flag] = True

    # ------------------------------------------------------------------------------------------
    # Incoherent atom bonds ---------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # Set if coherent bonds have to be checked
    coherent_bonds_flag = 'cohbonds'
    must_check_coherent_bonds = coherent_bonds_flag not in trust

    # If this analysis has been already passed then we skip the process
    if register.tests.get(coherent_bonds_flag, None) == True:
        must_check_stable_bonds = False

    # Run the coherent bonds analysis if necessary
    if must_check_coherent_bonds and structure.check_incoherent_bonds():
        print('FAIL: Uncoherent bonds were found')
        must_be_killed = coherent_bonds_flag not in mercy
        if must_be_killed:
            raise SystemExit('Failed to find coherent bonds')
        register.tests[coherent_bonds_flag] = False
        register.warnings.append('Bonds are not coherent. Some atoms may have less/more bonds than they should.')

    # Tag the test as succeed if we did not skip it
    if must_check_coherent_bonds and register.tests.get(coherent_bonds_flag, None) != False:
        register.tests[coherent_bonds_flag] = True

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
            raise SystemExit('ERROR: We are having splitted chains (e.g. ' + chain_name + ')')
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
            register.warnings.append('Atoms have been sorted to solve splitted residues')
            # Save the new order in the register
            # DANI: Dado que no reordenamos las topologías orignales (muchos formatos, mucho marrón) hay que guardar esto
            # DANI: Es para curarnos en salud, pero lo suyo sería poder exportar topologías de la API que ya tengan los datos bien
            register.cache['new_atom_order'] = structure.new_atom_order
            print('Sorting trajectory coordinates to fit the new structure atom sort...')
            structure.trajectory_atom_sorter(input_structure_file.path, input_trajectory_file.path, output_trajectory_file.path)

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


# Set a function to get the next letter from an input letter in alphabetic order
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
    'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
def get_new_letter(current_letters : list) -> str:
    new_letter = next((letter for letter in letters if letter not in current_letters), None)
    if not new_letter:
        raise SystemExit("There are no more letters")
    return new_letter

# Extract bonds from a source file
def get_bonds (bonds_source_file : 'File') -> list:
    if not bonds_source_file or not bonds_source_file.exists:
        return None
    # If we have the standard topology then get bonds from it
    if bonds_source_file.filename == TOPOLOGY_FILENAME:
        print('Bonds in the "' + bonds_source_file.filename + '" file will be used')
        standard_topology = None
        with open(bonds_source_file.path, 'r') as file:
            standard_topology = load(file)
        bonds = standard_topology.get('atom_bonds', None)
        if bonds:
            return bonds
    # In some ocasions, bonds may come inside a topology which can be parsed through pytraj
    if bonds_source_file.is_pytraj_supported:
        print('Bonds will be mined from "' + bonds_source_file.path + '"')
        pt_topology = pt.load_topology(filename=bonds_source_file.path)
        atom_bonds = [ [] for i in range(pt_topology.n_atoms) ]
        for bond in pt_topology.bonds:
            a,b = bond.indices
            atom_bonds[a].append(b)
            atom_bonds[b].append(a)
        # If there is any bonding data then return bonds
        if any(len(bonds) > 0 for bonds in atom_bonds):
            return atom_bonds
        # If all bonds are empty then it means the parsing failed or the pytraj topology has no bonds
        # We must guess them
        print(' Bonds could not be mined')
    # If we can not mine bonds then return None and they will be guessed further
    return None