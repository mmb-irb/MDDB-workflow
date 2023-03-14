from typing import Optional, List
from os import remove

# Import mdtoolbelt tools
from mdtoolbelt.structures import Structure

# Import local tools
from model_workflow.tools.get_safe_bonds import do_bonds_match, get_safe_bonds, get_safe_bonds_canonical_frame
from model_workflow.tools.get_pdb_frames import get_pdb_frame

# Analyze the topology looking for irregularities and then modify the topology to standarize the format
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

def topology_corrector (
    input_pdb_filename: str,
    output_topology_filename: str,
    input_trajectory_filename : Optional[str],
    output_trajectory_filename : Optional[str],
    register : dict,
    mercy : List[str],
    trust : List[str]
):

    print('----- Correcting topology -----')

    # Track if there has been any modification and then topology must be rewritten
    modified = False

    # Remove all bytes which can not be decoded by utf-8 codec
    # DANI: Esto era para que no fallase prody, pero ahora alomejor ya no es necesario
    #purgeNonUtf8(input_pdb_filename)

    # Import the pdb file and parse it to a structure object
    structure = Structure.from_pdb_file(input_pdb_filename)

    # ------------------------------------------------------------------------------------------
    # Missing/Non-standard elements ------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # It is important to fix elements before trying to fix bonds, since elements have an impact on bonds
    # VMD logic to find bonds relies in the atom element to set the covalent bond distance cutoff
    if structure.fix_atom_elements(trust=False):
        modified = True

    # ------------------------------------------------------------------------------------------
    # Unstable atom bonds ------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # Check and find if needed stable bonds
    must_check_stable_bonds = 'stabonds' not in trust
    if must_check_stable_bonds:

        # Using the trajectory, find the safe bonds (i.e. bonds stable along several frames)
        safe_bonds = get_safe_bonds(input_pdb_filename, input_trajectory_filename)
        # If the safe bonds do not match the structure bonds then we have to fix it
        if not do_bonds_match(structure.bonds, safe_bonds):
            modified = True
            print('WARNING: Default structure has wrong bonds')
            # Set the safe bonds as the structure bonds
            structure.bonds = safe_bonds
            # Find the first frame in the whole trajectory where safe bonds are respected
            safe_bonds_frame = get_safe_bonds_canonical_frame(input_pdb_filename, input_trajectory_filename, safe_bonds)
            # If there is no canonical frame then stop here since there must be a problem
            if safe_bonds_frame == None:
                print('There is no canonical frame for safe bonds. Is the trajectory not imaged?')
                must_be_killed = 'stabonds' not in mercy
                if must_be_killed:
                    raise RuntimeError()
            # Set also the safe bonds frame structure to mine its coordinates
            safe_bonds_frame_filename = get_pdb_frame(input_pdb_filename, input_trajectory_filename, safe_bonds_frame)
            safe_bonds_frame_structure = Structure.from_pdb_file(safe_bonds_frame_filename)
            # Set all coordinates in the main structure by copying the safe bonds frame coordinates
            for atom_1, atom_2 in zip(structure.atoms, safe_bonds_frame_structure.atoms):
                atom_1.coords = atom_2.coords
            # Remove the safe bonds frame since it is not requird anymore
            remove(safe_bonds_frame_filename)

    # ------------------------------------------------------------------------------------------
    # Incoherent atom bonds ---------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    must_check_coherent_bonds = 'cohbonds' not in trust
    if must_check_coherent_bonds and structure.check_incoherent_bonds():
        print('FAIL: Uncoherent bonds were found')
        must_be_killed = 'cohbonds' not in mercy
        if must_be_killed:
            raise SystemExit()

    # ------------------------------------------------------------------------------------------
    # Missing chains ---------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # Check if chains are missing. If so, create a new chainned topology and set it as the reference
    chains = structure.chains
    
    # In case there are not chains at all
    if len(chains) == 1 and ( chains[0].name == ' ' or chains[0].name == 'X' ):
        print('WARNING: chains are missing and they will be added')

        # Run the chainer
        structure.chainer()
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
        if input_trajectory_filename and structure.trajectory_atom_sorter:
            # Save a warning in the register
            register['warnings'].append('Atoms have been sorted to solve splitted residues')
            # Save the new order in the register
            # DANI: Dado que no reordenamos las topologías orignales (muchos formatos, mucho marrón) hay que guardar esto
            # DANI: Es para curarnos en salud, pero lo suyo sería poder exportar topologías de la API que ya tengan los datos bien
            register['new_atom_order'] = structure.new_atom_order
            print('Sorting trajectory coordinates to fit the new structure atom sort...')
            structure.trajectory_atom_sorter(input_pdb_filename, input_trajectory_filename, output_trajectory_filename)

    # ------------------------------------------------------------------------------------------
    # Repeated atoms ---------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    if structure.check_repeated_atoms(fix_atoms=True, display_summary=True):
        modified = True

    # ------------------------------------------------------------------------------------------
    # Final output -----------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # Write a new topology if any modification was done
    if modified:
        print(' The topology file has been modificated -> ' + output_topology_filename)
        structure.generate_pdb_file(output_topology_filename)
    else:
        print(' Everything is fine')

    print('-------------------------------')


# Set a function to get the next letter from an input letter in alphabetic order
letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
    'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
def get_new_letter(current_letters : list) -> str:
    new_letter = next((letter for letter in letters if letter not in current_letters), None)
    if not new_letter:
        raise SystemExit("There are no more letters")
    return new_letter

# Remove all bytes which can not be decoded by utf-8 codec
# This prevents the prody parsePDB function to return the following error:
# UnicodeDecodeError: 'utf-8' codec can't decode byte ...
def purgeNonUtf8 (filename : str):
    with open(filename, mode="r+", encoding="utf-8", errors= 'ignore') as file:
        lines = file.readlines()
        file.seek(0)
        for line in lines:
            file.write(line)
        file.truncate()