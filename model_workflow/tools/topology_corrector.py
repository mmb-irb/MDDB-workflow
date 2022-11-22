from typing import Optional

# Import mdtoolbelt tools
from mdtoolbelt.structures import Structure

# Import local tools
from model_workflow.tools.vmd_processor import vmd_chainer

# Analyze the topology looking for irregularities and then modify the topology to standarize the format
#
# Supported cases:
#
# * Missing chains -> Chains are added through VMD
# * Repeated chains -> Chains are renamed (e.g. A, G, B, G, C, G -> A, G, B, H, C, I)
# * Splitted residues -> Atoms are sorted together by residue. Trajectory coordinates are also sorted
# * Repeated residues -> Residues are renumerated (e.g. 1, 2, 3, 1, 2, 1, 2 -> 1, 2, 3, 4, 5, 6, 7)
# * Repeated atoms -> Atoms are renamed with their numeration (e.g. C, C, C, O, O -> C1, C2, C3, O1, O2)
# * Missing/Non-standard elements -> Atom elements are guessed when missing and standarized (e.g. ZN -> Zn)


def topology_corrector (
    input_pdb_filename: str,
    output_topology_filename: str,
    input_trajectory_filename : Optional[str] = None,
    output_trajectory_filename : Optional[str] = None):

    print('----- Correcting topology -----')

    # Track if there has been any modification and then topology must be rewritten
    modified = False

    # Remove all bytes which can not be decoded by utf-8 codec
    # DANI: Esto era para que no fallase prody, pero ahora alomejor ya no es necesario
    #purgeNonUtf8(input_pdb_filename)

    # Import the pdb file and parse it to a structure object
    structure = Structure.from_pdb_file(input_pdb_filename)

    # ------------------------------------------------------------------------------------------
    # Missing chains --------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # Check if chains are missing. If so, create a new chainned topology and set it as the reference
    chains = structure.chains
    
    # In case there are not chains at all
    if len(chains) == 1 and ( chains[0].name == ' ' or chains[0].name == 'X' ):
        print('WARNING: chains are missing and they will be added')

        # Use VMD to set chains according to fragments
        vmd_chainer(input_pdb_filename, input_pdb_filename)

        # Import the new pdb file and parse it again
        structure = Structure.from_pdb_file(input_pdb_filename)
        # Redefine chains as well
        chains = structure.chains

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
    # Repeated chains --------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    if structure.check_repeated_chains(fix_chains=True, display_summary=True):
        modified = True

    # ------------------------------------------------------------------------------------------
    # Repeated residues ------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    if structure.check_repeated_residues(fix_residues=True, display_summary=True):
        modified = True

        # Sort trajectory coordinates in case atoms were sorted
        if input_trajectory_filename and structure.trajectory_atom_sorter:
            print('Sorting trajectory coordinates to fit the new structure atom sort...')
            structure.trajectory_atom_sorter(input_pdb_filename, input_trajectory_filename, output_trajectory_filename)

    # ------------------------------------------------------------------------------------------
    # Repeated atoms --------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    if structure.check_repeated_atoms(fix_atoms=True, display_summary=True):
        modified = True

    # ------------------------------------------------------------------------------------------
    # Missing/Non-standard elements ------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    if structure.fix_atom_elements(trust=False):
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