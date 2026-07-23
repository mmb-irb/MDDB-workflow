from mddb_workflow.utils.auxiliar import warn, save_json, list_values_match, MISSING_CHARGES
from mddb_workflow.utils.auxiliar import MISSING_BONDS, JSON_SERIALIZABLE_MISSING_BONDS
from mddb_workflow.utils.auxiliar import round_to_hundredths, store_binary_data
from mddb_workflow.utils.constants import STANDARD_TOPOLOGY_FILENAME
from mddb_workflow.utils.type_hints import *

# Beyond this atom count the topology JSON may grow too large for MongoDB's 16 MB document limit
# Charges are then written as float32 binary files in GridFS and the topology stores only a flag
CHARGES_BIN_THRESHOLD = 300000

# Set a flag to state that the data is in GridFS
BIN_FLAG = 'gridfs'

def generate_topology (
    structure : 'Structure',
    md_charges : list[list[float]],
    residue_map : dict,
    pbc_selection : 'Selection',
    cg_selection : 'Selection',
    dummy_selection : 'Selection',
    output_directory : str,
):
    """ Prepare the standard topology file to be uploaded to the database. """
    # Set the main output filepath
    output_topology_filepath = f'{output_directory}/{STANDARD_TOPOLOGY_FILENAME}'

    # We assume that the structure will be coherent among MDs but note that this is actually checked

    # The structure will be a bunch of arrays
    # Atom data
    structure_atoms = structure.atoms
    atom_count = len(structure_atoms)
    atom_names = [ None ] * atom_count
    atom_elements = [ None ] * atom_count
    atom_residue_indices = [ None ] * atom_count
    for index, atom in enumerate(structure_atoms):
        atom_names[index] = atom.name
        atom_elements[index] = atom.element
        atom_residue_indices[index] = atom.residue.index

    # Set atom species: unique (name, element) pairs stored as compact [name, element] lists
    species_map = {}
    atom_species = []
    atom_species_indices = [None] * atom_count
    for i in range(atom_count):
        key = (atom_names[i], atom_elements[i])
        if key not in species_map:
            species_map[key] = len(atom_species)
            atom_species.append([atom_names[i], atom_elements[i]])
        atom_species_indices[i] = species_map[key]

    # Set the atom bonds
    # In order to make it more standard sort atom bonds by their indices
    # Also replace missing bonds exceptions by a json serializable flag
    # In order to save disk, do not list redundant bonds
    # e.g. if A-B is already listed do not list B-A
    # Atoms with more bonds must be the ones listing the most bonds to save disk
    # e.g. is more efficient 'D: A, B, C' than 'A: D, B: D, C: D'
    # To do so, we will have to this in multiple steps
    # Note that atoms will "loose" bonds to other atoms with more bonds, so the number may change

    # Make a copy of the original atom bonds
    unassigned_atom_bonds = [ 
        MISSING_BONDS if bonds == MISSING_BONDS else [ atom_index for atom_index in bonds ]
        for bonds in structure.bonds
    ]

    # Set the final optimized bonds
    # Serialize missing bonds exceptions
    final_atom_bonds = [
        JSON_SERIALIZABLE_MISSING_BONDS if bonds == MISSING_BONDS else []
        for bonds in unassigned_atom_bonds
    ]

    # Store atoms until we have not unassigned bonds left
    while True:
        # Counting the number of bonds per atom
        atom_bond_counts = [ 0 if bonds == MISSING_BONDS else len(bonds) for bonds in unassigned_atom_bonds ]
        # Find the atoms with most bonds and start assigning them final bonds
        highest_bond_count = max(atom_bond_counts)
        # If there are no more bonds to assign then we are done
        if highest_bond_count == 0: break
        # Iterate atom bonds looking for the atom with the maximum possible number of atoms
        for atom_index, atom_bonds in enumerate(unassigned_atom_bonds):
            if atom_bonds is MISSING_BONDS: continue
            # Repeat the count as it may have changed
            # Note that previous atoms may have claimed bonds from further atoms
            bonds_count = len(atom_bonds)
            if bonds_count < highest_bond_count: continue
            # If this atom has the maximum number of bonds then assign these bonds to it
            for bonded_atom_index in list(atom_bonds):
                final_atom_bonds[atom_index].append(bonded_atom_index)
                unassigned_atom_bonds[atom_index].remove(bonded_atom_index)
                unassigned_atom_bonds[bonded_atom_index].remove(atom_index)

    # Residue data
    structure_residues = structure.residues
    residue_count = len(structure_residues)
    residue_names = [ None ] * residue_count
    residue_numbers = [ None ] * residue_count
    # Icodes are saved as a dictionary since usually only a few residues have icodes (or none)
    # Resiude ids are used as keys and, when loaded to mongo, they will become strings
    # Saving icodes as an array would be unefficient since it will result in an array filled with nulls
    residue_icodes = {}
    residue_chain_indices = [ None ] * residue_count
    for index, residue in enumerate(structure_residues):
        residue_names[index] = residue.name
        residue_numbers[index] = residue.number
        if residue.icode:
            residue_icodes[index] = residue.icode
        residue_chain_indices[index] = residue.chain.index

    # In case there are not icodes at all set the icodes as None (i.e. null for mongo)
    if len(list(residue_icodes.keys())) == 0:
        residue_icodes = None

    # Chain data
    structure_chains = structure.chains
    chain_count = len(structure_chains)
    chain_names = [ None ] * chain_count
    for index, chain in enumerate(structure_chains):
        chain_names[index] = chain.name

    # Make a map of different possible values for each MD
    md_charges_map = []
    unique_md_charges = []
    for charges in md_charges:
        index_match = None
        # Check we have charges and, if not, set charges as None (i.e. null for mongo)
        has_charges = charges != MISSING_CHARGES and charges != None and len(charges) > 0
        # Normalize charges
        norm_charges = [ round_to_hundredths(c) for c in charges ] if has_charges else None
        # Check if these charges were found already
        for previous_index, previous_charges in enumerate(unique_md_charges):
            # We must check if charges match
            # If variable type is different then continue
            if type(norm_charges) != type(previous_charges): continue
            # Check if it is the same exact list
            # This may happen if lists come from the project
            # If not then comparte if charges match perfectly with previous values
            if norm_charges == previous_charges or (has_charges and list_values_match(norm_charges, previous_charges)):
                index_match = previous_index
                break
        # If there was no match then this is a new set of unique atom charges
        if index_match is None:
            index_match = len(unique_md_charges)
            unique_md_charges.append(norm_charges)
        # Add the current matched index to the map list
        md_charges_map.append(index_match)
    n_md_charges = len(unique_md_charges)
    # If there are no MD charges something went wrong
    if n_md_charges == 0: raise RuntimeError('No MD charges')
    # Large systems: charges are too heavy to inline in the topology JSON (MongoDB 16 MB limit)
    # Each unique charge set is written as a float32 binary file alongside the topology JSON
    # A .meta.json companion stores the set index and which MD indices share it,
    # so the API can resolve mdmap → file without re-reading the whole topology.
    if atom_count > CHARGES_BIN_THRESHOLD:
        # Group the indices of all MDs having the same charges
        md_indices_per_set = [[] for _ in range(n_md_charges)]
        for md_index, set_index in enumerate(md_charges_map):
            md_indices_per_set[set_index].append(md_index)
        # Track whenever at least one of the charges are exported
        is_bin : bool = False
        # Write every different MD charges to the binary format
        for set_index, (charges, md_indices) in enumerate(zip(unique_md_charges, md_indices_per_set)):
            # Handle missing charges
            if charges is None: continue
            is_bin = True
            # Write the binary file adding in the metadata the MDs it belongs to
            bin_path = f'{output_directory}/mdf.atom_charges_{set_index}.bin'
            store_binary_data(charges, 4, {
                'mds': md_indices,
                'x': {
                    'name': 'atoms',
                    'length': atom_count,
                },
                'bitsize': 32,
            }, bin_path)
        # Set the charges as the binary flag either if there is one or more MD charges
        final_atom_charges = BIN_FLAG if is_bin else None
    # Small systems: inline charges in the topology JSON as before
    else:
        # If there is only one set of charges then just set these charges as the value
        if n_md_charges == 1:
            final_atom_charges = unique_md_charges[0]
        # Otherwise set a more complex map
        else:
            final_atom_charges = { 'mdmap': md_charges_map, 'values': unique_md_charges }

    # Set custom atom selections
    selections = {}

    # Set a function to format selections in a database-efficient format
    # Parse selections to indices of chains, residues and atoms thus making the dict as light as possible
    # Alternatively the formatted selection may be the string 'all'
    def format_selection (selection : 'Selection') -> dict | str:
        if not selection: raise ValueError('Trying to format empty selection')
        if len(selection) == structure.atom_count: return 'all'
        atom_indices = set(selection.atom_indices)
        # Get indices of residues or chains the atoms of which are totally covered by the selection
        chain_indices = []
        residue_indices = []
        # Iterate structure chains to see if any of them fits entirely in the selection
        for chain in structure.chains:
            chain_atom_indices = set(chain.get_selection().atom_indices)
            # If the whole chain is included in the selection then add it
            if chain_atom_indices.issubset(atom_indices):
                chain_indices.append(chain.index)
                atom_indices -= chain_atom_indices
                continue
            # Iterate residues in this chain to see if any of them fits entirely in the selection
            for residue in chain.residues:
                residue_atom_indices = set(residue.get_selection().atom_indices)
                # If the whole residue is included in the selection then add it
                if residue_atom_indices.issubset(atom_indices):
                    residue_indices.append(residue.index)
                    atom_indices -= residue_atom_indices
        # Set the finally formatted selection
        formatted = { 'n': len(selection) }
        if len(chain_indices) > 0:
            formatted['ch'] = chain_indices
        if len(residue_indices) > 0:
            formatted['re'] = residue_indices
        if len(atom_indices) > 0:
            formatted['at'] = list(atom_indices).sort()
        return formatted
    
    # Add the main selections
    if pbc_selection: selections['pbc'] = format_selection(pbc_selection)
    if cg_selection: selections['cg'] = format_selection(cg_selection)
    if dummy_selection: selections['dummy'] = format_selection(dummy_selection)

    # Add residue classes as selections
    available_classes = set([ residue.classification for residue in structure.residues ])
    for classification in available_classes:
        selection = structure.select_by_classification(classification)
        selections[classification] = format_selection(selection)

    # Setup the final output
    topology = {
        'atom_species': atom_species,
        'atom_species_indices': atom_species_indices,
        'atom_charges': final_atom_charges,
        'atom_residue_indices': atom_residue_indices,
        'atom_bonds': final_atom_bonds,
        'residue_names': residue_names,
        'residue_numbers': residue_numbers,
        'residue_icodes': residue_icodes,
        'residue_chain_indices': residue_chain_indices,
        'chain_names': chain_names,
        # Residues map
        **residue_map,
        # Save also some residue indices lists here
        'selections': selections,
        'version': '0.1.0',
    }
    save_json(topology, output_topology_filepath)
