from model_workflow.utils.auxiliar import warn, save_json, MISSING_CHARGES
from model_workflow.utils.auxiliar import MISSING_BONDS, JSON_SERIALIZABLE_MISSING_BONDS
from model_workflow.utils.type_hints import *

def generate_topology (
    structure : 'Structure',
    charges : list[int],
    residue_map : dict,
    pbc_residues : list[int],
    cg_residues : list[int],
    output_filepath : str
):
    """ Prepare the standard topology file to be uploaded to the database. """

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

    # Set the atom bonds
    # In order to make it more standard sort atom bonds by their indices
    # Also replace missing bonds exceptions by a json serializable flag
    atom_bonds = []
    for atom_indices in structure.bonds:
        if atom_indices == MISSING_BONDS:
            atom_bonds.append(JSON_SERIALIZABLE_MISSING_BONDS)
            continue
        atom_bonds.append(sorted(atom_indices))

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

    # Check we have charges and, if not, set charges as None (i.e. null for mongo)
    has_charges = charges != MISSING_CHARGES and charges != None and len(charges) > 0
    atom_charges = charges if has_charges else None
    if not atom_charges:
        warn('Topology is missing atom charges')

    # Setup the final output
    topology = {
        'atom_names': atom_names,
        'atom_elements': atom_elements,
        'atom_charges': atom_charges,
        'atom_residue_indices': atom_residue_indices,
        'atom_bonds': atom_bonds,
        'residue_names': residue_names,
        'residue_numbers': residue_numbers,
        'residue_icodes': residue_icodes,
        'residue_chain_indices': residue_chain_indices,
        'chain_names': chain_names,
        # Residues map
        **residue_map,
        # Save also some residue indices lists here
        'pbc_residues': pbc_residues,
        'cg_residues': cg_residues,
    }
    save_json(topology, output_filepath)