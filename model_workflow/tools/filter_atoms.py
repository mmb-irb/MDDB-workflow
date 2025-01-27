from os import remove
from os.path import exists
from subprocess import run, PIPE, Popen
import json

import pytraj as pt

from model_workflow.utils.constants import STANDARD_TOPOLOGY_FILENAME, RAW_CHARGES_FILENAME, GREY_HEADER, COLOR_END
from model_workflow.utils.constants import STANDARD_SOLVENT_RESIDUE_NAMES, STANDARD_COUNTER_ION_ATOM_NAMES
from model_workflow.utils.constants import GROMACS_EXECUTABLE
from model_workflow.utils.structures import Structure
from model_workflow.utils.auxiliar import save_json
from model_workflow.utils.gmx_spells import get_tpr_atom_count
from model_workflow.utils.type_hints import *
from model_workflow.tools.get_charges import get_raw_charges

# Set the gromacs indices filename
index_filename = 'filter.ndx'
# Set the name for the group name in gromacs ndx file
filter_group_name = "not_water_or_counter_ions"

# Filter atoms of all input topologies by remvoing atoms and ions
# As an exception, some water and ions may be not removed if specified
# At the end, all topologies must match in atoms count
def filter_atoms (
    input_structure_file : 'File',
    input_trajectory_file : 'File',
    input_topology_file : 'File',
    output_structure_file : 'File',
    output_trajectory_file : 'File',
    output_topology_file : 'File',
    # Filter selection may be a custom selection or true
    # If true then we run a default filtering of water and counter ions
    filter_selection : Union[bool, str],
    filter_selection_syntax : str = 'vmd'
):

    # Handle missing filter selection
    if not filter_selection:
        return
    
    # Load the reference structure
    reference_structure = Structure.from_pdb_file(input_structure_file.path)

    # Parse the selection to be filtered
    # WARNING: Note that the structure is not corrected at this point and there may be limitations
    parsed_filter_selection = None
    if filter_selection == True: parsed_filter_selection = reference_structure.select_water_and_counter_ions() 
    else: parsed_filter_selection = reference_structure.select(filter_selection, syntax=filter_selection_syntax)

    # Invert the parsed selection to get the atoms to remain
    keep_selection = reference_structure.invert_selection(parsed_filter_selection)

    # Set the pytraj mask to filter the desired atoms from the structure
    filter_mask = keep_selection.to_pytraj()

    # Load the structure and trajectory
    trajectory = pt.iterload(input_trajectory_file.path, input_structure_file.path)
    pt_topology = trajectory.topology
    atoms_count = pt_topology.n_atoms

    # Set the filtered structure
    filtered_pt_topology = pt_topology[filter_mask]
    filtered_structure_atoms_count = filtered_pt_topology.n_atoms

    # Check if both the normal and the filtered topologies have the same number of atoms
    # If not, filter the whole trajectory and overwrite both topologies and trajectory
    print(f'Total number of atoms: {atoms_count}')
    print(f'Filtered number of atoms: {filtered_structure_atoms_count}')
    if filtered_structure_atoms_count < atoms_count:
        print('Filtering structure and trajectory...')
        # Filter both structure and trajectory using Gromacs, since is more efficient than pytraj
        # Set up an index file with all atom indices manually
        # As long as indices are for atoms and not residues there should never be any incompatibility
        keep_selection.to_ndx_file(output_filepath = index_filename)
        # Filter the trajectory
        xtc_filter(input_structure_file.path, input_trajectory_file.path, output_trajectory_file.path, index_filename)
        # Filter the structure
        pdb_filter(input_structure_file.path, output_structure_file.path, index_filename)

    # Filter topology according to the file format
    if input_topology_file and input_topology_file.exists:
        print('Filtering topology...')
        # Pytraj supported formats
        if input_topology_file.is_pytraj_supported():
            # Load the topology and count its atoms
            pt_topology = pt.load_topology(filename=input_topology_file.path)
            topology_atoms_count = pt_topology.n_atoms
            print(f'Topology atoms count: {topology_atoms_count}')
            # Filter the desired atoms using the mask and then count them
            filtered_pt_topology = pt_topology[filter_mask]
            filtered_topology_atoms_count = filtered_pt_topology.n_atoms
            # If there is a difference in atom counts then write the filtered topology
            if filtered_topology_atoms_count < topology_atoms_count:
                # WARNING: If the output topology is a symlink it will try to overwrite the origin
                # Remove it to avoid overwriting input data
                if output_topology_file.is_symlink(): output_topology_file.remove()
                # Now write the filtered topology
                pt.write_parm(
                    filename=output_topology_file.path,
                    top=filtered_pt_topology,
                    format=input_topology_file.get_pytraj_parm_format(),
                    overwrite=True
                )
        # Gromacs format format
        elif input_topology_file.format == 'tpr':
            # Get the input tpr atom count
            topology_atoms_count = get_tpr_atom_count(input_topology_file.path)
            print(f'Topology atoms count: {topology_atoms_count}')
            # If the number of atoms is greater than expected then filter the tpr file
            if topology_atoms_count > filtered_structure_atoms_count:
                if not exists(index_filename):
                    # In order to filter the tpr we need the filter.ndx file
                    # This must be generated from a pytraj supported topology that matches the number of atoms in the tpr file
                    raise ValueError('Topology atoms number does not match the structure atoms number and tpr files can not be filtered alone')
                tpr_filter(input_topology_file.path, output_topology_file.path, index_filename)
                # Get the output tpr atom count
                filtered_topology_atoms_count = get_tpr_atom_count(output_topology_file.path)
            else:
                filtered_topology_atoms_count = topology_atoms_count
        # Standard topology
        elif input_topology_file.filename == STANDARD_TOPOLOGY_FILENAME:
            standard_topology = None
            with open(input_topology_file.path, 'r') as file:
                standard_topology = json.load(file)
            topology_atoms_count = len(standard_topology['atom_names'])
            print(f'Topology atoms count: {topology_atoms_count}')
            # Make it match since there is no problem when these 2 do not match
            filtered_topology_atoms_count = filtered_structure_atoms_count
            # If the number of charges does not match the number of atoms then filter the topology
            if topology_atoms_count != filtered_structure_atoms_count:
                standard_topology_filter(input_topology_file, reference_structure, parsed_filter_selection, output_topology_file)
        # Raw charges
        elif input_topology_file.filename == RAW_CHARGES_FILENAME:
            charges = get_raw_charges(input_topology_file.path)
            # Nothing to do here. It better matches by defualt or we have a problem
            filtered_topology_atoms_count = len(charges)
            print(f'Topology atoms count: {filtered_topology_atoms_count}')
        else:
            raise ValueError(f'Topology file ({input_topology_file.filename}) is in a non supported format')

        print(f'Filtered topology atoms: {filtered_topology_atoms_count}')

        # Both filtered structure and topology must have the same number of atoms
        if filtered_structure_atoms_count != filtered_topology_atoms_count:
            print(f'Filtered structure atoms: {filtered_structure_atoms_count}')
            raise ValueError('Filtered atom counts in topology and charges does not match')

    # Remove the index file in case it was created
    if exists(index_filename):
        remove(index_filename)

    # Check if any of the output files does not exist
    # If so, then it means there was nothing to filter
    # However the output file is expected, so me make symlink
    if not output_structure_file.exists:
        output_structure_file.set_symlink_to(input_structure_file)
    if not output_trajectory_file.exists:
        output_trajectory_file.set_symlink_to(input_trajectory_file)
    if not output_topology_file.exists:
        output_topology_file.set_symlink_to(input_topology_file)

# Filter atoms in a pdb file
# This method conserves maximum resolution and chains
def pdb_filter (
    input_pdb_filename : str,
    output_pdb_filename : str,
    index_filename : str
):
    # Filter the trajectory
    p = Popen([
        "echo",
        filter_group_name,
    ], stdout=PIPE)
    logs = run([
        GROMACS_EXECUTABLE,
        "editconf",
        "-f",
        input_pdb_filename,
        '-o',
        output_pdb_filename,
        '-n',
        index_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE, stderr=PIPE).stdout.decode()
    p.stdout.close()

# Filter atoms in a xtc file
# Note that here we do not hide the stderr
# This is because it shows the progress
# Instead we color the output grey
def xtc_filter(
    structure_filename : str,
    input_trajectory_filename : str,
    output_trajectory_filename : str,
    index_filename : str
):
    print(GREY_HEADER)
    # Filter the trajectory
    p = Popen([
        "echo",
        filter_group_name,
    ], stdout=PIPE)
    logs = run([
        GROMACS_EXECUTABLE,
        "trjconv",
        "-s",
        structure_filename,
        "-f",
        input_trajectory_filename,
        '-o',
        output_trajectory_filename,
        '-n',
        index_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()
    print(COLOR_END)

# Filter atoms in both a pdb and a xtc file
def tpr_filter(
    input_tpr_filename : str,
    output_tpr_filename : str,
    index_filename : str
):
    # Filter the topology
    p = Popen([
        "echo",
        filter_group_name,
    ], stdout=PIPE)
    logs = run([
        GROMACS_EXECUTABLE,
        "convert-tpr",
        "-s",
        input_tpr_filename,
        '-o',
        output_tpr_filename,
        '-n',
        index_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE, stderr=PIPE).stdout.decode()
    p.stdout.close()

# Set a function to filter the standard topology file
# WARNING: This function has not been checked in depth to work properly
def standard_topology_filter (
    input_topology_file : 'File',
    reference_structure : 'Structure',
    parsed_filter_selection : 'Selection',
    output_topology_file : 'File'):

    # Load the topology
    topology = None
    with open(input_topology_file.path, 'r') as file:
        topology = json.load(file)

    # Get filtered atom, residues and chain indices
    atom_indices = parsed_filter_selection.atom_indices
    residue_indices = reference_structure.get_selection_residue_indices(parsed_filter_selection)
    chain_indices = reference_structure.get_selection_chain_indices(parsed_filter_selection)

    # Set backmapping
    atom_backmapping = { old_index: new_index for new_index, old_index in enumerate(atom_indices) }
    residue_backmapping = { old_index: new_index for new_index, old_index in enumerate(residue_indices) }
    chain_backmapping = { old_index: new_index for new_index, old_index in enumerate(chain_indices) }

    # Set a function to get substract specific values of a list given by its indices
    def filter_by_indices (values : list, indices : List[int]) -> list:
        if values == None:
            return None
        return [ values[i] for i in indices ]

    # Filter atomwise fields
    atom_names = filter_by_indices(topology['atom_names'], atom_indices)
    atom_elements = filter_by_indices(topology['atom_elements'], atom_indices)
    atom_charges = filter_by_indices(topology['atom_charges'], atom_indices)

    # Handle atom fields which require backmapping
    old_atom_residue_indices = filter_by_indices(topology['atom_residue_indices'], atom_indices)
    atom_residue_indices = [ residue_backmapping[index] for index in old_atom_residue_indices ]
    atom_bonds = None
    raw_atom_bonds = topology.get('atom_bonds', None)
    if raw_atom_bonds:
        old_atom_bonds = filter_by_indices(raw_atom_bonds, atom_indices)
        atom_bonds = [ [ atom_backmapping(bond) for bond in bonds ] for bonds in old_atom_bonds ]

    # Filter residuewise fields
    residue_names = filter_by_indices(topology['residue_names'], residue_indices)
    residue_numbers = filter_by_indices(topology['residue_numbers'], residue_indices)

    # Handle residue fields which require backmapping
    residue_indices_set = set(residue_indices)
    old_residue_chain_indices = filter_by_indices(topology['residue_chain_indices'], residue_indices)
    residue_chain_indices = [ chain_backmapping[index] for index in old_residue_chain_indices ]
    # Handle icodes if they exist
    residue_icodes = None
    raw_residue_icodes = topology['residue_icodes']
    if raw_residue_icodes:
        # Filter icodes
        old_residue_icodes = { index: icode for index, icode in raw_residue_icodes.items if index in residue_indices_set }
        # Backmap icodes
        residue_icodes = { residue_backmapping(index): icode for index, icode in old_residue_icodes.items() }
    # Handle PBC residues
    pbc_residues = [ residue_backmapping(index) for index in topology.get('pbc_residues', []) if index in residue_indices_set ]

    # Handle chainwise fields
    chain_names = filter_by_indices(topology['chain_names'], chain_indices)

    # Handle references
    references = None
    reference_types = None
    residue_reference_indices = None
    residue_reference_numbers = None
    # If they exist
    raw_references = topology['references']
    if raw_references:
        residue_reference_numbers = filter_by_indices(topology['residue_reference_numbers'], residue_indices)
        old_residue_reference_indices = filter_by_indices(topology['residue_reference_indices'], residue_indices)
        reference_indices = [ index for index in set(old_residue_reference_indices) if index != None ]
        references_backmapping = { old_index: new_index for new_index, old_index in enumerate(reference_indices) }
        references = [ raw_references[old_index] for old_index in references_backmapping.keys() ]
        raw_reference_types = topology.get('reference_types', None)
        if raw_reference_types:
            reference_types = [ raw_reference_types[old_index] for old_index in references_backmapping.keys() ]
        references_backmapping[None] = None # To residues with no reference
        residue_reference_indices = [ references_backmapping[index] for index in old_residue_reference_indices ]

    # Set the filtered topology
    output_topology = {
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
        'references': references,
        'reference_types': reference_types,
        'residue_reference_indices': residue_reference_indices,
        'residue_reference_numbers': residue_reference_numbers,
        'pbc_residues': pbc_residues,
    }
    # Wrtie the new topology
    save_json(output_topology, output_topology_file.path)