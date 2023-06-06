from os import remove
from os.path import exists
from subprocess import run, PIPE, Popen

import pytraj as pt

from model_workflow.tools.get_charges import get_raw_charges, raw_charges_filename, get_tpr_charges
from model_workflow.tools.formats import is_raw, is_pytraj_supported, is_tpr, get_pytraj_parm_format

from mdtoolbelt.selections import Selection

# Set the gromacs indices filename
index_filename = 'filter.ndx'
# Set the name for the group name in gromacs ndx file
filter_group_name = "not_water_or_counter_ions"

# Set the standard topology name
standard_topology_filename = 'topology.json'

# Filter atoms of all input topologies by remvoing atoms and ions
# As an exception, some water and ions may be not removed if specified
# At the end, all topologies must match in atoms count
def filter_atoms (
    topology_filename : str,
    trajectory_filename : str,
    charges_filename : str,
    filter_selection : str
):

    # Handle missing filter selection
    if not filter_selection:
        return

    # Load the topology and trajectory
    trajectory = pt.iterload(trajectory_filename, topology_filename)
    topology = trajectory.topology
    atoms_count = topology.n_atoms

    # Set the pytraj mask to filter the desired atoms from the topology
    filter_mask = get_filter_mask(topology_filename, filter_selection)

    # Set the filtered topology
    filtered_topology = topology[filter_mask]
    filtered_atoms_count = filtered_topology.n_atoms

    # Check if both the normal and the filtered topologies have the same number of atoms
    # If not, filter the whole trajectory and overwrite both topologies and trajectory
    print('Total number of atoms: ' + str(atoms_count))
    print('Filtered number of atoms: ' + str(filtered_atoms_count))
    if filtered_atoms_count < atoms_count:
        print('Filtering structure and trajectory...')
        # Filter both topology and trajectory using Gromacs, since is more efficient than pytraj
        # Set up an index file with all atom indices manually
        # As long as indices are for atoms and not residues there should never be any incompatibility
        pytraj_mask_2_gromacs_ndx(topology, { filter_group_name : filter_mask }, index_filename)
        # Filter the trajectory
        xtc_filter(topology_filename, trajectory_filename, trajectory_filename, index_filename)
        # Filter the topology
        pdb_filter(topology_filename, topology_filename, index_filename)

    # Filter charges according to the file format
    if charges_filename and exists(charges_filename):
        print('Processing charges')
        # Already parsed charges
        if charges_filename == standard_topology_filename:
            with open(standard_topology_filename, 'r') as file:
                standard_topology = load(file)
            charges = standard_topology['atom_charges']
            print('Charges count: ' + str(charges))
            # Make it match since there is no problem when these 2 do not match
            filtered_charges_atoms_count = filtered_atoms_count
            # If the number of charges does not match the number of atoms then just remove the standard topolog file
            # This way it will be generated again further but from the new filtered structure
            if charges != filtered_atoms_count:
                remove(standard_topology_filename)
        # Raw charges
        elif is_raw(charges_filename):
            charges = get_raw_charges(charges_filename)
            # Nothing to do here. It better matches by defualt or we have a problem
            filtered_charges_atoms_count = len(charges)
            print('Charges count: ' + str(filtered_charges_atoms_count))
        # Pytraj supported formats
        elif is_pytraj_supported(charges_filename):
            # Load the charges topology and count its atoms
            charges_topology = pt.load_topology(filename=charges_filename)
            charges_atoms_count = charges_topology.n_atoms
            print('Charges count: ' + str(charges_atoms_count))
            # Filter the desired atoms using the mask and then count them
            filter_mask = get_filter_mask(charges_filename, filter_selection)
            filtered_charges_topology = charges_topology[filter_mask]
            filtered_charges_atoms_count = filtered_charges_topology.n_atoms
            # If there is a difference in atom counts then write the filtered topology
            if filtered_charges_atoms_count < charges_atoms_count:
                pt.write_parm(
                    filename=charges_filename,
                    top=filtered_charges_topology,
                    format=get_pytraj_parm_format(charges_filename),
                    overwrite=True
                )
        # Gromacs format format
        elif is_tpr(charges_filename):
            # Extract charges from the tpr file and count them
            charges = get_tpr_charges(charges_filename)
            charges_atoms_count = len(charges)
            print('Charges count: ' + str(charges_atoms_count))
            # If the number of charges in greater than expected then filter the tpr file and extract charges again
            if charges_atoms_count > filtered_atoms_count:
                if not exists(index_filename):
                    # In order to filter the tpr we need the filter.ndx file
                    # This must be generated from a pytraj supported topology that matches the number of atoms in the tpr file
                    raise ValueError('Charges number does not match the structure atoms and tpr files can not be filtered alone')
                tpr_filter(charges_filename, charges_filename, index_filename)
                charges = get_tpr_charges(charges_filename)
            filtered_charges_atoms_count = len(charges)
        else:
            raise ValueError('Charges file (' + charges_filename + ') is in a non supported format')

        print('Filtered charges atoms: ' + str(filtered_charges_atoms_count))

        # Both filtered topology and charges must have the same number of atoms
        if filtered_atoms_count != filtered_charges_atoms_count:
            print('Filtered topology atoms: ' + str(filtered_atoms_count))
            raise ValueError('Filtered atom counts in topology and charges does not match')

    # Remove the index file in case it was created
    if exists(index_filename):
        remove(index_filename)

# Set the pytraj mask to filter the desired atoms from a specific topology
def get_filter_mask (source_topology_filename : str, filter_selection : str) -> str:
    # If the default filtering was rquested
    if filter_selection == True:
        return get_default_filter_mask(source_topology_filename)
    # The filter selection is meant to be in vmd selection format
    filter_mask = vmd_selection_2_pytraj_mask(source_topology_filename, filter_selection)
    return filter_mask

# Set the pytraj selection for not water
water_mask = '(:SOL,WAT,HOH,TIP)'
# Set the pytraj mask to filter default atoms
# i.e. water and counter ions with standard residue names
def get_default_filter_mask (source_topology_filename : str) -> str:
    filter_mask = '!' + water_mask
    counter_ions_mask = get_counter_ions_mask(source_topology_filename)
    if counter_ions_mask:
        filter_mask = '(!' + water_mask + '&!(' + counter_ions_mask + '))'
    return filter_mask

# Escape all vmd reserved characters from the selection string, which may use regular expression characters
vmd_reserved_characters = ['"','[',']']
def escape_vmd (selection : str) -> str:
    escaped_selection = selection
    for character in vmd_reserved_characters:
        escaped_selection = escaped_selection.replace(character, '\\' + character)
    return escaped_selection

# Set the vmd script filename
commands_filename = 'commands.vmd'
# Convert a vmd selection to a group of atom indices
def vmd_selection_2_pytraj_mask (source_topology_filename : str, filter_selection : str) -> dict:

    # Prepare a script for the VMD to geth a selection atom indices. This is Tcl lenguage
    with open(commands_filename, "w") as file:
        # Select the specified atoms and set the specified chain
        file.write('set atoms [atomselect top "' + escape_vmd(filter_selection) + '"]\n')
        file.write('$atoms list\n')
        file.write('exit\n')

    # Run VMD
    vmd_logs = run([
        "vmd",
        source_topology_filename,
        "-e",
        commands_filename,
        "-dispdev",
        "none"
    ], stdout=PIPE).stdout.decode()

    # Mine and parse the VMD logs into an atoms indices list
    atom_indices_line = None
    vmd_log_lines = vmd_logs.split('\n')
    for l, line in enumerate(vmd_log_lines):
        if line == 'atomselect0':
            atom_indices_line = vmd_log_lines[l + 1]
            break

    if atom_indices_line == None:
        print(vmd_logs)
        raise SystemExit('Something went wrong with VMD')

    # WARNING: Although pytraj atom indices goes from 0 to n-1, they go from 1 to n in the selection mask
    atom_indices = [ str(int(index) + 1) for index in atom_indices_line.split(' ') ]
    filter_mask = '@' + ','.join(atom_indices)
    return filter_mask

# Get a pytraj selection with all counter ions
counter_ions = ['K', 'NA', 'CL', 'CLA', 'SOD', 'POT']
def get_counter_ions_mask (topology_filename : str) -> str:
    pt_topology = pt.load_topology(filename=topology_filename)

    # Get all atoms from single atom residues
    single_atoms = []
    for residue in pt_topology.residues:
        if residue.n_atoms != 1:
            continue
        single_atoms.append(residue.first_atom_index)

    # Get a list with all topology atoms
    atoms = list(pt_topology.atoms)

    # Get atoms whose name matches any counter ion names list
    counter_ion_atoms = []
    for atom_index in single_atoms:
        atom = atoms[atom_index]
        atom_name = atom.name.upper()
        # Remove possible '+' and '-' signs by keeping only letters
        simple_atom_name = ''.join(filter(str.isalpha, atom_name))
        if simple_atom_name in counter_ions:
            # Atom indices go from 0 to n-1
            # Add +1 to the index since the mask selection syntax counts from 1 to n
            counter_ion_atoms.append(str(atom_index +1))
    
    # Return None in case there are no counter ion atoms
    if len(counter_ion_atoms) == 0:
        return None
    
    # Return atoms in a pytraj mask format
    counter_ions_mask = '@' + ','.join(counter_ion_atoms)
    return counter_ions_mask

# Use this function to create gromacs 'index.ndx' files to match pytraj mask selections in a given pytraj topology
# 'masks' is a dict where each entry key will be the gromacs group name and each entry value is the mask string
# e.g. { "my_solvent" : ":SOL,WAT,HOH", "my_ions" : "@17,18" }
def pytraj_mask_2_gromacs_ndx (topology : 'Topology', masks : dict, output_filename : str = 'index.ndx'):
    # Parse the masks object into the a dict with atom indices
    groups = pytraj_mask_2_atom_indices(topology, masks)
    # Now generate the .ndx file
    atom_indices_2_gromacs_ndx(groups, output_filename)

# Convert groups of pytraj masks to atom indices
# 'masks' is a dict where each entry key will be the gromacs group name and each entry value is the mask string
# e.g. { "my_solvent" : ":SOL,WAT,HOH", "my_ions" : "@17,18" }
def pytraj_mask_2_atom_indices (topology : 'Topology', masks : dict) -> dict:
    # Parse the masks object into a dict with atom indices
    groups = {}
    for group_name, mask in masks.items():
        # Find mask atom indices in the topology
        atom_indices = list(topology.atom_indices(mask))
        groups[group_name] = atom_indices
    return groups

# Use this function to create gromacs 'index.ndx' files from atom indexes
# 'groups' is a dict where each entry key will be the gromacs group name and each entry value is the atoms integer list
# e.g. { "my_solvent" : [1,2,3], "my_ions" : [4,5,6] }
def atom_indices_2_gromacs_ndx (groups : dict, output_filename : str = 'index.ndx'):
    with open(output_filename, "w") as file:
        for group_name, atom_indices in groups.items():
            # Add a header with the name for each group
            content = '[ ' + group_name + ' ]\n'
            count = 0
            for index in atom_indices:
                # Add a breakline each 15 indices
                count += 1
                if count == 15:
                    content += '\n'
                    count = 0
                # Add a space between indices
                # Atom indices go from 0 to n-1
                # Add +1 to the index since gromacs counts from 1 to n
                content += str(index + 1) + ' '
            file.write(content + '\n')

# Filter atoms in a pdb file
# This method conserves maximum resolution and chains
def pdb_filter(
    input_topology_filename : str,
    output_topology_filename : str,
    index_filename : str
):
    # Filter the trajectory
    p = Popen([
        "echo",
        filter_group_name,
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "editconf",
        "-f",
        input_topology_filename,
        '-o',
        output_topology_filename,
        '-n',
        index_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

# Filter atoms in a xtc file
def xtc_filter(
    topology_filename : str,
    input_trajectory_filename : str,
    output_trajectory_filename : str,
    index_filename : str
):
    # Filter the trajectory
    p = Popen([
        "echo",
        filter_group_name,
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "trjconv",
        "-s",
        topology_filename,
        "-f",
        input_trajectory_filename,
        '-o',
        output_trajectory_filename,
        '-n',
        index_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

# Filter atoms in both a pdb and a xtc file
def tpr_filter(
    input_topology_filename : str,
    output_topology_filename : str,
    index_filename : str
):
    # Filter the topology
    p = Popen([
        "echo",
        filter_group_name,
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "convert-tpr",
        "-s",
        input_topology_filename,
        '-o',
        output_topology_filename,
        '-n',
        index_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

# Filter a structure and a trajectory file together
# DANI: No se ha provado desde que se rehizo todo esto, cuidao
def get_unmembraned_trajectory (
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_topology_filename : str,
    output_trajectory_filename : str,
    structure : 'Structure',
    membranes : list
):
    # Get all membrane atom indices
    membrane_atom_indices = []
    for membrane in membranes:
        membrane_parsed_selection = structure.select(membrane['selection'], syntax='vmd')
        membrane_atom_indices += membrane_parsed_selection.atom_indices
    membrane_atom_indices = list(set(atom_indices))
    # Now get all non-membrane atom indices
    atom_indices = [ i for i, atom in enumerate(structure.atoms) if i not in membrane_atom_indices ]
    # Convert atom indices to a ndx file gromacs can read
    selection = Selection(atom_indices)
    selection_name = 'unmembraned_selection'
    ndx_atom_indices = selection.to_ndx(selection_name)
    with open(index_filename, 'w') as file:
        file.write(ndx_atom_indices)  
    # Filter the trajectory
    xtc_filter(input_topology_filename, input_trajectory_filename, output_trajectory_filename, index_filename)
    # Filter the topology
    pdb_filter(input_topology_filename, output_topology_filename, index_filename)