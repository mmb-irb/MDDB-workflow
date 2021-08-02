import os
from subprocess import run, PIPE, Popen

import pytraj as pt

from model_workflow.tools.get_charges import get_raw_charges, raw_charges_filename, get_tpr_charges
from model_workflow.tools.formats import is_raw, is_pytraj_supported, is_tpr, get_pytraj_parm_format

# Set the gromacs indices filename
index_filename = 'filter.ndx'
# Set the name for the group name in gromacs ndx file
filter_group_name = "not_water_or_counter_ions"

# Set the pytraj selection for not water
water_mask = '(:SOL,WAT,HOH)'

# Filter atoms of all input topologies by remvoing atoms and ions
# As an exception, some water and ions may be not removed if specified
# At the end, all topologies must match in atoms count
def filter_atoms (
    topology_filename : str,
    trajectory_filename : str,
    charges_filename : str,
    exceptions : list
):    

    # Handle missing exceptions
    if not exceptions:
        exceptions = []

    # Set the pytraj mask to filter the desired atoms from a specific topology
    def get_filter_mask (source_topology_filename : str) -> str:
        filter_mask = '!' + water_mask
        counter_ions_mask = get_counter_ions_mask(source_topology_filename)
        if counter_ions_mask:
            filter_mask = '(!' + water_mask + '&!(' + counter_ions_mask + '))'
        for exception in exceptions:
            filter_mask += '|' + exception['selection']
        return filter_mask

    # Load the topology and trajectory
    trajectory = pt.iterload(trajectory_filename, topology_filename)
    topology = trajectory.topology
    atoms_count = topology.n_atoms

    # Set the pytraj mask to filter the desired atoms from the topology
    filter_mask = get_filter_mask(topology_filename)

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
        gromacs_filter(topology_filename, trajectory_filename, topology_filename, trajectory_filename, index_filename)

    # Filter charges according to the file format
    if charges_filename and os.path.exists(charges_filename):
        print('Processing charges')
        # Raw charges
        if is_raw(charges_filename):
            charges = get_raw_charges(charges_filename)
            # Nothing to do here. It better matches by defualt or we have a problem
            filtered_charges_atoms_count = len(charges)
        # Pytraj supported formats
        elif is_pytraj_supported(charges_filename):
            # Load the charges topology and count its atoms
            charges_topology = pt.load_topology(filename=charges_filename)
            charges_atoms_count = charges_topology.n_atoms
            # Filter the desired atoms using the mask and then count them
            filter_mask = get_filter_mask(charges_filename)
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
            # If the number of charges in greater than expected then filter the tpr file and extract charges again
            if charges_atoms_count > filtered_atoms_count:
                if not os.path.exists(index_filename):
                    # In order to filter the tpr we need the filter.ndx file
                    # This must be generated from a pytraj supported topology that matches the number of atoms in the tpr file
                    raise ValueError('Charges number does not match the structure atoms and tpr files can not be filtered alone')
                gromacs_tpr_filter(charges_filename, charges_filename, index_filename)
                charges = get_tpr_charges(charges_filename)
            filtered_charges_atoms_count = len(charges)
        else:
            raise ValueError('Charges file is in a non supported format')

        print('Filtered charges atoms: ' + str(filtered_charges_atoms_count))

        # Both filtered topology and charges must have the same number of atoms
        if filtered_atoms_count != filtered_charges_atoms_count:
            print('Filtered topology atoms: ' + str(filtered_atoms_count))
            raise ValueError('Filtered atom counts in topology and charges does not match')

    # Remove the index file in case it was created
    if os.path.exists(index_filename):
        run([
            "rm",
            index_filename,
        ], stdout=PIPE).stdout.decode()

# Get a pytraj selection with all counter ions
counter_ions = ['K', 'NA', 'CL', 'CLA', 'SOD']
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
# Masks is a dict where each entry key will be the gromacs group name and each entry value is the mask string
# e.g. { "my_solvent" : ":SOL,WAT,HOH", "my_ions" : "@17,18" }
def pytraj_mask_2_gromacs_ndx (topology : 'Topology', masks : dict, output_filename : str = 'index.ndx'):
    with open(output_filename, "w") as file:
        for group_name, mask in masks.items():
            # Find mask atom indices in the topology
            atom_indices = list(topology.atom_indices(mask))
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

# Filter atoms in both a pdb and a xtc file
def gromacs_filter(
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_topology_filename : str,
    output_trajectory_filename : str,
    index_filename : str
):
    # First filter the base topology/structure (pdb) and then the trajectory (xtc)
    # Copy the original topology since we need it later to filter the trajectory
    # The original topology could be overwritten
    reference_topology = 'reference.topology.pdb'
    run([
        "cp",
        input_topology_filename,
        reference_topology
    ], stdout=PIPE).stdout.decode()

    # Filter the topology
    p = Popen([
        "echo",
        filter_group_name,
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "trjconv",
        "-s",
        reference_topology,
        "-f",
        input_trajectory_filename,
        '-o',
        output_topology_filename,
        '-n',
        index_filename,
        '-dump',
        '0',
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

    # Filter the trajectory
    p = Popen([
        "echo",
        filter_group_name,
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "trjconv",
        "-s",
        reference_topology,
        "-f",
        input_trajectory_filename,
        '-o',
        output_trajectory_filename,
        '-n',
        index_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

    # Remove the reference topology file
    run([
        "rm",
        reference_topology
    ], stdout=PIPE).stdout.decode()

# Filter atoms in both a pdb and a xtc file
def gromacs_tpr_filter(
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