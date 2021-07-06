from subprocess import run, PIPE, Popen

import pytraj as pt

from model_workflow.tools.get_charges import get_raw_charges, raw_charges_filename

# Set the name for the group name in gromacs ndx file
filter_group_name = "not_water_or_counter_ions"

# Set the pytraj selection for not water
not_water_mask = '!(:SOL,WAT,HOH)'

# Filter atoms of all input topologies by remvoing atoms and ions
# As an exception, some water and ions may be not removed if specified
# At the end, all topologies must match in atoms count
# DANI: Las exceptions no están del todo implementadas
# DANI: Falta implementarlas en el filtrado de topología y trayectoria
# DANI: Es decir, en el 'indexes'
def filter_atoms (
    topology_filename : str,
    trajectory_filename : str,
    charges_filename : str,
    exceptions : list
):    

    # Handle missing exceptions
    if not exceptions:
        exceptions = []

    # Load the topology and trajectory
    trajectory = pt.iterload(trajectory_filename, topology_filename)
    topology = trajectory.topology
    atoms_count = topology.n_atoms

    # Set if charges must be filtered
    # i.e. they are not a raw charges filename
    filtrable_charges = charges_filename and charges_filename != raw_charges_filename
    no_filtrable_charges = charges_filename and charges_filename == raw_charges_filename

    # Load the charges topology
    if filtrable_charges:
        charges_topology = pt.load_topology(filename=charges_filename)
        charges_atoms_count = charges_topology.n_atoms
    elif no_filtrable_charges:
        charges = get_raw_charges(charges_filename)
        charges_atoms_count = len(charges)

    # Set the pytraj mask to filter the desired atoms
    filter_string = not_water_mask
    counter_ions_filter = get_counter_ions_mask(topology_filename)
    if counter_ions_filter:
        filter_string = '(' + not_water_mask + '&!(' + counter_ions_filter + '))'
    for exception in exceptions:
        filter_string += '|' + exception['selection']

    # Set the filtered topology
    filtered_topology = topology[filter_string]
    filtered_atoms_count = filtered_topology.n_atoms

    # Set the filtered charges topology
    if filtrable_charges:
        filtered_charges_topology = charges_topology[filter_string]
        filtered_charges_atoms_count = filtered_charges_topology.n_atoms

        # Both filtered topologies must have the same number of atoms
        if filtered_atoms_count != filtered_charges_atoms_count:
            print('Base atoms: ' + str(filtered_atoms_count))
            print('Charges atoms: ' + str(filtered_charges_atoms_count))
            raise SystemExit('ERROR: Filtered atom counts in base and charges topologies does not match')

    # Check if both the normal and the filtered topologies have the same number of atoms
    # In not, filter the whole trajectory and overwrite both topologies and trajectory
    print('Total number of atoms: ' + str(atoms_count))
    print('Filtered number of atoms: ' + str(filtered_atoms_count))
    if filtered_atoms_count < atoms_count:
        print('Filtering structure and trajectory...')
        # Filter both topology and trajectory using Gromacs, since is more efficient than pytraj
        # Set up an 'index.ndx' file with all atom indices manually
        # As long as indices are for atoms and not residues there should never be any incompatibility
        indexes_filename = 'filter.ndx'
        pytraj_mask_2_gromacs_ndx(topology, { filter_group_name : filter_string }, indexes_filename)
        gromacs_filter(topology_filename, trajectory_filename, topology_filename, trajectory_filename, indexes_filename)
        #filtered_trajectory = trajectory[filter_string]
        #pt.write_traj(trajectory_filename, filtered_trajectory, overwrite=True)
        # pt.write_traj(
        #     filename=topology_filename,
        #     # DANI: No he encontrado otra manera de exportar a pdb con pytraj
        #     traj=filtered_trajectory[0:1],
        #     overwrite=True
        # )
    if charges_filename:
        print ('Total number of charges: ' + str(charges_atoms_count))
    if filtrable_charges and filtered_charges_atoms_count < charges_atoms_count:
        print('Filtering charges topology...')
        pt.write_parm(
            filename=charges_filename,
            top=filtered_charges_topology,
            format='amberparm',
            overwrite=True
        )
    # In case we have a raw charges file check the number of charges matches the number of fileterd atoms
    if no_filtrable_charges and filtered_charges_atoms_count < charges_atoms_count:
        if charges_count != filtered_atoms_count:
            raise SystemExit("Charges count in '" + raw_charges_filename + "' does not match the number of filtered atoms: " + str(charges_atoms_count))

# Get a pytraj selection with all counter ions
counter_ions = ['K', 'NA', 'CL']
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


# --------------------------------------------------------------------------------------------

# Filter atoms in both a pdb and a xtc file
def gromacs_filter(
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_topology_filename : str,
    output_trajectory_filename : str,
    indexes_filename : str
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
        indexes_filename,
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
        indexes_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

    # Remove the reference topology file
    run([
        "rm",
        reference_topology
    ], stdout=PIPE).stdout.decode()

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
