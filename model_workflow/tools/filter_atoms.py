from subprocess import run, PIPE, Popen

import pytraj as pt

not_water = '!(:SOL,WAT)'
not_ions = '!(@CL,Cl,Cl-,NA,Na,Na+,K,ZN,Zn)'
filter_base = not_water + '&' + not_ions

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

    # Load the charges topology
    charges_topology = pt.load_topology(filename=charges_filename)
    charges_atoms_count = charges_topology.n_atoms

    # Use pytraj to filter the charges topology
    # WARNING: It is saved in prmtop format by default
    filter_string = filter_base
    for exception in exceptions:
        filter_string += '|' + exception.selection

    # Set the filtered topology
    filtered_topology = topology[filter_string]
    filtered_atoms_count = filtered_topology.n_atoms

    # Set the filtered charges topology
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
        filtered_trajectory = trajectory[filter_string]
        pt.write_traj(trajectory_filename, filtered_trajectory)
        pt.write_parm(
            filename=topology_filename,
            top=filtered_topology,
            #format='infer',
            overwrite=True
        )
    if filtered_charges_atoms_count < charges_atoms_count:    
        pt.write_parm(
            filename=charges_filename,
            top=filtered_charges_topology,
            #format='infer',
            overwrite=True
        )


# --------------------------------------------------------------------------------------------

# DANI: No se usa
def gromacs_filter(
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_topology_filename : str,
    output_trajectory_filename : str
):
    # First filter the base topology/structure (pdb) and the trajectory (xtc)

    # Create indexes file to select only specific topology regions
    indexes = 'indexes.ndx'
    p = Popen([
        "echo",
        "!\"Water_and_ions\"\nq",
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "make_ndx",
        "-f",
        input_topology_filename,
        '-o',
        indexes,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

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
        "!Water_and_ions",
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
        indexes,
        '-dump',
        '0'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

    # Filter the trajectory
    p = Popen([
        "echo",
        "!Water_and_ions",
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
        indexes
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

    # Remove the reference topology file
    run([
        "rm",
        reference_topology
    ], stdout=PIPE).stdout.decode()