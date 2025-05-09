from model_workflow.utils.auxiliar import InputError, warn, CaptureOutput, load_json, MISSING_TOPOLOGY
from model_workflow.utils.constants import STANDARD_TOPOLOGY_FILENAME, GROMACS_EXECUTABLE
from model_workflow.utils.pyt_spells import find_first_corrupted_frame
from model_workflow.utils.gmx_spells import mine_system_atoms_count
from model_workflow.utils.vmd_spells import vmd_to_pdb
from model_workflow.utils.structures import Structure
from model_workflow.utils.file import File

from re import match, search
from typing import *
from subprocess import run, PIPE, Popen
from scipy.io import netcdf_file
import mdtraj as mdt
import pytraj as pyt

# Set some known message errors
NETCDF_DTYPE_ERROR = 'When changing to a larger dtype, its size must be a divisor of the total size in bytes of the last axis of the array.'
MDTRAJ_ATOM_MISMATCH_ERROR = r'xyz must be shape \(Any, ([0-9]*), 3\). You supplied  \(1, ([0-9]*), 3\)'
PYTRAJ_XTC_ATOM_MISMATCH_ERROR = r'Error: # atoms in XTC file \(([0-9]*)\) does not match # atoms in (topology|parm) [\w.-]* \(([0-9]*)\)'
GROMACS_ATOM_MISMATCH_ERROR = r'is larger than the number of atoms in the\ntrajectory file \(([0-9]*)\). There is a mismatch in the contents'

# List supported formats
TOPOLOGY_SUPPORTED_FORMATS = { 'tpr', 'top', 'prmtop', 'psf' }
TRAJECTORY_SUPPORTED_FORMATS = { 'xtc', 'trr', 'nc', 'dcd', 'crd', 'pdb', 'rst7' }
STRUCTURE_SUPPORTED_FORMATS = { *TOPOLOGY_SUPPORTED_FORMATS, 'pdb', 'gro' }
GROMACS_TRAJECTORY_SUPPORTED_FORMATS = { 'xtc', 'trr'}

# Auxiliar PDB file which may be generated to load non supported restart files
AUXILIAR_PDB_BILE = '.auxiliar.pdb'

# Check input files coherence and intergrity
# If there is any problem then raise an input error
def check_inputs (input_structure_file : 'File', input_trajectory_files : List['File'], input_topology_file : Union['File', Exception]):
    
    # Get a sample trajectory file and then check its format
    # All input trajectory files must have the same format
    trajectory_sample = input_trajectory_files[0]

    # Check input files are supported by the workflow
    if input_topology_file != MISSING_TOPOLOGY and input_topology_file.filename != STANDARD_TOPOLOGY_FILENAME and input_topology_file.format not in TOPOLOGY_SUPPORTED_FORMATS:
        if input_topology_file.format in { 'pdb', 'gro' }:
            raise InputError('A structure file is not supported as topology anymore. If there is no topology then use the argument "-top no"')
        raise InputError(f'Topology {input_topology_file.path} has a not supported format. Try one of these: {", ".join(TOPOLOGY_SUPPORTED_FORMATS)}')
    if trajectory_sample.format not in TRAJECTORY_SUPPORTED_FORMATS:
        raise InputError(f'Trajectory {trajectory_sample.path} has a not supported format. Try one of these: {", ".join(TRAJECTORY_SUPPORTED_FORMATS)}')
    if input_structure_file.format not in STRUCTURE_SUPPORTED_FORMATS:
        raise InputError(f'Structure {input_structure_file.path} has a not supported format. Try one of these: {", ".join(STRUCTURE_SUPPORTED_FORMATS)}')

    # Make sure the trajectory file is not corrupted

    # Check if reading the trajectory raises the following error
    # ValueError: When changing to a larger dtype, its size must be a divisor of the total size in bytes of the last axis of the array.
    # This error may happen with NetCDF files and it is a bit shady
    # Some tools may be able to read the first frames of the corrupted file: VMD and pytraj
    # Some other tools will instantly fail to read it: MDtraj and MDAnalysis
    if trajectory_sample.format == 'nc':
        try:
            # Iterate trajectory files
            for trajectory_file in input_trajectory_files:
                # This does not read the whole trajectory
                netcdf_file(trajectory_file.path, 'r')
        except Exception as error:
            # If the error message matches with a known error then report the problem
            error_message = str(error)
            if error_message == NETCDF_DTYPE_ERROR:
                warn(f'Corrupted trajectory file {trajectory_file.path}')
                pytraj_input_topology = input_topology_file if input_topology_file != MISSING_TOPOLOGY else input_structure_file
                first_corrupted_frame = find_first_corrupted_frame(pytraj_input_topology.path, trajectory_file.path)
                print(f' However some tools may be able to read the first {first_corrupted_frame} frames: VMD and PyTraj')
                raise InputError('Corrupted input trajectory file')
            # If we do not know the error then raise it as is
            else:
                raise error

    # Set a function to get atoms from a structure alone
    def get_structure_atoms (structure_file : 'File') -> int:
        # Get the number of atoms in the input structure
        structure = Structure.from_file(structure_file.path)
        return structure.atom_count

    # Set a function to get atoms from structure and trajectory together
    def get_structure_and_trajectory_atoms (structure_file : 'File', trajectory_file : 'File') -> Tuple[int, int]:
        # Note that declaring the iterator will not fail even when there is a mismatch
        trajectory = mdt.iterload(trajectory_file.path, top=structure_file.path, chunk=1)
        # We must consume the generator first value to make the error raise
        frame = next(trajectory)
        # Now obtain the number of atoms from the frame we just read
        trajectory_atom_count = frame.n_atoms
        # And still, it may happen that the topology has more atoms than the trajectory but it loads
        # MDtraj may silently load as many coordinates as possible and discard the rest of atoms in topology
        # This behaviour has been observed with a gromacs .top topology and a PDB used as trajectory
        # Two double check the match, load the topology alone with PyTraj
        topology = pyt.load_topology(structure_file.path)
        structure_atom_count = topology.n_atoms
        return structure_atom_count, trajectory_atom_count
    
    # Set a function to get atoms from topology and trajectory together
    def get_topology_and_trajectory_atoms (topology_file : 'File', trajectory_file : 'File') -> Tuple[int, int]:
        # To do so rely on different tools depending on the topology format
        # If there is no topology file then just compare strucutre and trajectory an exit
        if topology_file == MISSING_TOPOLOGY:
            # We do not have a topology atom count to return
            # Without a valid topology we can not count trajectory atoms either
            return None, None
        # If it is our standard topology then simply count the atom names
        # Get trajectory atoms using the structure instead
        if topology_file.filename == STANDARD_TOPOLOGY_FILENAME:
            # Parse the json and count atoms
            parsed_topology = load_json(topology_file.path)
            topology_atom_count = len(parsed_topology['atom_names'])
            # Without a valid topology we can not count trajectory atoms
            return topology_atom_count, None
        # For a TPR use Gromacs, which is its native tool
        if topology_file.format == 'tpr':
            # Make sure the trajectory is compatible with gromacs
            if trajectory_file.format not in GROMACS_TRAJECTORY_SUPPORTED_FORMATS:
                raise InputError('Why loading a TPR topology with a non-gromacs trajectory?')
            # Run Gromacs just to generate a structure using all atoms in the topology and coordinates in the first frame
            # If atoms do not match then we will see a specific error
            output_sample_file = File('.sample.gro')
            p = Popen([ "echo", "System" ], stdout=PIPE)
            process = run([
                GROMACS_EXECUTABLE,
                "trjconv",
                "-s",
                topology_file.path,
                "-f",
                trajectory_file.path,
                '-o',
                output_sample_file.path,
                "-dump",
                "0",
                '-quiet'
            ], stdin=p.stdout, stdout=PIPE, stderr=PIPE)
            logs = process.stdout.decode()
            p.stdout.close()
            # Always get error logs and mine topology atoms
            # Note that this logs include the output selection request from Gromacs
            # This log should be always there, even if there was a mismatch and then Gromacs failed
            error_logs = process.stderr.decode()
            topology_atom_count = mine_system_atoms_count(error_logs)
            # If the output does not exist at this point it means something went wrong with gromacs
            if not output_sample_file.exists:
                # Check if we know the error
                error_match = search(GROMACS_ATOM_MISMATCH_ERROR, error_logs)
                if error_match:
                    # Get the trajectory atom count
                    trajectory_atom_count = error_match[1]
                    return topology_atom_count, trajectory_atom_count
                # Otherwise just print the whole error logs and stop here anyway
                print(logs)
                print(error_logs)
                raise SystemExit('Something went wrong with GROMACS during the checking')
            # If we had an output then it means both topology and trajectory match in the number of atoms
            # Cleanup the file we just created and proceed
            output_sample_file.remove()
            # If there was no problem then it means the trajectory atom count matches the topology atom count
            trajectory_atom_count = topology_atom_count
            return topology_atom_count, trajectory_atom_count
        # For .top files we use PyTraj since MDtraj can not handle it
        if topology_file.format == 'top':
            # Note that calling ierload will print a error log when atoms do not match but will not raise a proper error
            # To capture the error log we must throw this command wrapped in a stdout redirect
            trajectory = None
            with CaptureOutput('stderr') as output:
                trajectory = pyt.iterload(trajectory_file.path, top=topology_file.path)
            logs = output.captured_text
            error_match = match(PYTRAJ_XTC_ATOM_MISMATCH_ERROR, logs)
            if error_match:
                topology_atom_count = error_match[3]
                trajectory_atom_count = error_match[1]
            # Now obtain the number of atoms from the frame we just read
            else:
                topology_atom_count = trajectory_atom_count = trajectory.n_atoms
            return topology_atom_count, trajectory_atom_count
        # At this point the topology should be supported by MDtraj
        # However, f the trajectory is a restart file MDtraj will not be able to read it
        # Make the conversion here, since restart files are single-frame trajectories this should be fast
        use_auxiliar_pdb = False
        if trajectory_file.format == 'rst7':
            # Generate the auxiliar PDB file
            vmd_to_pdb(topology_file.path, trajectory_file.path, AUXILIAR_PDB_BILE)
            use_auxiliar_pdb = True
        # For any other format use MDtraj
        try:
            # Note that declaring the iterator will not fail even when there is a mismatch
            trajectory_path = AUXILIAR_PDB_BILE if use_auxiliar_pdb else trajectory_file.path
            trajectory = mdt.iterload(trajectory_path, top=topology_file.path, chunk=1)
            # We must consume the generator first value to make the error raise
            frame = next(trajectory)
            # Now obtain the number of atoms from the frame we just read
            trajectory_atom_count = frame.n_atoms
            # And still, it may happen that the topology has more atoms than the trajectory but it loads
            # MDtraj may silently load as many coordinates as possible and discard the rest of atoms in topology
            # This behaviour has been observed with a gromacs .top topology and a PDB used as trajectory
            # Two double check the match, load the topology alone with PyTraj
            topology = pyt.load_topology(topology_file.path)
            topology_atom_count = topology.n_atoms
            return topology_atom_count, trajectory_atom_count
        except Exception as error:
            # If the error message matches with a known error then report the problem
            error_message = str(error)
            error_match = match(MDTRAJ_ATOM_MISMATCH_ERROR, error_message)
            if error_match:
                topology_atom_count = error_match[1]
                trajectory_atom_count = error_match[2]
                return topology_atom_count, trajectory_atom_count
            # If we do not know the error then raise it as is
            else:
                raise error
        
    # Get topology and trajectory atom counts
    topology_atom_count, trajectory_atom_count = get_topology_and_trajectory_atoms(input_topology_file, trajectory_sample)

    # If we have the trajectory atom count then it means we had a valid topology
    if trajectory_atom_count != None:

        # Make sure their atom counts match
        if topology_atom_count != trajectory_atom_count:
            raise InputError('Mismatch in the number of atoms between input files:\n' +
                f' Topology "{input_topology_file.path}" -> {topology_atom_count} atoms\n' +
                f' Trajectory "{trajectory_sample.path}" -> {trajectory_atom_count} atoms')
        
        # If the topology file is already the structure file then there is no need to check it
        if input_structure_file == input_topology_file:
            print(f'Topology and trajectory files match in number of atoms: {trajectory_atom_count}')
            return

        # If the counts match then also get the structure atom count and compare
        structure_atom_count = get_structure_atoms(input_structure_file)

        # Make sure it matches the topology and trajectory atom count
        if topology_atom_count != structure_atom_count:
            raise InputError('Mismatch in the structure input file number of atoms:\n'+
                f' Topology and trajectory -> {topology_atom_count} atoms\n' +
                f' Structure "{input_structure_file.path}" -> {structure_atom_count} atoms')
        
        # If we reached this point then it means everything is matching
        print(f'All input files match in number of atoms: {trajectory_atom_count}')
        return
        
    # Otherwise it means we had not a valid topology file
    # We must use the structure to find trajectory atoms
    structure_atom_count, trajectory_atom_count = get_structure_and_trajectory_atoms(input_structure_file, trajectory_sample)

    # Make sure their atom counts match
    if structure_atom_count != trajectory_atom_count:
        raise InputError('Mismatch in the number of atoms between input files:\n' +
            f' Structure "{input_structure_file.path}" -> {structure_atom_count} atoms\n' +
            f' Trajectory "{trajectory_sample.path}" -> {trajectory_atom_count} atoms')
            
    # If we have a number of topology atoms then make sure it matches the structure and trajectory atoms
    # This may happen if the topology is our standard topology file instead of a valid topology
    if topology_atom_count != None and topology_atom_count != trajectory_atom_count:
        raise InputError('Mismatch in the number of atoms between input files:\n' +
            f' Structure and trajectory -> {trajectory_atom_count} atoms\n' +
            f' Topology "{input_topology_file.path}" -> {topology_atom_count} atoms')
        
    # If we made it this far it means all checkings are good
    print(f'Input files match in number of atoms: {trajectory_atom_count}')