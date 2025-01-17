from model_workflow.utils.auxiliar import InputError, warn, CaptureOutput
from model_workflow.utils.constants import STANDARD_TOPOLOGY_FILENAME, GROMACS_EXECUTABLE
from model_workflow.utils.pyt_spells import find_first_corrupted_frame
from model_workflow.utils.structures import Structure
from model_workflow.utils.file import File

from re import match, search
from typing import List
from subprocess import run, PIPE, Popen
from scipy.io import netcdf_file
import mdtraj as mdt
import pytraj as pyt

# Set some known message errors
NETCDF_DTYPE_ERROR = 'When changing to a larger dtype, its size must be a divisor of the total size in bytes of the last axis of the array.'
MDTRAJ_ATOM_MISMATCH_ERROR = r'xyz must be shape \(Any, ([0-9]*), 3\). You supplied  \(1, ([0-9]*), 3\)'
PYTRAJ_XTC_ATOM_MISMATCH_ERROR = r'Error: # atoms in XTC file \(([0-9]*)\) does not match # atoms in (topology|parm) [\w.-]* \(([0-9]*)\)'
GROMACS_ATOM_MISMATCH_ERROR = r'is larger than the number of atoms in the\ntrajectory file \(([0-9]*)\). There is a mismatch in the contents'
GROMACS_SYSTEM_ATOMS = r'System\) has ([0-9]*) elements'

# List supported formats
TOPOLOGY_SUPPORTED_FORMATS = { 'tpr', 'top', 'prmtop', 'psf' }
TRAJECTORY_SUPPORTED_FORMATS = { 'xtc', 'trr', 'nc', 'dcd', 'mdcrd', 'pdb' }
STRUCTURE_SUPPORTED_FORMATS = { *TOPOLOGY_SUPPORTED_FORMATS, 'pdb', 'gro' }
GROMACS_TRAJECTORY_SUPPORTED_FORMATS = { 'xtc', 'trr'}

# Check input files coherence and intergrity
# If there is any problem then raise an input error
def check_inputs (input_structure_file : 'File', input_trajectory_files : List['File'], input_topology_file : 'File'):

    # Get a sample trajectory file and then check its format
    # All input trajectory files must have the same format
    trajectory_sample = input_trajectory_files[0]

    # Check input files are supported by the workflow
    if input_topology_file.filename != STANDARD_TOPOLOGY_FILENAME and input_topology_file.format not in TOPOLOGY_SUPPORTED_FORMATS:
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
                first_corrupted_frame = find_first_corrupted_frame(input_topology_file.path, trajectory_file.path)
                print(f' However some tools may be able to read the first {first_corrupted_frame} frames: VMD and PyTraj')
                raise InputError('Corrupted input trajectory file')
            # If we do not know the error then raise it as is
            else:
                raise error
            
    # Now make sure topology and trajectory match in number of atoms
    atom_count = None

    # To do so we will rely on MDtraj, which is able to read most formats
    # Also MDtraj would probably be the first library to read the input files further (conversion)
    # However MDtraj is not able to read TPR and of course it does not read our interal topology format
    if input_topology_file.filename == STANDARD_TOPOLOGY_FILENAME:
        print(f'We will skip the atom count matching check since we already have a standard {STANDARD_TOPOLOGY_FILENAME}')
        # DANI: Hay que hacer return aquí, porque sino luego el atom_count sigue siendo None y el checking del structure falla
        return
    elif input_topology_file.format == 'tpr':
        # Make sure the trajectory is compatible with gromacs
        if trajectory_sample.format not in GROMACS_TRAJECTORY_SUPPORTED_FORMATS:
            raise InputError('Why loading a TPR topology with a non-gromacs trajectory?')
        # Run Gromacs just to generate a structure using all atoms in the topology and coordinates in the first frame
        # If atoms do not match then we will see a specific error
        output_sample_file = File('.sample.gro')
        p = Popen([ "echo", "System" ], stdout=PIPE)
        process = run([
            GROMACS_EXECUTABLE,
            "trjconv",
            "-s",
            input_topology_file.path,
            "-f",
            trajectory_sample.path,
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
        system_atoms_match = search(GROMACS_SYSTEM_ATOMS, error_logs)
        if not system_atoms_match:
            raise ValueError('Failed to mine Gromacs error logs')
        atom_count = int(system_atoms_match[1])
        # If the output does not exist at this point it means something went wrong with gromacs
        if not output_sample_file.exists:
            # Check if we know the error
            error_match = search(GROMACS_ATOM_MISMATCH_ERROR, error_logs)
            if error_match:
                # Get the trajectory atom count
                trajectory_atom_count = error_match[1]
                raise InputError('Mismatch in the number of atoms between input files:\n' +
                    f' Topology "{input_topology_file.path}" -> {atom_count} atoms\n' +
                    f' Trajectory "{trajectory_sample.path}" -> {trajectory_atom_count} atoms')
            # Otherwise just print the whole error logs and stop here anyway
            print(logs)
            print(error_logs)
            raise SystemExit('Something went wrong with GROMACS during the checking')
        # If we had an output then it means both topology and trajectory match in the number of atoms
        # Cleanup the file we just created and proceed
        output_sample_file.remove()
    # For .top files we use PyTraj since MDtraj can not handle it
    elif input_topology_file.format == 'top':
        # Note that calling ierload will print a error log when atoms do not match but will not raise a proper error
        # To capture the error log we must throw this command wrapped in a stdout redirect
        trajectory = None
        with CaptureOutput('stderr') as output:
            trajectory = pyt.iterload(trajectory_sample.path, top=input_topology_file.path)
        logs = output.captured_text
        error_match = match(PYTRAJ_XTC_ATOM_MISMATCH_ERROR, logs)
        if error_match:
            topology_atoms = error_match[3]
            trajectory_atoms = error_match[1]
            raise InputError('Mismatch in the number of atoms between input files:\n' +
                f' Topology "{input_topology_file.path}" -> {topology_atoms} atoms\n' +
                f' Trajectory "{trajectory_sample.path}" -> {trajectory_atoms} atoms')
        # Now obtain the number of atoms from the frame we just read
        atom_count = trajectory.n_atoms
    # More any toher format use MDtraj
    else:
        try:
            # Note that declaring the iterator will not fail even when there is a mismatch
            trajectory = mdt.iterload(trajectory_sample.path, top=input_topology_file.path, chunk=1)
            # We must consume the generator first value to make the error raise
            frame = next(trajectory)
            # Now obtain the number of atoms from the frame we just read
            atom_count = frame.n_atoms
            # And still, it may happen that the topology has more atoms than the trajectory but it loads
            # MDtraj may silently load as many coordinates as possible and discard the rest of atoms in topology
            # This behaviour has been observed with a gromacs .top topology and a PDB used as trajectory
            # Two double check the match, load the topology alone with PyTraj
            topology = pyt.load_topology(input_topology_file.path)
            if topology.n_atoms != atom_count:
                raise InputError('Mismatch in the number of atoms between input files:\n' +
                    f' Topology "{input_topology_file.path}" -> {topology.n_atoms} atoms\n' +
                    f' Trajectory "{trajectory_sample.path}" -> {atom_count} atoms')
        except Exception as error:
            # If the error message matches with a known error then report the problem
            error_message = str(error)
            error_match = match(MDTRAJ_ATOM_MISMATCH_ERROR, error_message)
            if error_match:
                topology_atoms = error_match[1]
                trajectory_atoms = error_match[2]
                raise InputError('Mismatch in the number of atoms between input files:\n' +
                    f' Topology "{input_topology_file.path}" -> {topology_atoms} atoms\n' +
                    f' Trajectory "{trajectory_sample.path}" -> {trajectory_atoms} atoms')
            # If we do not know the error then raise it as is
            else:
                raise error
            
    # If we have an independent structure then check it also matches the number of atoms
    if input_structure_file != input_topology_file:
        # Get the number of atoms in the input structure
        structure = Structure.from_file(input_structure_file.path)
        # Make sure it matches
        if atom_count != structure.atom_count:
            raise InputError('Mismatch in the structure input file number of atoms:\n'+
                f' Topology and trajectory -> {atom_count} atoms\n' +
                f' Structure "{input_structure_file.path}" -> {structure.atom_count} atoms')
        
    # if we made it this far it means all checkings are good
    print(f'All input files match in number of atoms: {atom_count}')