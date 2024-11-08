from model_workflow.utils.auxiliar import InputError, warn
from model_workflow.utils.constants import TOPOLOGY_FILENAME
from model_workflow.utils.pyt_spells import find_first_corrupted_frame
from model_workflow.utils.structures import Structure

from re import match
from typing import List
from scipy.io import netcdf_file
import mdtraj as mdt

# Set some known message errors
NETCDF_DTYPE_ERROR = 'When changing to a larger dtype, its size must be a divisor of the total size in bytes of the last axis of the array.'
ATOM_MISMATCH_ERROR = r'xyz must be shape \(Any, ([0-9]*), 3\). You supplied  \(1, ([0-9]*), 3\)'

# Check input files coherence and intergrity
# If there is any problem then raise an input error
def check_inputs (input_structure_file : 'File', input_trajectory_files : List['File'], input_topology_file : 'File'):

    # Make sure the trajectory file is not corrupted

    # Check if reading the trajectory raises the following error
    # ValueError: When changing to a larger dtype, its size must be a divisor of the total size in bytes of the last axis of the array.
    # This error may happen with NetCDF files and it is a bit shady
    # Some tools may be able to read the first frames of the corrupted file: VMD and pytraj
    # Some other tools will instantly fail to read it: MDtraj and MDAnalysis

    # Get a sample trajectory file and then check its format
    # All input trajectory files must have the same format
    trajectory_sample = input_trajectory_files[0]
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
    if input_topology_file.filename == TOPOLOGY_FILENAME:
        print(f'We will skip the atom count matching check since we already have a standard {TOPOLOGY_FILENAME}')
        # DANI: Hay que hacer return aquí, porque sino luego el atom_count sigue siendo None y el checking del structure falla
        return
    elif input_topology_file.format == 'tpr':
        print('We will skip the atom count matching check since it is not yet implemented for TPR')
        # DANI: Hay que hacer return aquí, porque sino luego el atom_count sigue siendo None y el checking del structure falla
        return
    else:
        try:
            # Note that declaring the iterator will not fail even when there is a mismatch
            trajectory = mdt.iterload(trajectory_sample.path, top=input_topology_file.path, chunk=1)
            # We must consume the generator first value to make the error raise
            frame = next(trajectory)
            # Now obtain the number of atoms from the frame we just read
            atom_count = frame.n_atoms
        except Exception as error:
            # If the error message matches with a known error then report the problem
            error_message = str(error)
            error_match = match(ATOM_MISMATCH_ERROR, error_message)
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