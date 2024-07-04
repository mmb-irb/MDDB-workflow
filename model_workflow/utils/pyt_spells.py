import pytraj as pyt
from distutils.version import StrictVersion
from typing import Optional

from model_workflow.utils.auxiliar import InputError
from model_workflow.utils.selections import Selection

# Set pytraj supported formats
pytraj_supported_structure_formats = {'prmtop', 'pdb', 'parm7', 'mol2', 'psf', 'cif', 'top', 'sdf'}
pytraj_supported_trajectory_formats = {'xtc', 'trr', 'crd', 'mdcrd', 'nc', 'dcd'}

# Get the whole trajectory as a generator
def get_pytraj_trajectory (
    input_topology_filename : str,
    input_trajectory_filename : str,
    atom_selection : Optional['Selection'] = None):

    # Topology is mandatory to setup the pytraj trajectory
    if not input_topology_filename:
        raise SystemExit('Missing topology file to setup PyTraj trajectory')
    
    # Set the pytraj trayectory and get the number of frames
    # NEVER FORGET: The pytraj iterload does not accept a mask, but we can strip atomos later
    pyt_trajectory = pyt.iterload(input_trajectory_filename, input_topology_filename)
    
    # WARNING: This extra line prevents the error "Segment violation (core dumped)" in some pdbs
    # This happens with some random pdbs which pytraj considers to have 0 Mols
    # More info: https://github.com/Amber-MD/cpptraj/pull/820
    # DANI: Esto es útil en pytraj <= 2.0.5 pero hace fallar el código a partir de pytraj 2.0.6
    if StrictVersion(pyt.__version__) <= StrictVersion('2.0.5'):
        pyt_trajectory.top.start_new_mol()

    # Filter away atoms which ARE NOT in the atom selection
    if atom_selection:
        # Get atom indices for all atoms but the ones in the atom selection
        topology = pyt_trajectory.top
        all_atoms = set(range(topology.n_atoms))
        keep_atoms = set(atom_selection.atom_indices)
        strip_atoms = all_atoms - keep_atoms
        if len(strip_atoms) > 0:
            # Convert the strip atom indices to a pytraj mask string and strip the iterator
            mask = Selection(strip_atoms).to_pytraj()
            pyt_trajectory = pyt_trajectory.strip(mask)

    return pyt_trajectory

# Get the trajectory frames number
# LORE: This was tried also with mdtraj's iterload but pytraj was way faster
def get_frames_count (
    input_topology_file : 'File',
    input_trajectory_file : 'File') -> int:
    
    print('Counting number of frames...')

    if not input_trajectory_file.exists:
        raise InputError('Missing trajectroy file when counting frames: ' + input_trajectory_file.path)
    
    if not input_topology_file.exists:
        raise InputError('Missing topology file when counting frames: ' + input_topology_file.path)

    # Load the trajectory from pytraj
    pyt_trajectory = pyt.iterload(
        input_trajectory_file.path,
        input_topology_file.path)

    # Return the frames number
    frames = pyt_trajectory.n_frames
    print(' Frames: ' + str(frames))

    # If 0 frames were counted then there is something wrong with the file
    if frames == 0:
        raise InputError('Something went wrong when reading the trajectory')

    return frames
# Set function supported formats
get_frames_count.format_sets = [
    {
        'inputs': {
            'input_structure_filename': pytraj_supported_structure_formats,
            'input_trajectory_filename': pytraj_supported_trajectory_formats
        }
    }
]