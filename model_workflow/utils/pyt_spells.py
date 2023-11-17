import pytraj as pyt
from distutils.version import StrictVersion

# Set pytraj supported formats
pytraj_supported_structure_formats = {'prmtop', 'pdb', 'parm7', 'mol2', 'psf', 'cif', 'top', 'sdf'}
pytraj_supported_trajectory_formats = {'xtc', 'trr', 'crd', 'mdcrd', 'nc', 'dcd'}

# Get the whole trajectory as a generator
def get_pytraj_trajectory (
    input_topology_filename : str,
    input_trajectory_filename : str):

    # Topology is mandatory to setup the pytraj trajectory
    if not input_topology_filename:
        raise SystemExit('Missing topology file to setup PyTraj trajectory')
    
    # Set the pytraj trayectory and get the number of frames
    pyt_trajectory = pyt.iterload(input_trajectory_filename, input_topology_filename)
    
    # WARNING: This extra line prevents the error "Segment violation (core dumped)" in some pdbs
    # This happens with some random pdbs which pytraj considers to have 0 Mols
    # More info: https://github.com/Amber-MD/cpptraj/pull/820
    # DANI: Esto es útil en pytraj <= 2.0.5 pero hace fallar el código a partir de pytraj 2.0.6
    if StrictVersion(pyt.__version__) <= StrictVersion('2.0.5'):
        pyt_trajectory.top.start_new_mol()

    return pyt_trajectory

# Get the trajectory frames number using pytraj
def get_frames_count (
    input_topology_filename : str,
    input_trajectory_filename : str) -> int:

    # Load the trajectory from pytraj
    pyt_trajectory = get_pytraj_trajectory(
        input_topology_filename,
        input_trajectory_filename)

    # Return the frames number
    return pyt_trajectory.n_frames
# Set function supported formats
get_frames_count.format_sets = [
    {
        'inputs': {
            'input_structure_filename': pytraj_supported_structure_formats,
            'input_trajectory_filename': pytraj_supported_trajectory_formats
        }
    }
]