# This script is used to get standarized and corrected topology and trajectory files

import os

from model_workflow.tools.formats import is_pdb, is_psf, is_tpr, is_gro, is_xtc, is_dcd, is_netcdf, are_xtc, are_dcd, are_netcdf, is_raw, is_prmtop, is_top, is_psf, raw_charges_filename
from model_workflow.tools.gromacs_processor import topology_to_pdb, merge_xtc_files, get_first_frame
from model_workflow.tools.mdtraj_processor import merge_and_convert_traj
from model_workflow.tools.vmd_processor import vmd_processor, psf_to_pdb
from model_workflow.tools.filter_atoms import filter_atoms
from model_workflow.tools.image_and_fit import image_and_fit
from model_workflow.tools.topology_corrector import topology_corrector

# Process input files: topology, trajectory and charges
# These files must be processed all together since they depend on each other
# e.g. the topology file may require trajectory data to change the format
# e.g. charges are filtered according to atoms in topology
def process_input_files (
    input_topology_filename : str,
    input_trajectory_filenames : str,
    input_charges_filename : str,
    output_topology_filename : str,
    output_trajectory_filename : str,
    preprocess_protocol : int,
    exceptions : list
    ) -> None:

    # In case the topology is a 'pdb' file just rename it
    # This is done in first place since VMD does not handle tpr files
    if not os.path.exists(output_topology_filename) and is_pdb(input_topology_filename):
        os.rename(input_topology_filename, output_topology_filename)
        input_topology_filename = output_topology_filename

    # In case the topology is a 'tpr' or a 'gro' file convert it to pdb using gromacs
    # This is done in first place since VMD does not handle tpr files
    if not os.path.exists(output_topology_filename) and ( is_tpr(input_topology_filename) or is_gro(input_topology_filename)):
        topology_to_pdb(input_topology_filename, input_trajectory_filenames[0], output_topology_filename)
        input_topology_filename = output_topology_filename

    # In case the trajectory is a single 'xtc' file just rename it
    if (not os.path.exists(output_trajectory_filename)
        and len(input_trajectory_filenames) == 1
        and is_xtc(input_trajectory_filenames[0])
    ):
        os.rename(input_trajectory_filenames[0], output_trajectory_filename)
        input_trajectory_filenames = output_trajectory_filename

    # In case the trajectory files are all xtc merge them using Gromacs trjcat
    if not os.path.exists(output_trajectory_filename) and are_xtc(input_trajectory_filenames):
        merge_xtc_files(input_trajectory_filenames, output_trajectory_filename)

    # In case the trajectory files are all dcd merge them using MDtraj mdconvert
    if not os.path.exists(output_trajectory_filename) and are_dcd(input_trajectory_filenames):
        merge_and_convert_traj(input_trajectory_filenames, output_trajectory_filename)

    # In case the trajectory files are all netcdf merge them using MDtraj mdconvert
    if not os.path.exists(output_trajectory_filename) and are_netcdf(input_trajectory_filenames):
        merge_and_convert_traj(input_trajectory_filenames, output_trajectory_filename)

    # In case the topology is a psf file generate a pdb file from it using the first trajectory frame
    if not os.path.exists(output_topology_filename) and is_psf(input_topology_filename):
        if not os.path.exists(output_trajectory_filename):
            # DANI: Si tienes este problema habrá que implementar más cosas para resolverlo
            raise SystemExit('ERROR: I cannot convert psf to pdb if trajectories are not xtc or dcd')
        single_frame_filename = 'frame.xtc'
        get_first_frame(output_trajectory_filename, single_frame_filename)
        psf_to_pdb(input_topology_filename, single_frame_filename, output_topology_filename)

    # Process the topology and or trajectory files using VMD
    # Files are converted to supported formats and trajectory pieces are merged into a single file
    # In addition, some irregularities in the topology may be fixed by VMD
    # If the output topology and trajectory files already exists it is assumed they are already processed
    # WARNING: This is an emergency endpoint. VMD will "work" in many cases but it makes your RAM explode
    if not os.path.exists(output_topology_filename) or not os.path.exists(output_trajectory_filename):
        logs = vmd_processor(
            input_topology_filename,
            input_trajectory_filenames,
            output_topology_filename,
            output_trajectory_filename,
        )

    # Process charges file
    # Rename it according to standards
    output_charges_filename = None
    if input_charges_filename and os.path.exists(input_charges_filename):
        # If it is a 'charges.txt' then we do not have to rename it
        if is_raw(input_charges_filename):
            output_charges_filename = input_charges_filename
        # In other cases rename the charges file
        elif is_prmtop(input_charges_filename):
            output_charges_filename = 'topology.prmtop'
        elif is_top(input_charges_filename):
            output_charges_filename = 'topology.top'
        elif is_psf(input_charges_filename):
            output_charges_filename = 'topology.psf'
        elif is_tpr(input_charges_filename):
            output_charges_filename = 'topology.tpr'
        else:
            raise ValueError('Charges file is in a non supported format')
        
        if not os.path.exists(output_charges_filename):
            os.rename(input_charges_filename, output_charges_filename)

    # Filter atoms in bot topology and trajectory in order to remove solvent and counter ions
    filter_atoms(output_topology_filename, output_trajectory_filename, output_charges_filename, exceptions)

    # Image the trajectory if it is required
    # i.e. make the trajectory uniform avoiding atom jumps and making molecules to stay whole
    # Fit the trajectory by removing the translation and rotation if it is required
    if preprocess_protocol > 0:
        image_and_fit(output_topology_filename, output_trajectory_filename,
                      output_topology_filename, output_trajectory_filename, 
                      preprocess_protocol)

    # Examine and correct the topology file using ProDy
    topology_corrector(output_topology_filename, output_topology_filename)

# Find out if there is any file with a supported charges format in the current directory
# In that case, return the first match
# Raw energies have priority
def get_output_charges_filename():

    # Iterate over all files in the current directory
    for filename in os.listdir('.'):
        if is_raw(filename):
            return filename
        elif is_prmtop(filename):
            return 'topology.prmtop'
        elif is_top(filename):
            return 'topology.top'
        elif is_psf(filename):
            return 'topology.psf'
        elif is_tpr(filename):
            return 'topology.tpr'
    return None

# Find out if there is any of the standard topology filenames in the current directory
# In that case, return the first standard filename found
# Raw energies have priority
def find_charges_filename():

    # Set the standard filenames
    standard_filenames = [
        raw_charges_filename,
        'topology.prmtop',
        'topology.top',
        'topology.psf',
        'topology.tpr',
    ]

    # Return the first existing standard filename
    for standard_filename in standard_filenames:
        if os.path.exists(standard_filename):
            return standard_filename

    return None