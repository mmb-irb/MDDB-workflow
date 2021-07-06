# This script is used to get standarized and corrected topology and trajectory files

import os

from model_workflow.tools.gromacs_processor import tpr2pdb, merge_xtc_files, get_first_frame
from model_workflow.tools.mdtraj_processor import merge_and_convert_traj
from model_workflow.tools.vmd_processor import vmd_processor, psf_to_pdb
from model_workflow.tools.image_and_fit import image_and_fit
from model_workflow.tools.topology_corrector import topology_corrector

def is_pdb (filename : str) -> bool:
    return filename[-4:] == '.pdb'

def is_psf (filename : str) -> bool:
    return filename[-4:] == '.psf'

def is_tpr (filename : str) -> bool:
    return filename[-4:] == '.tpr'

def is_xtc (filename : str) -> bool:
    return filename[-4:] == '.xtc'

def is_dcd (filename : str) -> bool:
    return filename[-4:] == '.dcd'

def is_netcdf (filename : str) -> bool:
    return filename[-3:] == '.nc'

def are_xtc (filenames : list) -> bool:
    return all([ is_xtc(filename) for filename in filenames ])

def are_dcd (filenames : list) -> bool:
    return all([ is_dcd(filename) for filename in filenames ])

def are_netcdf (filenames : list) -> bool:
    return all([ is_netcdf(filename) for filename in filenames ])

def process_topology_and_trajectory (
    input_topology_filename : str,
    input_trajectory_filenames : str,
    output_topology_filename : str,
    output_trajectory_filename : str,
    preprocess_protocol : int
    ) -> None:

    # In case the topology is a 'pdb' file just rename it
    # This is done in first place since VMD does not handle tpr files
    if not os.path.exists(output_topology_filename) and is_pdb(input_topology_filename):
        os.rename(input_topology_filename, output_topology_filename)
        input_topology_filename = output_topology_filename

    # In case the topology is a 'tpr' file convert it to pdb using gromacs
    # This is done in first place since VMD does not handle tpr files
    if not os.path.exists(output_topology_filename) and is_tpr(input_topology_filename):
        pdb_filename = input_topology_filename[:-4] + '.pdb'
        tpr2pdb(input_topology_filename, input_trajectory_filenames[0], pdb_filename)
        input_topology_filename = pdb_filename

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

    # Image the trajectory if it is required
    # i.e. make the trajectory uniform avoiding atom jumps and making molecules to stay whole
    # Fit the trajectory by removing the translation and rotation if it is required
    if preprocess_protocol > 0:
        image_and_fit(output_topology_filename, output_trajectory_filename,
                      output_topology_filename, output_trajectory_filename, 
                      preprocess_protocol)

    # Examine and correct the topology file using ProDy
    topology_corrector(output_topology_filename, output_topology_filename)