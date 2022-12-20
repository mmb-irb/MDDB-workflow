# This script is used to get standarized and corrected topology and trajectory files

import os
from typing import List

from model_workflow.tools.formats import is_psf, is_tpr, is_raw, is_prmtop, is_top, is_psf, raw_charges_filename
from model_workflow.tools.filter_atoms import filter_atoms
from model_workflow.tools.image_and_fit import image_and_fit
from model_workflow.tools.topology_corrector import topology_corrector

from mdtoolbelt.conversions import convert

# Process input files: topology, trajectory and charges
# These files must be processed all together since they depend on each other
# e.g. the topology file may require trajectory data to change the format
# e.g. charges are filtered according to atoms in topology
def process_input_files (
    input_topology_filename : str,
    input_trajectory_filenames : List[str],
    input_charges_filename : str,
    output_topology_filename : str,
    output_trajectory_filename : str,
    preprocess_protocol : int,
    translation: list,
    filter_selection : str
    ) -> None:

    # Convert input topology and trajectories to output topology and trajectory
    convert(
        input_structure_filename=input_topology_filename,
        output_structure_filename=output_topology_filename,
        input_trajectory_filenames=input_trajectory_filenames,
        output_trajectory_filename=output_trajectory_filename,
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
    filter_atoms(output_topology_filename, output_trajectory_filename, output_charges_filename, filter_selection)

    # Image the trajectory if it is required
    # i.e. make the trajectory uniform avoiding atom jumps and making molecules to stay whole
    # Fit the trajectory by removing the translation and rotation if it is required
    if preprocess_protocol > 0:
        image_and_fit(output_topology_filename, output_trajectory_filename, output_charges_filename,
                      output_topology_filename, output_trajectory_filename, preprocess_protocol, translation)

    # Examine and correct the topology file
    topology_corrector(output_topology_filename, output_topology_filename, output_trajectory_filename, output_trajectory_filename)

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