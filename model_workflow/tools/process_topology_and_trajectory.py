# This script is used to get standarized and corrected topology and trajectory files

import os

from model_workflow.tools.vmd_processor import vmd_processor
from model_workflow.tools.image_and_fit import image_and_fit
from model_workflow.tools.topology_corrector import topology_corrector

def process_topology_and_trajectory (
    input_topology_filename : str,
    input_trajectory_filenames : str,
    output_topology_filename : str,
    output_trajectory_filename : str,
    preprocess_protocol : int
    ) -> None:

    # Process the topology and or trajectory files using VMD
    # Files are converted to supported formats and trajectory pieces are merged into a single file
    # In addition, some irregularities in the topology may be fixed by VMD
    # If the output topology and trajectory files already exists it is assumed they are already processed
    if not os.path.exists(input_topology_filename) or not os.path.exists(input_trajectory_filenames):
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
                      output_trajectory_filename, preprocess_protocol)

    # Examine and correct the topology file using ProDy
    topology_corrector(output_topology_filename, output_topology_filename)