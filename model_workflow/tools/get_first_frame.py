from model_workflow.utils.gmx_spells import get_first_frame as get_first_frame_gmx
from model_workflow.utils.type_hints import *

# DANI: No lo muevo a gmx spells porque allí ya hay un get_first_frame
# DANI: La otra función también tiene inputs con nombres estandares (para el converter)
def get_first_frame (
    structure_file : 'File',
    trajectory_file : 'File',
    output_filepath : str
    ):
    """Get the trajectory first frame in PDB format using Gromacs."""
    get_first_frame_gmx(structure_file.path, trajectory_file.path, output_filepath)