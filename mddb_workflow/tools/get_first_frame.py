from os import remove

from mddb_workflow.utils.auxiliar import get_auxiliar_filepath
from mddb_workflow.utils.gmx_spells import get_first_frame_structure
from mddb_workflow.utils.type_hints import *

# DANI: No lo muevo a gmx spells porque allí ya hay un get_first_frame
# DANI: La otra función también tiene inputs con nombres estandares (para el converter)
def get_first_frame (
    structure : 'Structure',
    trajectory_file : 'File',
    output_file : 'File',
    ):
    """Get the trajectory first frame in PDB format using Gromacs."""
    # Gromacs can not read mmCIF so we give it a PDB
    # Note that the first frame is used by tools which do not care about chains
    gromacs_structure_filepath = get_auxiliar_filepath('.first_frame_structure.pdb')
    structure.generate_pdb_file(gromacs_structure_filepath, show_warnings=False, accept_limitations=True)
    get_first_frame_structure(gromacs_structure_filepath, trajectory_file.path, output_file.path)
    remove(gromacs_structure_filepath)
