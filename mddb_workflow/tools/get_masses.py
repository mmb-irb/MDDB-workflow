import pytraj as pt

from mddb_workflow.utils.auxiliar import MISSING_TOPOLOGY
from mddb_workflow.utils.gmx_spells import read_and_parse_tpr
from mddb_workflow.utils.type_hints import *

# Topology formats which include atom masses and can be read through pytraj
# Note that PDB files are excluded since pytraj guesses their masses from the elements
PYTRAJ_MASSES_FORMATS = {'prmtop', 'top', 'psf'}


def get_masses(topology_file: Union['File', Exception]) -> Optional[list[float]]:
    """Get atom masses from the topology.
    Return None if there is no topology or its format has no masses (e.g. the standard topology).
    """
    if topology_file == MISSING_TOPOLOGY or not topology_file.exists: return None
    if topology_file.format in PYTRAJ_MASSES_FORMATS:
        topology = pt.load_topology(filename=topology_file.path)
        return [ float(mass) for mass in topology.mass ]
    if topology_file.format == 'tpr':
        parsed_tpr = read_and_parse_tpr(topology_file.path)
        return parsed_tpr['gromacs-topology']['mass']['value']
    return None
