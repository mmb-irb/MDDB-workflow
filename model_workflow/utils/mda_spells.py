import json
import numpy as np
import pickle

from model_workflow.utils.type_hints import *
from model_workflow.utils.constants import GREY_HEADER, COLOR_END
from model_workflow.utils.auxiliar import MISSING_CHARGES

from MDAnalysis.topology.TPRParser import TPRParser
#from MDAnalysis.topology.PDBParser import PDBParser # for class reference
from MDAnalysis.core.universe import Universe
from MDAnalysis.core.topology import Topology
from MDAnalysis.core.topologyattrs import (
    Atomnames,
    Bonds,
    ChainIDs,
    Charges,
    Elements,
    Resids,
    Resnames,
    Resnums,
    Segids
)

def to_MDAnalysis_topology(standard_topology_path : str) -> 'Topology':
    """
    Creates a MDAnalysis topology from a json topology file.

    The json file should contain the following keys:
    - atom_names
    - atom_elements
    - atom_charges
    - atom_bonds (list of lists of atom indices)
    - residue_names
    - residue_numbers
    - chain_names
    - atom_residue_indices
    - residue_chain_indices

    :param topology: path to the json file
    :returns: a MDAnalysis topology object
    """
    topology = json.load(open(standard_topology_path))

    # transform bond to non redundant tuples
    bonds = []
    for bond_from, bond_tos in enumerate(topology['atom_bonds']):
        for bond_to in bond_tos:
            bond = tuple(sorted([bond_from, bond_to]))
            bonds.append(bond)
    # sorted(set(bonds)) if you want to check for duplicates
    topology['bonds'] = set(bonds)
    atom_2_chain = np.array(topology['residue_chain_indices'])[topology['atom_residue_indices']]
    topology['chainIDs'] = np.array(topology['chain_names'])[atom_2_chain]


    attrs = [Atomnames (topology['atom_names']),
             Elements  (topology['atom_elements']),
             Charges   (topology['atom_charges']),
             Bonds     (topology['bonds']),
             Resnames  (topology['residue_names']),
             Resnums   (topology['residue_numbers']),
             ChainIDs  (topology['chainIDs']),
             Segids    (topology['chain_names']),
             Resids    (np.arange(len(topology['residue_names']))),
    ]

    mda_top = Topology(n_atoms=len(topology['atom_names']), 
                       n_res=len(topology['residue_names']), 
                       n_seg=len(topology['chain_names']),
                       attrs=attrs,
                       atom_resindex=topology['atom_residue_indices'],
                       residue_segindex=topology['residue_chain_indices']
    )
    return mda_top

def get_mda_universe_from_stopology (
        standard_topology_path : str, coordinates_file : str) -> 'Universe':
    """Create a MDAnalysis universe using data in the workflow.
    
    Args:
        standard_topology_path (str): Path to the standard topology file.
        coordinates_file (str): Path to the coordinates file (e.g., PDB, XTC)."""
    mda_topology = to_MDAnalysis_topology(standard_topology_path)
    # Create a MDAnalysis topology from the standard topology file
    return Universe(mda_topology, coordinates_file)

def get_mda_universe (structure_file : 'File',              # To load in MDAnalysis
                      trajectory_file : 'File',             # To load in MDAnalysis
                      reference_bonds : List[List[int]],    # To set the bonds
                      charges : List[float]) -> 'Universe': # To set the charges
    """Create a MDAnalysis universe using data in the workflow."""

    # Make MDAnalysis warnings and logs grey
    print(GREY_HEADER, end='\r')
    universe = Universe(structure_file.path, trajectory_file.path)
    # Set the atom bonds
    bonds = []
    for bond_from, bond_tos in enumerate(reference_bonds):
        for bond_to in bond_tos:
            bond = tuple(sorted([bond_from, bond_to]))
            bonds.append(bond)
    universe.add_TopologyAttr('bonds', set(bonds))

    # Set the charges
    if charges != MISSING_CHARGES and len(charges) > 0:
        universe.add_TopologyAttr('charges', np.array(charges, dtype=np.float32))
    print(COLOR_END, end='\r')
    return universe

# Get a cksum from a MDA universe for equality comparission
def get_mda_universe_cksum (universe) -> str:
    pickled = pickle.dumps(universe.atoms)
    return sum(pickled)

# Get TPR bonds using MDAnalysis
# WARNING: Sometimes this function takes additional constrains as actual bonds
# DANI: si miras los topology.bonds.values estos enlaces falsos van al final
# DANI: Lo veo porque los índices están en orden ascendente y veulven a empezar
# DANI: He pedido ayuda aquí https://github.com/MDAnalysis/mdanalysis/pull/463
def get_tpr_bonds_mdanalysis (tpr_filepath : str) -> List[ Tuple[int, int] ]:
    parser = TPRParser(tpr_filepath)
    topology = parser.parse()
    bonds = list(topology.bonds.values)
    return bonds