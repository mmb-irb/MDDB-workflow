import json
import numpy as np
import pickle

from model_workflow.utils.type_hints import *

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

def to_MDAnalysis_topology(standard_topology_file : 'File') -> 'Topology':
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
    topology = json.load(open(standard_topology_file.path))

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

def get_mda_universe (standard_topology_file : 'File', structure_file : 'File') -> 'Universe':
    """Create a MDAnalysis universe using data in the workflow."""
    mda_topology = to_MDAnalysis_topology(standard_topology_file)
    # Create a MDAnalysis topology from the standard topology file
    return Universe(mda_topology, structure_file.path)

# Get a cksum from a MDA universe for equality comparission
def get_mda_universe_cksum (universe) -> str:
    pickled = pickle.dumps(universe.atoms)
    return sum(pickled)