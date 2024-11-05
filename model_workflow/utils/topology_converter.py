import json
import numpy as np
#from MDAnalysis.topology.PDBParser import PDBParser # for class reference
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
def to_MDAnalysis_topology(top_js : str) -> 'Topology':
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

    :param top_js: path to the json file
    :returns: a MDAnalysis topology object
    """
    top_js = json.load(open(top_js))

    # transform bond to non redundant tuples
    bonds = []
    for bond_from, bond_tos in enumerate(top_js['atom_bonds']):
        for bond_to in bond_tos:
            bond = tuple(sorted([bond_from, bond_to]))
            bonds.append(bond)
    # sorted(set(bonds)) if you want to check for duplicates
    top_js['bonds'] = set(bonds)
    atom_2_chain = np.array(top_js['residue_chain_indices'])[top_js['atom_residue_indices']]
    top_js['chainIDs'] = np.array(top_js['chain_names'])[atom_2_chain]


    attrs = [Atomnames (top_js['atom_names']),
             Elements  (top_js['atom_elements']),
             Charges   (top_js['atom_charges']),
             Bonds     (top_js['bonds']),
             Resnames  (top_js['residue_names']),
             Resnums   (top_js['residue_numbers']),
             ChainIDs  (top_js['chainIDs']),
             Segids    (top_js['chain_names']),
             Resids    (np.arange(len(top_js['residue_names']))),
    ]

    mda_top = Topology(n_atoms=len(top_js['atom_names']), 
                       n_res=len(top_js['residue_names']), 
                       n_seg=len(top_js['chain_names']),
                       attrs=attrs,
                       atom_resindex=top_js['atom_residue_indices'],
                       residue_segindex=top_js['residue_chain_indices']
    )
    return mda_top