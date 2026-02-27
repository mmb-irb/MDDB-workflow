import re
import json
import numpy as np
import pickle

from mddb_workflow.utils.type_hints import *
from mddb_workflow.utils.constants import GREY_HEADER, COLOR_END
from mddb_workflow.utils.auxiliar import MISSING_CHARGES, MISSING_BONDS

import MDAnalysis
from MDAnalysis.topology.TPRParser import TPRParser
# from MDAnalysis.topology.PDBParser import PDBParser # for class reference
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


def to_MDAnalysis_topology(standard_topology_path: str) -> 'Topology':
    """Create a MDAnalysis topology from a json topology file.

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

    attrs = [Atomnames(topology['atom_names']),
             Elements(topology['atom_elements']),
             Bonds(topology['bonds']),
             Resnames(topology['residue_names']),
             Resnums(topology['residue_numbers']),
             ChainIDs(topology['chainIDs']),
             Segids(topology['chain_names']),
             Resids(np.arange(len(topology['residue_names']))),
    ]
    if topology['atom_charges']:
        attrs.append(Charges(topology['atom_charges']))

    mda_top = Topology(n_atoms=len(topology['atom_names']),
                       n_res=len(topology['residue_names']),
                       n_seg=len(topology['chain_names']),
                       attrs=attrs,
                       atom_resindex=topology['atom_residue_indices'],
                       residue_segindex=topology['residue_chain_indices']
    )
    return mda_top


def get_mda_universe_from_stopology(
        standard_topology_path: str, coordinates_file: str) -> 'Universe':
    """Create a MDAnalysis universe using data in the workflow.

    Args:
        standard_topology_path (str): Path to the standard topology file.
        coordinates_file (str): Path to the coordinates file (e.g., PDB, XTC).

    """
    mda_topology = to_MDAnalysis_topology(standard_topology_path)
    # Create a MDAnalysis topology from the standard topology file
    return Universe(mda_topology, coordinates_file)


def get_mda_universe(structure_file: 'File',               # To load in MDAnalysis
                     trajectory_file: 'File',              # To load in MDAnalysis
                     reference_bonds: list[list[int]],     # To set the bonds
                     charges: list[float]) -> 'Universe':  # To set the charges
    """Create a MDAnalysis universe using data in the workflow."""
    # Make MDAnalysis warnings and logs grey
    print(GREY_HEADER, end='\r')
    universe = Universe(structure_file.path, trajectory_file.path)
    # Set the atom bonds
    bonds = []
    for bond_from, bond_tos in enumerate(reference_bonds):
        if bond_tos == MISSING_BONDS: continue
        for bond_to in bond_tos:
            bond = tuple(sorted([bond_from, bond_to]))
            bonds.append(bond)
    universe.add_TopologyAttr('bonds', set(bonds))

    # Set the charges
    if charges and charges != MISSING_CHARGES and len(charges) > 0:
        universe.add_TopologyAttr('charges', np.array(charges, dtype=np.float32))

    # No elements can happen if we use the faith flag
    if not hasattr(universe.atoms, 'elements'):
        universe.guess_TopologyAttrs(to_guess=['elements'])
    print(COLOR_END, end='\r')
    return universe


def get_mda_universe_cksum(universe) -> str:
    """Get a cksum from a MDA universe for equality comparison."""
    pickled = pickle.dumps(universe.atoms)
    return sum(pickled)


# DANI: si miras los topology.bonds.values estos enlaces falsos van al final
# DANI: Lo veo porque los índices están en orden ascendente y veulven a empezar
# DANI: He pedido ayuda aquí https://github.com/MDAnalysis/mdanalysis/pull/463
def get_tpr_bonds_mdanalysis(tpr_filepath: str) -> list[tuple[int, int]]:
    """Get TPR bonds using MDAnalysis.
    WARNING: Sometimes this function takes additional constrains as actual bonds.
    """
    parser = TPRParser(tpr_filepath)
    topology = parser.parse()
    bonds = list(topology.bonds.values)
    return bonds


def get_all_acyl_chains(residue: 'MDAnalysis.Residue', min_length=6, removeHs=True) -> list[list[int]]:
    """Find all groups of connected Carbon atoms within a residue, including cyclic structures.

    Returns:
        list
            A list of lists, where each inner list contains atom indices of
            connected Carbon atoms forming a distinct group

    """
    def explore_carbon_group(start_atom, visited):
        """Explore a connected group of carbons."""
        to_visit = [start_atom]
        group = set()
        while to_visit:
            current_atom = to_visit.pop(0)
            if current_atom.index in visited:
                continue
            visited.add(current_atom.index)
            if current_atom.element == 'C':
                group.add(current_atom.index)
                # Add all bonded carbon atoms to the visit queue
                for bond in current_atom.bonds:
                    for bonded_atom in bond.atoms:
                        if (bonded_atom.element == 'C' and
                                bonded_atom.index not in visited):
                            to_visit.append(bonded_atom)
        return list(group)

    # Get all Carbon atoms in the residue
    carbon_atoms = residue.atoms.select_atoms('element C')
    visited = set()
    if removeHs:
        # Remove ester carbons (without hydrogens)
        for at in carbon_atoms:
            if 'H' not in at.bonded_atoms.elements:
                visited.add(at.index)
    carbon_groups = []
    # Find all distinct groups of connected carbons
    for carbon in carbon_atoms:
        if carbon.index not in visited:
            # Get all carbons connected to this one
            connected_carbons = explore_carbon_group(carbon, visited)
            if len(connected_carbons) >= min_length:  # Only add non-empty groups
                carbon_groups.append(sorted(connected_carbons))
    return carbon_groups


def get_acyl_chain_atom_names(
    residue: 'MDAnalysis.Residue',
    removeHs: bool = True,
    sort_key=None,
) -> list[list[str]]:
    """Get atom names for each acyl chain in a lipid residue.

    Returns:
        List of lists of atom name strings, one per chain.

    """
    chains = get_all_acyl_chains(residue, removeHs=removeHs)
    result = []
    for chain in chains:
        names = list(residue.universe.atoms[sorted(chain)].names)
        if sort_key is not None:
            names = sorted(names, key=sort_key)
        result.append(names)
    return result


def get_tail_atom_group(residue: 'MDAnalysis.Residue', include_hydrogens: bool = True) -> 'MDAnalysis.AtomGroup':
    """Get the tail (acyl chain) AtomGroup for a lipid residue.

    Finds connected carbon chains via `get_all_acyl_chains` and optionally
    includes bonded hydrogen atoms.

    Args:
        residue: An MDAnalysis Residue representing a single lipid molecule.
        include_hydrogens: If True, hydrogen atoms bonded to tail carbons are
            included in the returned AtomGroup.

    Returns:
        AtomGroup containing the tail atoms of the lipid.

    """
    chains = get_all_acyl_chains(residue)
    tail_idx = [idx for chain in chains for idx in chain]
    if not tail_idx:
        return residue.atoms[[]]  # empty AtomGroup
    tail_ag = residue.universe.atoms[tail_idx]
    if include_hydrogens:
        h_indices = [
            bonded.index
            for atom in tail_ag.atoms
            for bonded in atom.bonded_atoms
            if bonded.name.startswith('H')
        ]
        if h_indices:
            tail_ag = tail_ag | residue.universe.atoms[h_indices]
    return tail_ag


def get_head_tail_split(
    residue: 'MDAnalysis.Residue',
    include_hydrogens: bool = True,
) -> tuple['MDAnalysis.AtomGroup', 'MDAnalysis.AtomGroup']:
    """Split a lipid residue into headgroup and tail AtomGroups.

    Args:
        residue: An MDAnalysis Residue representing a single lipid molecule.
        include_hydrogens: Forwarded to `get_tail_atom_group`.

    Returns:
        (head_ag, tail_ag) tuple of AtomGroups.

    """
    tail_ag = get_tail_atom_group(residue, include_hydrogens=include_hydrogens)
    head_ag = residue.atoms - tail_ag
    return head_ag, tail_ag


# Regex matching CG tail bead names: C<digit(s)><chain letter> (e.g. C1A, C3B)
CG_TAIL_NAME_RE = re.compile(r'^C\d+[A-Z]$')


def get_cg_acyl_chains(residue: 'MDAnalysis.Residue') -> list[list[str]]:
    """Get CG acyl chain bead names grouped by chain letter.

    Identifies tail beads using `CG_TAIL_NAME_RE` and groups them by their
    trailing chain letter (e.g. C1A, C2A -> chain 'A').

    Args:
        residue: An MDAnalysis Residue representing a single CG lipid molecule.

    Returns:
        List of lists of bead name strings, one per chain (sorted by chain letter).

    """
    chains: dict[str, list[str]] = {}
    for atom in residue.atoms:
        m = CG_TAIL_NAME_RE.match(atom.name)
        if m:
            chain_letter = atom.name[-1]
            chains.setdefault(chain_letter, []).append(atom.name)
    return [names for _, names in sorted(chains.items())]


def get_cg_head_tail_split(
    residue: 'MDAnalysis.Residue',
) -> tuple['MDAnalysis.AtomGroup', 'MDAnalysis.AtomGroup']:
    """Split a coarse-grain lipid residue into headgroup and tail AtomGroups.

    Tail beads are identified via `get_cg_acyl_chains` using the CG naming
    convention: C<digit(s)><chain letter> (e.g. C1A, C3B).

    Args:
        residue: An MDAnalysis Residue representing a single CG lipid molecule.

    Returns:
        (head_ag, tail_ag) tuple of AtomGroups.

    """
    chains = get_cg_acyl_chains(residue)
    tail_names = {name for chain in chains for name in chain}
    mask = np.array([atom.name in tail_names for atom in residue.atoms])
    tail_ag = residue.atoms[mask]
    head_ag = residue.atoms[~mask]
    return head_ag, tail_ag
