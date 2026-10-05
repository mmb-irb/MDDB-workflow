import re
import numpy as np
import mdtraj as mdt
import pytraj as pt
from scipy.spatial import cKDTree

from mddb_workflow.tools.get_reduced_trajectory import get_reduced_trajectory
from mddb_workflow.utils.type_hints import *

# Set the distance cutoff in Ångstroms (Å) to consider two heavy atoms in contact
CONTACT_DISTANCE_CUTOFF = 5
# Set the minimum number of residues in contact (in any of both chains) to consider two chains in contact
# This prevents marginal contacts (e.g. a single side chain) to make a multimer
MIN_CONTACT_RESIDUES = 3
# Set the minimum number of different residue pairs with inter-chain hydrogen bonds between bases
# A single base pair may be casual, but several base pairs mean the strands are paired
MIN_BASE_PAIRS = 2

# Set the chain classifications we care about
PROTEIN_CLASSIFICATIONS = { 'protein' }
NUCLEIC_CLASSIFICATIONS = { 'dna', 'rna', 'nucleic' }

# Set the name of the multimer according to the number of chains
MULTIMER_NAMES = {
    1: 'monomer',
    2: 'dimer',
    3: 'trimer',
    4: 'tetramer',
    5: 'pentamer',
    6: 'hexamer',
    7: 'heptamer',
    8: 'octamer',
    9: 'nonamer',
    10: 'decamer',
    11: 'undecamer',
    12: 'dodecamer',
}
# Set the name for multimers with more chains than the ones above
MANY_CHAINS_MULTIMER_NAME = 'multimer'
# Set the prefixes for multimers with identical or different proteins
HOMO_PREFIX = 'homo'
HETERO_PREFIX = 'hetero'

# Set the nucleic strand keywords according to the number of strands
# The nucleic type is appended to them (e.g. 'double strand dna')
STRAND_NAMES = {
    1: 'single strand',
    2: 'double strand',
    3: 'triple strand',
}
# Set the name for groups with more strands than the ones above
MANY_STRANDS_NAME = 'multiple strand'
# Set the nucleic type for groups of strands with both dna and rna
HYBRID_NUCLEIC_TYPE = 'nucleic'

# Set the protein-nucleic complex keywords
PROTEIN_DNA_COMPLEX_KEYWORD = 'protein-dna complex'
PROTEIN_RNA_COMPLEX_KEYWORD = 'protein-rna complex'

# Set the names of nucleic atoms which are not part of the base
# Hydrogen bonds involving these atoms are not considered as base pairing
NUCLEIC_BACKBONE_ATOM_NAMES = { 'P', 'OP1', 'OP2', 'OP3', 'O1P', 'O2P', 'O3P' }

# Get the protein and nucleic chains in the structure
def get_polymer_chains (structure : 'Structure') -> tuple[list['Chain'], list['Chain']]:
    """Get the protein and nucleic chains in the structure.
    Chains with any coarse grain regions are ignored.
    """
    protein_chains = []
    nucleic_chains = []
    for chain in structure.chains:
        if chain.has_cg(): continue
        classification = chain.classification
        if classification in PROTEIN_CLASSIFICATIONS: protein_chains.append(chain)
        elif classification in NUCLEIC_CLASSIFICATIONS: nucleic_chains.append(chain)
    return protein_chains, nucleic_chains

# Get stable contacts along the trajectory between different chains
def get_chain_contacts (
    structure : 'Structure',
    structure_file : 'File',
    trajectory_file : 'File',
    snapshots : int,
    frames_limit : int = 100,
    # Minimum percent of frames (from 0 to 1) where a contact must happen to be considered
    min_frames_percent : float = 0.05,
) -> dict:
    """Find which polymer chains are in contact along the trajectory.
    For nucleic chains, find also which chains are base paired (i.e. hydrogen bonds between bases).
    Frames are taken along the whole trajectory (reduced trajectory).
    Only contacts and base pairings happening in a minimum percent of frames are kept.
    Note that hydrogen atoms are required to find hydrogen bonds.
    Chains are refered by their indices.
    """
    contacts = []
    base_pairings = []
    output = { 'frames': 0, 'contacts': contacts, 'base_pairings': base_pairings }
    # Get the chains we are interested in
    protein_chains, nucleic_chains = get_polymer_chains(structure)
    chains = protein_chains + nucleic_chains
    if len(chains) < 2: return output
    # Use a reduced trajectory in case the original trajectory has many frames
    reduced_trajectory_filepath, _, frames = get_reduced_trajectory(trajectory_file, snapshots, frames_limit)
    output['frames'] = frames

    # Count in how many frames every pair of chains is in contact according to their heavy atoms
    contact_counts = {}
    # Get heavy atom indices for every chain
    chain_atom_indices = {}
    for chain in chains:
        atom_indices = [ atom.index for atom in chain.atoms if atom.element != 'H' ]
        if len(atom_indices) == 0: continue
        chain_atom_indices[chain.index] = atom_indices
    chain_indices = list(chain_atom_indices.keys())
    # Load heavy atoms only
    # Note that atom indices must be sorted for mdtraj to keep their order
    loaded_atom_indices = np.array(sorted(sum(chain_atom_indices.values(), [])))
    trajectory = mdt.load(reduced_trajectory_filepath, top=structure_file.path, atom_indices=loaded_atom_indices)
    # Get the positions of every chain atom in the loaded atoms and their residue indices
    chain_positions = { index: np.searchsorted(loaded_atom_indices, atom_indices)
        for index, atom_indices in chain_atom_indices.items() }
    chain_residue_indices = { index: np.array([ structure.atoms[atom_index].residue_index for atom_index in atom_indices ])
        for index, atom_indices in chain_atom_indices.items() }
    # Iterate frames
    for frame_coordinates in trajectory.xyz:
        # Convert nanometers to Ångstroms
        frame_coordinates = frame_coordinates * 10
        chain_trees = { index: cKDTree(frame_coordinates[positions]) for index, positions in chain_positions.items() }
        for i, chain_index_1 in enumerate(chain_indices):
            tree_1 = chain_trees[chain_index_1]
            for chain_index_2 in chain_indices[i+1:]:
                tree_2 = chain_trees[chain_index_2]
                # Find atoms of chain 1 which are close to any atom of chain 2
                neighbours = tree_1.query_ball_tree(tree_2, CONTACT_DISTANCE_CUTOFF)
                atoms_1 = [ atom_1 for atom_1, atoms_2 in enumerate(neighbours) if len(atoms_2) > 0 ]
                if len(atoms_1) == 0: continue
                atoms_2 = list(set([ atom_2 for atoms_2 in neighbours for atom_2 in atoms_2 ]))
                residues_1 = set(chain_residue_indices[chain_index_1][atoms_1])
                residues_2 = set(chain_residue_indices[chain_index_2][atoms_2])
                # Check the contact to be relevant enough
                if max(len(residues_1), len(residues_2)) < MIN_CONTACT_RESIDUES: continue
                chain_pair = (chain_index_1, chain_index_2)
                contact_counts[chain_pair] = contact_counts.get(chain_pair, 0) + 1
    # Keep only contacts happening in enough frames
    for chain_pair, count in contact_counts.items():
        percent = count / frames
        if percent < min_frames_percent: continue
        contacts.append({ 'chains': list(chain_pair), 'percent': percent })

    # Count in how many frames every pair of nucleic chains in contact is base paired using pytraj hydrogen bonds
    # Two chains are base paired in a frame when they have hydrogen bonds between enough different base pairs
    contact_pairs = set([ tuple(contact['chains']) for contact in contacts ])
    nucleic_chain_indices = set([ chain.index for chain in nucleic_chains ])
    if not any(a in nucleic_chain_indices and b in nucleic_chain_indices for a, b in contact_pairs):
        return output
    pytraj_trajectory = pt.iterload(reduced_trajectory_filepath, structure_file.path)
    # Select nucleic atoms only
    nucleic_atom_indices = sum([ chain.atom_indices for chain in nucleic_chains ], [])
    nucleic_selection = structure.select_atom_indices(nucleic_atom_indices)
    hbonds = pt.hbond(pytraj_trajectory, mask=nucleic_selection.to_pytraj())
    # Set a function to get a structure atom from a pytraj residue number (1-based index) and an atom name
    # Double check the residue name and return None if it does not match
    def get_atom (pytraj_residue_number : int, residue_name : str, atom_name : str) -> Optional['Atom']:
        residue_index = pytraj_residue_number - 1
        if residue_index >= len(structure.residues): return None
        residue = structure.residues[residue_index]
        if residue.name[0:len(residue_name)] != residue_name: return None
        return residue.get_atom_by_name(atom_name)
    # For every frame and pair of chains, save the residue pairs with hydrogen bonds between bases
    frame_base_pairs = [ {} for frame in range(pytraj_trajectory.n_frames) ]
    # Use the hbond 'old keys' as it is done in the hydrogen bonds analysis
    # e.g. 'DG_23@O6-DC_4@N4-H41' (acceptor residue, acceptor atom, donor residue, donor atom, hydrogen)
    # Note that the first key is the total number of hydrogen bonds and it is skipped by the regex
    hbond_keys = hbonds._old_keys
    for i, hbond in enumerate(hbonds):
        match = re.match(r'(\w*)_(\d*)@(.*)-(\w*)_(\d*)@(.*)-(.*)', hbond_keys[i])
        if match is None: continue
        acceptor = get_atom(int(match.group(2)), match.group(1), match.group(3))
        donor = get_atom(int(match.group(5)), match.group(4), match.group(6))
        if acceptor is None or donor is None: continue
        # Skip hydrogen bonds within the same chain
        if acceptor.chain_index == donor.chain_index: continue
        # Skip hydrogen bonds involving atoms out of bases
        if not is_base_atom(acceptor) or not is_base_atom(donor): continue
        chain_pair = tuple(sorted([ acceptor.chain_index, donor.chain_index ]))
        residue_pair = tuple(sorted([ acceptor.residue_index, donor.residue_index ]))
        for frame, value in enumerate(hbond.values):
            if not value: continue
            frame_base_pairs[frame].setdefault(chain_pair, set()).add(residue_pair)
    # Count the frames where every pair of chains has enough base pairs
    pairing_counts = {}
    for base_pairs in frame_base_pairs:
        for chain_pair, residue_pairs in base_pairs.items():
            if len(residue_pairs) < MIN_BASE_PAIRS: continue
            pairing_counts[chain_pair] = pairing_counts.get(chain_pair, 0) + 1
    # Keep only base pairings between chains in contact happening in enough frames
    for chain_pair, count in pairing_counts.items():
        if chain_pair not in contact_pairs: continue
        percent = count / frames
        if percent < min_frames_percent: continue
        base_pairings.append({ 'chains': list(chain_pair), 'percent': percent })
    return output


def is_base_atom (atom : 'Atom') -> bool:
    """Check if a nucleic atom belongs to the base (i.e. it is not sugar or phosphate)."""
    if "'" in atom.name or '*' in atom.name: return False
    return atom.name not in NUCLEIC_BACKBONE_ATOM_NAMES


def get_multimerism_keywords (structure : 'Structure', protein_map : list[dict], chain_contacts : dict) -> list[str]:
    """Set multimerism keywords according to the chain contacts along the trajectory.
    Protein chains in contact are grouped in complexes which are labeled according to their number of chains
    (e.g. monomer, dimer, trimer...) and whether chains share the same UniProt reference (e.g. homodimer, heterodimer).
    Base paired nucleic chains are grouped to tell single/double strands.
    Protein chains in contact with nucleic chains are labeled as protein-dna/rna complexes.
    Coarse grain chains are ignored.
    """
    keywords = []
    # Get the chains we are interested in
    protein_chains, nucleic_chains = get_polymer_chains(structure)
    # Get pairs of chains in contact and base paired
    contact_pairs = set([ tuple(contact['chains']) for contact in chain_contacts['contacts'] ])
    paired_pairs = set([ tuple(pairing['chains']) for pairing in chain_contacts['base_pairings'] ])
    # Protein multimerism
    if len(protein_chains) > 0:
        keywords += get_protein_multimerism_keywords(protein_chains, contact_pairs, protein_map)
    # Nucleic strands
    if len(nucleic_chains) > 0:
        keywords += get_nucleic_strand_keywords(nucleic_chains, paired_pairs)
    # Protein-nucleic complexes
    if len(protein_chains) > 0 and len(nucleic_chains) > 0:
        keywords += get_protein_nucleic_complex_keywords(protein_chains, nucleic_chains, contact_pairs)
    return keywords


def get_connected_groups (chain_indices : list[int], pairs : set[tuple[int, int]]) -> list[set[int]]:
    """Given a list of chain indices and pairs of connected chains, get the groups of connected chains.
    Note that a chain may be connected to another chain through a third chain.
    """
    groups = []
    remaining = set(chain_indices)
    while remaining:
        group = set()
        pending = [ remaining.pop() ]
        while pending:
            current = pending.pop()
            group.add(current)
            for a, b in pairs:
                if a == current: neighbour = b
                elif b == current: neighbour = a
                else: continue
                if neighbour not in remaining: continue
                remaining.remove(neighbour)
                pending.append(neighbour)
        groups.append(group)
    return groups


def get_protein_multimerism_keywords (
    protein_chains : list['Chain'],
    contact_pairs : set[tuple[int, int]],
    protein_map : list[dict],
) -> list[str]:
    """Set the multimer keywords for every group of protein chains in contact.
    Use the UniProt reference of each chain to tell homo and hetero multimers.
    If a chain has no UniProt reference then use its sequence instead.
    """
    keywords = []
    # Get the UniProt accession of every chain from the protein map
    # Note that chain names are expected to be unique at this point
    chain_uniprots = {}
    for chain_data in protein_map:
        match = chain_data.get('match', None)
        if not match: continue
        reference = match.get('ref', None)
        # Reference may also be a flag (e.g. no referable, not found) or None
        if type(reference) != dict: continue
        chain_uniprots[chain_data['name']] = reference['uniprot']
    # Set a function to get a chain identity, which is used to compare chains
    def get_chain_identity (chain : 'Chain') -> str:
        uniprot = chain_uniprots.get(chain.name, None)
        if uniprot: return uniprot
        return chain.get_sequence()
    # Get groups of protein chains in contact
    chains_by_index = { chain.index: chain for chain in protein_chains }
    groups = get_connected_groups(list(chains_by_index.keys()), contact_pairs)
    for group in groups:
        multimer_name = MULTIMER_NAMES.get(len(group), MANY_CHAINS_MULTIMER_NAME)
        keywords.append(multimer_name)
        # Monomers have no homo/hetero labels
        if len(group) == 1: continue
        identities = set([ get_chain_identity(chains_by_index[chain_index]) for chain_index in group ])
        prefix = HOMO_PREFIX if len(identities) == 1 else HETERO_PREFIX
        keywords.append(prefix + multimer_name)
    # Remove duplicates while keeping the order
    return list(dict.fromkeys(keywords))


def get_nucleic_strand_keywords (
    nucleic_chains : list['Chain'],
    paired_pairs : set[tuple[int, int]],
) -> list[str]:
    """Set strand keywords for every group of base paired nucleic chains.
    Keywords include the nucleic type of the group (e.g. 'single strand rna', 'double strand dna').
    Groups with both dna and rna (or chains with both) are labeled as nucleic (e.g. 'double strand nucleic').
    """
    keywords = []
    chains_by_index = { chain.index: chain for chain in nucleic_chains }
    groups = get_connected_groups(list(chains_by_index.keys()), paired_pairs)
    for group in groups:
        strand_name = STRAND_NAMES.get(len(group), MANY_STRANDS_NAME)
        classifications = set([ chains_by_index[chain_index].classification for chain_index in group ])
        nucleic_type = classifications.pop() if classifications in [{ 'dna' }, { 'rna' }] else HYBRID_NUCLEIC_TYPE
        keywords.append(f'{strand_name} {nucleic_type}')
    # Remove duplicates while keeping the order
    return list(dict.fromkeys(keywords))


def get_protein_nucleic_complex_keywords (
    protein_chains : list['Chain'],
    nucleic_chains : list['Chain'],
    contact_pairs : set[tuple[int, int]],
) -> list[str]:
    """Set protein-dna/rna complex keywords if any protein chain is in contact with a nucleic chain."""
    keywords = []
    protein_chain_indices = set([ chain.index for chain in protein_chains ])
    for nucleic_chain in nucleic_chains:
        if not any((a == nucleic_chain.index and b in protein_chain_indices)
            or (b == nucleic_chain.index and a in protein_chain_indices) for a, b in contact_pairs): continue
        classification = nucleic_chain.classification
        if classification in { 'dna', 'nucleic' }: keywords.append(PROTEIN_DNA_COMPLEX_KEYWORD)
        if classification in { 'rna', 'nucleic' }: keywords.append(PROTEIN_RNA_COMPLEX_KEYWORD)
    # Remove duplicates while keeping the order
    return list(dict.fromkeys(keywords))
