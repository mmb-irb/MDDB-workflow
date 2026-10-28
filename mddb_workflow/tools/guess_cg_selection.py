import os
import numpy as np
from collections import Counter
from scipy.spatial import cKDTree

from mddb_workflow.tools.get_bonds import mine_topology_bonds
from mddb_workflow.utils.auxiliar import InputError, MISSING_BONDS, get_auxiliar_filepath
from mddb_workflow.utils.constants import CG_ATOM_ELEMENT, STANDARD_SOLVENT_RESIDUE_NAMES
from mddb_workflow.utils.constants import PROTEIN_RESIDUE_NAME_LETTERS, NUCLEIC_RESIDUE_NAME_LETTERS
from mddb_workflow.utils.selections import Selection
from mddb_workflow.utils.gmx_spells import count_xtc_frames
from mddb_workflow.utils.pyt_spells import get_reduced_pytraj_trajectory
from mddb_workflow.utils.type_hints import *

# Words in the force field names which tell us it is a coarse grain force field
CG_FORCE_FIELD_KEYWORDS = ['martini', 'sirah', 'elnedyn', 'oxdna', '3spn']

# Number of heavy atoms in every residue inside a chain
# Atomistic force fields have at least these atoms (hydrogens apart)
PROTEIN_HEAVY_ATOMS = {
    'A': 5, 'R': 11, 'N': 8, 'D': 8, 'C': 6, 'Q': 9, 'E': 9, 'G': 4, 'H': 10, 'I': 8,
    'L': 8, 'K': 9, 'M': 8, 'F': 11, 'P': 7, 'S': 6, 'T': 7, 'W': 14, 'Y': 12, 'V': 7,
}
NUCLEIC_HEAVY_ATOMS = { 'A': 21, 'C': 19, 'G': 22, 'T': 20, 'U': 20 }

# Atom names which every atomistic residue has
PROTEIN_MAIN_ATOM_NAMES = [{'N', 'CA', 'C'}]
NUCLEIC_MAIN_ATOM_NAMES = [{"C3'", "C4'"}, {'C3*', 'C4*'}]

# If a group of residues has less atoms than this ratio of the expected heavy atoms then it is coarse grain
CG_ATOM_COUNT_RATIO = 0.6
# If atoms in a group of residues have less than this ratio of residues with the main atom names then it is coarse grain
CG_ATOM_NAMES_RATIO = 0.5
# Atomistic bonds are always shorter than this distance (Å) while coarse grain bonds are longer
CG_BOND_LENGTH = 2.2
# Number of frames along the trajectory where distances are measured
# Taking the median along frames we skip punctual clashes and bonds broken by periodic boundary conditions
FRAMES_LIMIT = 10

CG, AA = 'CG', 'AA'


def has_only_cg_force_fields(input_force_fields: Optional[list[str] | str]) -> bool:
    """Check if all input force fields are coarse grain force fields."""
    if not input_force_fields: return False
    force_fields = [input_force_fields] if type(input_force_fields) == str else input_force_fields
    return all(any(keyword in force_field.lower() for keyword in CG_FORCE_FIELD_KEYWORDS)
        for force_field in force_fields)


def get_residue_groups(structure: 'Structure') -> dict[str, list['Residue']]:
    """Split residues in groups which are expected to be either all coarse grain or all atomistic.
    Protein and nucleic residues are grouped by chain.
    The rest of residues are grouped by residue name and number of atoms, so for instance all lipids of a kind go together.
    """
    groups = {}
    for residue in structure.residues:
        if residue.name in PROTEIN_RESIDUE_NAME_LETTERS:
            group_name = f'protein in chain "{residue.chain.name}"'
        elif residue.name in NUCLEIC_RESIDUE_NAME_LETTERS:
            group_name = f'nucleic in chain "{residue.chain.name}"'
        else:
            group_name = f'residue {residue.name} ({residue.atom_count} atoms)'
        groups.setdefault(group_name, []).append(residue)
    return groups


def is_polymer(residues: list['Residue']) -> bool:
    return residues[0].name in PROTEIN_RESIDUE_NAME_LETTERS or residues[0].name in NUCLEIC_RESIDUE_NAME_LETTERS


def vote_by_atom_names(residues: list['Residue']) -> Optional[str]:
    """Check if residues have the atom names which every atomistic residue has.
    This only works for protein and nucleic residues.
    """
    if not is_polymer(residues): return None
    main_atom_names = PROTEIN_MAIN_ATOM_NAMES if residues[0].name in PROTEIN_RESIDUE_NAME_LETTERS else NUCLEIC_MAIN_ATOM_NAMES
    matching_residues = 0
    for residue in residues:
        atom_names = set(atom.name for atom in residue.atoms)
        if any(names <= atom_names for names in main_atom_names): matching_residues += 1
    return CG if matching_residues / len(residues) < CG_ATOM_NAMES_RATIO else AA


def vote_by_atom_count(residues: list['Residue']) -> Optional[str]:
    """Compare the number of atoms with the number of heavy atoms we would expect.
    This only works for protein and nucleic residues.
    """
    if not is_polymer(residues): return None
    is_protein = residues[0].name in PROTEIN_RESIDUE_NAME_LETTERS
    expected_atom_count = 0
    for residue in residues:
        if is_protein: expected_atom_count += PROTEIN_HEAVY_ATOMS.get(PROTEIN_RESIDUE_NAME_LETTERS[residue.name], 0)
        else: expected_atom_count += NUCLEIC_HEAVY_ATOMS.get(NUCLEIC_RESIDUE_NAME_LETTERS[residue.name], 0)
    # Residues with an unusual letter (e.g. X) have no expected number of atoms
    if expected_atom_count == 0: return None
    atom_count = sum(residue.atom_count for residue in residues)
    return CG if atom_count / expected_atom_count < CG_ATOM_COUNT_RATIO else AA


def get_frames_coordinates(structure: 'Structure', trajectory_file: 'File') -> np.ndarray:
    """Get atom coordinates in a few frames along the trajectory."""
    # Write the structure in a mmCIF file to be used as topology for pytraj
    auxiliar_mmcif_filepath = get_auxiliar_filepath('.cg_guess.cif')
    structure.generate_mmcif_file(auxiliar_mmcif_filepath)
    snapshots = count_xtc_frames(trajectory_file)
    reduced_trajectory, _, _ = get_reduced_pytraj_trajectory(auxiliar_mmcif_filepath, trajectory_file.path, snapshots, FRAMES_LIMIT)
    coordinates = np.array([ frame.xyz for frame in reduced_trajectory.iterframe() ])
    os.remove(auxiliar_mmcif_filepath)
    return coordinates


def vote_by_topology_bonds(residues: list['Residue'], coordinates: np.ndarray, topology_bonds: Optional[list]) -> Optional[str]:
    """Measure the length of the topology bonds along frames."""
    if not topology_bonds: return None
    atom_indices = set(atom_index for residue in residues for atom_index in residue.atom_indices)
    bonds = [ (atom_index, bonded_index) for atom_index in atom_indices
        if topology_bonds[atom_index] != MISSING_BONDS
        for bonded_index in topology_bonds[atom_index]
        if atom_index < bonded_index and bonded_index in atom_indices ]
    if not bonds: return None
    bonds = np.array(bonds)
    lengths = np.linalg.norm(coordinates[:, bonds[:, 0]] - coordinates[:, bonds[:, 1]], axis=2)
    return CG if np.median(np.median(lengths, axis=0)) > CG_BOND_LENGTH else AA


def vote_by_closest_atoms(residues: list['Residue'], coordinates: np.ndarray) -> Optional[str]:
    """Measure the distance from every atom to its closest atom, as a guess of its bond length.
    Single atom residues are skipped unless they are part of a chain (e.g. one bead per nucleotide).
    Otherwise the closest atom of an ion would be another ion far away.
    """
    if not is_polymer(residues): residues = [ residue for residue in residues if residue.atom_count > 1 ]
    atom_indices = [ atom_index for residue in residues for atom_index in residue.atom_indices ]
    if len(atom_indices) < 2: return None
    distances = [ cKDTree(frame[atom_indices]).query(frame[atom_indices], k=2)[0][:, 1] for frame in coordinates ]
    return CG if np.median(np.median(distances, axis=0)) > CG_BOND_LENGTH else AA


def guess_cg_selection(
    structure: 'Structure',
    trajectory_file: 'File',
    topology_file: Union['File', Exception],
    input_force_fields: Optional[list[str] | str],
    verbose: bool = False,
) -> 'Selection':
    """Guess which atoms are coarse grain beads.

    If all input force fields are coarse grain then all atoms are coarse grain.
    Otherwise residues are split in groups (see get_residue_groups) and every group is checked in 4 ways:
    1. Topology bonds: coarse grain bonds are much longer than atomistic bonds
    2. Closest atoms: in atomistic residues every atom has a bonded atom very close
    3. Atom count: coarse grain residues have much less atoms than heavy atoms in the real residue
    4. Atom names: atomistic residues have standard names (e.g. CA in proteins)
    Distances are measured in a few frames along the trajectory and we keep the median.
    The most voted option wins. In case of tie, votes are prioritized in the same order.
    Atoms which are already flagged with the coarse grain element make their group coarse grain.
    """
    if has_only_cg_force_fields(input_force_fields):
        if verbose: print(f' All force fields are coarse grain ({input_force_fields}) -> All atoms are CG')
        return structure.select_all()
    # Get topology bonds, if any
    topology_bonds = mine_topology_bonds(topology_file)
    if topology_bonds and len(topology_bonds) != structure.atom_count:
        raise InputError(f'Topology atom count ({len(topology_bonds)}) does not match the structure atom count ({structure.atom_count})')
    coordinates = get_frames_coordinates(structure, trajectory_file)
    # Decide every group of residues
    groups = get_residue_groups(structure)
    decisions = {}
    undecided_groups = {}
    for group_name, residues in groups.items():
        if any(atom.element == CG_ATOM_ELEMENT for residue in residues for atom in residue.atoms):
            decisions[group_name] = CG
            if verbose: print(f'  {group_name}: atoms are already flagged as CG')
            continue
        votes = [
            vote_by_topology_bonds(residues, coordinates, topology_bonds),
            vote_by_closest_atoms(residues, coordinates),
            vote_by_atom_count(residues),
            vote_by_atom_names(residues),
        ]
        valid_votes = [ vote for vote in votes if vote ]
        if not valid_votes:
            undecided_groups[group_name] = residues
            continue
        counts = Counter(valid_votes)
        decision = valid_votes[0] if counts[CG] == counts[AA] else counts.most_common(1)[0][0]
        decisions[group_name] = decision
        if verbose: print(f'  {group_name}: {decision} (votes: bonds {votes[0]}, closest atoms {votes[1]}, atom count {votes[2]}, atom names {votes[3]})')
    # Groups with no votes left are single atom residues (e.g. ions) so their atoms tell us nothing
    # If the rest of the system agrees then use the same decision
    # In hybrid systems single atom water is coarse grain, since atomistic water has always 3 atoms or more
    # The rest (e.g. ions) are considered atomistic
    if undecided_groups:
        decided_values = set(decisions.values())
        if not decided_values:
            raise InputError('Could not guess if residues are coarse grain or not since they all have a single atom.\n'
                ' Please provide the coarse grain atom selection (VMD syntax) instead of "auto".\n'
                ' Use the "-cg" argument or set the "cg_selection" field in the inputs file.\n'
                ' e.g. -cg "resname W ION"')
        for group_name, residues in undecided_groups.items():
            if len(decided_values) == 1:
                decisions[group_name] = next(iter(decided_values))
                reason = 'like the rest of the system'
            elif residues[0].name in STANDARD_SOLVENT_RESIDUE_NAMES:
                decisions[group_name] = CG
                reason = 'single atom water in a hybrid system'
            else:
                decisions[group_name] = AA
                reason = 'single atom residue in a hybrid system'
            if verbose: print(f'  {group_name}: {decisions[group_name]} ({reason})')
    # Join all coarse grain atoms in a single selection
    atom_indices = [ atom_index
        for group_name, residues in groups.items() if decisions[group_name] == CG
        for residue in residues for atom_index in residue.atom_indices ]
    return Selection(sorted(atom_indices))


# Atoms closer than this distance (Å) are considered bonded when looking for united atoms
UA_BOND_LENGTH = 1.9
# Hydrogens closer than this distance (Å) to a carbon are considered bonded to it
UA_HYDROGEN_BOND_LENGTH = 1.3
# If less than this ratio of carbons which should have hydrogens actually have them then it is united atom
UA_HYDROGENS_RATIO = 0.5


def vote_ua_by_hydrogens(atoms: list['Atom'], coordinates: np.ndarray) -> Optional[bool]:
    """Check if carbons are missing their hydrogens.
    Carbons bonded to 1 or 2 heavy atoms should always have hydrogens in atomistic force fields.
    Note that other carbons may have no hydrogens (e.g. C=O carbons) so they are not checked.
    """
    elements = [ atom.element or atom.guess_element() for atom in atoms ]
    carbon_indices = [ atom.index for atom, element in zip(atoms, elements) if element == 'C' ]
    if not carbon_indices: return None
    heavy_indices = [ atom.index for atom, element in zip(atoms, elements) if element != 'H' ]
    hydrogen_indices = [ atom.index for atom, element in zip(atoms, elements) if element == 'H' ]
    # Count the bonded heavy atoms of every carbon along frames (the carbon itself is also counted so we subtract it)
    heavy_counts = [ [ len(neighbours) - 1 for neighbours in cKDTree(frame[heavy_indices]).query_ball_point(frame[carbon_indices], r=UA_BOND_LENGTH) ]
        for frame in coordinates ]
    must_have_hydrogens = np.median(heavy_counts, axis=0) <= 2
    if not np.any(must_have_hydrogens): return None
    if not hydrogen_indices: return True
    hydrogen_distances = [ cKDTree(frame[hydrogen_indices]).query(frame[carbon_indices])[0] for frame in coordinates ]
    have_hydrogens = np.median(hydrogen_distances, axis=0) < UA_HYDROGEN_BOND_LENGTH
    return np.mean(have_hydrogens[must_have_hydrogens]) < UA_HYDROGENS_RATIO


def guess_ua_selection(
    structure: 'Structure',
    trajectory_file: 'File',
    cg_selection: 'Selection',
    verbose: bool = False,
) -> 'Selection':
    """Guess which atoms are united atoms: carbons with their hydrogens absorbed.
    Residues are split in the same groups used to guess coarse grain (see get_residue_groups).
    Groups including coarse grain atoms are skipped.
    Note that the whole group is selected, including atoms which are not united atoms (e.g. oxygens).
    """
    coordinates = get_frames_coordinates(structure, trajectory_file)
    cg_atom_indices = set(cg_selection.atom_indices)
    atom_indices = []
    for group_name, residues in get_residue_groups(structure).items():
        atoms = [ atom for residue in residues for atom in residue.atoms ]
        if any(atom.index in cg_atom_indices for atom in atoms): continue
        is_united_atom = vote_ua_by_hydrogens(atoms, coordinates)
        if is_united_atom is None: continue
        if verbose: print(f'  {group_name}: {"UA" if is_united_atom else "AA"}')
        if is_united_atom: atom_indices += [ atom.index for atom in atoms ]
    return Selection(sorted(atom_indices))
