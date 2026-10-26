"""Recognise coarse grained lipids and their head and tail beads from residue and bead names.

Bead names follow the Martini lipid conventions (Martini 2 and 3 lipidomes, insane, TS2CG...), so lipids are found
without bead types or bonds. The rule is benchmarked against the parameter sets of Martini and other coarse grained
force fields in test/benchmarks/cg_lipid_naming.
"""
import re
import numpy as np
from mddb_workflow.utils.constants import LIPIDS_RESIDUE_NAMES
from mddb_workflow.utils.type_hints import *

# Bead names of coarse grained (Martini) lipid tails (see coarse_grain_lipid_tails)
# Bead number and chain letter (C1A, D2B, T1A), with the acyl group number in cardiolipins (C1A1)
CG_TAIL_BEAD = re.compile(r'^([CDT])(\d{1,2})([A-F])(\d?)$')
# Chain letter and bead number in branched tails of ionizable lipids (CA1, CB3)
CG_BRANCHED_TAIL_BEAD = re.compile(r'^([CD])([A-D])(\d{1,2})$')
# Bead number without chain letter in single tail lipids, like fatty acids and detergents (C1, C2, D3)
CG_SINGLE_TAIL_BEAD = re.compile(r'^([CD])(\d{1,2})$')
# Head beads of single tail lipids, to tell them from other molecules with C1, C2... beads
CG_SINGLE_TAIL_HEADS = {'COO', 'COOH', 'NC3', 'NH3', 'SO3', 'SO4', 'TAP', 'PEO1', 'B1'}
# Sterol ring beads, and sterol head beads: hydroxyl ROH or R0 when esterified or glycosylated
CG_RING_BEAD = re.compile(r'^R\d{1,2}$')
CG_STEROL_HEADS = ('ROH', 'R0')


def coarse_grain_lipids(universe: 'Universe', cg_residues: list[int], forced_lipids: Optional[list[str]] = None
        ) -> tuple[Optional['AtomGroup'], Optional['AtomGroup']]:
    """Find coarse grained lipids and get their head beads (one per lipid) and their tail beads.

    Lipids are residues with a known or forced lipid residue name, and coarse grained residues whose bead names follow
    the Martini lipid conventions (see coarse_grain_lipid_tails), which also give the tail beads.
    """
    lipid_names = list(set(LIPIDS_RESIDUE_NAMES) | set(forced_lipids or []))
    named = np.isin(universe.residues.resnames, lipid_names)
    coarse_grained = np.zeros(len(universe.residues), dtype=bool)
    coarse_grained[np.asarray(cg_residues, dtype=int)] = True
    # Lipids have at least a head and two tail beads
    n_atoms = np.bincount(universe.atoms.resindices, minlength=len(universe.residues))
    candidates = np.flatnonzero((named | coarse_grained) & (n_atoms >= 3))
    # Atom indices of every residue
    residue_atoms = np.argsort(universe.atoms.resindices, kind='stable')
    starts = np.concatenate([[0], np.cumsum(n_atoms)[:-1]])
    atom_names = universe.atoms.names
    known_tails = {}
    heads, tails = [], []
    n_by_beads = 0
    for resindex in candidates:
        atoms = residue_atoms[starts[resindex]:starts[resindex] + n_atoms[resindex]]
        names = tuple(atom_names[atoms])
        if names not in known_tails: known_tails[names] = coarse_grain_lipid_tails(names)
        tail_mask = known_tails[names]
        if not named[resindex]:
            if tail_mask is None: continue
            n_by_beads += 1
        heads.append(atoms[coarse_grain_head(names, tail_mask)])
        if tail_mask is not None: tails.extend(atoms[tail_mask])
    if len(heads) == 0:
        print("No coarse-grained lipid atoms found in the structure.")
        return None, None
    print(f"Found {len(heads)} lipid residues: {len(heads) - n_by_beads} by residue name and {n_by_beads} by bead names.")
    return universe.atoms[heads], universe.atoms[tails]


def tail_bead(name: str) -> Optional[tuple[str, int]]:
    """Get the chain and the bead number of a bead named like a lipid tail chain bead (C1A, D2B, C1A1, CA1), or None."""
    match = CG_TAIL_BEAD.match(name)
    if match: return match.group(3) + match.group(4), int(match.group(2))
    match = CG_BRANCHED_TAIL_BEAD.match(name)
    if match: return match.group(2), int(match.group(3))
    return None


def coarse_grain_lipid_tails(names: list[str]) -> Optional[np.ndarray]:
    """Get the tail beads of a coarse grained lipid from its bead names, as a mask, or None if it is not a lipid.

    Bead names follow the Martini lipid conventions (Martini 2 and 3 lipidomes, insane, TS2CG...):
    - Fatty lipids have a head or linker bead and tail beads named by number and chain (C1A, D2A... C1B...).
      A chain A or B must start at bead 1, with at least 3 beads in a chain or 2 beads in at least two chains
      (C8-C10 tails), but in no more than 3 chains (the pyrrole rings A-D of chlorophyll are not tails).
    - Sterols have an ROH (or R0) head bead and at least 3 ring beads (R1, R2...). Their tails are the ring and C beads.
    - Single tail lipids (fatty acids and detergents) have a known head bead and a tail C1, C2, C3... without chain letter.
    """
    names = list(names)
    beads = [tail_bead(name) for name in names]
    tails = np.array([bead is not None for bead in beads], dtype=bool)
    chains = {}
    for bead in beads:
        if bead: chains.setdefault(bead[0], set()).add(bead[1])
    # Only chains starting at bead 1 count, and there must be a chain A or B
    started = [(chain, len(numbers)) for chain, numbers in chains.items() if 1 in numbers]
    if not tails.all() and any(chain[0] in 'AB' and length >= 2 for chain, length in started):
        lengths = sorted((length for _, length in started), reverse=True)
        if lengths[0] >= 3 or (len(lengths) >= 2 and lengths[1] >= 2 and len(chains) <= 3): return tails
    rings = np.array([bool(CG_RING_BEAD.match(name)) for name in names])
    singles = np.array([bool(CG_SINGLE_TAIL_BEAD.match(name)) for name in names])
    sterol_rings = rings & ~np.isin(names, CG_STEROL_HEADS)
    if any(head in names for head in CG_STEROL_HEADS) and sterol_rings.sum() >= 3: return sterol_rings | singles
    numbers = [int(CG_SINGLE_TAIL_BEAD.match(name).group(2)) for name in np.array(names)[singles]]
    if len(numbers) >= 3 and min(numbers) == 1 and CG_SINGLE_TAIL_HEADS & set(names) and not rings.any() and 'BB' not in names:
        return singles
    return None


def coarse_grain_head(names: list[str], tails: Optional[np.ndarray]) -> int:
    """Get the position of the head bead of a coarse grained lipid.

    It is the phosphate bead PO4, the sterol head bead ROH or R0, the only bead named P* (other force fields)
    or the first bead which is not a tail bead, since Martini lipid beads go from the head to the tails.
    """
    names = list(names)
    if names.count('PO4') == 1: return names.index('PO4')
    for head in CG_STEROL_HEADS:
        if head in names: return names.index(head)
    phosphates = [i for i, name in enumerate(names) if name.startswith('P')]
    if len(phosphates) == 1: return phosphates[0]
    if tails is not None and not tails.all(): return int(np.argmin(tails))
    return 0


def coarse_grain_tail_chains(names: list[str]) -> list[list[str]]:
    """Get the tail bead names of a coarse grained lipid by chain, sorted by chain (e.g. [['C1A', 'D2A', 'C3A'], ['C1B', 'C2B']]).

    Chains are the beads named by number and chain (see tail_bead), so sterol rings and single tails without chain
    letter have no chains.
    """
    chains = {}
    for name in names:
        bead = tail_bead(name)
        if bead: chains.setdefault(bead[0], []).append(name)
    return [chains[chain] for chain in sorted(chains)]


def coarse_grain_head_tail_split(residue: 'Residue') -> tuple['AtomGroup', 'AtomGroup']:
    """Split a coarse grained lipid residue into its head and tail beads (see coarse_grain_lipid_tails).

    Lipids whose bead names do not follow the conventions (e.g. lipids found by residue name with too short tails)
    have as tails the beads named like tail chain beads, if any.
    """
    names = residue.atoms.names
    tails = coarse_grain_lipid_tails(names)
    if tails is None: tails = np.array([tail_bead(name) is not None for name in names], dtype=bool)
    return residue.atoms[~tails], residue.atoms[tails]
