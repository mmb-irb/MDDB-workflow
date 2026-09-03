from mddb_workflow.utils.auxiliar import ranger
from mddb_workflow.utils.type_hints import *

# References with these ids carry no meaningful biological entity and must be excluded
FORBIDDEN_REFERENCES = {'noref', 'notfound'}

# LORE: This was done originally on the fly by the API for the "pointers" endpoint
# LORE: However when we eliminated the storage of global topologies this had to be stored in project data
# LORE: But pre-calculating this and stroing it in the project is a win win situtation
# LORE: Now the API's responses are faster and we can query projects based on these values
def get_coverage_and_presence(residue_map: dict, reference_lengths: dict) -> tuple[dict, dict, dict, dict]:
    """Calculate per-reference coverage, presence, and their residue range strings.

    coverage (COVER): fraction of a protein reference's residues present in the system.
        e.g. if a 500-residue UniProt entry has 300 residues in the simulation -> 0.6
    presence (PRESE): fraction of system residues belonging to a given reference.
        e.g. if 200 out of 1000 system residues belong to protein P12345 -> 0.2
    covered residues (COVRES): range notation of the UniProt residue positions covered.
        e.g. '1-150,200-300'
    present residues (PRESRES): range notation of the system residue indices belonging to this reference.
        e.g. '0-149,300-449'

    All metrics replicate the calculations formerly done on-the-fly by the API using
    the topology collection, which is no longer available globally.  Pre-computing them
    here and storing them in metadata.json makes them queryable and avoids runtime cost.
    """
    map_references = residue_map.get('references')
    map_reference_types = residue_map.get('reference_types')
    map_residue_indices = residue_map.get('residue_reference_indices')
    map_residue_numbers = residue_map.get('residue_reference_numbers')

    cover: dict = {}
    prese: dict = {}
    covres: dict = {}
    presres: dict = {}

    # Nothing to compute if the system has no reference mapping
    if not map_references or not map_residue_indices:
        return cover, prese, covres, presres

    total_residues = len(map_residue_indices)
    if total_residues == 0:
        return cover, prese, covres, presres

    # Collect system residue indices per reference (used for presence and present_residues)
    ref_system_indices: dict[str, list] = {ref: [] for ref in map_references}
    # Collect unique UniProt residue numbers per protein reference (used for coverage and covered_residues)
    ref_covered_numbers: dict[str, set] = {
        ref: set()
        for ref, rtype in zip(map_references, map_reference_types)
        if rtype == 'protein'
    }

    # Single pass over all system residues to populate both accumulators
    for residue_idx, ref_idx in enumerate(map_residue_indices):
        if ref_idx is None:
            continue
        ref_id = map_references[ref_idx]
        ref_system_indices[ref_id].append(residue_idx)
        # residue_reference_numbers gives the 1-based UniProt position within the reference sequence
        if ref_id in ref_covered_numbers and map_residue_numbers:
            num = map_residue_numbers[residue_idx]
            if num is not None:
                ref_covered_numbers[ref_id].add(num)

    for ref_id, ref_type in zip(map_references, map_reference_types):
        if ref_id in FORBIDDEN_REFERENCES:
            continue
        # Presence and present_residues: both proteins and inchikeys
        if ref_type in ('protein', 'inchikey'):
            system_indices = ref_system_indices[ref_id]  # already sorted (iterated in index order)
            prese[ref_id] = round(len(system_indices) / total_residues, 4)
            presres[ref_id] = ranger(system_indices)
        # Coverage and covered_residues: only proteins have a meaningful full-length reference sequence
        if ref_type == 'protein' and ref_id in reference_lengths:
            seq_len = reference_lengths[ref_id]
            if seq_len > 0:
                covered = list(ref_covered_numbers[ref_id])
                cover[ref_id] = round(len(covered) / seq_len, 4)
                covres[ref_id] = ranger(covered)

    return cover, prese, covres, presres
