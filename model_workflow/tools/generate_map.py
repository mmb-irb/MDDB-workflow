import json

from pathlib import Path

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo

# from Bio import Align
# from Bio.Align import substitution_matrices

# # Set the aligner
# aligner = Align.PairwiseAligner(mode='local')
# aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
# aligner.open_gap_score = -10
# aligner.extend_gap_score = -0.5

# Load the references
# Load the reference for residue names and letters
resources = str(Path(__file__).parent.parent / "utils" / "resources")
references_source = resources + '/references.json'
residues_source = resources + '/residues.json'
residues_reference = None
with open(residues_source, 'r') as file:
    residues_reference = json.load(file)
stored_references = None
with open(references_source, 'r') as file:
    stored_references = json.load(file)

# For each stored topology reference, align the reference sequence with the topology sequence
def generate_map (structure : 'Structure') -> dict:
    # For each input topology reference, align the reference sequence with the topology sequence
    # Each 'reference' must include the reference name and the align map
    references = []
    for reference in stored_references:
        # Find the corresponding stored topology reference and get its residues sequence
        name = reference['name']
        reference_sequence = reference['sequence']
        reference_map = map_sequence(reference_sequence, structure)
        if all(residue == None for residue in reference_map):
            continue
        references.append({
            'name': name,
            'map': reference_map,
        })

    # Reformat data according to a new system (introduced later)
    residues_count = len(structure.residues)
    reference_names = []
    residue_reference_indices = [ None ] * residues_count
    residue_reference_numbers = [ None ] * residues_count
    for index, reference in enumerate(references):
        reference_names.append(reference['name'])
        for r, residue in enumerate(reference['map']):
            if residue == None:
                continue
            residue_reference_indices[r] = index
            residue_reference_numbers[r] = residue
    # If there are not references at the end then set all feilds as None, in order to save space
    if len(reference_names) == 0:
        reference_names = None
        residue_reference_indices = None
        residue_reference_numbers = None

    residues_map = {
        'references': reference_names,
        'residue_reference_indices': residue_reference_indices,
        'residue_reference_numbers': residue_reference_numbers,
    }
    return residues_map

# Align a reference aminoacid sequence with each chain sequence in a topology
# NEVER FORGET: This system relies on the fact that topology chains are not repeated
def map_sequence (ref_sequence : str, structure : 'Structure') -> list:
    sequences = get_chain_sequences(structure)
    mapping = []
    for s in sequences:
        sequence = sequences[s]
        sequence_map = align(ref_sequence, sequence)
        mapping += sequence_map
    return mapping

# Get each chain name and aminoacids sequence in a topology
# Output format example: {'A':'VNLTT', 'B':'SVASQ', ...}
# WARNING: Prody sometimes splits chains with no reason
# All chains with the same name are joined together to avoid this problem, even when the split is real
# WARNING: Prody's 'getSequence' tool sometimes adds wrong aminoacid characters
# e.g. Glycan residues '2MA' and 'YMA' are taken as adenine and glycine respectively
# This is due to the prody custom map in 'mod_res_map.dat'
# For this reason, another tool is used to obtain the residues sequence
def get_chain_sequences (structure : 'Structure') -> dict:
    sequences = {}
    chains = structure.chains
    for chain in chains:
        name = chain.name
        sequence = ''
        for residue in chain.residues:
            resname = residue.name
            letter = resname_2_letter(resname)
            sequence += letter
        # Save sequences by chain name (i.e. chain id or chain letter)
        if sequences.get(name, False):
            sequences[name] += sequence
        else:
            sequences[name] = sequence
    return sequences

# Set a function to transform residue names in residue letters
# If the residue name is not recognized then return 'X'
# e.g. 'ARG' -> 'R', 'WTF' -> 'X'
def resname_2_letter(resname : str) -> str:
    ref = residues_reference.get(resname, False)
    return ref if ref else 'X'

# Align two aminoacid sequences
# Return a list with the reference residue indexes (values)
# which match each new sequence residues indexes (indexes)
def align (ref_sequence : str, new_sequence : str) -> list:

    #print('- REFERENCE\n' + ref_sequence + '\n- NEW\n' + new_sequence)

    # If the new sequence is all 'X' stop here, since this would make the alignment infinite
    # Then an array filled with None is returned
    if all([ letter == 'X' for letter in new_sequence ]):
        return [ None for letter in new_sequence ]

    # Return the new sequence as best aligned as possible with the reference sequence
    alignments = pairwise2.align.localds(ref_sequence, new_sequence, MatrixInfo.blosum62, -10, -0.5)
    # DANI: Habría que hacerlo de esta otra forma según el deprecation warning (arriba hay más código)
    # DANI: El problema es que el output lo tiene todo menos la sequencia en formato alienada
    # DANI: i.e. formato '----VNLTT', que es justo el que necesito
    #alignments = aligner.align(ref_sequence, new_sequence)

    # In case there are no alignments it means the current chain has nothing to do with this reference
    # Then an array filled with None is returned
    if len(alignments) == 0:
        return [ None for letter in new_sequence ]

    # Several alignments may be returned, specially when it is a difficult or impossible alignment
    # Output format example: '----VNLTT'
    best_alignment = alignments[0]
    aligned_sequence = best_alignment[1]
    print(format_alignment(*alignments[0]))
    score = alignments[0][2]
    # WARNING: Do not use 'aligned_sequence' length here since it has the total sequence length
    normalized_score = score / len(new_sequence)
    print('Normalized score: ' + str(normalized_score))

    # If the normalized score does not reaches the minimum we consider the alignment is not valid
    # It may happen when the reference goes for a specific chain but we must map all chains
    # This 1 has been found experimentally
    # Non maching sequence may return a 0.1-0.3 normalized score
    # Matching sequence may return >4 normalized score
    if normalized_score < 1:
        print('Not valid alignment')
        return [ None for letter in new_sequence ]

    # Match each residue
    aligned_mapping = []
    aligned_index = 0
    for l, letter in enumerate(aligned_sequence):
        # Guions are skipped
        if letter == '-':
            continue
        # Get the current residue of the new sequence
        equivalent_letter = new_sequence[aligned_index]
        if not letter == equivalent_letter:
            raise SystemExit('Something was wrong :S')
        # 'X' residues cannot be mapped since reference sequences should never have any 'X'
        if letter == 'X':
            aligned_mapping.append(None)
            aligned_index += 1
            continue
        # Otherwise add the equivalent aligned index to the mapping
        # WARNING: Add +1 since uniprot residue counts start at 1, not 0
        aligned_mapping.append(l + 1)
        aligned_index += 1

    return aligned_mapping