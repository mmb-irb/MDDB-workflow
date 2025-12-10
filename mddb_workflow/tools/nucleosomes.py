from mddb_workflow.tools.generate_map import align
from mddb_workflow.utils.type_hints import *

# Use human nucleosome as reference
# Note that histones are strongly conserved between species
# Thus no matter which species we use as reference the alignment should work
REFERENCE_HISTONE_SEQUENCES = {
    # UniProt Q6FI13
    'H2A': 'MSGRGKQGGKARAKAKSRSSRAGLQFPVGRVHRLLRKGNYAERVGAGAPVYMAAVLEYLTAEILELAGNAARDNKKTRIIPRHLQLAIRNDEELNKLLGKVTIAQGGVLPNIQAVLLPKKTESHHKAKGK',
    # UniProt P33778
    'H2B': 'MPEPSKSAPAPKKGSKKAITKAQKKDGKKRKRSRKESYSIYVYKVLKQVHPDTGISSKAMGIMNSFVNDIFERIAGEASRLAHYNKRSTITSREIQTAVRLLLPGELAKHAVSEGTKAVTKYTSSK',
    # UniProt P84243
    'H3': 'MARTKQTARKSTGGKAPRKQLATKAARKSAPSTGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSAAIGALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA',
    # UniProt P62805
    'H4': 'MSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG',
}

# Set a cutoff for the DNA length in number of nucleotides
DNA_LENGTH_CUTOFF = 100

# Set if the system has a nucleosome
# To do so we compare protein sequences
# Also make sure DNA is long enought to wrap around the histone octamer
def has_nucleosome (structure : 'Structure') -> bool:
    histone_counts = {}
    # Align protein sequences and make sure we have 2 sequences of each histone type
    protein_sequences = structure.get_sequences('protein')
    for sequence in protein_sequences:
        for histone_name, reference_sequence in REFERENCE_HISTONE_SEQUENCES.items():
            if align(sequence, reference_sequence):
                current_count = histone_counts.get(histone_name, 0)
                histone_counts[histone_name] = current_count + 1
                break
    # If at least one of the histones was not found then we are done already
    if len(histone_counts) < 4: return False
    # If at least one of the histones did not have 2 matches then there are no nucleosomes
    for count in histone_counts.values():
        if count < 2: return False
    # So far we assume we have an octamer
    # Now make sure there is DNA around it
    dna_total_sequence_length = 0
    dna_sequences = structure.get_sequences('dna')
    for sequence in dna_sequences:
        dna_total_sequence_length += len(sequence)
    if dna_total_sequence_length < DNA_LENGTH_CUTOFF * 2:
        return False
    # If we made it this far then there may be a nucleosome
    return True