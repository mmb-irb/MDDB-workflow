import prody
import json

from pathlib import Path

from Bio import pairwise2
from Bio.SubsMat import MatrixInfo

# Load the reference for residue names and letters
resources = str(Path(__file__).parent.parent / "utils" / "resources")
residues_source = resources + '/residues.json'
residues_reference = None
with open(residues_source, 'r') as file:
    residues_reference = json.load(file)

# Align a reference aminoacid sequence with each chain sequence in a topology
# NEVER FORGET: This system relies on the fact that topology chains are not repeated
def generate_map (ref_sequence : str, topology_filename : str) -> list:
    sequences = get_chain_sequences(topology_filename)
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
def get_chain_sequences (topology_filename : str) -> dict:
    topology = prody.parsePDB(topology_filename)
    sequences = {}
    chains = topology.iterChains()
    for chain in chains:
        name = chain.getChid()
        # Error prone prody alternative
        # sequence = chain.getSequence(allres=True)
        sequence = ''
        for residue in chain.iterResidues():
            resname = residue.getResname()
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

    # Devuelve la segunda sequencia alineada lo mejor que puede y ya
    alignments = pairwise2.align.localds(ref_sequence, new_sequence, MatrixInfo.blosum62, -10, -0.5)

    # In case there are no alignments it means the current chain has nothing to do with this reference
    # Then an array filled with None is returned
    if len(alignments) == 0:
        return [ None for letter in new_sequence ]

    # DANI: En principio siempre hay solo 1, pero está dentro de una array. hay que hacer más pruebas
    # Output format example: '----VNLTT'
    aligned_sequence = alignments[0][1]
    print(aligned_sequence)

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
        aligned_mapping.append(l)
        aligned_index += 1

    return aligned_mapping