from mddb_workflow.tools.nucleosomes import has_nucleosome
from mddb_workflow.utils.type_hints import *

# Set the keywords as constants
# Basic keywords
PROTEIN_KEYWORD = 'protein'
DNA_KEYWORD = 'dna'
RNA_KEYWORD = 'rna'
LIPID_KEYWORD = 'lipid'
CARBOHYDRATE_KEYWORD = 'carbohydrate'
COUNTER_ION_KEYWORD = 'counter ion'
NON_COUNTER_ION_KEYWORD = 'non-counter ion'
SOLVENT_KEYWORD = 'solvent'
OTHER_KEYWORD = 'other'
# Only keywords
PROTEIN_ONLY_KEYWORD = 'protein only'
DNA_ONLY_KEYWORD = 'dna only'
RNA_ONLY_KEYWORD = 'rna only'
LIPID_ONLY_KEYWORD = 'lipid only'
CARBOHYDRATE_ONLY_KEYWORD = 'carbohydrate only'
SOLVENT_ONLY_KEYWORD = 'solvent only'
LIGAND_ONLY_KEYWORD = 'ligand only'
# Mapped features
MEMBRANE_KEYWORD = 'membrane'
LIGAND_KEYWORD = 'ligand'
# More complex keywords
NUCLEOSOME_KEYWORD = 'nucleosome'

# System keywords are useful to find different types of systems later in the browser
# Try to assign as many keywords as possible to the current system
def get_system_keywords (
    structure : 'Structure',
    ligand_count : int,
    membrane_count : int,
) -> list[str]:
    keywords = []

    # Basic system keywords inferred from the nature of the residues
    # Protein
    protein_selection = structure.select_protein()
    has_protein = len(protein_selection) > 0
    if has_protein:
        keywords.append(PROTEIN_KEYWORD)
    # DNA
    dna_selection = structure.select_by_classification('dna')
    has_dna = len(dna_selection) > 0
    if has_dna:
        keywords.append(DNA_KEYWORD)
    # RNA
    rna_selection = structure.select_by_classification('rna')
    has_rna = len(rna_selection) > 0
    if has_rna:
        keywords.append(RNA_KEYWORD)
    # Lipids
    lipids_selection = structure.select_lipids()
    has_lipids = len(lipids_selection) > 0
    if has_lipids:
        keywords.append(LIPID_KEYWORD)
    # Carbohydrates
    carbohydrates_selection = structure.select_carbohydrates()
    has_carbohydrates = len(carbohydrates_selection) > 0
    if has_carbohydrates:
        keywords.append(CARBOHYDRATE_KEYWORD)
    # Counter ions
    counter_ions_selection = structure.select_counter_ions()
    has_counter_ions = len(counter_ions_selection) > 0
    if has_counter_ions:
        keywords.append(COUNTER_ION_KEYWORD)
    # Non counter ions
    ions_selection = structure.select_ions()
    non_counter_ions_selection = ions_selection - counter_ions_selection
    has_non_counter_ions = len(non_counter_ions_selection) > 0
    if has_non_counter_ions:
        keywords.append(NON_COUNTER_ION_KEYWORD)
    # Solvent
    solvent_selection = structure.select_water()
    has_solvent = len(solvent_selection) > 0
    if has_solvent:
        keywords.append(SOLVENT_KEYWORD)
    # Other
    other_selection = structure.select_all() - protein_selection - dna_selection - rna_selection \
        - lipids_selection - carbohydrates_selection - ions_selection
    has_other = len(other_selection) > 0
    if has_other:
        keywords.append(OTHER_KEYWORD)

    # Set the "only" keywords
    
    # Set the solvent-only conditions first
    solvent_part_keywords = set([SOLVENT_KEYWORD, COUNTER_ION_KEYWORD, NON_COUNTER_ION_KEYWORD])
    if all( keyword in solvent_part_keywords for keyword in keywords ):
        keywords.append(SOLVENT_ONLY_KEYWORD)

    # The rest of only-keywords do not care about solvent
    relevant_only_keywords = set(keywords) - solvent_part_keywords
    # Thus there must be only one category of the basic types
    # Note that any small molecule falling in 'other' will prevent protein/dna/rna of being "only"
    if len(relevant_only_keywords) == 1:
        if PROTEIN_KEYWORD in relevant_only_keywords:
            keywords.append(PROTEIN_ONLY_KEYWORD)
        elif DNA_KEYWORD in relevant_only_keywords:
            keywords.append(DNA_ONLY_KEYWORD)
        elif RNA_KEYWORD in relevant_only_keywords:
            keywords.append(RNA_ONLY_KEYWORD)
        elif LIPID_KEYWORD in relevant_only_keywords:
            keywords.append(LIPID_ONLY_KEYWORD)
        elif CARBOHYDRATE_KEYWORD in relevant_only_keywords:
            keywords.append(CARBOHYDRATE_ONLY_KEYWORD)
        # If there is only other then we issue no additional keyword

    # Set if a membrane was found
    if membrane_count > 0:
        keywords.append(MEMBRANE_KEYWORD)

    # Set if ligands were found
    if ligand_count > 0:
        keywords.append(LIGAND_KEYWORD)
        # If there is only ligands
        if len(relevant_only_keywords) == 1 and OTHER_KEYWORD in relevant_only_keywords:
            keywords.append(LIGAND_ONLY_KEYWORD)

    # More complicated context-dependent keywords

    # Find out if we have a nucleosome
    if has_protein and has_dna and has_nucleosome(structure):
        keywords.append(NUCLEOSOME_KEYWORD)

    # Display issued system keywords
    print('Issued system keywords: ' + ', '.join(keywords))

    return keywords