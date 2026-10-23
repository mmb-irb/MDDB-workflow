from mddb_workflow.tools.nucleosomes import has_nucleosome
from mddb_workflow.tools.multimerism import get_multimerism_keywords
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
# Atomic vs United Atom vs Coarse Grain labels
# Abreviations are added in the keyword to in case a user makes a query with them
AA_KEYWORD = 'all-atom (aa)'
UA_KEYWORD = 'united atom (ua)'
CG_KEYWORD = 'coarse grained (cg)'
# Hybrid labels saying which resolutions are combined
# DANI: "multiscale" me gusta mucho pero puede ser ambiguo
HYBRID_KEYWORDS = {
    frozenset(['aa', 'ua']): 'hybrid resolution aa/ua',
    frozenset(['aa', 'cg']): 'hybrid resolution aa/cg',
    frozenset(['ua', 'cg']): 'hybrid resolution ua/cg',
    frozenset(['aa', 'ua', 'cg']): 'hybrid resolution aa/ua/cg',
}


# Changes in the function can be tested with: pytest -k "TestRunAll and pmeta" -s
def get_system_keywords(
    structure: 'Structure',
    inchikey_map: list[dict],
    lipids_selection: 'Selection',
    membrane_map: dict,
    cg_selection: 'Selection',
    ua_selection: 'Selection',
    protein_map: list[dict],
    chain_contacts: dict,
) -> list[str]:
    """System keywords are useful to find different types of systems later in the browser
    Try to assign as many keywords as possible to the current system.
    """
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
    # Ligands
    ligands_selection = structure.select_ligands(inchikey_map)
    has_ligand = bool(ligands_selection)
    if has_ligand:
        keywords.append(LIGAND_KEYWORD)
    # Lipids
    # WARNING: Do not rely in the structure logic for lipid identification, it is precarious
    # WARNING: Instead we use selection made out ut lipidmaps references
    has_lipids = len(lipids_selection) > 0
    if has_lipids:
        keywords.append(LIPID_KEYWORD)
    # Carbohydrates
    carbohydrates_selection = structure.select_carbohydrates()
    # If there is only one carbohydrate, treat it as a ligand
    if len(structure.get_selection_residues(carbohydrates_selection)) == 1:
        carbohydrates_selection = carbohydrates_selection - ligands_selection
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
        - lipids_selection - carbohydrates_selection - ions_selection - ligands_selection - solvent_selection
    has_other = len(other_selection) > 0
    if has_other:
        keywords.append(OTHER_KEYWORD)

    # Set the "only" keywords
    # Set the solvent-only conditions first
    solvent_part_keywords = set([SOLVENT_KEYWORD, COUNTER_ION_KEYWORD, NON_COUNTER_ION_KEYWORD])
    if all(keyword in solvent_part_keywords for keyword in keywords):
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
        elif LIGAND_KEYWORD in relevant_only_keywords:
            keywords.append(LIGAND_ONLY_KEYWORD)
        # If there is only other then we issue no additional keyword

    # Set if a membrane was found
    membrane_count = membrane_map['n_mems']
    has_membrane = membrane_count > 0
    if has_membrane:
        keywords.append(MEMBRANE_KEYWORD)

    # More complicated context-dependent keywords

    # Find out if we have a nucleosome
    if has_protein and has_dna and has_nucleosome(structure):
        keywords.append(NUCLEOSOME_KEYWORD)

    # Find out protein multimerism (e.g. homodimer), nucleic strands (e.g. double strand)
    # and protein-nucleic complexes according to chain contacts along the trajectory
    keywords += get_multimerism_keywords(structure, protein_map, chain_contacts)

    # Set atomic vs united atom vs coarse-grained keywords
    # Atomistic water and ions are not taken in count, since they are also used in united atom systems
    aa_selection = structure.select_all() - cg_selection - ua_selection - solvent_selection - ions_selection
    resolutions = set()
    if aa_selection: resolutions.add('aa')
    if ua_selection: resolutions.add('ua')
    if cg_selection: resolutions.add('cg')
    if len(resolutions) > 1: keywords.append(HYBRID_KEYWORDS[frozenset(resolutions)])
    elif 'cg' in resolutions: keywords.append(CG_KEYWORD)
    elif 'ua' in resolutions: keywords.append(UA_KEYWORD)
    else: keywords.append(AA_KEYWORD)

    # Display issued system keywords
    print('Issued system keywords: ' + ', '.join(keywords))

    return keywords
