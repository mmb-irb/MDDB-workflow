# CONSTANTS ---------------------------------------------------------------------------

# Database
DEFAULT_API_URL = 'https://mmb-dev.mddbr.eu/api'

# Selections
# Set a standard selection for protein and nucleic acid backbones in vmd syntax
PROTEIN_AND_NUCLEIC = 'protein or nucleic'
PROTEIN_AND_NUCLEIC_BACKBONE = "(protein and name N CA C) or (nucleic and name P O5' O3' C5' C4' C3')"

# Inputs file
DEFAULT_INPUTS_FILENAME = 'inputs.yaml'
ACCEPTED_INPUT_FILENAMES = [
    DEFAULT_INPUTS_FILENAME, # The default
    'inputs.yml', # Another extension of yaml files
    'inputs.json' # Legacy inputs file
]

# Markov State Model input filenames
DEFAULT_POPULATIONS_FILENAME = 'populations.json'
DEFAULT_TRANSITIONS_FILENAME = 'transitions.json'

# Input files processing intermediate steps
# We name differenlty every intermediate file and we never rename/overwrite any input or intermediate file
# This allows us to know where we were in case the process was interrupted and not repeat steps on reset
# Intermediate files are removed at the end of the process if it was successful

INCOMPLETE_PREFIX = '.incomplete_'

CONVERTED_STRUCTURE = 'converted.pdb'
CONVERTED_TRAJECTORY = 'converted.xtc'

FILTERED_STRUCTURE = 'filtered.pdb'
FILTERED_TRAJECTORY = 'filtered.xtc'

IMAGED_STRUCTURE = 'imaged.pdb'
IMAGED_TRAJECTORY = 'imaged.xtc'

# Input and output core files
TOPOLOGY_FILENAME = 'topology.json'
STRUCTURE_FILENAME = 'structure.pdb'
TRAJECTORY_FILENAME = 'trajectory.xtc'

# Intermediate filenames
REGISTER_FILENAME = '.register.json'

# Files saving resorted bonds and charges when we have to resort atoms
# Note that these files have priority when loading both bonds and charges
RESORTED_CHARGES_FILENAME = 'resorted_charges.json'
RESORTED_BONDS_FILENAME = 'resorted_bonds.json'

# An old system for when original topology is very wrong and charges must be provided manually
RAW_CHARGES_FILENAME = 'charges.txt'
# Accepted topology formats for atomic charges mining
ACCEPTED_TOPOLOGY_FORMATS = ['tpr', 'top', 'psf', 'prmtop', 'prm7']

# Set generated file names
FIRST_FRAME_FILENAME = 'first_frame.pdb'
AVERAGE_STRUCTURE_FILENAME = 'average.pdb'

# Set the reference labels according to the reference file used
REFERENCE_LABELS = {
    FIRST_FRAME_FILENAME: 'firstframe',
    AVERAGE_STRUCTURE_FILENAME: 'average'
}

# Set output files generated to be uploaded to the database

# Set the references filename
OUTPUT_REFERENCES_FILENAME = 'references.json'
# Set the ligands filename
OUTPUT_LIGANDS_FILENAME = 'ligands.json'

# Set the metadata filename
OUTPUT_METADATA_FILENAME = 'metadata.json'

# Set the screenshot filename
OUTPUT_SCREENSHOT_FILENAME = 'mdf.screenshot.jpg'

# Set analyses files to be generated
OUTPUT_INTERACTIONS_FILENAME = 'mda.interactions.json'
OUTPUT_RMSDS_FILENAME = 'mda.rmsds.json'
OUTPUT_TMSCORES_FILENAME = 'mda.tmscores.json'
OUTPUT_RMSF_FILENAME = 'mda.rmsf.json'
OUTPUT_RGYR_FILENAME = 'mda.rgyr.json'
OUTPUT_PCA_FILENAME = 'mda.pca.json'
OUTPUT_PCA_PROJECTION_PREFIX = 'mdt.pca_trajectory'
OUTPUT_PCA_CONTACTS_FILENAME = 'mda.pca_contacts.json'
OUTPUT_RMSD_PERRES_FILENAME = 'mda.rmsd_perres.json'
OUTPUT_RMSD_PAIRWISE_FILENAME = 'mda.rmsd_pairwise.json'
OUTPUT_DIST_PERRES_FILENAME = 'mda.dist_perres.json'
OUTPUT_HBONDS_FILENAME = 'mda.hbonds.json'
OUTPUT_SASA_FILENAME = 'mda.sasa.json'
OUTPUT_ENERGIES_FILENAME = 'mda.energies.json'
OUTPUT_POCKETS_FILENAME = 'mda.pockets.json'
OUTPUT_POCKET_STRUCTURES_PREFIX = 'mdf.pocket' # WARNING: If this is changed then the pockets function must be updated as well
OUTPUT_HELICAL_PARAMETERS_FILENAME = 'mda.helical.json'
OUTPUT_MARKOV_FILENAME = 'mda.markov.json'

# Set folder names for some analyses which generate a lot of intermediate step files
ENERGIES_FOLDER = 'energies'
POCKETS_FOLDER = 'mdpocket'

# Set problematic signs for directory/folder names
FORBIDEN_DIRECTORY_CHARACTERS = ['.', ',', ';', ':']

# Default parameters
DEFAULT_RMSD_CUTOFF = 9
DEFAULT_INTERACTION_CUTOFF = 0.1

# Set register cache flags
SNAPSHOTS_FLAG = 'snapshots'

# Set the different test flags
STABLE_BONDS_FLAG = 'stabonds'
COHERENT_BONDS_FLAG = 'cohbonds'
TRAJECTORY_INTEGRITY_FLAG = 'intrajrity'
CORRECT_ELEMENTS = 'elements'
REFERENCE_SEQUENCE_FLAG = 'refseq'
STABLE_INTERACTIONS_FLAG = 'interact'

# State all the available checkings, which may be trusted
AVAILABLE_CHECKINGS = [ STABLE_BONDS_FLAG, COHERENT_BONDS_FLAG, TRAJECTORY_INTEGRITY_FLAG ]
# State all critical process failures, which are to be lethal for the workflow unless mercy is given
AVAILABLE_FAILURES = AVAILABLE_CHECKINGS + [ CORRECT_ELEMENTS, REFERENCE_SEQUENCE_FLAG, STABLE_INTERACTIONS_FLAG ]

# Set which tests are to be run when some input files are modified
STRUCTURE_TESTS = [STABLE_BONDS_FLAG, COHERENT_BONDS_FLAG]
TRAJECTORY_TESTS = [STABLE_BONDS_FLAG, TRAJECTORY_INTEGRITY_FLAG]
TOPOLOGY_TESTS = [STABLE_BONDS_FLAG, COHERENT_BONDS_FLAG]

# Terminal colors
# https://stackoverflow.com/questions/287871/how-do-i-print-colored-text-to-the-terminal
GREEN_HEADER = '\033[92m'
CYAN_HEADER = '\033[96m'
YELLOW_HEADER = '\033[93m'
RED_HEADER = '\033[91m'
GREY_HEADER = '\033[90m'
COLOR_END = '\033[0m'

# Set a dictionary to parse an internal raw name to a pretty human firendly name
NICE_NAMES = {
    STABLE_BONDS_FLAG: 'Stable bonds test',
    COHERENT_BONDS_FLAG: 'Coherent bonds test',
    TRAJECTORY_INTEGRITY_FLAG: 'Trajectory integrity test',
    CORRECT_ELEMENTS: 'Correct elements',
    REFERENCE_SEQUENCE_FLAG: 'Reference sequence match',
    STABLE_INTERACTIONS_FLAG: 'Interactions are stable'
}

# Set the "standard" file format of every possible file extension
# Note that some formats have different possible extension (e.g. nc, cdf, netcdf)
EXTENSION_FORMATS = {
    # Topologies
    'tpr': 'tpr',
    'top': 'top',
    'psf': 'psf',
    'prmtop': 'prmtop',
    'prm7': 'prmtop',
    'txt': 'txt', # charges.txt
    # Structures
    'pdb': 'pdb',
    'gro': 'gro',
    # Trajectories
    'xtc': 'xtc',
    'trr': 'trr',
    'dcd': 'dcd',
    'nc': 'nc',
    'cdf': 'nc',
    'netcdf': 'nc',
    'crd': 'crd',
    'mdcrd': 'crd',
    'trj': 'crd',
    # Other
    'json': 'json',
    'yaml': 'yaml',
    'yml': 'yaml',
    'npy': 'npy',
    'in': 'txt'
}

# Topology and trajectory file formats supported by PyTraj
PYTRAJ_SUPPORTED_FORMATS = set([
    # Topologies
    'prmtop', 'top', 'psf', 'pdb'
    # Trajectories
    'nc', 'crd', 'dcd', 'trr', 'xtc'
])

# From GitHub:
# ParmFormatDict = {
#     "AMBERPARM": AMBERPARM,
#     "PDBFILE": PDBFILEPARM,
#     "MOL2FILE": MOL2FILEPARM,
#     "CHARMMPSF": CHARMMPSF,
#     "CIFFILE": CIFFILE,
#     "GMXTOP": GMXTOP,
#     "SDFFILE": SDFFILE,
#     "TINKER": TINKERPARM,
#     "UNKNOWN_PARM": UNKNOWN_PARM,
# }

# Set some flags requeired to write files with pytraj
PYTRAJ_PARM_FORMAT = {
    'prmtop': 'AMBERPARM',
    'psf': 'CHARMMPSF',
    'top': 'GMXTOP',
    'pdb': 'PDBFILE'
}

# Elements supported while correcting atom elements
SUPPORTED_POLYMER_ELEMENTS = set([ 'C', 'N', 'O', 'H', 'P', 'S' ])
SUPPORTED_ION_ELEMENTS = set([ 'K', 'F', 'Cl', 'Na', 'Zn', 'Mg', 'Fe', 'Br', 'Mn', 'I', 'Ca' ])
SUPPORTED_ELEMENTS = SUPPORTED_POLYMER_ELEMENTS.union(SUPPORTED_ION_ELEMENTS)

# Set a dictionary with all residue names and their equivalent letters
RESIDUE_NAME_LETTERS = {
    # Amino acids
    'ALA':'A',
    'ALAN':'A',
    'ALAC':'A',
    'ARG':'R',
    'ARGN':'R',
    'ARGC':'R',
    'ASN':'N',
    'ASNN':'N',
    'ASNC':'N',
    'ASP':'D',
    'ASPN':'D',
    'ASPC':'D',
    'CYS':'C',
    'CYSN':'C',
    'CYSC':'C',
    'CYH':'C',
    'CSH':'C',
    'CSS':'C',
    'CYX':'C',
    'CYP':'C',
    'GLN':'Q',
    'GLNN':'Q',
    'GLNC':'Q',
    'GLU':'E',
    'GLUN':'E',
    'GLUC':'E',
    'GLY':'G',
    'GLYN':'G',
    'GLYC':'G',
    'HIS':'H',
    'HISN':'H',
    'HISC':'H',
    'HID':'H',
    'HIE':'H',
    'HIP':'H',
    'HSD':'H',
    'HSE':'H',
    'ILE':'I',
    'ILEN':'I',
    'ILEC':'I',
    'ILU':'I',
    'LEU':'L',
    'LEUN':'L',
    'LEUC':'L',
    'LYS':'K',
    'LYSN':'K',
    'LYSC':'K',
    'MET':'M',
    'METN':'M',
    'METC':'M',
    'PHE':'F',
    'PHEN':'F',
    'PHEC':'F',
    'PRO':'P',
    'PRON':'P',
    'PROC':'P',
    'PRØ':'P',
    'PR0':'P',
    'PRZ':'P',
    'SER':'S',
    'SERN':'S',
    'SERC':'S',
    'THR':'T',
    'THRN':'T',
    'THRC':'R',
    'TRP':'W',
    'TRPN':'W',
    'TRPC':'W',
    'TRY':'W',
    'TYR':'Y',
    'TYRN':'Y',
    'TYRC':'Y',
    'VAL':'V',
    'VALN':'V',
    'VALC':'V',
    # Nucleotides
    'A': 'A',
    'A3': 'A',
    'A5': 'A',
    'DA': 'A',
    'RA': 'A',
    'C': 'C',
    'C3': 'C',
    'C5': 'C',
    'DC': 'C',
    'RC': 'C',
    'T': 'T',
    'T3': 'T',
    'T5': 'T',
    'DT': 'T',
    'G': 'G',
    'G3': 'G',
    'G5': 'G',
    'DG': 'G',
    'RG': 'G',
    'U': 'U',
    'U3': 'U',
    'U5': 'U',
    'RU': 'U',
    'DA3': 'A',
    'DA5': 'A',
    'DT3': 'T',
    'DT5': 'T',
    'DC3': 'C',
    'DC5': 'C',
    'DG3': 'G',
    'DG5': 'G',
}

# Set typical residue names to guess what residues are
STANDARD_SOLVENT_RESIDUE_NAMES = {'SOL', 'WAT', 'HOH', 'TIP', 'TP3'}
STANDARD_COUNTER_ION_ATOM_NAMES = {'K', 'NA', 'CL', 'CLA', 'SOD', 'POT'}

# Topology flags

# Set a flag to represent a protein which is not referable (e.g. antibodies, synthetic constructs)
NO_REFERABLE_FLAG = 'noref'

# Set a flag to represent a not found reference
NOT_FOUND_FLAG = 'notfound'