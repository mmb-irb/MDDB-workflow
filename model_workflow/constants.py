# CONSTANTS ---------------------------------------------------------------------------

# Database
DEFAULT_API_URL = 'https://mdposit-dev.mddbr.eu/api'

# Selections
# Set a standard selection for protein and nucleic acid backbones in vmd syntax
PROTEIN_AND_NUCLEIC = 'protein or nucleic'
PROTEIN_AND_NUCLEIC_BACKBONE = "(protein and name N CA C) or (nucleic and name P O5' O3' C5' C4' C3')"

# Purely input filenames
DEFAULT_INPUTS_FILENAME = 'inputs.json'
DEFAULT_POPULATIONS_FILENAME = 'populations.json'
DEFAULT_TRANSITIONS_FILENAME = 'transitions.json'

# Input files processing intermediate steps
# We name differenlty every intermediate file and we never rename/overwrite any input or intermediate file
# This allows us to know where we were in case the process was interrupted and not repeat steps on reset
# Intermediate files are removed at the end of the process if it was successful

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

# An old system for when original topology is very wrong and charges must be provided manually
RAW_CHARGES_FILENAME = 'charges.txt'
# Accepted topology formats for atomic charges mining
ACCEPTED_TOPOLOGY_FORMATS = ['top', 'psf', 'prmtop', 'prm7']

# Set generated file names
FIRST_FRAME_FILENAME = 'first_frame.pdb'
AVERAGE_STRUCTURE_FILENAME = 'average.pdb'

# Set output files generated to be uploaded to the database

# Set the output metadata file
OUTPUT_METADATA_FILENAME = 'metadata.json'

# Set the output screenshot filename
OUTPUT_SCREENSHOT_FILENAME = 'screenshot.jpg'

# Set analyses files to be generated
OUTPUT_INTERACTIONS_FILENAME = 'md.interactions.json'
OUTPUT_RMSDS_FILENAME = 'md.rmsds.json'
OUTPUT_TMSCORES_FILENAME = 'md.tmscores.json'
OUTPUT_RMSF_FILENAME = 'md.rmsf.json'
OUTPUT_RGYR_FILENAME = 'md.rgyr.json'
OUTPUT_PCA_FILENAME = 'md.pca.json'
OUTPUT_PCA_PROJECTION_PREFIX = 'pca.trajectory'
OUTPUT_PCA_CONTACTS_FILENAME = 'md.pca_contacts.json'
OUTPUT_RMSD_PERRES_FILENAME = 'md.rmsd.perres.json'
OUTPUT_RMSD_PAIRWISE_FILENAME = 'md.rmsd.pairwise.json'
OUTPUT_DIST_PERRES_FILENAME = 'md.dist.perres.json'
OUTPUT_HBONDS_FILENAME = 'md.hbonds.json'
OUTPUT_SASA_FILENAME = 'md.sasa.json'
OUTPUT_ENERGIES_FILENAME = 'md.energies.json'
OUTPUT_POCKETS_FILENAME = 'md.pockets.json'
OUTPUT_HELICAL_PARAMETERS_FILENAME = 'md.helical.json'
OUTPUT_MARKOV_FILENAME = 'md.markov.json'

# Default parameters
DEFAULT_RMSD_CUTOFF = 9
DEFAULT_INTERACTION_CUTOFF = 0.1

# State all the available checkings, which may be trusted
AVAILABLE_CHECKINGS = [ 'stabonds', 'cohbonds', 'intrajrity' ]
# State all critical process failures, which are to be lethal for the workflow unless mercy is given
AVAILABLE_FAILURES = AVAILABLE_CHECKINGS + [ 'refseq', 'interact' ]