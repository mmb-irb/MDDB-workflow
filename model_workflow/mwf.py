#!/usr/bin/env python

# This is the starter script

# Import python libraries
from argparse import ArgumentParser, RawTextHelpFormatter
import os
import sys
import io
import math
from pathlib import Path
import urllib.request
import json

# Import local tools
from model_workflow.tools.vmd_processor import vmd_processor
from model_workflow.tools.filter_atoms import filter_atoms
from model_workflow.tools.image_and_fit import image_and_fit
from model_workflow.tools.topology_corrector import topology_corrector
from model_workflow.tools.process_input_files import process_input_files, find_charges_filename, get_output_charges_filename
from model_workflow.tools.topology_manager import setup_structure
from model_workflow.tools.get_pytraj_trajectory import get_pytraj_trajectory, get_reduced_pytraj_trajectory
from model_workflow.tools.get_first_frame import get_first_frame
from model_workflow.tools.get_average import get_average
from model_workflow.tools.process_interactions import process_interactions
from model_workflow.tools.generate_metadata import generate_metadata
from model_workflow.tools.generate_map import generate_map_online
from model_workflow.tools.generate_topology import generate_topology
from model_workflow.tools.get_summarized_trajectory import get_summarized_trajectory
from model_workflow.tools.get_frames_count import get_frames_count
from model_workflow.tools.get_charges import get_charges
from model_workflow.tools.remove_trash import remove_trash

# Import local analyses
from model_workflow.analyses.rmsds import rmsds
from model_workflow.analyses.tmscores import tmscores
from model_workflow.analyses.rmsf import rmsf
from model_workflow.analyses.rgyr import rgyr
from model_workflow.analyses.pca import pca
from model_workflow.analyses.pca_contacts import pca_contacts
from model_workflow.analyses.rmsd_per_residue import rmsd_per_residue
from model_workflow.analyses.rmsd_pairwise import rmsd_pairwise
from model_workflow.analyses.distance_per_residue import distance_per_residue
from model_workflow.analyses.hydrogen_bonds import hydrogen_bonds
from model_workflow.analyses.sasa import sasa
from model_workflow.analyses.energies import energies
from model_workflow.analyses.pockets import pockets

# Make the system output stream to not be buffered
# This is useful to make prints work on time in Slurm
# Otherwise, output logs are written after the script has fully run
# Note that this fix affects all modules and built-ins
unbuffered = io.TextIOWrapper(open(sys.stdout.fileno(), 'wb', 0), write_through=True)
sys.stdout = unbuffered

# Set a standard selection for protein and nucleic acid backbones in vmd syntax
protein_and_nucleic_backbone = "(protein and name N CA C) or (nucleic and name P O5' O3' C5' C4' C3')"

# Provisional fix for SSL bypass
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

# CLASSES -------------------------------------------------------------------------

class Dependency:
    _depend = True
    # The 'func' is the function to get the value
    # The 'args' are the arguments for the 'func'
    # The 'filename' is the name of the file which is expected to be required for this dependency
    # When a filename is passed the dependecy 'value' is replaced by this 'filename' when exists
    def __init__ (self, func, args = {}, alias = ''):
        self.func = func
        self.args = args
        self.alias = alias
        self._value = None

    # For each 'Dependency' instance in args call its dependency function
    def parse_args(self):

        # For a given non Dependency instance, return the value as it is
        # For a given Dependency instance, get the real value of it
        # Getting the value means also calling its dependency function
        def parse(value):
            # When it is a Dependency or an inheritor class
            if hasattr(value, '_depend'):
                parsed_value = value.value
                return parsed_value
            return value
        
        # Now parse each arg
        parsed_args = {}
        for attr, value in self.args.items():
            # If it is a list we must check each field
            if type(value) == list:
                parsed_args[attr] = [ parse(element) for element in value ]
            # Otherwise, simply parse the value
            else:
                parsed_args[attr] = parse(value)
        return parsed_args

    # Get the dependency value
    # Execute the corresponding function if we do not have the value yet
    # Once the value has been calculated save it as an internal variabel
    def get_value(self):
        value = self._value
        if value:
            return value
        parsed_args = self.parse_args()
        if self.alias:
            sys.stdout.write('Running "' + self.alias + '"\n')
        value = self.func(**parsed_args)
        self._value = value # Save the result for possible further use
        return value

    value = property(get_value, None, None, "The dependency value")

class File(Dependency):
    # The 'func' is the function to generate the file
    # If there is no func it means the file will be always there (i.e. it is an input file)
    # The 'args' are the arguments for the 'func'
    def __init__ (self, filename : str, func, args = {}, alias = ''):
        self.filename = filename
        super().__init__(func, args, alias)

    def exists (self) -> bool:
        if not self.filename:
            return False
        return os.path.exists(self.filename)

    # Get the dependency value, which in file cases is the filename
    # Execute the corresponding function if we do not have the value yet
    # Once the value has been calculated save it as an internal variabel
    def get_value(self):
        if not self.filename:
            return None
        if self.exists():
            return self.filename
        parsed_args = self.parse_args()
        sys.stdout.write('Generating "' + self.filename + '" file\n')
        self.func(**parsed_args)
        return self.filename

    value = property(get_value, None, None, "The dependency value")

# SETUP ---------------------------------------------------------------------------

# Set output filenames
# WARNING: Output filenames can not be changed easily since the loader uses these names
# WARNING: Thus, these are also the filenames in the database, API, and client
OUTPUT_pdb_filename = 'md.imaged.rot.dry.pdb'
OUTPUT_trajectory_filename = 'md.imaged.rot.xtc'
OUTPUT_charges_filename = get_output_charges_filename()
OUTPUT_first_frame_filename = 'firstFrame.pdb'
OUTPUT_average_structure_filename = 'average.pdb'
OUTPUT_average_frame_filename = 'average.xtc'
OUTPUT_metadata_filename = 'metadata.json'
OUTPUT_topology_filename = 'topology.json'
OUTPUT_interactions_filename = 'md.interactions.json'
OUTPUT_rmsds_filename = 'md.rmsds.json'
OUTPUT_tmscores_filename = 'md.tmscores.json'
OUTPUT_rmsf_filename = 'md.rmsf.xvg'
OUTPUT_rgyr_filename = 'md.rgyr.xvg'
OUTPUT_pca_filename = 'md.pca.json'
OUTPUT_pca_trajectory_projection_prefix = 'pca.trajectory'
OUTPUT_rmsdperres_filename = 'md.rmsd.perres.json'
OUTPUT_rmsdpairwise_filename = 'md.rmsd.pairwise.json'
OUTPUT_distperres_filename = 'md.dist.perres.json'
OUTPUT_hbonds_filename = 'md.hbonds.json'
OUTPUT_sasa_filename = 'md.sasa.json'
OUTPUT_energies_filename = 'md.energies.json'
OUTPUT_pockets_filename = 'md.pockets.json'

# Define all dependencies
# Dependencies are tools and files that are required by some analyses
# They are only run/generated when some analysis requires them

# First of all define the inputs dict
# This dict will be updated as soon as the workflow starts with the input filenames
inputs = {}
missing_input_exception = Exception('Missing input')
# Set a function to retrieve 'inputs' values and handle missing keys
def get_input (name : str):
    value = inputs.get(name, missing_input_exception)
    if value == missing_input_exception:
        raise SystemExit('ERROR: Missing input "' + name + '"')
    return value
# Set a function to load the inputs file
def load_inputs (inputs_filename : str):
    # Check if the inputs file exists
    # If it does not, then try to download it
    if not os.path.exists(inputs_filename):
        # Check we have the inputs required to download the file
        server_url = get_input('url')
        project = get_input('project')
        # If not, then there is nothing to do
        if not server_url or not project:
            raise SystemExit('ERROR: Missing inputs file "' + inputs_filename + '"')
        # Download the inputs json file if it does not exists
        sys.stdout.write('Downloading inputs (' + inputs_filename + ')\n')
        project_url = server_url + '/api/rest/current/projects/' + project
        inputs_url = project_url + '/inputs/'
        urllib.request.urlretrieve(inputs_url, inputs_filename)
        # Write the inputs file in a pretty formatted way
        with open(inputs_filename, 'r+') as file:
            file_content = json.load(file)
            file.seek(0)
            json.dump(file_content, file, indent=4)
            file.truncate()
    # Update the inputs object with the contents of the inputs json file
    with open(inputs_filename, 'r') as file:
        inputs_data = json.load(file)
        inputs.update(inputs_data)

# Set functions to get input files names and, in case they are missing, try to download them
# DANI: Que la descarga de cada input file sea independiente debería ser útil para ahorrar descargas innecesarias
# DANI: Sin embargo es poco probable que ahorremos nada ahora mismo, porque el 'process_input_files' lo necesita todo
# DANI: Hay que darle una vuelta -> Usar el mdtb convert tal vez sería una solución

# Get the input pdb filename from the inputs
# If the file is not found try to download it
def get_input_pdb_filename ():
    # Get the pdb filename from the inputs
    original_pdb_filename = get_input('input_topology_filename')
    # Check if the file exists
    if os.path.exists(original_pdb_filename):
        return original_pdb_filename
    # If not, try to download it
    # Check we have the inputs required to download the file
    server_url = get_input('url')
    project = get_input('project')
    # If not, then there is nothing to do
    if not server_url or not project:
        raise SystemExit('ERROR: Missing input pdb file "' + original_pdb_filename + '"')
    # Download the file
    project_url = server_url + '/api/rest/current/projects/' + project
    sys.stdout.write('Downloading structure (' + original_pdb_filename + ')\n')
    topology_url = project_url + '/files/' + original_pdb_filename
    try:
        urllib.request.urlretrieve(topology_url, original_pdb_filename)
    except urllib.error.HTTPError as error:
        if error.code == 404:
            raise SystemExit('ERROR: Missing input pdb file "' + original_pdb_filename + '"')
        else:
            raise ValueError('Something went wrong with the MDposit request: ' + topology_url)
    return original_pdb_filename

# Get the input trajectory filename(s) from the inputs
# If file(s) are not found try to download it
def get_input_trajectory_filenames (sample : bool = False):
    # Get the trajectory filename(s) from the inputs
    original_trajectory_filenames = get_input('input_trajectory_filenames')
    # Check if the file exists
    if all([ os.path.exists(filename) for filename in original_trajectory_filenames ]):
        return original_trajectory_filenames
    # If not, try to download them
    # Check we have the inputs required to download the files
    server_url = get_input('url')
    project = get_input('project')
    # If not, then there is nothing to do
    if not server_url or not project:
        raise SystemExit('ERROR: Missing input trajectory files "' + ', '.join(original_trajectory_filenames) + '"')
    # Download each trajectory file (ususally it will be just one)
    project_url = server_url + '/api/rest/current/projects/' + project
    # In case only a sample is requested we use the trajectory endpoint to get the first 10 frames only
    if sample:
        original_trajectory_filename = original_trajectory_filenames[0]
        sys.stdout.write('Downloading trajectory (' + original_trajectory_filename + ')\n')
        trajectory_url = project_url + '/files/trajectory?format=xtc&frames=1:10:1'
        try:
            urllib.request.urlretrieve(trajectory_url, original_trajectory_filename)
        except urllib.error.HTTPError as error:
            if error.code == 404:
                raise SystemExit('ERROR: Missing input trajectory file "' + original_trajectory_filename + '"')
            else:
                raise ValueError('Something went wrong with the MDposit request: ' + trajectory_url)
        return [ original_trajectory_filename ]
    # Otherwise, we download the requested
    for original_trajectory_filename in original_trajectory_filenames:
        # Skip already existing files
        if os.path.exists(original_trajectory_filename):
            continue
        sys.stdout.write('Downloading trajectory (' + original_trajectory_filename + ')\n')
        trajectory_url = project_url + '/files/' + original_trajectory_filename
        try:
            urllib.request.urlretrieve(trajectory_url, original_trajectory_filename)
        except urllib.error.HTTPError as error:
            if error.code == 404:
                raise SystemExit('ERROR: Missing input trajectory file "' + original_trajectory_filename + '"')
            else:
                raise ValueError('Something went wrong with the MDposit request: ' + trajectory_url)
    return original_trajectory_filenames

# Get the input charges filename from the inputs
# If the file is not found try to download it
def get_input_charges_filename ():
    # Get the charges filename from the inputs
    original_charges_filename = get_input('input_charges_filename')
    # Check if the file exists
    if original_charges_filename and os.path.exists(original_charges_filename):
        return original_charges_filename
    # If not, try to download it
    # Check we have the inputs required to download the file
    server_url = get_input('url')
    project = get_input('project')
    # If not, then there is nothing to do
    if not server_url or not project:
        raise SystemExit('ERROR: Missing input charges file "' + original_charges_filename + '"')
    # Check which files are available for this project
    project_url = server_url + '/api/rest/current/projects/' + project
    files_url = project_url + '/files'
    response = urllib.request.urlopen(files_url)
    files = json.loads(response.read())
    # In case there is no specified charges we must find out which is the charges filename
    # It may be a topology with charges or it may be a raw charges text file
    if not original_charges_filename:
        # Check if there is any topology file (e.g. topology.prmtop, topology.top, ...)
        charges_file = next((f for f in files if f['filename'][0:9] == 'topology.'), None)
        # In case there is no topology, check if there is raw charges (i.e. 'charges.txt')
        if not charges_file:
            charges_file = next((f for f in files if f['filename'] == 'charges.txt'), None)
        # In case there is neither raw charges we send a warning and stop here
        if not charges_file:
            print('WARNING: There are no charges in this project')
        else:
            # In case we found a charges file set the input charges filename as this
            original_charges_filename = charges_file['filename']
    # Download the charges file
    if original_charges_filename and not os.path.exists(original_charges_filename):
        sys.stdout.write('Downloading charges (' + original_charges_filename + ')\n')
        charges_url = project_url + '/files/' + original_charges_filename
        try:
            urllib.request.urlretrieve(charges_url, original_charges_filename)
        except urllib.error.HTTPError as error:
            if error.code == 404:
                raise SystemExit('ERROR: Missing input topology file "' + original_charges_filename + '"')
            else:
                raise ValueError('Something went wrong with the MDposit request: ' + charges_url)
        # In the special case that the topology is from Gromacs (i.e. 'topology.top')...
        # We must also download all itp files, if any
        if original_charges_filename == 'topology.top':
            itp_filenames = [f['filename'] for f in files if f['filename'][-4:] == '.itp']
            for itp_filename in itp_filenames:
                sys.stdout.write('Downloading itp file (' + itp_filename + ')\n')
                itp_url = project_url + '/files/' + itp_filename
                urllib.request.urlretrieve(itp_url, itp_filename)
        # Set the final charges filename as the downloaded charges filename
        #charges_filename.filename = original_charges_filename
    return original_charges_filename

# Get some input values which are passed through command line instead of the inputs file
preprocess_protocol = Dependency(get_input, {'name': 'preprocess_protocol'})
translation = Dependency(get_input, {'name': 'translation'})
filter_selection = Dependency(get_input, {'name': 'filter_selection'})
pca_fit_selection = Dependency(get_input, {'name': 'pca_fit_selection'})
pca_selection = Dependency(get_input, {'name': 'pca_selection'})
skip_checkings = Dependency(get_input, {'name': 'skip_checkings'})
sample_trajectory = Dependency(get_input, {'name': 'sample_trajectory'})

# Get input filenames from the already updated inputs

# The original topology and trajectory filenames are the input filenames
original_topology_filename = Dependency(get_input_pdb_filename, {})
original_trajectory_filenames = Dependency(get_input_trajectory_filenames, {'sample': sample_trajectory})
original_charges_filename = Dependency(get_input_charges_filename, {})
# The 'inputs file' is a json file which must contain all parameters to run the workflow
inputs_filename = Dependency(get_input, {'name': 'inputs_filename'})

# Extract some additional input values from the inputs json file
input_interactions = Dependency(get_input, {'name': 'interactions'})
ligands = Dependency(get_input, {'name': 'ligands'})
membranes = Dependency(get_input, {'name': 'membranes'})
forced_references = Dependency(get_input, {'name': 'forced_references'})

# Define intermediate tools and files

# Process the topology and or trajectory files using VMD
# Files are converted to supported formats and trajectory pieces are merged into a single file
# In addition, some irregularities in the topology may be fixed by VMD
# If the output topology and trajectory files already exists it is assumed they are already processed
vmd = Dependency(vmd_processor, {
    'input_topology_filename': original_topology_filename,
    'input_trajectory_filenames': original_trajectory_filenames,
    'output_topology_filename': OUTPUT_pdb_filename,
    'output_trajectory_filename': OUTPUT_trajectory_filename,
}, 'vmd')

# Preprocessing
# This is equivalent to running 'vmd', 'imaging' and 'corrector' together
process_input_files = Dependency(process_input_files, {
    'input_topology_filename': original_topology_filename,
    'input_trajectory_filenames': original_trajectory_filenames,
    'input_charges_filename': original_charges_filename,
    'output_topology_filename': OUTPUT_pdb_filename,
    'output_trajectory_filename': OUTPUT_trajectory_filename,
    'preprocess_protocol': preprocess_protocol,
    'translation': translation,
    'filter_selection' : filter_selection,
})

# Main topology and trajectory files
# These may be created from the original topology and trajectory files after some processing
# These files are then widely used along the workflow
pdb_filename = File(OUTPUT_pdb_filename,
    process_input_files.func,
    process_input_files.args)
trajectory_filename = File(OUTPUT_trajectory_filename,
    process_input_files.func,
    process_input_files.args)

# Charges filename
# This may be created from the original charges filename
# This file is used for the energies analysis
charges_filename = File(OUTPUT_charges_filename,
    process_input_files.func,
    process_input_files.args)

# Filter atoms to remove water and ions
# As an exception, some water and ions may be not removed if specified
# WARNING: This is the independent call for this function
# WARNING: In the canonical workflow, this function is called inside 'process_input_files'
filtering = Dependency(filter_atoms, {
    'topology_filename' : original_topology_filename,
    'trajectory_filename' : original_trajectory_filenames,
    'charges_filename' : original_charges_filename,
    'filter_selection': filter_selection,
}, 'filter')

# Image the trajectory if it is required
# i.e. make the trajectory uniform avoiding atom jumps and making molecules to stay whole
# Fit the trajectory by removing the translation and rotation if it is required
# WARNING: This is the independent call for this function
# WARNING: In the canonical workflow, this function is called inside 'process_input_files'
imaging = Dependency(image_and_fit, {
    'input_topology_filename': original_topology_filename,
    'input_trajectory_filename': original_trajectory_filenames,
    'input_tpr_filename' : original_charges_filename,
    'output_topology_filename': original_topology_filename,
    'output_trajectory_filename': original_trajectory_filenames,
    'preprocess_protocol': preprocess_protocol,
    'translation': translation,
}, 'imaging')

# Examine and correct the topology file
# WARNING: This is the independent call for this function
# WARNING: In the canonical workflow, this function is called inside 'process_input_files'
corrector = Dependency(topology_corrector, {
    'input_pdb_filename': original_topology_filename,
    'output_topology_filename': original_topology_filename,
    'input_trajectory_filename': original_trajectory_filenames,
    'output_trajectory_filename': original_trajectory_filenames,
}, 'corrector')

# Set a parsed structure/topology with useful features
# This object also include functions to convert residue numeration from one format to another
structure = Dependency(setup_structure, {
    'pdb_filename': pdb_filename
})

# Set the pytraj trayectory
pt_trajectory = Dependency(get_pytraj_trajectory, {
    'input_topology_filename': pdb_filename,
    'input_trajectory_filename': trajectory_filename
})

# Get the first trajectory frame
first_frame_filename = File(OUTPUT_first_frame_filename, get_first_frame, {
    'input_topology_filename' : pdb_filename,
    'input_trajectory_filename' : trajectory_filename,
    'first_frame_filename' : OUTPUT_first_frame_filename
})
# Get the average structure in pdb format
average_structure_filename = File(OUTPUT_average_structure_filename, get_average, {
    'pytraj_trajectory': pt_trajectory,
    'output_average_filename': OUTPUT_average_structure_filename
})
# Get the average structure in frame format
# It is further loaded to database and used to represent pockets
# DANI: Ahora mismo no se usan, pero algún día deberían volver a usarse
average_frame_filename = File(OUTPUT_average_frame_filename, get_average, {
    'pytraj_trajectory': pt_trajectory,
    'output_average_filename': OUTPUT_average_frame_filename
})

# Get additional metadata usedin the workflow which is not in the inputs

# Count the number of snapshots
snapshots = Dependency(get_frames_count, {
    'input_topology_filename': pdb_filename,
    'input_trajectory_filename': trajectory_filename
}, 'snapshots')

# Find out residues and interface residues for each interaction
interactions = Dependency(process_interactions, {
    'topology_filename': pdb_filename,
    'trajectory_filename': trajectory_filename,
    'interactions': input_interactions,
    'structure': structure,
    'interactions_file': OUTPUT_interactions_filename,
}, 'interactions')

# Find out residues and interface residues for each interaction
charges = Dependency(get_charges, {
    'charges_source_filename': charges_filename,
}, 'charges')

# Prepare the metadata output file
residues_map = Dependency(generate_map_online, {
    'structure': structure,
    'forced_references': forced_references
}, 'map')

# Prepare the metadata output file
metadata_filename = File(OUTPUT_metadata_filename, generate_metadata, {
    'input_topology_filename': pdb_filename,
    'input_trajectory_filename': trajectory_filename,
    'inputs_filename': inputs_filename,
    'snapshots': snapshots,
    'residues_map': residues_map,
    'output_metadata_filename': OUTPUT_metadata_filename
}, 'metadata')

# Prepare the topology output file
topology_filename = File(OUTPUT_topology_filename, generate_topology, {
    'structure': structure,
    'charges': charges,
    'residues_map': residues_map,
    'output_topology_filename': OUTPUT_topology_filename
}, 'topology')

# Pack up all tools which may be called directly from the console
tools = [
    vmd,
    filtering,
    imaging,
    corrector,
    interactions,
    snapshots,
    charges,
]

# Define all analyses
analyses = [
    # WARNING: This analysis is fast enought to use the full trajectory
    # WARNING: However, the output file size depends on the trajectory size
    # WARNING: In very long trajectories the number of points may make the client go slow when loading data
    File(OUTPUT_rmsds_filename, rmsds, {
        "input_trajectory_filename": trajectory_filename,
        "output_analysis_filename": OUTPUT_rmsds_filename,
        'frames_limit': 5000,
        'first_frame_filename': first_frame_filename,
        'average_structure_filename': average_structure_filename,
        'structure': structure,
        'skip_checkings': skip_checkings,
    }, 'rmsds'),
    # Here we set a small frames limit since this anlaysis is a bit slow
    File(OUTPUT_tmscores_filename, tmscores, {
        "input_topology_filename": pdb_filename,
        "input_trajectory_filename": trajectory_filename,
        "output_analysis_filename": OUTPUT_tmscores_filename,
        'first_frame_filename': first_frame_filename,
        'average_structure_filename': average_structure_filename,
        'structure' : structure,
        'frames_limit': 200,
    }, 'tmscores'),
    # This analysis is fast and the output size depends on the number of atoms only
    # For this reason here it is used the whole trajectory with no frames limit
    File(OUTPUT_rmsf_filename, rmsf, {
        "input_topology_filename": pdb_filename,
        "input_trajectory_filename": trajectory_filename,
        "output_analysis_filename": OUTPUT_rmsf_filename
    }, 'rmsf'),
    # WARNING: This analysis is fast enought to use the full trajectory instead of the reduced one
    # WARNING: However, the output file size depends on the trajectory size
    # WARNING: In very long trajectories the number of points may make the client go slow when loading data
    File(OUTPUT_rgyr_filename, rgyr, {
        "input_topology_filename": pdb_filename,
        "input_trajectory_filename": trajectory_filename,
        "output_analysis_filename": OUTPUT_rgyr_filename,
        'frames_limit': 5000,
    }, 'rgyr'),
    # WARNING: This analysis will generate several output files
    # File 'pca.average.pdb' is generated by the PCA and used by the client
    # File 'covar.log' is generated by the PCA but never used
    File(OUTPUT_pca_filename, pca, {
        "input_topology_filename": pdb_filename,
        "input_trajectory_filename": trajectory_filename,
        "output_analysis_filename": OUTPUT_pca_filename,
        "output_trajectory_projections_prefix": OUTPUT_pca_trajectory_projection_prefix,
        'frames_limit': 2000,
        'structure': structure,
        'fit_selection': pca_fit_selection,
        'analysis_selection': pca_selection,
    }, 'pca'),
    # DANI: Intenta usar mucha memoria, hay que revisar
    # DANI: Puede saltar un error de imposible alojar tanta memoria
    # DANI: Puede comerse toda la ram y que al final salte un error de 'Terminado (killed)'
    # DANI: De momento me lo salto
    # Dependency(pca_contacts, {
    #     "trajectory": trajectory_filename,
    #     "topology": pdb_filename,
    #     "interactions": interactions,
    #     "output_analysis_filename": "md.pca.contacts.json"
    # }, 'pcacons'),
    # WARNING: This analysis is fast enought to use the full trajectory instead of the reduced one
    # WARNING: However, the output file size depends on the trajectory size. It may be pretty big
    File(OUTPUT_rmsdperres_filename, rmsd_per_residue, {
        'input_topology_filename': pdb_filename,
        'input_trajectory_filename': trajectory_filename,
        "output_analysis_filename": OUTPUT_rmsdperres_filename,
        'structure': structure,
        'membranes': membranes,
        'frames_limit': 100,
    }, 'rmsdperres'),
    # WARNING: This analysis is fast enought to use the full trajectory instead of the reduced one
    # WARNING: However, the output file size depends on the trajectory size exponentially. It may be huge
    File(OUTPUT_rmsdpairwise_filename, rmsd_pairwise, {
        'input_topology_filename': pdb_filename,
        'input_trajectory_filename': trajectory_filename,
        "output_analysis_filename": OUTPUT_rmsdpairwise_filename,
        "interactions": interactions,
        'frames_limit': 200,
    }, 'rmsdpairwise'),
    # WARNING: This analysis is not fast enought to use the full trajectory. It would take a while
    File(OUTPUT_distperres_filename, distance_per_residue, {
        'input_topology_filename': pdb_filename,
        'input_trajectory_filename': trajectory_filename,
        "output_analysis_filename": OUTPUT_distperres_filename,
        "interactions": interactions,
        'frames_limit': 200,
    }, 'distperres'),
    # WARNING: This analysis is fast enought to use the full trajectory instead of the reduced one
    # WARNING: However, the output file size depends on the trajectory
    # WARNING: Files have no limit, but analyses must be no heavier than 16Mb in BSON format
    # WARNING: In case of large surface interaction the output analysis may be larger than the limit
    # DANI: Esto no puede quedar así
    # DANI: Me sabe muy mal perder resolución con este análisis, porque en cáculo es muy rápido
    # DANI: Hay que crear un sistema de carga en mongo alternativo para análisis pesados
    File(OUTPUT_hbonds_filename, hydrogen_bonds, {
        'input_topology_filename': pdb_filename,
        'input_trajectory_filename': trajectory_filename,
        "output_analysis_filename": OUTPUT_hbonds_filename,
        'structure': structure,
        "interactions": interactions,
        'frames_limit': 200,
    }, 'hbonds'),
    File(OUTPUT_sasa_filename, sasa, {
        "input_topology_filename": pdb_filename,
        "input_trajectory_filename": trajectory_filename,
        "output_analysis_filename": OUTPUT_sasa_filename,
        'structure': structure,
        'membranes': membranes,
        'frames_limit': 100,
    }, 'sasa'),
    File(OUTPUT_energies_filename, energies, {
        "input_topology_filename": pdb_filename,
        "input_trajectory_filename": trajectory_filename,
        "output_analysis_filename": OUTPUT_energies_filename,
        'structure': structure,
        "interactions": interactions,
        'charges': charges,
        'frames_limit': 100,
    }, 'energies'),
    File(OUTPUT_pockets_filename, pockets, {
        "input_topology_filename": pdb_filename,
        "input_trajectory_filename": trajectory_filename,
        "output_analysis_filename": OUTPUT_pockets_filename,
        'structure': structure,
        'membranes': membranes,
        'frames_limit': 100,
    }, 'pockets'),
]

# Set a list with all dependencies to be required if the whole workflow is run
workflow = [ metadata_filename, topology_filename, *analyses ]

# Set a list with all dependencies which may be requested independently
requestables = [ *analyses, *tools, metadata_filename, topology_filename ]

# Main ---------------------------------------------------------------------------------            

# Main function
def main():

    # Parse input arguments from the console
    args = parser.parse_args()
    input_topology_filename = args.input_topology_filename
    input_trajectory_filenames = args.input_trajectory_filenames
    inputs_filename =  args.inputs_filename
    input_charges_filename = args.input_charges_filename

    # Set the command line inputs
    # The vars function converts the args object to a dictionary
    inputs.update(vars(args))

    # Load the inputs file
    load_inputs(inputs_filename)

    # If download is passed as True then exit as soon as the setup is finished
    if args.download:
        original_topology_filename.value
        original_trajectory_filenames.value
        original_charges_filename.value
        return

    # Run tools which must be run always
    # They better be fast
    process_input_files.value

    # If setup is passed as True then exit as soon as the setup is finished
    if args.setup:
        return

    # Run the requested analyses
    if args.include and len(args.include) > 0:
        sys.stdout.write(f"Executing specific dependencies: " + str(args.include) + '\n')
        # Include only the specified dependencies
        requested_dependencies = [ dep for dep in requestables if dep.alias in args.include ]
        for dependency in requested_dependencies:
            dependency.value
    # Run all analyses if none was specified
    else:
        # Include all non excluded dependencies
        requested_dependencies = workflow
        if args.exclude and len(args.exclude) > 0:
            requested_dependencies = [ dep for dep in workflow if dep.alias not in args.exclude ]
        for dependency in requested_dependencies:
            dependency.value

    # Remove gromacs backups
    remove_trash()

    sys.stdout.write("Done!\n")


# Define console arguments to call the workflow
parser = ArgumentParser(description="MoDEL Workflow", formatter_class=RawTextHelpFormatter)

# Get current directory
current_directory = Path.cwd()

# Set default input filenames
# They may be modified through console command arguments
DEFAULT_input_topology_filename = 'md.imaged.rot.dry.pdb'
DEFAULT_input_trajectory_filenames = ['md.imaged.rot.xtc']
DEFAULT_inputs_filename = 'inputs.json'
DEFAULT_input_charges_filename = find_charges_filename()

# Set optional arguments
parser.add_argument(
    "-dir", "--working_dir",
    default=current_directory,
    help="Directory where to perform analysis. "
    "If empty, will use current directory.")

parser.add_argument(
    "-proj", "--project",
    default=None,
    help="If given a project name, trajectory and "
    "topology files will be downloaded from remote server.")

parser.add_argument(
    "-url",
    default="https://mdposit-dev.bsc.es",
    help="URL from where to download project")

parser.add_argument(
    "-top", "--input_topology_filename",
    default=DEFAULT_input_topology_filename,
    help="Path to topology filename")

parser.add_argument(
    "-traj", "--input_trajectory_filenames",
    #type=argparse.FileType('r'),
    nargs='*',
    default=DEFAULT_input_trajectory_filenames,
    help="Path to trajectory filename")

parser.add_argument(
    "-in", "--inputs_filename",
    default=DEFAULT_inputs_filename,
    help="Path to inputs filename")

parser.add_argument(
    "-char", "--input_charges_filename",
    default=DEFAULT_input_charges_filename, # There is no default since many formats may be possible
    help="Path to charges topology filename")

parser.add_argument(
    "-pr", "--preprocess_protocol",
    type=int,
    default=0,
    help=("Set how the trajectory must be imaged and fitted (i.e. centered, without translation or rotation)\n"
        "These protocolos may help in some situations, but the imaging step can not be fully automatized\n"
        "If protocols do not work, the gromacs parameters must be modified manually\n"
        "Available protocols:\n"
        "0. Do nothing (default) -> The trajectory is already imaged and fitted\n"
        "1. No imaging, only fitting -> The trajectory is already imaged but not fitted\n"
        "2. Basic imaging -> Atoms are centered automatically\n"
        "   Recommended for single molecules only\n"
        "3. Translated imaging -> Manually translate everything before imaging\n"
        "   Recommended for interacting molecules\n"
        "4. Allowed jump imaging -> Residues are centered automatically. The 'nojump' and the fitting steps are skipped\n"
        "   Recommended for proteins inside membranes\n"
        "   * Note that a .tpr topology is required in order to run protocol 4"))

parser.add_argument(
    "-trans", "--translation",
    nargs='*',
    default=[0,0,0],
    help=("Set the x y z translation for the imaging process (only protocol 3)\n"
        "e.g. -trans 0.5 -1 0"))

parser.add_argument(
    "-filt", "--filter_selection",
    nargs='?',
    default=False,
    const=True,
    help=("Atoms selection to be filtered in VMD format\n"
        "If the argument is passed alone (i.e. with no selection) then water and counter ions are filtered"))

parser.add_argument(
    "-pcafit", "--pca_fit_selection",
    default=protein_and_nucleic_backbone,
    help="Atom selection for the pca fitting in vmd syntax")

parser.add_argument(
    "-pcasel", "--pca_selection",
    default=protein_and_nucleic_backbone,
    help="Atom selection for pca analysis in vmd syntax")

parser.add_argument(
    "-d", "--download",
    action='store_true',
    help="If passed, only download required files. Then exits.")

parser.add_argument(
    "-s", "--setup",
    action='store_true',
    help="If passed, only download required files and run mandatory dependencies. Then exits.")

parser.add_argument(
    "-sck", "--skip_checkings",
    action='store_true',
    help="If passed, skip RMSD checking and proceed with the workflow")

parser.add_argument(
    "-smp", "--sample_trajectory",
    action='store_true',
    help="If passed, download just the 10 first frames of the trajectory instead of it all")

# Set a list with the alias of all requestable dependencies
choices = [ dependency.alias for dependency in requestables ]

parser.add_argument(
    "-i", "--include",
    nargs='*',
    choices=choices,
    help="Set the unique analyses or tools to be run. All other steps will be skipped")

parser.add_argument(
    "-e", "--exclude",
    nargs='*',
    choices=choices,
    help=("Set the unique analyses or tools to be skipped. All other steps will be run.\n"
        "If the 'include' argument is passed the 'exclude' argument will be ignored.\n"
        "WARNING: If an excluded dependecy is required by others then it will be run anyway"))
