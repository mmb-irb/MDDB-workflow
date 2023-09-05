#!/usr/bin/env python

# This is the starter script

# Import python libraries
from argparse import ArgumentParser, RawTextHelpFormatter, Action
import os
from os.path import exists
import sys
import io
import math
from pathlib import Path
import urllib.request
import json
from typing import Optional, Union, List

# Import local tools
from model_workflow.tools.topology_corrector import topology_corrector
from model_workflow.tools.process_input_files import process_input_files, find_charges_filename
from model_workflow.tools.topology_manager import setup_structure
from model_workflow.tools.get_pytraj_trajectory import get_pytraj_trajectory
from model_workflow.tools.get_first_frame import get_first_frame
from model_workflow.tools.get_average import get_average
from model_workflow.tools.process_interactions import process_interactions
from model_workflow.tools.get_pbc_residues import get_pbc_residues
from model_workflow.tools.generate_metadata import generate_metadata
from model_workflow.tools.generate_map import generate_map_online
from model_workflow.tools.generate_topology import generate_topology
from model_workflow.tools.get_summarized_trajectory import get_summarized_trajectory
from model_workflow.tools.get_frames_count import get_frames_count
from model_workflow.tools.get_charges import get_charges
from model_workflow.tools.remove_trash import remove_trash
from model_workflow.tools.register import start_register, save_register
from model_workflow.tools.get_screenshot import get_screenshot


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
#from model_workflow.analyses.hydrogen_bonds_2 import hydrogen_bonds
from model_workflow.analyses.hydrogen_bonds import hydrogen_bonds
from model_workflow.analyses.sasa import sasa
from model_workflow.analyses.energies import energies
from model_workflow.analyses.pockets import pockets
from model_workflow.analyses.rmsd_check import check_sudden_jumps
from model_workflow.analyses.helical_parameters import helical_parameters
from model_workflow.analyses.markov import markov

# Make the system output stream to not be buffered
# This is useful to make prints work on time in Slurm
# Otherwise, output logs are written after the script has fully run
# Note that this fix affects all modules and built-ins
unbuffered = io.TextIOWrapper(open(sys.stdout.fileno(), 'wb', 0), write_through=True)
sys.stdout = unbuffered

# Set a standard selection for protein and nucleic acid backbones in vmd syntax
protein_and_nucleic = 'protein or nucleic'
protein_and_nucleic_backbone = "(protein and name N CA C) or (nucleic and name P O5' O3' C5' C4' C3')"

# CLASSES -------------------------------------------------------------------------

dependency_unname = 'unnamed dependency'
class Dependency:
    _depend = True
    # The 'func' is the function to get the value
    # The 'args' are the arguments for the 'func'
    # The 'filename' is the name of the file which is expected to be required for this dependency
    # When a filename is passed the dependency 'value' is replaced by this 'filename' when exists
    def __init__ (self, func, args = {}, alias = dependency_unname):
        self.func = func
        self.args = args
        self.alias = alias
        self._value = None
        # Set an internal flag to know if a a dependency has been calculated
        # Checking if _value is None is not enough, since the result of a calculation may be None
        self._done = False

    def __str__ (self):
        return self.alias

    def __repr__ (self):
        return self.alias

    # For each 'Dependency' instance in args call its dependency function
    def parse_args (self):

        # For a given non Dependency instance, return the value as it is
        # For a given Dependency instance, get the real value of it
        # Getting the value means also calling its dependency function
        def parse (value):
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
    def get_value (self):
        value = self._value
        if self._done:
            return value
        parsed_args = self.parse_args()
        if self.alias and self.alias != dependency_unname:
            sys.stdout.write('Running "' + self.alias + '"\n')
        value = self.func(**parsed_args)
        self._value = value # Save the result for possible further use
        self._done = True
        return value

    value = property(get_value, None, None, "The dependency value")

    # Reset the dependency for it to be considered as not run
    def reset (self):
        self._value = None
        self._done = False

class File(Dependency):
    # The 'func' is the function to generate the file
    # If there is no func it means the file will be always there (i.e. it is an input file)
    # The 'args' are the arguments for the 'func'
    # The 'must_remake' flag is used to generate the file and thus overwrite the previous file once per workflow run
    def __init__ (self, filename : str, func, args = {}, alias = '', must_remake : bool = False):
        self.filename = filename
        self.must_remake = must_remake
        super().__init__(func, args, alias)

    def __str__ (self):
        return self.filename

    def __repr__ (self):
        return self.filename

    def exists (self) -> bool:
        if (self.must_remake and not self._done) or not self.filename:
            return False
        return exists(self.filename)

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
        self._done = True
        return self.filename

    value = property(get_value, None, None, "The dependency value")

# SETUP ---------------------------------------------------------------------------

# Set output filenames
# WARNING: Output filenames can not be changed easily since the loader uses these names
# WARNING: Thus, these are also the filenames in the database, API, and client
OUTPUT_provisional_pdb_filename = 'provisional.pdb'
OUTPUT_pdb_filename = 'md.imaged.rot.dry.pdb'
OUTPUT_trajectory_filename = 'md.imaged.rot.xtc'
#OUTPUT_pdb_filename = 'structure.pdb'
#OUTPUT_trajectory_filename = 'trajectory.xtc'
OUTPUT_first_frame_filename = 'firstFrame.pdb'
OUTPUT_average_structure_filename = 'average.pdb'
OUTPUT_average_frame_filename = 'average.xtc'
OUTPUT_metadata_filename = 'metadata.json'
OUTPUT_topology_filename = 'topology.json'
OUTPUT_interactions_filename = 'md.interactions.json'
OUTPUT_rmsds_filename = 'md.rmsds.json'
OUTPUT_tmscores_filename = 'md.tmscores.json'
OUTPUT_rmsf_filename = 'md.rmsf.json'
OUTPUT_rgyr_filename = 'md.rgyr.json'
OUTPUT_pca_filename = 'md.pca.json'
OUTPUT_pca_trajectory_projection_prefix = 'pca.trajectory'
OUTPUT_rmsdperres_filename = 'md.rmsd.perres.json'
OUTPUT_rmsdpairwise_filename = 'md.rmsd.pairwise.json'
OUTPUT_distperres_filename = 'md.dist.perres.json'
OUTPUT_hbonds_filename = 'md.hbonds.json'
OUTPUT_sasa_filename = 'md.sasa.json'
OUTPUT_energies_filename = 'md.energies.json'
OUTPUT_pockets_filename = 'md.pockets.json'
OUTPUT_helical_parameters_filename = 'md.helical.parameters.json'
OUTPUT_screenshot_filename = 'screenshot.jpg'
OUTPUT_markov_filename = 'md.markov.json'

# State all the available checkings, which may be trusted
available_checkings = [ 'stabonds', 'cohbonds', 'intrajrity' ]

# State all critical process failures, which are to be lethal for the workflow unless mercy is given
available_failures = available_checkings + [ 'refseq', 'interact' ]

# Define all dependencies
# Dependencies are tools and files that are required by some analyses
# They are only run/generated when some analysis requires them

# First of all define the inputs dict
# This dict will be updated as soon as the workflow starts with the input filenames
inputs = {}
missing_input_exception = Exception('Missing input')
# Set a function to retrieve 'inputs' values and handle missing keys
def get_input (name : str, missing_input_callback = missing_input_exception):
    value = inputs.get(name, missing_input_callback)
    if value == missing_input_exception:
        raise SystemExit('ERROR: Missing input "' + name + '"')
    return value
# Set a function to load the inputs file
def load_inputs (inputs_filename : str):
    # Check if the inputs file exists
    # If it does not, then try to download it
    if not exists(inputs_filename):
        # Check we have the inputs required to download the file
        server_url = get_input('database_url')
        project = get_input('project')
        # If not, then there is nothing to do
        if not server_url or not project:
            #raise SystemExit('ERROR: Missing inputs file "' + inputs_filename + '"')
            return
        # Download the inputs json file if it does not exists
        sys.stdout.write('Downloading inputs (' + inputs_filename + ')\n')
        project_url = server_url + '/rest/current/projects/' + project
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
def get_input_pdb_filename () -> str:
    # Get the pdb filename from the inputs
    original_pdb_filename = get_input('input_topology_filename')
    # Check if the file exists
    if exists(original_pdb_filename):
        return original_pdb_filename
    # If not, try to download it
    # Check we have the inputs required to download the file
    server_url = get_input('database_url')
    project = get_input('project')
    # If not, then there is nothing to do
    if not server_url or not project:
        raise SystemExit('ERROR: Missing input pdb file "' + original_pdb_filename + '"')
    # Download the file
    project_url = server_url + '/rest/current/projects/' + project
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
def get_input_trajectory_filenames (sample : bool = False) -> str:
    # Get the trajectory filename(s) from the inputs
    original_trajectory_filenames = get_input('input_trajectory_filenames')
    # Check if the file exists
    if all([ exists(filename) for filename in original_trajectory_filenames ]):
        return original_trajectory_filenames
    # If not, try to download them
    # Check we have the inputs required to download the files
    server_url = get_input('database_url')
    project = get_input('project')
    # If not, then there is nothing to do
    if not server_url or not project:
        raise SystemExit('ERROR: Missing input trajectory files "' + ', '.join(original_trajectory_filenames) + '"')
    # Download each trajectory file (ususally it will be just one)
    project_url = server_url + '/rest/current/projects/' + project
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
        if exists(original_trajectory_filename):
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
def get_input_charges_filename () -> str:
    # Get the charges filename from the inputs
    original_charges_filename = get_input('input_charges_filename')
    # Check if the file exists
    if original_charges_filename and exists(original_charges_filename):
        return original_charges_filename
    # If not, try to download it
    # Check we have the inputs required to download the file
    server_url = get_input('database_url')
    project = get_input('project')
    # If not, then there is nothing to do
    if original_charges_filename and (not server_url or not project):
        raise SystemExit('ERROR: Missing input charges file "' + original_charges_filename + '"')
    # Check which files are available for this project
    if project:
        project_url = server_url + '/rest/current/projects/' + project
        # Check if the project has a topology and download it in json format if so
        topology_url = project_url + '/topology'
        topology = None
        try:
            sys.stdout.write('Downloading charges (topology.json)\n')
            response = urllib.request.urlopen(topology_url)
            topology = json.loads(response.read())
        except:
            pass
        # If we have the standard topology then export it and stop here
        if topology:
            with open(OUTPUT_topology_filename, 'w') as file:
                json.dump(topology, file, indent=4)
            return OUTPUT_topology_filename
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
        if original_charges_filename and not exists(original_charges_filename):
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
    return original_charges_filename

# Get the input populations filename from the inputs
# If the file is not found and we have download inputs then try to download it
# If the file is not found and we don't have download inputs then its okay, this is an optional input
def get_populations (populations_filename : str) -> Optional[ List[float] ]:
    # Get the populations filename from the inputs
    input_populations_filename = get_input('input_populations_filename')
    # If the file already exists then read and parse it and return its content
    if exists(input_populations_filename):
        with open(input_populations_filename, 'r') as file:
            populations = json.load(file)
        return populations
    # Try to download the file
    # Check we have the inputs required to download the file
    server_url = get_input('database_url')
    project = get_input('project')
    # If we have not download inputs then there is nothing to do
    if not server_url or not project:
        return None
    # Get a list with the names of files in this project
    project_url = server_url + '/rest/current/projects/' + project
    files_url = project_url + '/files'
    files_response = None
    try:
        files_response = urllib.request.urlopen(files_url)
    except:
        raise SystemExit('Something went wrong with the MDposit request: ' + files_url)
    files = json.loads(files_response.read())
    filenames = [ f['filename'] for f in files ]
    # If the populations file is not in the list then there is nothing to do
    if not input_populations_filename in filenames:
        return None
    # Download the file
    sys.stdout.write('Downloading populations (' + input_populations_filename + ')\n')
    populations_url = files_url + '/' + input_populations_filename
    try:
        urllib.request.urlretrieve(populations_url, input_populations_filename)
    except urllib.error.HTTPError as error:
        if error.code == 404:
            raise SystemExit('ERROR: Missing input populations file "' + input_populations_filename + '"')
        else:
            raise ValueError('Something went wrong with the MDposit request: ' + populations_url)
    # Read and parse the donwloaded file and return its content
    with open(input_populations_filename, 'r') as file:
        populations = json.load(file)
    return populations

# Get some input values which are passed through command line instead of the inputs file
image = Dependency(get_input, {'name': 'image'})
fit = Dependency(get_input, {'name': 'fit'})
translation = Dependency(get_input, {'name': 'translation'})
filter_selection = Dependency(get_input, {'name': 'filter_selection'})
pca_fit_selection = Dependency(get_input, {'name': 'pca_fit_selection'})
pca_selection = Dependency(get_input, {'name': 'pca_selection'})
mercy = Dependency(get_input, {'name': 'mercy'})
trust = Dependency(get_input, {'name': 'trust'})
sample_trajectory = Dependency(get_input, {'name': 'sample_trajectory'})
input_type = Dependency(get_input, {'name': 'type'})
# DANI: Esto es temporal, creo. El input 'type' debió llamarse 'is_time_dependend/related' y ser un boolean
# DANI: Algun día lo cambiaré, pero esto implica cambio en db, api y cliente también
def check_is_time_dependent (input_type : str) -> bool:
    if input_type == 'trajectory':
        return True
    elif input_type == 'ensemble':
        return False
    raise SystemExit('Not supported input type value: ' + input_type)
is_time_dependend = Dependency(check_is_time_dependent, {'input_type': input_type})

# Get input filenames from the already updated inputs

# The original topology and trajectory filenames are the input filenames
original_topology_filename = Dependency(get_input_pdb_filename, {})
original_trajectory_filenames = Dependency(get_input_trajectory_filenames, {'sample': sample_trajectory})
original_charges_filename = Dependency(get_input_charges_filename, {})
# The 'inputs file' is a json file which must contain all parameters to run the workflow
inputs_filename = Dependency(get_input, {'name': 'inputs_filename'})
populations_filename = Dependency(get_input, {'name': 'input_populations_filename'})

# Extract some additional input values from the inputs json file
input_interactions = Dependency(get_input, {'name': 'interactions'})
input_pbc_selection = Dependency(get_input, {'name': 'pbc_selection', 'missing_input_callback': None})
forced_references = Dependency(get_input, {'name': 'forced_references'})
pdb_ids = Dependency(get_input, {'name': 'pdbIds'})
time_length = Dependency(get_input, {'name': 'length'})

# Set the register here
# Note that register is called after input have been fully defined
register = Dependency(start_register, { 'inputs': inputs })

# Define intermediate tools and files

# Preprocessing
# This is equivalent to running 'vmd', 'imaging' and 'corrector' together
process_input_files = Dependency(process_input_files, {
    'input_topology_filename': original_topology_filename,
    'input_trajectory_filenames': original_trajectory_filenames,
    'input_charges_filename': original_charges_filename,
    'output_topology_filename': OUTPUT_provisional_pdb_filename,
    'output_trajectory_filename': OUTPUT_trajectory_filename,
    'image': image,
    'fit': fit,
    'translation': translation,
    'filter_selection' : filter_selection,
    'pbc_selection' : input_pbc_selection,
})

# Main topology and trajectory files
# These may be created from the original topology and trajectory files after some processing
# These files are then widely used along the workflow
provisional_pdb_filename = File(OUTPUT_provisional_pdb_filename,
    process_input_files.func,
    process_input_files.args)
trajectory_filename = File(OUTPUT_trajectory_filename,
    process_input_files.func,
    process_input_files.args)

# Charges filename
# This may be created from the original charges filename
# Note that there is not a defined name since it may change
# Atom charges are extracted from this file
# DANI: This is the name of a file which is generated through the 'process_input_files' function
# DANI: Althought this should be a File and not a Dependency there are two reasons to set it like this
# DANI: 1 - The name of the file is variable
# DANI: 2 - This file is never requested alone
charges_filename = Dependency(find_charges_filename, {})

# Count the number of snapshots
# This value is used widely along the workflow instead of counting frames again to be more efficient
# Note that some trajectory formats require to read the whole file to count frames (e.g. xtc)
# The logic to count frames is powered by pytraj and thus it handles many different trajectory formats
snapshots = Dependency(get_frames_count, {
    'input_topology_filename': provisional_pdb_filename,
    'input_trajectory_filename': trajectory_filename
}, 'snapshots')

# Examine and correct the structure file
corrector = Dependency(topology_corrector, {
    'input_pdb_filename': OUTPUT_provisional_pdb_filename,
    'output_topology_filename': OUTPUT_pdb_filename,
    'input_trajectory_filename': trajectory_filename,
    'output_trajectory_filename': trajectory_filename,
    'input_charges_filename': charges_filename,
    'snapshots' : snapshots,
    'register' : register,
    'mercy' : mercy,
    'trust' : trust,
}, 'corrector')

# Main topology and trajectory files
# These may be created from the original topology and trajectory files after some processing
# These files are then widely used along the workflow
pdb_filename = File(OUTPUT_pdb_filename,
    corrector.func,
    corrector.args,
    must_remake=True)

# Set a parsed structure/topology with useful features
# IMPORTANT: Note that the pdb file at this point is already corrected
# This object also include functions to convert residue numeration from one format to another
structure = Dependency(setup_structure, {
    'pdb_filename': pdb_filename,
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

# Get additional metadata used in the workflow which is not in the inputs

# Get the cutoff for the interactions below
interaction_cutoff = Dependency(get_input, {'name': 'interaction_cutoff'})
# Find out residues and interface residues for each interaction
interactions = Dependency(process_interactions, {
    'input_interactions': input_interactions,
    'topology_filename': pdb_filename,
    'trajectory_filename': trajectory_filename,
    'structure': structure,
    'snapshots' : snapshots,
    'interactions_file': OUTPUT_interactions_filename,
    'mercy' : mercy,
    'frames_limit': 1000,
    'interaction_cutoff': interaction_cutoff
}, 'interactions')

# Find the PBC residues
pbc_residues = Dependency(get_pbc_residues, {
    'structure': structure,
    'input_pbc_selection': input_pbc_selection
}, 'pbc_residues')

# Find out residues and interface residues for each interaction
charges = Dependency(get_charges, {
    'charges_source_filename': charges_filename,
}, 'charges')

# Load the populations
populations = Dependency(get_populations, { 'populations_filename': populations_filename })

# Prepare the metadata output file
residues_map = Dependency(generate_map_online, {
    'structure': structure,
    'register': register,
    'mercy': mercy,
    'forced_references': forced_references,
    'pdb_ids': pdb_ids,
}, 'map')

# Prepare the metadata output file
# It is cheap to remake and this may solve problems when changing inputs if metadata is already done
metadata_filename = File(OUTPUT_metadata_filename, generate_metadata, {
    'input_topology_filename': pdb_filename,
    'input_trajectory_filename': trajectory_filename,
    'inputs_filename': inputs_filename,
    'snapshots': snapshots,
    'residues_map': residues_map,
    'interactions': interactions,
    'output_metadata_filename': OUTPUT_metadata_filename,
    'register': register
}, 'metadata', must_remake=True)

# Prepare the topology output file
topology_filename = File(OUTPUT_topology_filename, generate_topology, {
    'structure': structure,
    'charges': charges,
    'residues_map': residues_map,
    'pbc_residues': pbc_residues,
    'output_topology_filename': OUTPUT_topology_filename
}, 'topology')

# Prepare screenshot file 
screenshot = File(OUTPUT_screenshot_filename, get_screenshot, {
    'input_structure_filename' : OUTPUT_pdb_filename,
    'output_screenshot_filename' : OUTPUT_screenshot_filename,
}, 'screenshot')

# Get the cutoff for the test below
rmsd_cutoff = Dependency(get_input, {'name': 'rmsd_cutoff'})
# Set a test to check trajectory integrity
sudden_jumps = Dependency(check_sudden_jumps, {
    'input_structure_filename': pdb_filename,
    'input_trajectory_filename': trajectory_filename,
    'structure': structure,
    'pbc_residues': pbc_residues,
    'register': register,
    #'time_length': time_length,
    'check_selection': protein_and_nucleic,
    'standard_deviations_cutoff': rmsd_cutoff,
})

# Pack up all tools which may be called directly from the console
tools = [
    corrector,
    interactions,
    snapshots,
    charges,
    residues_map,
    screenshot,
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
        'snapshots': snapshots,
        'structure': structure,
        'pbc_residues': pbc_residues,
    }, 'rmsds'),
    # Here we set a small frames limit since this anlaysis is a bit slow
    File(OUTPUT_tmscores_filename, tmscores, {
        "input_topology_filename": pdb_filename,
        "input_trajectory_filename": trajectory_filename,
        "output_analysis_filename": OUTPUT_tmscores_filename,
        'first_frame_filename': first_frame_filename,
        'average_structure_filename': average_structure_filename,
        'structure' : structure,
        'pbc_residues': pbc_residues,
        'snapshots': snapshots,
        'frames_limit': 200,
    }, 'tmscores'),
    # This analysis is fast and the output size depends on the number of atoms only
    # For this reason here it is used the whole trajectory with no frames limit
    File(OUTPUT_rmsf_filename, rmsf, {
        "input_topology_filename": pdb_filename,
        "input_trajectory_filename": trajectory_filename,
        "output_analysis_filename": OUTPUT_rmsf_filename,
        'structure': structure,
        'pbc_residues': pbc_residues,
    }, 'rmsf'),
    # WARNING: This analysis is fast enought to use the full trajectory instead of the reduced one
    # WARNING: However, the output file size depends on the trajectory size
    # WARNING: In very long trajectories the number of points may make the client go slow when loading data
    File(OUTPUT_rgyr_filename, rgyr, {
        "input_topology_filename": pdb_filename,
        "input_trajectory_filename": trajectory_filename,
        "output_analysis_filename": OUTPUT_rgyr_filename,
        'snapshots': snapshots,
        'frames_limit': 5000,
        'structure': structure,
        'pbc_residues': pbc_residues,
    }, 'rgyr'),
    # WARNING: This analysis will generate several output files
    # File 'pca.average.pdb' is generated by the PCA and used by the client
    # File 'covar.log' is generated by the PCA but never used
    File(OUTPUT_pca_filename, pca, {
        "input_topology_filename": pdb_filename,
        "input_trajectory_filename": trajectory_filename,
        "output_analysis_filename": OUTPUT_pca_filename,
        "output_trajectory_projections_prefix": OUTPUT_pca_trajectory_projection_prefix,
        'snapshots' : snapshots,
        'frames_limit': 2000,
        'structure': structure,
        'fit_selection': pca_fit_selection,
        'analysis_selection': pca_selection,
        'pbc_residues': pbc_residues,
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
        'pbc_residues': pbc_residues,
        'snapshots': snapshots,
        'frames_limit': 100,
    }, 'rmsdperres'),
    # WARNING: This analysis is fast enought to use the full trajectory instead of the reduced one
    # WARNING: However, the output file size depends on the trajectory size exponentially. It may be huge
    File(OUTPUT_rmsdpairwise_filename, rmsd_pairwise, {
        'input_topology_filename': pdb_filename,
        'input_trajectory_filename': trajectory_filename,
        "output_analysis_filename": OUTPUT_rmsdpairwise_filename,
        "interactions": interactions,
        'structure': structure,
        'pbc_residues': pbc_residues,
        'snapshots': snapshots,
        'frames_limit': 200,
        'overall_selection': "name CA or name C5"
    }, 'rmsdpairwise'),
    # WARNING: This analysis is not fast enought to use the full trajectory. It would take a while
    File(OUTPUT_distperres_filename, distance_per_residue, {
        'input_topology_filename': pdb_filename,
        'input_trajectory_filename': trajectory_filename,
        "output_analysis_filename": OUTPUT_distperres_filename,
        "interactions": interactions,
        'snapshots': snapshots,
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
        'snapshots': snapshots,
        'frames_limit': 200,
        # 'is_time_dependend': is_time_dependend,
        # 'time_splits': 100,
        # 'populations': populations
    }, 'hbonds'),
    File(OUTPUT_sasa_filename, sasa, {
        "input_topology_filename": pdb_filename,
        "input_trajectory_filename": trajectory_filename,
        "output_analysis_filename": OUTPUT_sasa_filename,
        'structure': structure,
        'pbc_residues': pbc_residues,
        'snapshots': snapshots,
        'frames_limit': 100,
    }, 'sasa'),
    File(OUTPUT_energies_filename, energies, {
        "input_topology_filename": pdb_filename,
        "input_trajectory_filename": trajectory_filename,
        "output_analysis_filename": OUTPUT_energies_filename,
        'structure': structure,
        "interactions": interactions,
        'charges': charges,
        'snapshots': snapshots,
        'frames_limit': 100,
    }, 'energies'),
    File(OUTPUT_pockets_filename, pockets, {
        "input_topology_filename": pdb_filename,
        "input_trajectory_filename": trajectory_filename,
        "output_analysis_filename": OUTPUT_pockets_filename,
        'structure': structure,
        'pbc_residues': pbc_residues,
        'snapshots': snapshots,
        'frames_limit': 100,
    }, 'pockets'),
    File(OUTPUT_helical_parameters_filename, helical_parameters, {
        'input_topology_filename': pdb_filename,
        'input_trajectory_filename': trajectory_filename,
        "output_analysis_filename": OUTPUT_helical_parameters_filename,
        'structure': structure,
        'frames_limit': 1000,
    }, 'helical'),
    File(OUTPUT_markov_filename, markov, {
        "input_topology_filename": pdb_filename,
        "input_trajectory_filename": trajectory_filename,
        "output_analysis_filename": OUTPUT_markov_filename,
        'structure': structure,
        'populations': populations,
        'rmsd_selection': protein_and_nucleic,
    }, 'markov'),
]

# Set a list with all dependencies to be required if the whole workflow is run
# Metadata is the last dependency since it contains the regiter warnings
basic_dependencies = [ topology_filename, *analyses, metadata_filename ]

# Set a list with all dependencies which may be requested independently
requestables = [ *analyses, *tools, metadata_filename, topology_filename ]

# Main ---------------------------------------------------------------------------------            

# Function called through argparse
def main ():

    # Parse input arguments from the console
    args = parser.parse_args()
    # Set the command line inputs
    # The vars function converts the args object to a dictionary
    workflow_args = vars(args)
    # Call the actual main function
    workflow(**workflow_args)

# Set default input filenames
# They may be modified through console command arguments
#DEFAULT_working_directory = str(Path.cwd())
DEFAULT_input_topology_filename = OUTPUT_pdb_filename
DEFAULT_input_trajectory_filenames = [OUTPUT_trajectory_filename]
DEFAULT_inputs_filename = 'inputs.json'
DEFAULT_input_charges_filename = find_charges_filename()
DEFAULT_input_populations_filename = 'populations.json'
DEFAULT_database_url = 'https://mdposit-dev.mddbr.eu/api'
DEFAULT_rmsd_cutoff = 9
DEFAULT_interaction_cutoff = 0.1

# The actual main function
def workflow (
    #working_dir : str = DEFAULT_working_directory,
    input_topology_filename : str = DEFAULT_input_topology_filename,
    input_trajectory_filenames : List[str] = DEFAULT_input_trajectory_filenames,
    inputs_filename : str = DEFAULT_inputs_filename,
    input_charges_filename : str = DEFAULT_input_charges_filename,
    input_populations_filename : str = DEFAULT_input_populations_filename,
    database_url : str = DEFAULT_database_url,
    project : Optional[str] = None,
    sample_trajectory : bool = False,
    download : bool = False,
    setup : bool = False,
    include : Optional[List[str]] = None,
    exclude : Optional[List[str]] = None,
    filter_selection : Union[bool, str] = False,
    image : bool = False, fit : bool = False,
    translation : List[float] = [0, 0, 0],
    mercy : Union[ List[str], bool ] = [],
    trust : Union[ List[str], bool ] = [],
    pca_selection : str = protein_and_nucleic_backbone,
    pca_fit_selection : str = protein_and_nucleic_backbone,
    rmsd_cutoff : float = DEFAULT_rmsd_cutoff,
    interaction_cutoff : float = DEFAULT_interaction_cutoff
):

    # Fix the input_trajectory_filenames argument: in case it is a string convert it to a list
    if type(input_trajectory_filenames) == str:
        input_trajectory_filenames = [input_trajectory_filenames]

    # Fix the mercy input, if needed
    # If a boolean is passed instead of a list we set its corresponding value
    if type(mercy) == bool:
        if mercy:
            mercy = available_failures
        else:
            mercy = []

    # Fix the trust input, if needed
    # If a boolean is passed instead of a list we set its corresponding value
    if type(trust) == bool:
        if trust:
            trust = available_checkings
        else:
            trust = []
        

    # Update the inputs variable with all current function arguments
    # WARNING: Do not declare any variable over this line or it will be included in the inputs and thus in the register
    inputs.update(locals())

    # Reset all dependencies
    # This is useful in case we run this function more than once without reseting the module
    # Note that this script was originally called only from argparse so this was not necessary
    for var_value in dict(globals()).values():
        if isinstance(var_value, Dependency):
            var_value.reset()

    # Load the inputs file
    load_inputs(inputs_filename)

    # If download is passed as True then download all data we may need further and exit
    if download:
        original_topology_filename.value
        original_trajectory_filenames.value
        original_charges_filename.value
        populations.value
        # Note that there is no need to save the register if we just downloaded data
        return

    # Run tools which must be run always
    # They better be fast
    process_input_files.value

    # Run some checkings
    must_check_trajectory_integrity = 'intrajrity' not in trust
    if must_check_trajectory_integrity and sudden_jumps.value:
        must_be_killed = 'intrajrity' not in mercy
        if must_be_killed:
            raise SystemExit('Failed RMSD checking')

    # If setup is passed as True then exit as soon as the setup is finished
    if setup:
        pdb_filename.value
        trajectory_filename.value
        charges_filename.value
        # Save the register
        save_register(register.value)
        return

    # Run the requested analyses
    if include and len(include) > 0:
        sys.stdout.write(f"Executing specific dependencies: " + str(include) + '\n')
        # Include only the specified dependencies
        requested_dependencies = [ dep for dep in requestables if dep.alias in include ]
        for dependency in requested_dependencies:
            dependency.value
    # Run all analyses if none was specified
    else:
        # Include all non excluded dependencies
        requested_dependencies = basic_dependencies
        if exclude and len(exclude) > 0:
            requested_dependencies = [ dep for dep in basic_dependencies if dep.alias not in exclude ]
        for dependency in requested_dependencies:
            dependency.value

    # Remove gromacs backups
    remove_trash()

    # Save the register
    save_register(register.value)

    sys.stdout.write("Done!\n")


# Define console arguments to call the workflow
parser = ArgumentParser(description="MoDEL Workflow", formatter_class=RawTextHelpFormatter)

# Set optional arguments
# parser.add_argument(
#     "-dir", "--working_dir",
#     default=DEFAULT_working_directory,
#     help="Directory where to perform analysis. "
#     "If empty, will use current directory.")

parser.add_argument(
    "-proj", "--project",
    default=None,
    help="If given a project name, trajectory and "
    "topology files will be downloaded from remote server.")

parser.add_argument(
    "-url", "--database_url",
    default=DEFAULT_database_url,
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
    "-pop", "--input_populations_filename",
    default=DEFAULT_input_populations_filename,
    help="Path to populations filename")

parser.add_argument(
    "-img", "--image",
    action='store_true',
    help="Set if the trajectory is to be imaged")

parser.add_argument(
    "-fit", "--fit",
    action='store_true',
    help="Set if the trajectory is to be fitted (both rotation and translation)")

parser.add_argument(
    "-trans", "--translation",
    nargs='*',
    default=[0,0,0],
    help=("Set the x y z translation for the imaging process\n"
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

# Set a custom argparse action to handle the following 2 arguments
# This is done becuase it is not possible to combine nargs='*' with const
# https://stackoverflow.com/questions/72803090/argparse-how-to-create-the-equivalent-of-const-with-nargs
class custom (Action):
    # If argument is not passed -> default
    # If argument is passed empty -> const
    # If argument is passed with values -> values
    def __call__(self, parser, namespace, values, option_string=None):
        if values:
            setattr(namespace, self.dest, values)
        else:
            setattr(namespace, self.dest, self.const)

parser.add_argument(
    "-t", "--trust",
    type=str,
    nargs='*',
    default=[],
    action=custom,
    const=available_checkings,
    choices=available_checkings,
    help=("If passed, do not run the specified checking. Note that all checkings are skipped if passed alone.\n"
        "Available protocols:\n"
        "- stabonds - Stable bonds\n"
        "- cohbonds - Coherent bonds\n"
        "- intrajrity - Trajectory integrity")
)

parser.add_argument(
    "-m", "--mercy",
    type=str,
    nargs='*',
    default=[],
    action=custom,
    const=available_failures,
    choices=available_failures,
    help=("If passed, do not kill the process when any of the specfied checkings fail and proceed with the workflow.\n"
        "Note that all checkings are allowed to fail if the argument is passed alone.\n"
        "Available protocols:\n"
        "- stabonds - Stable bonds\n"
        "- cohbonds - Coherent bonds\n"
        "- intrajrity - Trajectory integrity\n"
        "- refseq - Reference sequence")
)

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
    "-rcut", "--rmsd_cutoff",
    type=float,
    default=DEFAULT_rmsd_cutoff,
    help=("Set the cutoff for the RMSD sudden jumps analysis to fail (default " + str(DEFAULT_rmsd_cutoff) + ").\n"
        "This cutoff stands for the number of standard deviations away from the mean an RMSD value is to be.\n"))

parser.add_argument(
    "-icut", "--interaction_cutoff",
    type=float,
    default=DEFAULT_interaction_cutoff,
    help=("Set the cutoff for the interactions analysis to fail (default " + str(DEFAULT_interaction_cutoff) + ").\n"
        "This cutoff stands for percent of the trajectory where the interaction happens (from 0 to 1).\n"))
