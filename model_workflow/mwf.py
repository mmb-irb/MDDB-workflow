#!/usr/bin/env python

# This is the starter script

# Import python libraries
from os import remove
from os.path import exists
import sys
import io
import math
from pathlib import Path
import urllib.request
import json
import numpy
from datetime import datetime
from typing import Optional, Union, List

# Constants
from model_workflow.constants import *

# Import local tools
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
from model_workflow.tools.get_screenshot import get_screenshot
from model_workflow.tools.formats import get_file_standard_format
from model_workflow.tools.filter_atoms import filter_atoms
from model_workflow.tools.image_and_fit import image_and_fit
from model_workflow.tools.structure_corrector import structure_corrector
from model_workflow.tools.httpsf import mount

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

# Import mdtoolbelt tools
from mdtoolbelt.conversions import convert
from mdtoolbelt.structures import Structure

# Make the system output stream to not be buffered
# This is useful to make prints work on time in Slurm
# Otherwise, output logs are written after the script has fully run
# Note that this fix affects all modules and built-ins
unbuffered = io.TextIOWrapper(open(sys.stdout.fileno(), 'wb', 0), write_through=True)
sys.stdout = unbuffered

# Set an special exception for input errors
missing_input_exception = Exception('Missing input')

# The project is the main handler
class Project:
    def __init__ (self,
        # The local directory where the project takes place
        directory : str = '.',
        # Accession of the project in the database, given that this project is already uploaded
        accession : Optional[str] = None,
        # URL to query for missing files when an accession is provided
        database_url : str = DEFAULT_API_URL,
        # A file containing a lof of inputs related to metadata, MD simulation parameters and analysis configurations
        inputs_filename : str = DEFAULT_INPUTS_FILENAME,
        # The input topology filename
        # Multiple formats are accepted but the default is our own parsed json topology
        input_topology_filename : str = None,
        # Input structure and trajectory filenames
        input_structure_filename : str = STRUCTURE_FILENAME,
        input_trajectory_filenames : Union[ str, List[str] ] = [ TRAJECTORY_FILENAME ],
        populations_filename : str = DEFAULT_POPULATIONS_FILENAME,
        transitions_filename : str = DEFAULT_TRANSITIONS_FILENAME,
        # DANI: A partir de aquí no estoy muy convencido de que estos argumentos vayan dentro del project
        filter_selection : Union[bool, str] = False,
        image : bool = False,
        fit : bool = False,
        translation : List[float] = [0, 0, 0],
        mercy : Union[ List[str], bool ] = [],
        trust : Union[ List[str], bool ] = [],
        pca_selection : str = PROTEIN_AND_NUCLEIC_BACKBONE,
        pca_fit_selection : str = PROTEIN_AND_NUCLEIC_BACKBONE,
        rmsd_cutoff : float = DEFAULT_RMSD_CUTOFF,
        interaction_cutoff : float = DEFAULT_INTERACTION_CUTOFF,
        # Set it we must download just a few frames instead of the whole trajectory
        sample_trajectory : bool = False
    ):
        # Save input parameters
        self.directory = directory
        self.accession = accession
        self.database_url = database_url
        # Set the project URL in case we have the required data
        self.project_url = None
        if database_url and accession:
            self.project_url = database_url + '/rest/current/projects/' + accession

        # Set internal variables for the input filenames
        self._inputs_filename = inputs_filename
        self._input_topology_filename = input_topology_filename
        if not input_topology_filename:
            self._input_topology_filename = self.guess_input_topology_filename()
        self._input_structure_filename = input_structure_filename
        self._input_trajectory_filenames = input_trajectory_filenames
        # Fix the input_trajectory_filenames argument: in case it is a string convert it to a list
        if type(input_trajectory_filenames) == list:
            self._input_trajectory_filenames = input_trajectory_filenames
        elif type(input_trajectory_filenames) == str:
            self._input_trajectory_filenames = [input_trajectory_filenames]
        else:
            raise Exception('Input trajectory filenames must be a list of strings or a string')
        self._populations_filename = populations_filename
        self._transitions_filename = transitions_filename

        self.filter_selection = filter_selection
        self.image = image
        self.fit = fit
        self.translation = translation
        self.mercy = mercy
        # Fix the mercy input, if needed
        # If a boolean is passed instead of a list we set its corresponding value
        if type(mercy) == bool:
            if mercy:
                self.mercy = AVAILABLE_FAILURES
            else:
                self.mercy = []
        self.trust = trust
        # Fix the trust input, if needed
        # If a boolean is passed instead of a list we set its corresponding value
        if type(trust) == bool:
            if trust:
                self.trust = AVAILABLE_CHECKINGS
            else:
                self.trust = []
        self.pca_selection = pca_selection
        self.pca_fit_selection = pca_fit_selection
        self.rmsd_cutoff = rmsd_cutoff
        self.interaction_cutoff = interaction_cutoff
        self.sample_trajectory = sample_trajectory

        # Set the inputs, where values from the inputs file will be stored
        self._inputs = None

        # Set the processed topology filename, which depends on the input topology filename
        # Note that this file is different from the standard topology, although it may be standard as well
        self._topology_filename = self.inherit_topology_filename()

        # Other values which may be found/calculated on demand
        self._available_files = None
        self._sudden_jumps = None
        self._snapshots = None
        self._structure = None
        self._pytraj_topology = None
        self._interactions = None
        self._pbc_residues = None
        self._charges = None
        self._populations = None
        self._transitions = None
        self._residues_map = None

        # Set a new entry for the register
        # This is useful to track previous workflow runs and problems
        self.register = {
            'date': datetime.today().strftime('%d-%m-%Y %H:%M:%S'),
            'inputs': {
                'directory': directory,
                'accession': accession,
                'database_url': database_url,
                'inputs_filename': inputs_filename,
                'input_topology_filename': input_topology_filename,
                'input_structure_filename': input_structure_filename,
                'input_trajectory_filenames': input_trajectory_filenames,
                'populations_filename': populations_filename,
                'transitions_filename': transitions_filename,
                'filter_selection': filter_selection,
                'image': image,
                'fit': fit,
                'translation': translation,
                'mercy': mercy,
                'trust': trust,
                'pca_selection': pca_selection,
                'pca_fit_selection': pca_fit_selection,
                'rmsd_cutoff': rmsd_cutoff,
                'interaction_cutoff': interaction_cutoff,
                'sample_trajectory': sample_trajectory,
            },
            'warnings': [],
        }

    # Check input files exist when their filenames are read
    # If they do not exist then try to download them
    # If the download is not possible then raise an error

    # Inputs filename ------------

    # Set a function to load the inputs file
    def get_inputs_filename (self) -> str:
        # There must be an inputs filename
        if not self._inputs_filename:
            raise Exception('Not defined inputs filename')
        # If the file already exists then we are done
        if exists(self._inputs_filename):
            return self._inputs_filename
        # Try to download it
        # If we do not have the required parameters to download it then we surrender here
        if not self.project_url:
            raise Exception('Missing inputs file "' + self._inputs_filename + '"')
        # Download the inputs json file if it does not exists
        sys.stdout.write('Downloading inputs (' + self._inputs_filename + ')\n')
        inputs_url = self.project_url + '/inputs/'
        urllib.request.urlretrieve(inputs_url, self._inputs_filename)
        # Rewrite the inputs file in a pretty formatted way
        with open(self._inputs_filename, 'r+') as file:
            file_content = json.load(file)
            file.seek(0)
            json.dump(file_content, file, indent=4)
            file.truncate()
        return self._inputs_filename
    inputs_filename = property(get_inputs_filename, None, None, "Inputs filename (read only)")

    # Topology filename ------------

    # If there is not input topology filename, which is possible, we must guess the topology filename
    # Note that if we can download then the name will always be the remote topology name (topology.json)
    def guess_input_topology_filename (self) -> Optional[str]:
        # If we can download then the we use the remote topology filename (topology.json)
        if self.project_url:
            return TOPOLOGY_FILENAME
        # Otherwise we must guess among the local files
        # If the default topology filename exists then use it
        default_topology_path = self.directory + '/' + TOPOLOGY_FILENAME
        if exists(default_topology_path):
            return default_topology_path
        # Find if the raw charges file is present
        raw_charges_path = self.directory + '/' + RAW_CHARGES_FILENAME
        if exists(raw_charges_path):
            return raw_charges_path
        # Otherwise, find all possible accepted topology formats
        for topology_format in ACCEPTED_TOPOLOGY_FORMATS:
            topology_path = self.directory + '/topology.' + topology_format
            if exists(topology_path):
                return topology_path
        return None        

    # Get the input charges filename from the inputs
    # If the file is not found try to download it
    def get_input_topology_filename (self) -> str:
        # There must be a topology filename
        if not self._input_topology_filename:
            return None
        # If the file already exists then we are done
        if exists(self._input_topology_filename):
            return self._input_topology_filename
        # Try to download it
        # If we do not have the required parameters to download it then we surrender here
        if not self.project_url:
            raise Exception('Missing input topology file "' + inputs_filename + '"')
        # If the input topology name is the standard then download it from the topology endpoint
        if self._input_topology_filename == TOPOLOGY_FILENAME:
            # Check if the project has a topology and download it in json format if so
            topology_url = self.project_url + '/topology'
            try:
                sys.stdout.write('Downloading topology (' + self._input_topology_filename + ')\n')
                urllib.request.urlretrieve(topology_url, self._input_topology_filename)
                # Rewrite the topology file in a pretty formatted way
                with open(self._input_topology_filename, 'r+') as file:
                    file_content = json.load(file)
                    file.seek(0)
                    json.dump(file_content, file, indent=4)
                    file.truncate()
            except:
                raise Exception('Something where wrong while downloading the topology')
            return self._input_topology_filename
        # Otherwise, try to download it using the files endpoint
        # Note that this is not usually required
        topology_url = self.project_url + '/files/' + self._input_topology_filename
        try:
            sys.stdout.write('Downloading topology (' + self._input_topology_filename + ')\n')
            urllib.request.urlretrieve(topology_url, self._input_topology_filename)
        except:
            raise Exception('Something where wrong while downloading the topology')
        # Before we finish, in case the topology is a '.top' file, we may need to download the itp files as well
        if self._input_topology_filename == 'topology.top':
            # Find available .itp files and download each of them
            itp_filenames = [filename for filename in self.available_files if filename[-4:] == '.itp']
            for itp_filename in itp_filenames:
                sys.stdout.write('Downloading itp file (' + itp_filename + ')\n')
                itp_url = self.project_url + '/files/' + itp_filename
                urllib.request.urlretrieve(itp_url, itp_filename)
        return self._input_topology_filename
    input_topology_filename = property(get_input_topology_filename, None, None, "Input topology filename (read only)")

    # Input structure filename ------------

    # Get the input pdb filename from the inputs
    # If the file is not found try to download it
    def get_input_structure_filename (self) -> str:
        # There must be an input structure filename
        if not self._input_structure_filename:
            raise Exception('Not defined input structure filename')
        # If the file already exists then we are done
        if exists(self._input_structure_filename):
            return self._input_structure_filename
        # Try to download it
        # If we do not have the required parameters to download it then we surrender here
        if not self.project_url:
            raise Exception('Missing inputs file "' + self._input_structure_filename + '"')
        sys.stdout.write('Downloading structure (' + self._input_structure_filename + ')\n')
        # Set the download URL
        structure_url = self.project_url + '/files/' + self._input_structure_filename
        # If the structure filename is the standard structure filename then use the structure endpoint instead
        if self._input_structure_filename == STRUCTURE_FILENAME:
            structure_url = self.project_url + '/structure'
        # Download the file
        try:
            urllib.request.urlretrieve(structure_url, self._input_structure_filename)
        except urllib.error.HTTPError as error:
            if error.code == 404:
                raise Exception('ERROR: Missing input pdb file "' + self._input_structure_filename + '"')
            else:
                raise Exception('Something went wrong with the MDposit request: ' + structure_url)
        return self._input_structure_filename
    input_structure_filename = property(get_input_structure_filename, None, None, "Input structure filename (read only)")

    # Input trajectory filename ------------

    # Get the input trajectory filename(s) from the inputs
    # If file(s) are not found try to download it
    def get_input_trajectory_filenames (self) -> str:
        # There must be an input structure filename
        if not self._input_trajectory_filenames or len(self._input_trajectory_filenames) == 0:
            raise Exception('Not defined input trajectory filenames')
        # If all files already exists then we are done
        missing_input_trajectory_files = [ filename for filename in self._input_trajectory_filenames if not exists(filename) ]
        if len(missing_input_trajectory_files) == 0:
            return self._input_trajectory_filenames
        # Try to download the missing files
        # If we do not have the required parameters to download it then we surrender here
        if not self.project_url:
            raise Exception('Missing input trajectory files: ' + ', '.join(missing_input_trajectory_files))
        # Download each trajectory file (ususally it will be just one)
        for trajectory_filename in self._input_trajectory_filenames:
            sys.stdout.write('Downloading trajectory (' + trajectory_filename + ')\n')
            # Set the trajectory URL, which may depend in different parameters
            trajectory_url = self.project_url + '/files/' + trajectory_filename
            if trajectory_filename == TRAJECTORY_FILENAME:
                trajectory_url = self.project_url + '/trajectory?format=xtc'
                if self.sample_trajectory:
                    trajectory_url += '&frames=1:10:1'
            # Download the file
            try:
                urllib.request.urlretrieve(trajectory_url, trajectory_filename)
            except urllib.error.HTTPError as error:
                if error.code == 404:
                    raise Exception('Missing input trajectory file "' + trajectory_filename + '"')
                raise Exception('Something went wrong with the MDposit request: ' + trajectory_url)
        return self._input_trajectory_filenames
    input_trajectory_filenames = property(get_input_trajectory_filenames, None, None, "Input trajectory filenames (read only)")

    # Populations filename ------------

    def get_populations_filename (self) -> str:
        self.get_file(self._populations_filename)
        return self._populations_filename
    populations_filename = property(get_populations_filename, None, None, "MSM equilibrium populations filename (read only)")

    # Transitions filename ------------

    def get_transitions_filename (self) -> str:
        self.get_file(self._transitions_filename)
        return self._transitions_filename
    transitions_filename = property(get_transitions_filename, None, None, "MSM transition probabilities filename (read only)")

    # ---------------------------------

    # Remote available files
    def get_available_files (self) -> List[str]:
        # If we already have a stored value then return it
        if self._available_files != None:
            return self._available_files
        # If we have not remote access then there are not available files
        if not self.project_url:
            return []
        # Get the available files
        files_url = self.project_url + '/files'
        response = urllib.request.urlopen(files_url)
        filenames = json.loads(response.read())
        self._available_files = filenames
        return self._available_files
    available_files = property(get_available_files, None, None, "Remote available files (read only)")

    # Check if a file exists
    # If not, try to download it from the database
    # If the file is not found in the database it is fine, we do not even warn the user
    # Note that this function is used to get populations and transitions files, which are not common
    def get_file (self, filename : str) -> bool:
        # If it exists we are done
        if exists(filename):
            return True
        # Try to download the missing file
        # If we do not have the required parameters to download it then we surrender here
        if not self.project_url:
            return False
        # Check if the file is among the available remote files
        # If it is no then stop here
        if filename not in self.available_files:
            return False
        sys.stdout.write('Downloading file (' + filename + ')\n')
        # Set the download URL
        file_url = self.project_url + '/files/' + filename
        # Download the file
        try:
            urllib.request.urlretrieve(file_url, filename)
        except urllib.error.HTTPError as error:
            if error.code == 404:
                print('Missing file "' + filename + '"')
                return False
            else:
                raise Exception('Something went wrong with the MDposit request: ' + file_url)
        return True

    # Input file values -----------------------------------------

    # First of all set input themselves

    # Get inputs
    def get_inputs (self) -> dict:
        # If inputs are already loaded then return them
        if self._inputs:
            return self._inputs
        # Otherwise, load inputs from the inputs file
        with open(self.inputs_filename, 'r') as file:
            inputs_data = json.load(file)
            self._inputs = inputs_data
        # Finally return the updated inputs
        return self._inputs
    inputs = property(get_inputs, None, None, "Inputs from the inputs file (read only)")

    # Then set getters for every value in the inputs file

    # Set a function to get 'inputs' values and handle missing keys
    def input_getter (name : str, missing_input_callback = missing_input_exception):
        def getter (self):
            value = self.inputs.get(name, missing_input_callback)
            if value == missing_input_exception:
                raise Exception('Missing input "' + name + '"')
            return value
        return getter

    # Assign the getters
    input_interactions = property(input_getter('interactions'), None, None, "Interactions to be analyzed (read only)")
    pbc_selection = property(input_getter('pbc_selection'), None, None, "Selection of atoms which are still in periodic boundary conditions (read only)")
    forced_references = property(input_getter('forced_references'), None, None, "Uniprot IDs to be used first when aligning protein sequences (read only)")
    pdb_ids = property(input_getter('type'), None, None, "Protein Data Bank IDs used for the setup of the system (read only)")
    input_type = property(input_getter('type'), None, None, "Set if its a trajectory or an ensemble (read only)")

    # Set additional values infered from input values

    # Set if trajectory is time dependent
    def check_is_time_dependent (self) -> bool:
        if self.input_type == 'trajectory':
            return True
        elif self.input_type == 'ensemble':
            return False
        raise SystemExit('Not supported input type value: ' + self.input_type)
    is_time_dependend = property(check_is_time_dependent, None, None, "Check if trajectory frames are time dependent (read only)")

    # Register ---------------

    # Save the current register to a file
    # In case there is a previous register, read it and append data to it
    def save_register (self):
        register_data = []
        # Check if there is previous register data
        if exists(REGISTER_FILENAME):
            with open(REGISTER_FILENAME, 'r') as file:
                register_data = json.load(file)
        # Add current register data
        register_data.append(self.register)
        # Write it to a json file
        with open(REGISTER_FILENAME, 'w') as file:
            json.dump(register_data, file, indent=4)

    # Processed files ----------------------------------------------------

    # Set the expected output topology filename given the input topology filename
    # Note that topology formats are conserved
    def inherit_topology_filename (self) -> str:
        if not self._input_topology_filename:
            return None
        if self._input_topology_filename == TOPOLOGY_FILENAME:
            return self._input_topology_filename
        if self._input_topology_filename == TOPOLOGY_FILENAME:
            return self._input_topology_filename
        standard_format = get_file_standard_format(self._input_topology_filename)
        return 'topology.' + standard_format        

    # Process input files to generate the processed files
    # This process corrects and standarizes the topology, the trajectory and the structure
    def process_input_files (self):

        # Set the output filenames
        output_structure_filename = STRUCTURE_FILENAME
        output_trajectory_filename = TRAJECTORY_FILENAME
        output_topology_filename = self._topology_filename

        print('--- Processing input files ---')

        # --- CONVERTING AND MERGING ------------------------------------------------------------

        # Set output filenames for the already converted structure
        input_structure_format = get_file_standard_format(self.input_structure_filename)
        output_structure_format = get_file_standard_format(output_structure_filename)
        converted_structure_filename = CONVERTED_STRUCTURE
        # If input structure already matches the output format then avoid the renaming
        if input_structure_format == output_structure_format:
            onverted_structure_filename = self.input_structure_filename
        # Set output filenames for the already converted structure
        # Input trajectories should have all the same format
        input_trajectories_formats = set([ get_file_standard_format(trajectory) for trajectory in self.input_trajectory_filenames ])
        if len(input_trajectories_formats) > 1:
            raise Exception('All input trajectory files must have the same format')
        input_trajectories_format = input_trajectories_formats[0]
        output_trajectory_format = get_file_standard_format(output_trajectory_filename)
        converted_trajectory_filename = CONVERTED_TRAJECTORY
        # If input trajectory already matches the output format and is unique then avoid the renaming
        if input_trajectory_format == output_trajectory_format and len(self.input_trajectory_filenames) == 1:
            converted_trajectory_filename = self.input_trajectory_filenames[0]

        # Convert input structure and trajectories to output structure and trajectory
        if not exists(converted_structure_filename) or not exists(converted_trajectory_filename):
            print(' - Converting and merging')
            convert(
                input_structure_filename = self.input_structure_filename,
                output_structure_filename = converted_structure_filename,
                input_trajectory_filenames = self.input_trajectory_filenames,
                output_trajectory_filename = converted_trajectory_filename,
            )

        # Topologies are never converted, but they are kept in their original format

        # --- FILTERING ATOMS ------------------------------------------------------------

        # Find out if we need to filter
        # i.e. check if there is a selection filter and it matches some atoms
        must_filter = bool(self.filter_selection and Structure.from_pdb_file(converted_structure_filename).select(self.filter_selection))

        # Set output filenames for the already filtered structure and trajectory
        filtered_structure = FILTERED_STRUCTURE if must_filter else converted_structure_filename
        filtered_trajectory = FILTERED_TRAJECTORY if must_filter else converted_trajectory_filename
        
        # Filter atoms in structure, trajectory and topology if required and not done yet
        # Note that this is the only step affecting topology and thus here we output the definitive topology
        if must_filter and (not exists(filtered_structure) or not exists(filtered_trajectory)):
            print(' - Filtering atoms')
            filter_atoms(
                input_structure_filename = converted_structure_filename,
                input_trajectory_filename = converted_trajectory_filename,
                input_topology_filename = self.input_topology_filename, # We use input topology
                output_structure_filename = filtered_structure,
                output_trajectory_filename = filtered_trajectory,
                output_topology_filename = output_topology_filename, # We genereate the definitive topology
                filter_selection = self.filter_selection
            )

        # --- IMAGING AND FITTING ------------------------------------------------------------

        # There is no logical way to know if the trajectory is already imaged or it must be imaged
        # We rely exclusively in input flags
        must_image = self.image or self.fit

        # Set output filenames for the already filtered structure and trajectory
        imaged_structure = IMAGED_STRUCTURE if must_image else filtered_structure
        imaged_trajectory = IMAGED_TRAJECTORY if must_image else filtered_trajectory

        # Image the trajectory if it is required
        # i.e. make the trajectory uniform avoiding atom jumps and making molecules to stay whole
        # Fit the trajectory by removing the translation and rotation if it is required
        if must_image and (not exists(imaged_structure) or not exists(imaged_trajectory)):
            print(' - Imaging and fitting')
            image_and_fit(filtered_structure, filtered_trajectory, output_topology_filename,
                imaged_structure, imaged_trajectory, self.image, self.fit, self.translation, self.pbc_selection)

        # --- CORRECTING STRUCTURE ------------------------------------------------------------

        # Note that this step, although it is foucsed in the structure, requires also the trajectory
        # Also the trajectory may be altered in very rare cases where coordinates must be resorted

        # There is no possible reason to not correct the structure
        # This is the last step so the output files will be named as the output files of the whole processing

        # WARNING:
        # For the correcting function we need the number of snapshots and at this point it should not be defined
        # Snapshots are calculated by default from the already processed structure and trajectory
        # For this reason we can not rely on the public snaphsots getter
        # We must calculate snapshots here using last step structure and trajectory
        self._snapshots = get_frames_count(imaged_structure, imaged_trajectory)

        # Correct the structure
        structure_corrector(
            input_structure_filename = imaged_structure,
            input_trajectory_filename = imaged_trajectory,
            input_topology_filename = output_topology_filename,
            output_structure_filename = output_structure_filename,
            output_trajectory_filename = output_trajectory_filename,
            snapshots = self._snapshots,
            register = self.register,
            mercy = self.mercy,
            trust = self.trust
        )

        # --- Cleanup intermediate files

        intermediate_filenames = [
            CONVERTED_STRUCTURE, CONVERTED_TRAJECTORY,
            FILTERED_STRUCTURE, FILTERED_TRAJECTORY,
            IMAGED_STRUCTURE, IMAGED_TRAJECTORY
        ]
        for filename in intermediate_files:
            if exists(filename):
                remove(filename)

    # Get the processed topology
    def get_topology_filename (self) -> str:
        # If the file already exists then we are done
        if exists(self._topology_filename):
            return self._topology_filename
        # Otherwise, process input files to generate the processed topology
        self.process_input_files()
        return self._topology_filename
    topology_filename = property(get_topology_filename, None, None, "Topology filename (read only)")

    # Get the processed structure
    def get_structure_filename (self) -> str:
        # If the file already exists then we are done
        if exists(STRUCTURE_FILENAME):
            return STRUCTURE_FILENAME
        # Otherwise, process input files to generate the processed structure
        self.process_input_files()
        return STRUCTURE_FILENAME
    structure_filename = property(get_structure_filename, None, None, "Structure filename (read only)")

    # Get the processed trajectory
    def get_trajectory_filename (self) -> str:
        # If the file already exists then we are done
        if exists(TRAJECTORY_FILENAME):
            return TRAJECTORY_FILENAME
        # Otherwise, process input files to generate the processed trajectory
        self.process_input_files()
        return TRAJECTORY_FILENAME
    trajectory_filename = property(get_trajectory_filename, None, None, "Trajectory filename (read only)")

    # ---------------------------------------------------------------------------------
    # Others values which may be found/calculated and files to be generated on demand
    # ---------------------------------------------------------------------------------

    # Trajectory snapshots
    def get_snapshots (self) -> str:
        # If we already have a stored value then return it
        # Note that checking if the trajectory filename exists triggers all the processing logic
        # The processing logic is able to set the internal snapshots value as well
        if self.trajectory_filename and self._snapshots != None:
            return self._snapshots
        # Otherwise we must find the value
        # This happens when the input files are already porcessed and thus we did not yet count the frames
        self._snapshots = get_frames_count(self.structure_filename, self.trajectory_filename)
        return self._snapshots
    snapshots = property(get_snapshots, None, None, "Trajectory snapshots (read only)")

    # Parsed structure
    def get_structure (self) -> 'Structure':
        # If we already have a stored value then return it
        if self._structure:
            return self._structure
        # Otherwise we must set the structure
        # Note that this is not only the mdtoolbelt structure, but it also contains additional logic
        self._structure = setup_structure(self.structure_filename)
        return self._structure
    structure = property(get_structure, None, None, "Parsed structure (read only)")

    # Pytraj trajectory
    def get_pytraj_trajectory (self) -> 'TrajectoryIterator':
        # If we already have a stored value then return it
        if self._pytraj_topology:
            return self._pytraj_topology
        # Otherwise we must set the pytarj trajectory
        self._pytraj_topology = get_pytraj_trajectory(self.structure_filename, self.trajectory_filename)
        return self._pytraj_topology
    pytraj_trajectory = property(get_pytraj_trajectory, None, None, "Pytraj trajectory (read only)")

    # First frame filename
    def get_first_frame_filename (self) -> str:
        # If the file already exists then send it
        if exists(FIRST_FRAME_FILENAME):
            return FIRST_FRAME_FILENAME
        # Otherwise, generate it
        get_first_frame(self.structure_filename, self.trajectory_filename, FIRST_FRAME_FILENAME)
        return FIRST_FRAME_FILENAME
    first_frame_filename = property(get_first_frame_filename, None, None, "First frame (read only)")

    # Average structure filename
    def get_average_structure_filename (self) -> str:
        # If the file already exists then send it
        if exists(AVERAGE_STRUCTURE_FILENAME):
            return AVERAGE_STRUCTURE_FILENAME
        # Otherwise, generate it
        get_average(self.pytraj_trajectory, AVERAGE_STRUCTURE_FILENAME)
        return AVERAGE_STRUCTURE_FILENAME
    average_structure_filename = property(get_average_structure_filename, None, None, "Average structure filename (read only)")

    # The interactions
    # This is a bit exceptional since it is a value to be used and an analysis file to be generated
    def get_interactions (self) -> List[dict]:
        # If we already have a stored value then return it
        if self._interactions != None:
            return self._interactions
        print('-> Processing interactions')
        # Otherwise, process interactions
        self._interactions = process_interactions(
            input_interactions = self.input_interactions,
            topology_filename = self.structure_filename,
            trajectory_filename = self.trajectory_filename,
            structure = self.structure,
            snapshots = self.snapshots,
            interactions_file = OUTPUT_INTERACTIONS_FILENAME,
            mercy = self.mercy,
            frames_limit = 1000,
            interaction_cutoff = self.interaction_cutoff
        )
        return self._interactions
    interactions = property(get_interactions, None, None, "Interactions (read only)")

    # Indices of residues in periodic boundary conditions
    def get_pbc_residues (self) -> List[int]:
        # If we already have a stored value then return it
        if self._pbc_residues:
            return self._pbc_residues
        # Otherwise we must find the value
        self._pbc_residues = get_pbc_residues(self.structure, self.pbc_selection)
        return self._pbc_residues
    pbc_residues = property(get_pbc_residues, None, None, "Indices of residues in periodic boundary conditions (read only)")

    # Atom charges
    def get_charges (self) -> List[float]:
        # If we already have a stored value then return it
        if self._charges:
            return self._charges
        # Otherwise we must find the value
        self._charges = get_charges(self.topology_filename)
        return self._charges
    charges = property(get_charges, None, None, "Atom charges (read only)")

    # Equilibrium populations from a MSM
    def get_populations (self) -> List[float]:
        # If we already have a stored value then return it
        if self._populations:
            return self._populations
        # Otherwise we must find the value
        self._populations = read_file(self.populations_filename)
        return self._populations
    populations = property(get_populations, None, None, "Equilibrium populations from a MSM (read only)")

    # Transition probabilities from a MSM
    def get_transitions (self) -> List[List[float]]:
        # If we already have a stored value then return it
        if self._transitions:
            return self._transitions
        # Otherwise we must find the value
        self._transitions = read_file(self.transitions_filename)
        return self._transitions
    transitions = property(get_transitions, None, None, "Transition probabilities from a MSM (read only)")

    # Residues mapping
    def get_residues_map (self) -> dict:
        # If we already have a stored value then return it
        if self._residues_map:
            return self._residues_map
        print('-> Mapping')
        # Otherwise we must find the value
        self._residues_map = generate_map_online(
            structure = self.structure,
            register = self.register,
            mercy = self.mercy,
            forced_references = self.forced_references,
            pdb_ids = self.pdb_ids,
        )
        return self._residues_map
    residues_map = property(get_residues_map, None, None, "Residues mapping (read only)")

    # Metadata filename
    def get_metadata_filename (self) -> str:
        # If the file already exists then send it
        if exists(OUTPUT_METADATA_FILENAME):
            return OUTPUT_METADATA_FILENAME
        print('-> Generating metadata')
        # Otherwise, generate it
        generate_metadata(
            input_topology_filename = self.structure_filename,
            input_trajectory_filename = self.trajectory_filename,
            inputs_filename = self.inputs_filename,
            structure = self.structure,
            snapshots = self.snapshots,
            residues_map = self.residues_map,
            interactions = self.interactions,
            register = self.register,
            output_metadata_filename = OUTPUT_METADATA_FILENAME,
        )
        return OUTPUT_METADATA_FILENAME
    metadata_filename = property(get_metadata_filename, None, None, "Metadata filename (read only)")

    # Standard topology filename
    def get_standard_topology_filename (self) -> str:
        # If the file already exists then send it
        if exists(TOPOLOGY_FILENAME):
            return TOPOLOGY_FILENAME
        print('-> Generating topology')
        # Otherwise, generate it
        generate_topology(
            structure = self.structure,
            charges = self.charges,
            residues_map = self.residues_map,
            pbc_residues = self.pbc_residues,
            output_topology_filename = TOPOLOGY_FILENAME
        )
        return TOPOLOGY_FILENAME
    standard_topology_filename = property(get_standard_topology_filename, None, None, "Standard topology filename (read only)")

    # Screenshot filename
    def get_screenshot_filename (self) -> str:
        # If the file already exists then send it
        if exists(OUTPUT_SCREENSHOT_FILENAME):
            return OUTPUT_SCREENSHOT_FILENAME
        print('-> Generating screenshot')
        # Otherwise, generate it
        get_screenshot(
            input_structure_filename = self.structure_filename,
            output_screenshot_filename = OUTPUT_SCREENSHOT_FILENAME,
        )
        return OUTPUT_SCREENSHOT_FILENAME
    screenshot_filename = property(get_screenshot_filename, None, None, "Screenshot filename (read only)")

    # Sudden jumps test
    def get_sudden_jumps (self) -> bool:
        # If we already have a stored value then return it
        if self._sudden_jumps != None:
            return self._sudden_jumps
        # Otherwise we must find the value
        self._sudden_jumps = check_sudden_jumps(
            input_structure_filename = self.structure_filename,
            input_trajectory_filename = self.trajectory_filename,
            structure = self.structure,
            pbc_residues = self.pbc_residues,
            mercy = self.mercy,
            trust = self.trust,
            register = self.register,
            # time_length = self.time_length,
            check_selection = PROTEIN_AND_NUCLEIC,
            standard_deviations_cutoff = self.rmsd_cutoff,
        )
        return self._sudden_jumps
    sudden_jumps = property(get_sudden_jumps, None, None, "Sudden jumps test (read only)")

    # ---------------------------------------------------------------------------------
    # Analyses
    # ---------------------------------------------------------------------------------

    # RMSDs
    def generate_rmsds_analysis (self):
        # Do not run the analysis if the output file already exists
        if exists(OUTPUT_RMSDS_FILENAME):
            return
        print('-> Running RMSDs analysis')
        # WARNING: This analysis is fast enought to use the full trajectory
        # WARNING: However, the output file size depends on the trajectory size
        # WARNING: In very long trajectories the number of points may make the client go slow when loading data
        rmsds(
            input_trajectory_filename = self.trajectory_filename,
            output_analysis_filename = OUTPUT_RMSDS_FILENAME,
            frames_limit = 5000,
            first_frame_filename = self.first_frame_filename,
            average_structure_filename = self.average_structure_filename,
            snapshots = self.snapshots,
            structure = self.structure,
            pbc_residues = self.pbc_residues,
        )

    # TM scores
    def generate_tmscores_analysis (self):
        # Do not run the analysis if the output file already exists
        if exists(OUTPUT_TMSCORES_FILENAME):
            return
        print('-> Running TM scores analysis')
        # Here we set a small frames limit since this anlaysis is a bit slow
        tmscores(
            input_topology_filename = self.structure_filename,
            input_trajectory_filename = self.trajectory_filename,
            output_analysis_filename = OUTPUT_TMSCORES_FILENAME,
            first_frame_filename = self.first_frame_filename,
            average_structure_filename = self.average_structure_filename,
            structure = self.structure,
            pbc_residues = self.pbc_residues,
            snapshots = self.snapshots,
            frames_limit = 200,
        )

    # RMSF, atom fluctuation
    def generate_rmsf_analysis (self):
        # Do not run the analysis if the output file already exists
        if exists(OUTPUT_RMSF_FILENAME):
            return
        print('-> Running RMSF analysis')
        # This analysis is fast and the output size depends on the number of atoms only
        # For this reason here it is used the whole trajectory with no frames limit
        rmsf(
            input_topology_filename = self.structure_filename,
            input_trajectory_filename = self.trajectory_filename,
            output_analysis_filename = OUTPUT_RMSF_FILENAME,
            structure = self.structure,
            pbc_residues = self.pbc_residues,
        )

    # RGYR, radius of gyration
    def generate_rgyr_analysis (self):
        # Do not run the analysis if the output file already exists
        if exists(OUTPUT_RGYR_FILENAME):
            return
        print('-> Running RGYR analysis')
        # WARNING: This analysis is fast enought to use the full trajectory instead of the reduced one
        # WARNING: However, the output file size depends on the trajectory size
        # WARNING: In very long trajectories the number of points may make the client go slow when loading data
        rgyr(
            input_topology_filename = self.structure_filename,
            input_trajectory_filename = self.trajectory_filename,
            output_analysis_filename = OUTPUT_RGYR_FILENAME,
            snapshots = self.snapshots,
            frames_limit = 5000,
            structure = self.structure,
            pbc_residues = self.pbc_residues,
        )

    # PCA, principal component analysis
    def generate_pca_analysis (self):
        # Do not run the analysis if the output file already exists
        if exists(OUTPUT_PCA_FILENAME):
            return
        print('-> Running PCA analysis')
        # WARNING: This analysis will generate several output files
        # File 'pca.average.pdb' is generated by the PCA and it was used by the client but not anymore
        # File 'covar.log' is generated by the PCA but never used
        pca(
            input_topology_filename = self.structure_filename,
            input_trajectory_filename = self.trajectory_filename,
            output_analysis_filename = OUTPUT_PCA_FILENAME,
            output_trajectory_projections_prefix = OUTPUT_PCA_PROJECTION_PREFIX,
            snapshots = self.snapshots,
            frames_limit = 2000,
            structure = self.structure,
            fit_selection = self.pca_fit_selection,
            analysis_selection = self.pca_selection,
            pbc_residues = self.pbc_residues,
        )

    # PCA contacts
    # DANI: Intenta usar mucha memoria, hay que revisar
    # DANI: Puede saltar un error de imposible alojar tanta memoria
    # DANI: Puede comerse toda la ram y que al final salte un error de 'Terminado (killed)'
    # DANI: De momento me lo salto
    # def generate_pca_contacts (self):
    #     # Do not run the analysis if the output file already exists
    #     if exists("md.pca.contacts.json"):
    #         return
    #     print('-> Running PCA contacts analysis')
    #     pca_contacts(
    #         trajectory = self.trajectory_filename,
    #         topology = self.pdb_filename,
    #         interactions = self.interactions,
    #         output_analysis_filename = "md.pca.contacts.json"
    #     )

    # RMSD per residue
    def genereate_rmsd_perres_analysis (self):
        # Do not run the analysis if the output file already exists
        if exists(OUTPUT_RMSD_PERRES_FILENAME):
            return
        print('-> Running RMSD per residue analysis')
        # WARNING: This analysis is fast enought to use the full trajectory instead of the reduced one
        # WARNING: However, the output file size depends on the trajectory size. It may be pretty big
        rmsd_per_residue(
            input_topology_filename = self.structure_filename,
            input_trajectory_filename = self.trajectory_filename,
            output_analysis_filename = OUTPUT_RMSD_PERRES_FILENAME,
            structure = self.structure,
            pbc_residues = self.pbc_residues,
            snapshots = self.snapshots,
            frames_limit = 100,
        )

    # RMSD pairwise
    def genereate_rmsd_pairwise_analysis (self):
        # Do not run the analysis if the output file already exists
        if exists(OUTPUT_RMSD_PAIRWISE_FILENAME):
            return
        print('-> Running RMSD pairwise analysis')
        # WARNING: This analysis is fast enought to use the full trajectory instead of the reduced one
        # WARNING: However, the output file size depends on the trajectory size exponentially. It may be huge
        rmsd_pairwise(
            input_topology_filename = self.structure_filename,
            input_trajectory_filename = self.trajectory_filename,
            output_analysis_filename = OUTPUT_RMSD_PAIRWISE_FILENAME,
            interactions = self.interactions,
            structure = self.structure,
            pbc_residues = self.pbc_residues,
            snapshots = self.snapshots,
            frames_limit = 200,
            overall_selection = "name CA or name C5"
        )

    # Distance per residue
    def generate_dist_perres_analysis (self):
        # Do not run the analysis if the output file already exists
        if exists(OUTPUT_DIST_PERRES_FILENAME):
            return
        print('-> Running distance per residue analysis')
        # WARNING: This analysis is not fast enought to use the full trajectory. It would take a while
        distance_per_residue(
            input_topology_filename = self.structure_filename,
            input_trajectory_filename = self.trajectory_filename,
            output_analysis_filename = OUTPUT_DIST_PERRES_FILENAME,
            interactions = self.interactions,
            snapshots = self.snapshots,
            frames_limit = 200,
        )

    # Hydrogen bonds
    def generate_hbonds_analysis (self):
        # Do not run the analysis if the output file already exists
        if exists(OUTPUT_HBONDS_FILENAME):
            return
        print('-> Running hydrogen bonds analysis')
        # WARNING: This analysis is fast enought to use the full trajectory instead of the reduced one
        # WARNING: However, the output file size depends on the trajectory
        # WARNING: Files have no limit, but analyses must be no heavier than 16Mb in BSON format
        # WARNING: In case of large surface interaction the output analysis may be larger than the limit
        # DANI: Esto no puede quedar así
        # DANI: Me sabe muy mal perder resolución con este análisis, porque en cáculo es muy rápido
        # DANI: Hay que crear un sistema de carga en mongo alternativo para análisis pesados
        hydrogen_bonds(
            input_topology_filename = self.structure_filename,
            input_trajectory_filename = self.trajectory_filename,
            output_analysis_filename = OUTPUT_HBONDS_FILENAME,
            structure = self.structure,
            interactions = self.interactions,
            snapshots = self.snapshots,
            frames_limit = 200,
            # is_time_dependend = self.is_time_dependend,
            # time_splits = 100,
            # populations = self.populations
        )

    # SASA, solvent accessible surfave analysis
    def generate_sas_analysis (self):
        # Do not run the analysis if the output file already exists
        if exists(OUTPUT_SASA_FILENAME):
            return
        print('-> Running SAS analysis')
        # Run the analysis
        sasa(
            input_topology_filename = self.structure_filename,
            input_trajectory_filename = self.trajectory_filename,
            output_analysis_filename = OUTPUT_SASA_FILENAME,
            structure = self.structure,
            pbc_residues = self.pbc_residues,
            snapshots = self.snapshots,
            frames_limit = 100,
        )

    # Energies
    def generate_energies_analysis (self):
        # Do not run the analysis if the output file already exists
        if exists(OUTPUT_ENERGIES_FILENAME):
            return
        print('-> Running energies analysis')
        # Run the analysis
        energies(
            input_topology_filename = self.structure_filename,
            input_trajectory_filename = self.trajectory_filename,
            output_analysis_filename = OUTPUT_ENERGIES_FILENAME,
            structure = self.structure,
            interactions = self.interactions,
            charges = self.charges,
            snapshots = self.snapshots,
            frames_limit = 100,
        )

    # Pockets
    def generate_pockets_analysis (self):
        # Do not run the analysis if the output file already exists
        if exists(OUTPUT_POCKETS_FILENAME):
            return
        print('-> Running pockets analysis')
        # Run the analysis
        pockets(
            input_topology_filename = self.structure_filename,
            input_trajectory_filename = self.trajectory_filename,
            output_analysis_filename = self.OUTPUT_POCKETS_FILENAME,
            structure = self.structure,
            pbc_residues = self.pbc_residues,
            snapshots = self.snapshots,
            frames_limit = 100,
        )

    # Helical parameters
    # DANI: Al final lo reimplementará Subamoy (en python) osea que esto no lo hacemos de momento
    # def generate_helical_analysis (self):
    #     # Do not run the analysis if the output file already exists
    #     if exists(OUTPUT_HELICAL_PARAMETERS_FILENAME):
    #         return
    #     print('-> Running helical analysis')
    #     # Run the analysis
    #     helical_parameters(
    #         input_topology_filename = self.structure_filename,
    #         input_trajectory_filename = self.trajectory_filename,
    #         output_analysis_filename = OUTPUT_HELICAL_PARAMETERS_FILENAME,
    #         structure = self.structure,
    #         frames_limit = 1000,
    #     )

    # Markov
    def generate_markov_analysis (self):
        # Do not run the analysis if the output file already exists
        if exists(OUTPUT_MARKOV_FILENAME):
            return
        print('-> Running Markov analysis')
        # Run the analysis
        markov(
            input_topology_filename = self.structure_filename,
            input_trajectory_filename = self.trajectory_filename,
            output_analysis_filename = OUTPUT_MARKOV_FILENAME,
            structure = self.structure,
            populations = self.populations,
            #transitions = self.transitions,
            rmsd_selection = PROTEIN_AND_NUCLEIC,
        )

    # ---------------------------------------------------------------------------------
    # List of requestables for the console
    # ---------------------------------------------------------------------------------

    # Input files
    # They may be requested to download
    input_files = {
        'istructure': get_input_structure_filename,
        'itrajectory': get_input_trajectory_filenames,
        'itopology': get_input_topology_filename,
        'inputs': get_inputs_filename,
        'populations': get_populations_filename,
        'transitions': get_transitions_filename
    }

    processed_files = {
        'structure': get_structure_filename,
        'trajectory': get_trajectory_filename,
        'topology': get_topology_filename,
    }

    # List all the available analyses
    analyses = {
        'dist': generate_dist_perres_analysis,
        'energies': generate_energies_analysis,
        'hbonds': generate_hbonds_analysis,
        #'helical': generate_helical_analysis,
        'markov': generate_markov_analysis,
        'pca': generate_pca_analysis,
        #'pcacons': generate_pca_contacts,
        'pockets': generate_pockets_analysis,
        'rgyr': generate_rgyr_analysis,
        'rmsds': generate_rmsds_analysis,
        'perres': genereate_rmsd_pairwise_analysis,
        'pairwise': genereate_rmsd_perres_analysis,
        'rmsf': generate_rmsf_analysis,
        'sas': generate_sas_analysis,
        'tmscore': generate_tmscores_analysis,
    }

    # Sort the available processes and analyses by sections
    requestables = {
        **input_files,
        **processed_files,
        **analyses,
        'interactions': get_interactions,
        'snapshots': get_snapshots,
        'charges': get_charges,
        'mapping': get_residues_map,
        'screenshot': get_screenshot_filename,
        'stopology': get_standard_topology_filename,
        'metadata': get_metadata_filename
    }

    

# AUXILIAR FUNCTIONS ---------------------------------------------------------------------------

# Set a function to read a file which may be in differen formats
# DANI: En cuanto se concrete el formato de los markov esta función no hará falta
def read_file (filename : str) -> dict:
    # Get the file format
    file_format = filename.split('.')[-1]
    # Read numpy files
    if file_format == 'npy':
        return numpy.load(filename)
    # Read JSON files
    if file_format == 'json':
        with open(filename, 'r') as file:
            return json.load(file)


# The actual main function
def workflow (
    # Project parameters
    project_parameters : dict = {},
    # The actual workflow parameters
    # Download only
    download : bool = False,
    # Download and correct only
    setup : bool = False,
    # Run only specific analyses/processes
    include : Optional[List[str]] = None,
    # Run everything but specific analyses/processes
    exclude : Optional[List[str]] = None,
):

    # Initiate the project handler
    handler = Project(**project_parameters)

    # Check input files are present and download them if they are missing and it is possible
    for getter in handler.input_files.values():
        getter(handler)

    # If download is passed as True then exit here
    # Note that there is no need to save the register if we just downloaded data
    if download:
        return

    # Process input files if needed
    for getter in handler.processed_files.values():
        getter(handler)

    # Now that we have the processed files run any additional tests here

    # Check the trajectory has not sudden jumps
    # Calling the value is enought to trigger the logic
    # The logic will warn us if something is wrong
    # DANI: Estría bien anotar en el registro o así los tests que han pasado, para no repetirlos cada vez
    handler.sudden_jumps

    # If setup is passed as True then exit here
    if setup:
        # Save the register and exit here
        handler.save_register()
        return

    # Run the requested analyses
    if include and len(include) > 0:
        sys.stdout.write(f"Executing specific dependencies: " + ', '.join(include) + '\n')
        # Include only the specified dependencies
        requested = [ getter for name, getter in handler.requestables.items() if name in include ]
        for getter in requested:
            getter(handler)
        # Save the register and exit here
        handler.save_register()
        return

    # Set the default requests, when there are not specific requests
    # Request all the analysis, the metadata, the standard topology and the screenshot
    requests = [
        *handler.analyses.keys(),
        'metadata',
        'stopology',
        'screenshot'
    ]

    # If the exclude parameter was passed then remove excluded requests from the default requests
    if exclude and len(exclude) > 0:
        sys.stdout.write(f"Excluding specific dependencies: " + ', '.join(exclude) + '\n')
        requests = [ name for name in requests if name not in exclude ]
        
    # Run the requests
    for request in requests:
        getter = handler.requestables[request]
        getter(handler)

    # Remove gromacs backups
    # DANI: Esto iría mejor en otro sitio
    remove_trash()

    # Save the register
    handler.save_register()

    sys.stdout.write("Done!\n")