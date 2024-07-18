#!/usr/bin/env python

# This is the starter script

# Import python libraries
from os import chdir, symlink, rename, walk, mkdir, getcwd
from os.path import exists, isabs
import sys
import io
import math
from pathlib import Path
import urllib.request
import json
import numpy
from glob import glob
from inspect import getfullargspec
from typing import Optional, Union, List, Callable

# Constants
from model_workflow.utils.constants import *

# Import local tools
from model_workflow.tools.topology_manager import setup_structure
from model_workflow.tools.get_pytraj_trajectory import get_pytraj_trajectory
from model_workflow.tools.get_first_frame import get_first_frame
from model_workflow.tools.get_average import get_average
from model_workflow.tools.get_bonds import get_safe_bonds, get_bonds_canonical_frame
from model_workflow.tools.process_interactions import process_interactions
from model_workflow.tools.get_pbc_residues import get_pbc_residues
from model_workflow.tools.generate_metadata import generate_project_metadata, generate_md_metadata
from model_workflow.tools.generate_ligands_desc import generate_ligand_mapping
from model_workflow.tools.residue_mapping import generate_residue_mapping
from model_workflow.tools.generate_map import generate_protein_mapping, get_sequence_metadata
from model_workflow.tools.generate_topology import generate_topology
from model_workflow.tools.get_summarized_trajectory import get_summarized_trajectory
from model_workflow.tools.get_charges import get_charges
from model_workflow.tools.remove_trash import remove_trash
from model_workflow.tools.get_screenshot import get_screenshot
from model_workflow.tools.filter_atoms import filter_atoms
from model_workflow.tools.image_and_fit import image_and_fit
from model_workflow.tools.structure_corrector import structure_corrector
from model_workflow.tools.fix_gromacs_masses import fix_gromacs_masses

# Import local utils
#from model_workflow.utils.httpsf import mount
from model_workflow.utils.auxiliar import InputError, warn, load_json, save_json, load_yaml
from model_workflow.utils.register import Register
from model_workflow.utils.conversions import convert
from model_workflow.utils.structures import Structure
from model_workflow.utils.file import File
from model_workflow.utils.pyt_spells import get_frames_count

# Import local analyses
from model_workflow.analyses.rmsds import rmsds
from model_workflow.analyses.tmscores import tmscores
from model_workflow.analyses.rmsf import rmsf
from model_workflow.analyses.rgyr import rgyr
from model_workflow.analyses.pca import pca
#from model_workflow.analyses.pca_contacts import pca_contacts
from model_workflow.analyses.rmsd_per_residue import rmsd_per_residue
from model_workflow.analyses.rmsd_pairwise import rmsd_pairwise
from model_workflow.analyses.clusters import clusters_analysis
from model_workflow.analyses.distance_per_residue import distance_per_residue
#from model_workflow.analyses.hydrogen_bonds_2 import hydrogen_bonds
from model_workflow.analyses.hydrogen_bonds import hydrogen_bonds
from model_workflow.analyses.sasa import sasa
from model_workflow.analyses.energies import energies
from model_workflow.analyses.pockets import pockets
from model_workflow.analyses.rmsd_check import check_trajectory_integrity
#from model_workflow.analyses.helical_parameters import helical_parameters
from model_workflow.analyses.markov import markov

# Make the system output stream to not be buffered
# This is useful to make prints work on time in Slurm
# Otherwise, output logs are written after the script has fully run
# Note that this fix affects all modules and built-ins
unbuffered = io.TextIOWrapper(open(sys.stdout.fileno(), 'wb', 0), write_through=True)
sys.stdout = unbuffered

# Set an special exception for input errors
missing_input_exception = Exception('Missing input')

# Run a fix in gromacs if not done before
# Note that this is run always at the moment the code is read, no matter the command or calling origin
fix_gromacs_masses()

# A Molecular Dynamics (MD) is the union of a structure and a trajectory
# Having this data several analyses are possible
# Note that an MD is always defined inside of a Project and thus it has additional topology and metadata
class MD:
    def __init__ (self,
        # The parent project this MD belongs to
        project : 'Project',
        # The number of the MD according to its accession
        number : int,
        # The local directory where the MD takes place
        directory : str,
        # Input structure and trajectory files
        input_structure_filepath : str,
        input_trajectory_filepaths : List[str],
    ):
        # Save the inputs
        self.project = project
        if not project:
            raise Exception('Project is mandatory to instantiate a new MD')
        # Save the MD number
        self.number = number
        # Set the MD accession and request URL
        self.accession = None
        self.url = None
        if self.project.database_url and self.project.accession:
            self.accession = self.project.accession + '.' + str(self.number)
            self.url = self.project.database_url + '/rest/current/projects/' + self.accession
        # Save the directory
        self.directory = remove_final_slash(directory)
        # If the directory does not exists then we may have 2 different scenarios
        if not exists(self.directory):
            # If we have an URL to donwload then it means we must download input files
            # Thus we simpy create the missing directory and it be filled further
            if self.url:
                mkdir(self.directory)
            # Otherwise we are supposed to find input files locally
            # If the directory does not exist we have nothing to do so we raise an error
            else:
                raise InputError(f'MD directory {self.directory} does not exist')
        # Input structure filepath
        # They may be relative to the project directory (unique) or relative to the MD directory (one per MD)
        # If the path is absolute then it is considered unique
        # If the file does not exist and it is to be downloaded then it is downloaded for each MD
        # Priorize the MD directory over the project directory
        project_structure = File(input_structure_filepath)
        md_structure = File(self.md_pathify(input_structure_filepath))
        if isabs(input_structure_filepath) or (project_structure.exists and not md_structure.exists):
            self._input_structure_file = project_structure
        else:
            self._input_structure_file = md_structure
        # Input trajectory filepaths
        # They are always relative
        self.input_trajectory_filepaths = [ self.md_pathify(path) for path in input_trajectory_filepaths ]
        # Now parse the filepaths to actual files
        # In case there are glob characters we must parse the paths
        self._input_trajectory_files = []
        for path in self.input_trajectory_filepaths:
            has_spread_syntax = '*' in path
            if has_spread_syntax:
                if self.url:
                    raise InputError('Spread syntax in trajectory input filepaths is not supported when downloading remote files')
                for parsed_path in glob(path):
                    trajectory_file = File(parsed_path)
                    self._input_trajectory_files.append(trajectory_file)
            else:
                trajectory_file = File(path)
                self._input_trajectory_files.append(trajectory_file)
        # Check we successfully defined some trajectory file
        if len(self._input_trajectory_files) == 0:
            raise InputError('No trajectory file was reached in paths ' + ', '.join(self.input_trajectory_filepaths))

        # Processed structure and trajectory files
        self._structure_file = None
        self._trajectory_file = None

        # Other values which may be found/calculated on demand
        self._md_inputs = None
        self._available_files = None
        self._snapshots = None
        self._structure = None
        self._pytraj_topology = None
        self._interactions = None
        self._reference_frame = None

        # Tests
        self._trajectory_integrity = None

        # Set a new MD specific register
        # In case the directory is the project directory itself, use the project register
        register_file = File(self.md_pathify(REGISTER_FILENAME))
        if register_file.path == self.project.register.file.path:
            self.register = self.project.register
        else:
            self.register = Register(register_file)

    def __repr__ (self):
        return '<MD (' + str(len(self.structure.atoms)) + ' atoms)>'

    # Given a filename or relative path, add the MD directory path at the beginning
    def md_pathify (self, filename_or_relative_path : str) -> str:
        return self.directory + '/' + filename_or_relative_path

    # Input structure filename ------------

    # Get the input pdb filename from the inputs
    # If the file is not found try to download it
    def get_input_structure_file (self) -> str:
        # There must be an input structure filename
        if not self._input_structure_file:
            raise InputError('Not defined input structure filename')
        # If the file already exists then we are done
        if self._input_structure_file.exists:
            return self._input_structure_file
        # Try to download it
        # If we do not have the required parameters to download it then we surrender here
        if not self.url:
            raise InputError('Missing input structure file "' + self._input_structure_file.filename + '"')
        sys.stdout.write('Downloading structure (' + self._input_structure_file.filename + ')\n')
        # Set the download URL
        structure_url = self.url + '/files/' + self._input_structure_file.filename
        # If the structure filename is the standard structure filename then use the structure endpoint instead
        if self._input_structure_file.filename == STRUCTURE_FILENAME:
            structure_url = self.url + '/structure'
        # Download the file
        try:
            urllib.request.urlretrieve(structure_url, self._input_structure_file.path)
        except urllib.error.HTTPError as error:
            if error.code == 404:
                raise Exception('Missing input pdb file "' + self._input_structure_file.filename + '"')
            else:
                raise Exception('Something went wrong with the MDposit request: ' + structure_url)
        return self._input_structure_file
    input_structure_file = property(get_input_structure_file, None, None, "Input structure filename (read only)")

    # Input trajectory filename ------------

    # Get the input trajectory filename(s) from the inputs
    # If file(s) are not found try to download it
    def get_input_trajectory_files (self) -> str:
        # There must be an input structure filename
        if not self._input_trajectory_files or len(self._input_trajectory_files) == 0:
            print(self._input_trajectory_files)
            raise InputError('Not defined input trajectory filenames')
        # Find missing trajectory files
        missing_input_trajectory_files = []
        for trajectory_file in self._input_trajectory_files:
            if not trajectory_file.exists:
                missing_input_trajectory_files.append(trajectory_file)
        # If all files already exists then we are done
        if len(missing_input_trajectory_files) == 0:
            return self._input_trajectory_files
        # Try to download the missing files
        # If we do not have the required parameters to download it then we surrender here
        if not self.url:
            missing_filepaths = [ trajectory_file.path for trajectory_file in missing_input_trajectory_files ]
            raise InputError('Missing input trajectory files: ' + ', '.join(missing_filepaths))
        # Download each trajectory file (ususally it will be just one)
        for trajectory_file in self._input_trajectory_files:
            sys.stdout.write('Downloading trajectory (' + trajectory_file.filename + ')\n')
            # Set the trajectory URL, which may depend in different parameters
            trajectory_url = self.url + '/files/' + trajectory_file.filename
            if trajectory_file.filename == TRAJECTORY_FILENAME:
                trajectory_url = self.url + '/trajectory?format=xtc'
                if self.project.sample_trajectory:
                    trajectory_url += '&frames=1:10:1'
            # Download the file
            try:
                urllib.request.urlretrieve(trajectory_url, trajectory_file.path)
            except urllib.error.HTTPError as error:
                if error.code == 404:
                    raise Exception('Missing input trajectory file "' + trajectory_file.filename + '"')
                raise Exception('Something went wrong with the MDposit request: ' + trajectory_url)
        return self._input_trajectory_files
    input_trajectory_files = property(get_input_trajectory_files, None, None, "Input trajectory filenames (read only)")

    # MD specific inputs
    def get_md_inputs (self) -> dict:
        # If we already have a value stored then return it
        if self._md_inputs:
            return self._md_inputs
        # Otherwise we must find its value
        if self.project.input_mds:
            for md in self.project.input_mds:
                name = md['name']
                directory = name_2_directory(name)
                if directory == self.directory:
                    self._md_inputs = md
                    return self._md_inputs
        # If this MD directory has not associated inputs then it means it was forced through command line
        # We set a provisional MD inputs for it
        provisional_name = directory_2_name(self.directory)
        self._md_inputs = { 'name': provisional_name }
        return self._md_inputs

    md_inputs = property(get_md_inputs, None, None, "MD specific inputs (read only)")

    # ---------------------------------

    # Remote available files
    def get_available_files (self) -> List[str]:
        # If we already have a stored value then return it
        if self._available_files != None:
            return self._available_files
        # If we have not remote access then there are not available files
        if not self.url:
            return []
        # Get the available files
        files_url = self.url + '/files'
        response = urllib.request.urlopen(files_url)
        filenames = json.loads(response.read())
        self._available_files = filenames
        return self._available_files
    available_files = property(get_available_files, None, None, "Remote available files (read only)")

    # Check if a file exists
    # If not, try to download it from the database
    # If the file is not found in the database it is fine, we do not even warn the user
    # Note that this function is used to get populations and transitions files, which are not common
    def get_file (self, target_file : File) -> bool:
        # If it exists we are done
        if target_file.exists:
            return True
        # Try to download the missing file
        # If we do not have the required parameters to download it then we surrender here
        if not self.url:
            return False
        # Check if the file is among the available remote files
        # If it is no then stop here
        if target_file.filename not in self.available_files:
            return False
        sys.stdout.write('Downloading file (' + target_file.filename + ')\n')
        # Set the download URL
        file_url = self.url + '/files/' + target_file.filename
        # Download the file
        try:
            urllib.request.urlretrieve(file_url, target_file.path)
        except urllib.error.HTTPError as error:
            if error.code == 404:
                print('Missing file "' + target_file.filename + '"')
                return False
            else:
                raise Exception('Something went wrong with the MDposit request: ' + file_url)
        return True

    # Processed files ----------------------------------------------------      

    # Process input files to generate the processed files
    # This process corrects and standarizes the topology, the trajectory and the structure
    def process_input_files (self):

        # Set the input filepaths
        input_structure_file = self.input_structure_file
        input_trajectory_files = self.input_trajectory_files
        input_topology_file = self.project.input_topology_file

        # Set the output filepaths
        output_structure_file = File(self.md_pathify(STRUCTURE_FILENAME))
        output_trajectory_file = File(self.md_pathify(TRAJECTORY_FILENAME))
        output_topology_file = File(self.project.topology_filepath)

        # If all output files already exist we may skip the processing
        outputs_exist = output_structure_file.exists and output_trajectory_file.exists and output_topology_file.exists

        # Check which tests are to be run
        required_tests = set()

        # If there is no structure then we must run some tests
        if not output_structure_file.exists:
            required_tests.update(STRUCTURE_TESTS)
        # If the structure was modified since the last time then we must run these tests as well
        elif self.register.is_file_modified(output_structure_file):
            message = 'Structure was modified since the last processing or is new'
            warn(message)
            required_tests.update(STRUCTURE_TESTS)

        # If there is no trajectory then we must run some tests
        if not output_trajectory_file.exists:
            required_tests.update(TRAJECTORY_TESTS)
        # If the trajectory was modified since the last time then we must run these tests as well
        elif self.register.is_file_modified(output_trajectory_file):
            message = 'Trajectory was modified since the last processing or is new'
            warn(message)
            required_tests.update(TRAJECTORY_TESTS)

        # If there is no topology then we must run some tests
        if not output_topology_file.exists:
            required_tests.update(TOPOLOGY_TESTS)
        # If the topology was modified since the last time then we must run these tests as well
        elif self.project.register.is_file_modified(output_topology_file):
            message = 'Topology was modified since the last processing or is new'
            warn(message)
            required_tests.update(TOPOLOGY_TESTS)

        # If any of the required tests was already passed then reset its value and warn the user
        repeated_tests = [ test for test in required_tests if self.register.tests.get(test, None) == True ]
        if len(repeated_tests) > 0:
            print('  The following tests will be run again: ' + ', '.join(repeated_tests))
            for test in repeated_tests:
                self.register.update_test(test, None)

        # Check if the processing parameters (filter, image, etc.) have changed since the last time
        # If so, then we must reset all tests and rerun the processing
        previous_processed_parameters = self.register.cache.get(PROCESSED, None)
        current_processed_parameters = {
            'filter': self.project.filter_selection,
            'image': self.project.image,
            'fit': self.project.fit,
        }
        # Compare current and previous values parameter by parameters
        same_processed_paramaters = previous_processed_parameters == current_processed_parameters
        if previous_processed_parameters and not same_processed_paramaters:
            # Warn the user that there has been a change in input parameters
            message = 'There is a change in the processing parameters'
            print(YELLOW_HEADER + 'WARNING: ' + COLOR_END + message)
            print(' Processed files will be remade')

        # If output files already exist and not test is to be run then we skip the processing
        # Check also if all availables tests were actually passed in the last run
        # They may be skipped or allowed to fail
        # Also make sure processing parameters are the same that the last time
        if outputs_exist and len(required_tests) == 0 and self.all_tests_succeeded() and same_processed_paramaters:
            return

        print('-> Processing input files')

        # --- CONVERTING AND MERGING ------------------------------------------------------------

        # Set the output format for the already converted structure
        input_structure_format = self.input_structure_file.format
        output_structure_format = output_structure_file.format
        converted_structure_filepath = self.md_pathify(CONVERTED_STRUCTURE)
        # If input structure already matches the output format then avoid the renaming
        if input_structure_format == output_structure_format:
            converted_structure_filepath = input_structure_file.path
        # Set the output file for the already converted structure
        converted_structure_file = File(converted_structure_filepath)
        # Input trajectories should have all the same format
        input_trajectory_formats = set([ trajectory_file.format for trajectory_file in input_trajectory_files ])
        if len(input_trajectory_formats) > 1:
            raise InputError('All input trajectory files must have the same format')
        # Set the output format for the already converted trajectory
        input_trajectories_format = list(input_trajectory_formats)[0]
        output_trajectory_format = output_trajectory_file.format
        # Set the output file for the already converted trajectory
        converted_trajectory_filepath = self.md_pathify(CONVERTED_TRAJECTORY)
        # If input trajectory already matches the output format and is unique then avoid the renaming
        if input_trajectories_format == output_trajectory_format and len(input_trajectory_files) == 1:
            converted_trajectory_filepath = input_trajectory_files[0].path
        converted_trajectory_file = File(converted_trajectory_filepath)
        # Join all input trajectory paths
        input_trajectory_paths = [ trajectory_file.path for trajectory_file in input_trajectory_files ]

        # Set an intermeidate file for the trajectory while it is being converted
        # This prevents using an incomplete trajectory in case the workflow is suddenly interrupted while converting
        incompleted_converted_trajectory_filepath = self.md_pathify(INCOMPLETE_PREFIX + CONVERTED_TRAJECTORY)
        incompleted_converted_trajectory_file = File(incompleted_converted_trajectory_filepath)
        # If there is an incomplete trajectory then remove it
        if incompleted_converted_trajectory_file.exists:
            incompleted_converted_trajectory_file.remove()

        # Convert input structure and trajectories to output structure and trajectory
        if not converted_structure_file.exists or not converted_trajectory_file.exists:
            print(' * Converting and merging')
            convert(
                input_structure_filepath = input_structure_file.path,
                output_structure_filepath = converted_structure_file.path,
                input_trajectory_filepaths = input_trajectory_paths,
                output_trajectory_filepath = incompleted_converted_trajectory_file.path,
            )
            # Once converted, rename the trajectory file as completed
            rename(incompleted_converted_trajectory_file.path, converted_trajectory_file.path)

        # Topologies are never converted, but they are kept in their original format

        # --- FILTERING ATOMS ------------------------------------------------------------

        # Find out if we need to filter
        # i.e. check if there is a selection filter and it matches some atoms
        must_filter = bool(self.project.filter_selection)

        # Set output filenames for the already filtered structure and trajectory
        # Note that this is the only step affecting topology and thus here we output the definitive topology
        filtered_structure_file = File(self.md_pathify(FILTERED_STRUCTURE)) if must_filter else converted_structure_file
        filtered_trajectory_file = File(self.md_pathify(FILTERED_TRAJECTORY)) if must_filter else converted_trajectory_file
        filtered_topology_file = output_topology_file if must_filter else input_topology_file

        # Set an intermeidate file for the trajectory while it is being filtered
        # This prevents using an incomplete trajectory in case the workflow is suddenly interrupted while filtering
        incompleted_filtered_trajectory_filepath = self.md_pathify(INCOMPLETE_PREFIX + FILTERED_TRAJECTORY)
        incompleted_filtered_trajectory_file = File(incompleted_filtered_trajectory_filepath)
        # If there is an incomplete trajectory then remove it
        if incompleted_filtered_trajectory_file.exists:
            incompleted_filtered_trajectory_file.remove()

        # Check if any output file is missing
        missing_filter_output = not filtered_structure_file.exists or not filtered_trajectory_file.exists

        # Check if parameters have changed
        # Note that for this specific step only filtering is important
        previous_filtered_parameters = self.register.cache.get(FILTERED, None)
        current_filtered_parameters = { 'filter': self.project.filter_selection }
        same_filtered_parameters = previous_filtered_parameters == current_filtered_parameters
        
        # Filter atoms in structure, trajectory and topology if required and not done yet
        if must_filter and (missing_filter_output or not same_filtered_parameters):
            print(' * Filtering atoms')
            filter_atoms(
                input_structure_file = converted_structure_file,
                input_trajectory_file = converted_trajectory_file,
                input_topology_file = input_topology_file, # We use input topology
                output_structure_file = filtered_structure_file,
                output_trajectory_file = incompleted_filtered_trajectory_file,
                output_topology_file = filtered_topology_file, # We genereate the definitive topology
                filter_selection = self.project.filter_selection
            )
            # Once filetered, rename the trajectory file as completed
            rename(incompleted_filtered_trajectory_file.path, filtered_trajectory_file.path)
            # Update the cache
            self.register.update_cache(FILTERED, current_filtered_parameters)

        # --- IMAGING AND FITTING ------------------------------------------------------------

        # There is no logical way to know if the trajectory is already imaged or it must be imaged
        # We rely exclusively in input flags
        must_image = self.project.image or self.project.fit

        # Set output filenames for the already filtered structure and trajectory
        imaged_structure_file = File(self.md_pathify(IMAGED_STRUCTURE)) if must_image else filtered_structure_file
        imaged_trajectory_file = File(self.md_pathify(IMAGED_TRAJECTORY)) if must_image else filtered_trajectory_file

        # Set an intermeidate file for the trajectory while it is being imaged
        # This prevents using an incomplete trajectory in case the workflow is suddenly interrupted while imaging
        incompleted_imaged_trajectory_filepath = self.md_pathify(INCOMPLETE_PREFIX + IMAGED_TRAJECTORY)
        incompleted_imaged_trajectory_file = File(incompleted_imaged_trajectory_filepath)
        # If there is an incomplete trajectory then remove it
        if incompleted_imaged_trajectory_file.exists:
            incompleted_imaged_trajectory_file.remove()

        # Check if any output file is missing
        missing_imaged_output = not imaged_structure_file.exists or not imaged_trajectory_file.exists

        # Check if parameters have changed
        # Note that for this step the filter parameters is also important
        previous_imaged_parameters = self.register.cache.get(IMAGED, None)
        current_imaged_parameters = {
            'filter': self.project.filter_selection,
            'image': self.project.image,
            'fit': self.project.fit,
        }
        same_imaged_parameters = previous_imaged_parameters == current_imaged_parameters

        # Image the trajectory if it is required
        # i.e. make the trajectory uniform avoiding atom jumps and making molecules to stay whole
        # Fit the trajectory by removing the translation and rotation if it is required
        if must_image and (missing_imaged_output or not same_imaged_parameters):
            print(' * Imaging and fitting')
            image_and_fit(
                input_structure_file = filtered_structure_file,
                input_trajectory_file = filtered_trajectory_file,
                input_topology_file = filtered_topology_file, # This is optional if there are no PBC residues
                output_structure_file = imaged_structure_file,
                output_trajectory_file = incompleted_imaged_trajectory_file,
                image = self.project.image,
                fit = self.project.fit,
                translation = self.project.translation,
                pbc_selection = self.project.pbc_selection
            )
            # Once imaged, rename the trajectory file as completed
            rename(incompleted_imaged_trajectory_file.path, imaged_trajectory_file.path)
            # Update the cache
            self.register.update_cache(IMAGED, current_imaged_parameters)

        # --- CORRECTING STRUCTURE ------------------------------------------------------------

        # Note that this step, although it is foucsed in the structure, requires also the trajectory
        # Also the trajectory may be altered in very rare cases where coordinates must be resorted

        # There is no possible reason to not correct the structure
        # This is the last step so the output files will be named as the output files of the whole processing

        # WARNING:
        # For the correcting function we need the number of snapshots and at this point it should not be defined
        # Snapshots are calculated by default from the already processed structure and trajectory
        # For this reason we can not rely on the public snapshots getter
        # We must calculate snapshots here using last step structure and trajectory
        # If we already have a value in the register cache then use it
        cached_snapshots = self.register.cache.get(SNAPSHOTS_FLAG, None)
        if cached_snapshots != None:
            self._snapshots = cached_snapshots
        # Othwerise count the number of snaphsots
        else:
            self._snapshots = get_frames_count(imaged_structure_file, imaged_trajectory_file)
            # Save the snapshots value in the register cache as well
            self.register.update_cache(SNAPSHOTS_FLAG, self._snapshots)

        # WARNING:
        # We may need to resort atoms in the structure corrector function
        # In such case, bonds and charges must be resorted as well and saved apart to keep values coherent
        # Bonds are calculated during the structure corrector but atom charges must be extracted no
        self.project._charges = get_charges(filtered_topology_file)

        print(' * Correcting structure')

        # Correct the structure
        # This function reads and or modifies the following MD variables:
        #   snapshots, safe_bonds, register, mercy, trust
        structure_corrector(
            input_structure_file = imaged_structure_file,
            input_trajectory_file = imaged_trajectory_file,
            input_topology_file = filtered_topology_file,
            output_structure_file = output_structure_file,
            output_trajectory_file = output_trajectory_file,
            MD = self
        )

        # Now we must rename files to match the output file in case there is any missmatch
        # Some processed files may remain with some intermediate filename
        input_and_output_files = [
            (input_structure_file, imaged_structure_file, output_structure_file),
            (input_trajectory_files[0], imaged_trajectory_file, output_trajectory_file),
            (input_topology_file, filtered_topology_file, output_topology_file)
        ]
        for input_file, processed_file, output_file in input_and_output_files:
            if output_file.exists:
                continue
            # There is also a chance that the input files have not been modified
            # This means the input format has already the output format and it is not to be imaged, fitted or corrected
            # However we need the output files to exist and we dont want to rename the original ones to conserve them
            # In order to not duplicate data, we will setup a symbolic link to the input files with the output filepaths
            if processed_file == input_file:
                output_file.set_symlink_to(input_file)
            # Some processed files may remain with some intermediate filename
            else:
                rename(processed_file.path, output_file.path)

        # Save the internal variables
        self._structure_file = output_structure_file
        self._trajectory_file = output_trajectory_file
        self.project._topology_file = output_topology_file

        # Register the last modification times of the recently processed files
        # This way we know if they have been modifed in the future and checkings need to be rerun
        self.register.update_mtime(output_structure_file)
        self.register.update_mtime(output_trajectory_file)
        self.project.register.update_mtime(output_topology_file)

        # Update the parameters used to get the last processed structure and trajectory files
        self.register.update_cache(PROCESSED, current_processed_parameters)

        # --- Cleanup intermediate files

        # Set a list of intermediate files
        intermediate_files = set([
            converted_structure_file, converted_trajectory_file,
            filtered_structure_file, filtered_trajectory_file,
            imaged_structure_file, imaged_trajectory_file,
        ])
        # Set also a list of input files
        inputs_files = set([ input_structure_file, *input_trajectory_files, input_topology_file ])
        # We must make sure an intermediate file is not actually an input file before deleting it
        removable_files = intermediate_files - inputs_files
        # Now delete every removable file
        for removable_file in removable_files:
            # Note that a broken symlink does not 'exists'
            if removable_file.exists or removable_file.is_symlink():
                removable_file.remove()

        # --- RUNNING FINAL TESTS ------------------------------------------------------------

        # Note that some tests have been run already
        # e.g. stable bonds is run in the structure corrector function

        # Note that tests here do not modify any file

        # Check the trajectory has not sudden jumps
        self.is_trajectory_integral()

        # Make a final summary
        print('Tests summary:')
        for test_name in AVAILABLE_CHECKINGS:
            test_result = self.register.tests.get(test_name, None)
            # Print things pretty
            test_nice_name = NICE_NAMES[test_name]
            test_nice_result = None
            if test_result == None:
                test_nice_result = YELLOW_HEADER + 'Not run' + COLOR_END
            elif test_result == False:
                test_nice_result = RED_HEADER + 'Failed' + COLOR_END
            elif test_result == True:
                test_nice_result = GREEN_HEADER + 'Passed' + COLOR_END
            else:
                raise ValueError()
            
            print(' - ' + test_nice_name + ' -> ' + test_nice_result)

        # Issue some warnings if failed or never run tests are skipped
        for test_name in AVAILABLE_CHECKINGS:
            # If test was not skipped then proceed
            if test_name not in self.project.trust:
                continue
            # If test passed in a previous run the proceed
            test_result = self.register.tests.get(test_name)
            if test_result == True:
                continue
            # If test failed in a previous run we can also proceed
            # The failing warning must be among the inherited warnings, so there is no need to add more warnings here
            elif test_result == False:
                continue
            # If the test has been always skipped then issue a warning
            elif test_result == None:
                # Set the test skip flag and remove previous warnings
                test_skip_flag = 'skip_' + test_name
                self.register.remove_warnings(test_skip_flag)
                # Get test pretty name
                test_nice_name = NICE_NAMES[test_name]
                # Issue the corresponding warning            
                self.register.add_warning(test_skip_flag, test_nice_name + ' was skipped and never run before')
            else:
                raise ValueError('Test value is not supported')

    # Check all tests to be succeeded in the last run
    def all_tests_succeeded (self) -> bool:
        for checking in AVAILABLE_CHECKINGS:
            test_result = self.register.tests.get(checking, None)
            if test_result != True:
                return False
        return True

    # Get the processed structure
    def get_structure_file (self) -> str:
        # If we have a stored value then return it
        # This means we already found or generated this file
        if self._structure_file:
            return self._structure_file
        # Set the file
        structure_filepath = self.md_pathify(STRUCTURE_FILENAME)
        self._structure_file = File(structure_filepath)
        # Run the processing logic
        self.process_input_files()
        # Now that the file is sure to exist we return it
        return self._structure_file
    structure_file = property(get_structure_file, None, None, "Structure filename (read only)")

    # Get the processed trajectory
    def get_trajectory_file (self) -> str:
        # If we have a stored value then return it
        # This means we already found or generated this file
        if self._trajectory_file:
            return self._trajectory_file
        # If the file already exists then we are done
        trajectory_filepath = self.md_pathify(TRAJECTORY_FILENAME)
        self._trajectory_file = File(trajectory_filepath)
        # Run the processing logic
        self.process_input_files()
        # Now that the file is sure to exist we return it
        return self._trajectory_file
    trajectory_file = property(get_trajectory_file, None, None, "Trajectory filename (read only)")

    # Get the processed topology from the project
    def get_topology_file (self) -> str:
        return self.project.topology_file
    topology_file = property(get_topology_file, None, None, "Topology filename from the project (read only)")

    # ---------------------------------------------------------------------------------
    # Others values which may be found/calculated and files to be generated on demand
    # ---------------------------------------------------------------------------------

    # Trajectory snapshots
    def get_snapshots (self) -> str:
        # If we already have a stored value then return it
        # WARNING: Do not remove the self.trajectory_file checking
        # Note that checking if the trajectory file exists triggers all the processing logic
        # The processing logic is able to set the internal snapshots value as well so this avoid repeating the process
        if self.trajectory_file and self._snapshots != None:
            return self._snapshots
        # If we already have a value in the register cache then use it
        cached_value = self.register.cache.get(SNAPSHOTS_FLAG, None)
        if cached_value != None:
            return cached_value
        print('-> Counting snapshots')
        # Otherwise we must find the value
        # This happens when the input files are already porcessed and thus we did not yet count the frames
        self._snapshots = get_frames_count(self.structure_file, self.trajectory_file)
        # Save the snapshots value in the register cache as well
        self.register.update_cache(SNAPSHOTS_FLAG, self._snapshots)
        return self._snapshots
    snapshots = property(get_snapshots, None, None, "Trajectory snapshots (read only)")

    # Safe bonds
    def get_safe_bonds (self) -> List[List[int]]:
        return self.project.safe_bonds
    safe_bonds = property(get_safe_bonds, None, None, "Atom bonds to be trusted (read only)")

    # Parsed structure
    def get_structure (self) -> 'Structure':
        # If we already have a stored value then return it
        if self._structure:
            return self._structure
        # Otherwise we must set the structure
        # Note that this is not only the structure class, but it also contains additional logic
        self._structure = setup_structure(self.structure_file.path)
        # If the stable bonds test failed and we had mercy then it is sure our structure will have wrong bonds
        # In order to make it coherent with the topology we will mine topology bonds from here and force them in the structure
        # If we fail to get bonds from topology then just go along with the default structure bonds
        if self.register.tests.get(STABLE_BONDS_FLAG, None) == False:
            self._structure.bonds = self.safe_bonds
        return self._structure
    structure = property(get_structure, None, None, "Parsed structure (read only)")

    # Pytraj trajectory
    def get_pytraj_trajectory (self) -> 'TrajectoryIterator':
        # If we already have a stored value then return it
        if self._pytraj_topology:
            return self._pytraj_topology
        # Otherwise we must set the pytarj trajectory
        self._pytraj_topology = get_pytraj_trajectory(
            input_topology_filename = self.structure_file.path,
            input_trajectory_filename = self.trajectory_file.path
        )
        return self._pytraj_topology
    pytraj_trajectory = property(get_pytraj_trajectory, None, None, "Pytraj trajectory (read only)")

    # First frame filename
    def get_first_frame_file (self) -> str:
        # If the file already exists then send it
        first_frame_filepath = self.md_pathify(FIRST_FRAME_FILENAME)
        first_frame_file = File(first_frame_filepath)
        if first_frame_file.exists:
            return first_frame_file
        # Otherwise, generate it
        get_first_frame(
            input_structure_filename = self.structure_file.path,
            input_trajectory_filename = self.trajectory_file.path,
            first_frame_filename = first_frame_file.path
        )
        return first_frame_file
    first_frame_file = property(get_first_frame_file, None, None, "First frame (read only)")

    # Average structure filename
    def get_average_structure_file (self) -> str:
        # If the file already exists then send it
        average_structure_filepath = self.md_pathify(AVERAGE_STRUCTURE_FILENAME)
        average_structure_file = File(average_structure_filepath)
        if average_structure_file.exists:
            return average_structure_file
        # Otherwise, generate it
        get_average(
            pytraj_trajectory = self.pytraj_trajectory,
            output_average_filename = average_structure_file.path
        )
        return average_structure_file
    average_structure_file = property(get_average_structure_file, None, None, "Average structure filename (read only)")

    # MD metadata filename
    def get_metadata_file (self, overwrite : bool = False) -> File:
        # Set the metadata file
        metadata_filepath = self.md_pathify(OUTPUT_METADATA_FILENAME)
        metadata_file = File(metadata_filepath)
        # If the file already exists then send it
        if metadata_file.exists and not overwrite:
            return metadata_file
        print('-> Generating MD metadata')
        # Otherwise, generate it
        generate_md_metadata(
            md_inputs = self.md_inputs, # DANI: No sería mejor pasarle los inputs?
            structure = self.structure,
            snapshots = self.snapshots,
            reference_frame = self.reference_frame,
            register = self.register,
            output_metadata_filename = metadata_file.path,
        )
        return metadata_file
    metadata_file = property(get_metadata_file, None, None, "Project metadata filename (read only)")

    # The interactions
    # This is a bit exceptional since it is a value to be used and an analysis file to be generated
    def get_interactions (self, overwrite : bool = False) -> List[dict]:
        # If we already have a stored value then return it
        if self._interactions != None:
            return self._interactions
        # Set the interactions file
        interactions_filepath = self.md_pathify(OUTPUT_INTERACTIONS_FILENAME)
        interactions_file = File(interactions_filepath)
        # If the file already exists then interactions will be read from it
        # If the overwrite argument is passed we must delete it here
        if interactions_file.exists and overwrite:
            interactions_file.remove()
        print('-> Processing interactions')
        # Otherwise, process interactions
        self._interactions = process_interactions(
            input_interactions = self.project.input_interactions,
            structure_file = self.structure_file,
            trajectory_file = self.trajectory_file,
            structure = self.structure,
            snapshots = self.snapshots,
            interactions_file = interactions_file,
            mercy = self.project.mercy,
            frames_limit = 1000,
            interaction_cutoff = self.project.interaction_cutoff
        )
        return self._interactions
    interactions = property(get_interactions, None, None, "Interactions (read only)")

    # Indices of residues in periodic boundary conditions
    # WARNING: Do not inherit project pbc residues
    # WARNING: It may trigger all the processing logic of the reference MD when there is no need
    def get_pbc_residues (self) -> List[int]:
        # If we already have a stored value then return it
        if self.project._pbc_residues:
            return self.project._pbc_residues
        # Otherwise we must find the value
        self.project._pbc_residues = get_pbc_residues(self.structure, self.project.pbc_selection)
        return self.project._pbc_residues
    pbc_residues = property(get_pbc_residues, None, None, "Indices of residues in periodic boundary conditions (read only)")

    # Atom charges
    # Inherited from project
    def get_charges (self) -> List[float]:
        return self.project.charges
    charges = property(get_charges, None, None, "Atom charges (read only)")

    # Equilibrium populations from a MSM
    # Inherited from project
    def get_populations (self) -> List[float]:
        return self.project.populations
    populations = property(get_populations, None, None, "Equilibrium populations from a MSM (read only)")

    # Transition probabilities from a MSM
    # Inherited from project
    def get_transitions (self) -> List[List[float]]:
        return self.project.transitions
    transitions = property(get_transitions, None, None, "Transition probabilities from a MSM (read only)")

    # Residues mapping
    # Inherited from project
    def get_protein_map (self) -> dict:
        return self.project.protein_map
    protein_map = property(get_protein_map, None, None, "Residues mapping (read only)")

    # Reference frame
    # Frame to be used when representing the MD
    def get_reference_frame (self) -> dict:
        # If we already have a stored value then return it
        # Note that this value is usually assigned at the structure_corrector
        if self._reference_frame:
            return self._reference_frame
        # Otherwise we must find the value
        # Get some input values
        structure_filepath = self.structure_file.path
        trajectory_filepath = self.trajectory_file.path
        # If the reference frame was not found because input files were not yet processed then now it should be available
        if self._reference_frame:
            return self._reference_frame
        # Otherwise it means we have input files which are not to be processed
        # So we must calculate from here the reference frame
        # Get the rest of inputs
        atom_elements = [ atom.element for atom in self.structure.atoms ]
        # Find the first frame in the whole trajectory where safe bonds are respected
        self._reference_frame = get_bonds_canonical_frame(
            structure_filepath = structure_filepath,
            trajectory_filepath = trajectory_filepath,
            snapshots = self.snapshots,
            reference_bonds = self.safe_bonds,
            atom_elements = atom_elements
        )
        return self._reference_frame
    reference_frame = property(get_reference_frame, None, None, "Reference frame to be used to represent the MD (read only)")

    # ---------------------------------------------------------------------------------
    # Tests
    # ---------------------------------------------------------------------------------

    # Sudden jumps test
    def is_trajectory_integral (self) -> Optional[bool]:
        # If we already have a stored value then return it
        if self._trajectory_integrity != None:
            return self._trajectory_integrity
        # Otherwise we must find the value
        self._trajectory_integrity = check_trajectory_integrity(
            input_structure_filename = self.structure_file.path,
            input_trajectory_filename = self.trajectory_file.path,
            structure = self.structure,
            pbc_residues = self.pbc_residues,
            mercy = self.project.mercy,
            trust = self.project.trust,
            register = self.register,
            # time_length = self.time_length,
            check_selection = ALL_ATOMS,
            standard_deviations_cutoff = self.project.rmsd_cutoff,
        )
        return self._trajectory_integrity

    # ---------------------------------------------------------------------------------
    # Analyses
    # ---------------------------------------------------------------------------------

    # RMSDs
    def run_rmsds_analysis (self, overwrite : bool = False):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_RMSDS_FILENAME)
        if exists(output_analysis_filepath) and not overwrite:
            return
        print('-> Running RMSDs analysis')
        # WARNING: This analysis is fast enought to use the full trajectory
        # WARNING: However, the output file size depends on the trajectory size
        # WARNING: In very long trajectories the number of points may make the client go slow when loading data
        rmsds(
            trajectory_file = self.trajectory_file,
            first_frame_file = self.first_frame_file,
            average_structure_file = self.average_structure_file,
            output_analysis_filepath = output_analysis_filepath,
            frames_limit = 5000,
            snapshots = self.snapshots,
            structure = self.structure,
            pbc_residues = self.pbc_residues,
            ligand_map = self.project.ligand_map,
        )

    # TM scores
    def run_tmscores_analysis (self, overwrite : bool = False):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_TMSCORES_FILENAME)
        if exists(output_analysis_filepath) and not overwrite:
            return
        print('-> Running TM scores analysis')
        # Here we set a small frames limit since this anlaysis is a bit slow
        tmscores(
            input_trajectory_file = self.trajectory_file,
            output_analysis_filename = output_analysis_filepath,
            first_frame_file = self.first_frame_file,
            average_structure_file = self.average_structure_file,
            structure = self.structure,
            pbc_residues = self.pbc_residues,
            snapshots = self.snapshots,
            frames_limit = 200,
        )

    # RMSF, atom fluctuation
    def run_rmsf_analysis (self, overwrite : bool = False):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_RMSF_FILENAME)
        if exists(output_analysis_filepath) and not overwrite:
            return
        print('-> Running RMSF analysis')
        # This analysis is fast and the output size depends on the number of atoms only
        # For this reason here it is used the whole trajectory with no frames limit
        rmsf(
            input_topology_filename = self.structure_file.path,
            input_trajectory_filename = self.trajectory_file.path,
            output_analysis_filename = output_analysis_filepath,
            structure = self.structure,
            pbc_residues = self.pbc_residues,
        )

    # RGYR, radius of gyration
    def run_rgyr_analysis (self, overwrite : bool = False):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_RGYR_FILENAME)
        if exists(output_analysis_filepath):
            return
        print('-> Running RGYR analysis')
        # WARNING: This analysis is fast enought to use the full trajectory instead of the reduced one
        # WARNING: However, the output file size depends on the trajectory size
        # WARNING: In very long trajectories the number of points may make the client go slow when loading data
        rgyr(
            input_topology_file = self.structure_file,
            input_trajectory_file = self.trajectory_file,
            output_analysis_filepath = output_analysis_filepath,
            snapshots = self.snapshots,
            frames_limit = 5000,
            structure = self.structure,
            pbc_residues = self.pbc_residues,
        )

    # PCA, principal component analysis
    def run_pca_analysis (self, overwrite : bool = False):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_PCA_FILENAME)
        if exists(output_analysis_filepath) and not overwrite:
            return
        print('-> Running PCA analysis')
        # WARNING: This analysis will generate several output files
        # File 'pca.average.pdb' is generated by the PCA and it was used by the client but not anymore
        # File 'covar.log' is generated by the PCA but never used
        pca(
            input_topology_file = self.structure_file,
            input_trajectory_file = self.trajectory_file,
            output_analysis_filepath = output_analysis_filepath,
            output_trajectory_projections_prefix = self.md_pathify(OUTPUT_PCA_PROJECTION_PREFIX),
            snapshots = self.snapshots,
            frames_limit = 2000,
            structure = self.structure,
            fit_selection = self.project.pca_fit_selection,
            analysis_selection = self.project.pca_selection,
            pbc_residues = self.pbc_residues,
        )

    # PCA contacts
    # DANI: Intenta usar mucha memoria, hay que revisar
    # DANI: Puede saltar un error de imposible alojar tanta memoria
    # DANI: Puede comerse toda la ram y que al final salte un error de 'Terminado (killed)'
    # DANI: De momento me lo salto
    # def run_pca_contacts (self, overwrite : bool = False):
    #     # Do not run the analysis if the output file already exists
    #     output_analysis_filepath = self.md_pathify(OUTPUT_PCA_CONTACTS_FILENAME)
    #     if exists(output_analysis_filepath) and not overwrite:
    #         return
    #     print('-> Running PCA contacts analysis')
    #     pca_contacts(
    #         trajectory = self.trajectory_file.path,
    #         topology = self.pdb_filename,
    #         interactions = self.interactions,
    #         output_analysis_filename = output_analysis_filepath
    #     )

    # RMSD per residue
    def run_rmsd_perres_analysis (self, overwrite : bool = False):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_RMSD_PERRES_FILENAME)
        if exists(output_analysis_filepath) and not overwrite:
            return
        print('-> Running RMSD per residue analysis')
        # WARNING: This analysis is fast enought to use the full trajectory instead of the reduced one
        # WARNING: However, the output file size depends on the trajectory size. It may be pretty big
        rmsd_per_residue(
            input_topology_filename = self.structure_file.path,
            input_trajectory_filename = self.trajectory_file.path,
            output_analysis_filename = output_analysis_filepath,
            structure = self.structure,
            pbc_residues = self.pbc_residues,
            snapshots = self.snapshots,
            frames_limit = 100,
        )

    # RMSD pairwise
    def run_rmsd_pairwise_analysis (self, overwrite : bool = False):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_RMSD_PAIRWISE_FILENAME)
        if exists(output_analysis_filepath) and not overwrite:
            return
        print('-> Running RMSD pairwise analysis')
        # WARNING: This analysis is fast enought to use the full trajectory instead of the reduced one
        # WARNING: However, the output file size depends on the trajectory size exponentially. It may be huge
        rmsd_pairwise(
            input_topology_filename = self.structure_file.path,
            input_trajectory_filename = self.trajectory_file.path,
            output_analysis_filename = output_analysis_filepath,
            interactions = self.interactions,
            structure = self.structure,
            pbc_residues = self.pbc_residues,
            snapshots = self.snapshots,
            frames_limit = 200,
            overall_selection = "name CA or name C5"
        )

    # Clusters
    def run_clusters_analysis (self, overwrite : bool = False):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_CLUSTERS_FILENAME)
        if exists(output_analysis_filepath) and not overwrite:
            return
        output_screenshot_filepath = self.md_pathify(OUTPUT_CLUSTER_SCREENSHOT_FILENAMES)
        print('-> Running clusters analysis')
        clusters_analysis(
            input_structure_filename = self.structure_file.path,
            input_trajectory_filename = self.trajectory_file.path,
            interactions = self.interactions,
            structure = self.structure,
            output_analysis_filename = output_analysis_filepath,
            output_screenshots_filename = output_screenshot_filepath,
        )

    # Distance per residue
    def run_dist_perres_analysis (self, overwrite : bool = False):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_DIST_PERRES_FILENAME)
        if exists(output_analysis_filepath) and not overwrite:
            return
        print('-> Running distance per residue analysis')
        # WARNING: This analysis is not fast enought to use the full trajectory. It would take a while
        distance_per_residue(
            input_topology_filename = self.structure_file.path,
            input_trajectory_filename = self.trajectory_file.path,
            output_analysis_filename = output_analysis_filepath,
            interactions = self.interactions,
            snapshots = self.snapshots,
            frames_limit = 200,
        )

    # Hydrogen bonds
    def run_hbonds_analysis (self, overwrite : bool = False):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_HBONDS_FILENAME)
        if exists(output_analysis_filepath) and not overwrite:
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
            input_topology_filename = self.structure_file.path,
            input_trajectory_filename = self.trajectory_file.path,
            output_analysis_filename = output_analysis_filepath,
            structure = self.structure,
            interactions = self.interactions,
            snapshots = self.snapshots,
            frames_limit = 200,
            # is_time_dependend = self.is_time_dependend,
            # time_splits = 100,
            # populations = self.populations
        )

    # SASA, solvent accessible surfave analysis
    def run_sas_analysis (self, overwrite : bool = False):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_SASA_FILENAME)
        if exists(output_analysis_filepath) and not overwrite:
            return
        print('-> Running SAS analysis')
        # Run the analysis
        sasa(
            input_topology_filename = self.structure_file.path,
            input_trajectory_filename = self.trajectory_file.path,
            output_analysis_filename = output_analysis_filepath,
            structure = self.structure,
            pbc_residues = self.pbc_residues,
            snapshots = self.snapshots,
            frames_limit = 100,
        )

    # Energies
    def run_energies_analysis (self, overwrite : bool = False):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_ENERGIES_FILENAME)
        if exists(output_analysis_filepath) and not overwrite:
            return
        print('-> Running energies analysis')
        # Run the analysis
        energies(
            input_trajectory_file = self.trajectory_file,
            output_analysis_filename = output_analysis_filepath,
            energies_folder = self.md_pathify(ENERGIES_FOLDER),
            structure = self.structure,
            interactions = self.interactions,
            charges = self.charges,
            snapshots = self.snapshots,
            frames_limit = 100,
        )

    # Pockets
    def run_pockets_analysis (self, overwrite : bool = False):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_POCKETS_FILENAME)
        if exists(output_analysis_filepath) and not overwrite:
            return
        print('-> Running pockets analysis')
        # Run the analysis
        pockets(
            structure_file = self.structure_file,
            trajectory_file = self.trajectory_file,
            pockets_prefix = self.md_pathify(OUTPUT_POCKET_STRUCTURES_PREFIX),
            output_analysis_filepath = output_analysis_filepath,
            mdpocket_folder = self.md_pathify(POCKETS_FOLDER),
            pbc_residues = self.pbc_residues,
            snapshots = self.snapshots,
            frames_limit = 100,
        )

    # Helical parameters
    # DANI: Al final lo reimplementará Subamoy (en python) osea que esto no lo hacemos de momento
    # def run_helical_analysis (self, overwrite : bool = False):
    #     # Do not run the analysis if the output file already exists
    #     output_analysis_filepath = self.md_pathify(OUTPUT_HELICAL_PARAMETERS_FILENAME)
    #     if exists(output_analysis_filepath) and not overwrite:
    #         return
    #     print('-> Running helical analysis')
    #     # Run the analysis
    #     helical_parameters(
    #         input_topology_filename = self.structure_file.path,
    #         input_trajectory_filename = self.trajectory_file.path,
    #         output_analysis_filename = output_analysis_filepath,
    #         structure = self.structure,
    #         frames_limit = 1000,
    #     )

    # Markov
    def run_markov_analysis (self, overwrite : bool = False):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_MARKOV_FILENAME)
        if exists(output_analysis_filepath) and not overwrite:
            return
        print('-> Running Markov analysis')
        # Run the analysis
        markov(
            input_topology_filename = self.structure_file.path,
            input_trajectory_filename = self.trajectory_file.path,
            output_analysis_filename = output_analysis_filepath,
            structure = self.structure,
            populations = self.populations,
            #transitions = self.transitions,
            rmsd_selection = PROTEIN_AND_NUCLEIC,
        )

# The project is the main project
# A project is a set of related MDs
# These MDs share all or most topology and metadata
class Project:
    def __init__ (self,
        # The local directory where the project takes place
        directory : str = '.',
        # Accession of the project in the database, given that this project is already uploaded
        accession : Optional[str] = None,
        # URL to query for missing files when an accession is provided
        database_url : str = DEFAULT_API_URL,
        # A file containing a lof of inputs related to metadata, MD simulation parameters and analysis configurations
        inputs_filepath : str = None,
        # The input topology filename
        # Multiple formats are accepted but the default is our own parsed json topology
        input_topology_filepath : Optional[str] = None,
        # Input structure filepath
        # It may be both relative to the project directory or to every MD directory
        input_structure_filepath : Optional[str] = None,
        # Input trajectory filepaths
        # These files are searched in every MD directory so the path MUST be relative
        input_trajectory_filepaths : Optional[str] = None,
        # Set the different MD directories to be run
        # Each MD directory must contain a trajectory and may contain a structure
        md_directories : Optional[List[str]] = None,
        # Reference MD directory
        # Project functions which require structure or trajectory will use the ones from the reference MD
        # If no reference is passed then the first directory is used
        reference_md_index : Optional[int] = None,        
        # Input populations and transitions (MSM only)
        populations_filepath : str = DEFAULT_POPULATIONS_FILENAME,
        transitions_filepath : str = DEFAULT_TRANSITIONS_FILENAME,
        # Processing and analysis instructions
        filter_selection : Union[bool, str] = False,
        pbc_selection : Optional[str] = None,
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
        sample_trajectory : bool = False,
    ):
        # Save input parameters
        self.directory = remove_final_slash(directory)
        self.accession = accession
        self.database_url = database_url
        # Set the project URL in case we have the required data
        self.url = None
        if database_url and accession:
            self.url = database_url + '/rest/current/projects/' + accession

        # Save inputs for the register, even if they are not used in the class

        # Set the inputs file
        # Set the expected default name in case there is no inputs file since it may be downloaded
        self._inputs_file = File(DEFAULT_INPUTS_FILENAME)
        # If there is an input filepath then use it
        if inputs_filepath:
            self._inputs_file = File(inputs_filepath)
        # Otherwise guess the inputs file using the accepted filenames
        else:
            for filename in ACCEPTED_INPUT_FILENAMES:
                inputs_file = File(filename)
                if inputs_file.exists:
                    self._inputs_file = inputs_file
                    break
        # Set the input topology file
        self._input_topology_filepath = input_topology_filepath
        self._input_topology_file = None
        # Input structure and trajectory filepaths
        # Do not parse them to files yet, let this to the MD class
        self._input_structure_filepath = input_structure_filepath
        self._input_trajectory_filepaths = input_trajectory_filepaths
        if self._input_trajectory_filepaths:
            self.check_input_trajectory_filepaths()
        # Input populations and transitions for MSM
        self.populations_filepath = populations_filepath
        self._populations_file = File(self.populations_filepath)
        self.transitions_filepath = transitions_filepath
        self._transitions_file = File(self.transitions_filepath)

        # Set the processed topology filepath, which depends on the input topology filename
        # Note that this file is different from the standard topology, although it may be standard as well
        self._topology_filepath = None
        self._topology_file = None

        # Set the standard topology file
        self._standard_topology_file = None

        # Set the MD directories
        self._md_directories = md_directories
        if self._md_directories:
            self.check_md_directories()

        # Set the reference MD
        self._reference_md = None
        self._reference_md_index = reference_md_index

        # Set the rest of inputs
        # Note that the filter selection variable is not handled here at all
        # This is just pased to the filtering function which knows how to handle the default
        self.filter_selection = filter_selection
        # PBC selection may come from the console or from the inputs
        self._pbc_selection = pbc_selection
        self.image = image
        self.fit = fit
        self.translation = translation
        self.mercy = mercy
        # Fix the mercy input, if needed
        # If a boolean is passed instead of a list then we set its corresponding value
        if type(mercy) == bool:
            if mercy:
                self.mercy = AVAILABLE_FAILURES
            else:
                self.mercy = []
        self.trust = trust
        # Fix the trust input, if needed
        # If a boolean is passed instead of a list then we set its corresponding value
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

        # Other values which may be found/calculated on demand
        self._available_files = None
        self._pbc_residues = None
        self._safe_bonds = None
        self._charges = None
        self._populations = None
        self._transitions = None
        self._protein_map = None
        self._ligand_map = None
        self.pubchem_name_list = None
        self._residue_map = None
        self._mds = None

        # Force a couple of extraordinary files which is generated if atoms are resorted
        self.resorted_bonds_file = File(RESORTED_BONDS_FILENAME)
        self.resorted_charges_file = File(RESORTED_CHARGES_FILENAME)

        # Set a new entry for the register
        # This is useful to track previous workflow runs and problems
        register_file = File(REGISTER_FILENAME)
        self.register = Register(register_file)

    # Given a filename or relative path, add the project directory path at the beginning
    def project_pathify (self, filename_or_relative_path : str) -> str:
        return self.directory + '/' + filename_or_relative_path

    # Check input trajectory file paths to be right
    # If there is any problem then fix it or directly raise an input error
    def check_input_trajectory_filepaths (self):
        # WARNING: Input trajectory filepaths may be both a list or a single string
        # However we must keep a list
        if type(self._input_trajectory_filepaths) == list:
            pass 
        elif type(self._input_trajectory_filepaths) == str:
            self._input_trajectory_filepaths = [ self._input_trajectory_filepaths ]
        else:
            raise InputError('Input trajectory filepaths must be a list of strings or a string')
        # Check no trajectory path is absolute
        for path in self._input_trajectory_filepaths:
            if path[0] == '/':
                raise InputError('Trajectory paths MUST be relative, not absolute (' + path + ')')

    # Check MD directories to be right
    # If there is any problem then directly raise an input error
    def check_md_directories (self):
        # Check there is at least one MD
        if len(self._md_directories) < 1:
            raise InputError('There must be at least one MD')
        # Check there are not duplicated MD directories
        if len(set(self._md_directories)) != len(self._md_directories):
            raise InputError('There are duplicated MD directories')

    # Set a function to get input structure file path
    def get_input_structure_filepath (self) -> str:
        # If we already have a value then return it
        # If this value was passed through command line then it would be set as the internal value already
        if self._input_structure_filepath:
            return self._input_structure_filepath
        # Otherwise we must find it
        # Check if the inputs file has the value
        if self.is_inputs_file_available():
            # Get the input value, whose key must exist
            inputs_value = self.inputs.get('input_structure_filepath', missing_input_exception)
            if inputs_value == missing_input_exception:
                raise InputError('Missing input "input_structure_filepath"')
            # If there is a valid input then use it
            if inputs_value:
                self._input_structure_filepath = inputs_value
                return self._input_structure_filepath
        # If there is not any input then return the default
        return STRUCTURE_FILENAME
    input_structure_filepath = property(get_input_structure_filepath, None, None, "Input structure file path (read only)")

    # Set a function to get input trajectory file paths
    def get_input_trajectory_filepaths (self) -> str:
        # If we already have a value then return it
        # If this value was passed through command line then it would be set as the internal value already
        if self._input_trajectory_filepaths:
            return self._input_trajectory_filepaths
        # Otherwise we must find it
        # Check if the inputs file has the value
        if self.is_inputs_file_available():
            # Get the input value, whose key must exist
            inputs_value = self.inputs.get('input_trajectory_filepaths', missing_input_exception)
            if inputs_value == missing_input_exception:
                raise InputError('Missing input "input_trajectory_filepaths"')
            # If there is a valid input then use it
            if inputs_value:
                self._input_trajectory_filepaths = inputs_value
                self.check_input_trajectory_filepaths()
                return self._input_trajectory_filepaths
        # If there is not any input then return the default
        return [ TRAJECTORY_FILENAME ]
    input_trajectory_filepaths = property(get_input_trajectory_filepaths, None, None, "Input trajectory file paths (read only)")

    # Set a function to get MD directories
    def get_md_directories (self) -> list:
        # If MD directories are already declared then return them
        if self._md_directories:
            return self._md_directories
        # Otherwise use the default MDs
        self._md_directories = []
        # Use the MDs from the inputs file when available
        if self.is_inputs_file_available() and self.input_mds:
            for input_md in self.input_mds:
                # Set the directory from the MD name
                name = input_md['name']
                directory = name_2_directory(name)
                self._md_directories.append(directory)
        # Otherwise, guess MD directories by checking which directories include a register file
        else:
            available_directories = sorted(next(walk(self.directory))[1])
            for directory in available_directories:
                if exists(directory + '/' + REGISTER_FILENAME):
                    self._md_directories.append(directory)
            # If we found no MD directory then it means MDs were never declared before
            if len(self._md_directories) == 0:
                raise InputError('Impossible to know which are the MD directories. '
                    'You can either declare them using the "-mdir" option or by providing an inputs file')
        self.check_md_directories()
        return self._md_directories
    md_directories = property(get_md_directories, None, None, "MD directories (read only)")

    # Set the reference MD index
    def get_reference_md_index (self) -> int:
        # If we are already have a value then return it
        if self._reference_md_index:
            return self._reference_md_index
        # Otherwise we must find the reference MD index
        # If the inputs file is available then it must declare the reference MD index
        if self.is_inputs_file_available():
            self._reference_md_index = self.input_reference_md_index
        # Otherwise we simply set the first MD as the reference and warn the user about this
        else:
            #print('WARNING: No reference MD was specified. The first MD will be used as reference.')
            self._reference_md_index = 0
        return self._reference_md_index
    reference_md_index = property(get_reference_md_index, None, None, "Reference MD index (read only)")

    # Set the reference MD
    def get_reference_md (self) -> int:
        # If we are already have a value then return it
        if self._reference_md:
            return self._reference_md
        # Otherwise we must find the reference MD
        self._reference_md = self.mds[self.reference_md_index]
        return self._reference_md
    reference_md = property(get_reference_md, None, None, "Reference MD (read only)")

    # Setup the MDs
    def get_mds (self) -> list:
        # If MDs are already declared then return them
        if self._mds:
            return self._mds
        # Now instantiate a new MD for each declared MD and save the reference MD
        self._mds = []
        for n, md_directory in enumerate(self.md_directories, 1):
            md = MD(
                project = self, number = n, directory = md_directory,
                input_structure_filepath = self.input_structure_filepath,
                input_trajectory_filepaths = self.input_trajectory_filepaths,
            )
            self._mds.append(md)
        return self._mds
    mds = property(get_mds, None, None, "Available MDs (read only)")

    # Check input files exist when their filenames are read
    # If they do not exist then try to download them
    # If the download is not possible then raise an error

    # Inputs filename ------------

    # Set a function to check if inputs file is available
    # Note that asking for it when it is not available will lead to raising an input error
    def is_inputs_file_available (self) -> bool:
        # If name is not declared then it is impossible to reach it
        if not self._inputs_file:
            return False
        # If the file already exists then it is available
        if self._inputs_file.exists:
            return True
        # If it does not exist but it may be downloaded then it is available
        if self.url:
            return True
        return False

    # Set a function to load the inputs file
    def get_inputs_file (self) -> File:
        # There must be an inputs filename
        if not self._inputs_file:
            raise InputError('Not defined inputs filename')
        # If the file already exists then we are done
        if self._inputs_file.exists:
            return self._inputs_file
        # Try to download it
        # If we do not have the required parameters to download it then we surrender here
        if not self.url:
            raise InputError('Missing inputs file "' + self._inputs_file.filename + '"')
        # Download the inputs json file if it does not exists
        sys.stdout.write('Downloading inputs (' + self._inputs_file.filename + ')\n')
        inputs_url = self.url + '/inputs'
        # In case this is a json file we must specify the format in the query
        is_json = self._inputs_file.format == 'json'
        if is_json:
            inputs_url += '?format=json'
        urllib.request.urlretrieve(inputs_url, self._inputs_file.path)
        # Rewrite the inputs file in a pretty formatted way in case it is a JSON file
        if is_json:
            file_content = load_json(self._inputs_file.path)
            save_json(file_content, self._inputs_file.path, indent = 4)
        return self._inputs_file
    inputs_file = property(get_inputs_file, None, None, "Inputs filename (read only)")

    # Topology filename ------------

    # If there is not input topology filename, which is possible, we must guess the topology filename
    # Note that if we can download then the name will always be the remote topology name (topology.json)
    def guess_input_topology_filename (self) -> Optional[str]:
        # If we can download then the we use the remote topology filename (topology.json)
        if self.url:
            return TOPOLOGY_FILENAME
        # Otherwise we must guess among the local files
        # If the default topology filename exists then use it
        if exists(TOPOLOGY_FILENAME):
            return TOPOLOGY_FILENAME
        # Find if the raw charges file is present
        if exists(RAW_CHARGES_FILENAME):
            return RAW_CHARGES_FILENAME
        # Otherwise, find all possible accepted topology formats
        for topology_format in ACCEPTED_TOPOLOGY_FORMATS:
            topology_filename = 'topology.' + topology_format
            if exists(topology_filename):
                return topology_filename
        return None        

    # Get the input topology filepath from the inputs
    def get_input_topology_filepath (self) -> File:
        # If we already have a value then return it
        # If this value was passed through command line then it would be set as the internal value already
        if self._input_topology_filepath:
            return self._input_topology_filepath
        # Otherwise we must find it
        # Check if the inputs file has the value
        if self.is_inputs_file_available():
            # Get the input value, whose key must exist
            inputs_value = self.inputs.get('input_topology_filepath', missing_input_exception)
            if inputs_value == missing_input_exception:
                raise InputError('Missing input "input_topology_filepath"')
            # If there is a valid input then use it
            if inputs_value:
                self._input_topology_filepath = inputs_value
                return self._input_topology_filepath
        # Otherwise we must guess which is the topology file
        guess = self.guess_input_topology_filename()
        if guess:
            self._input_topology_filepath = guess
            return self._input_topology_filepath
        raise InputError('Missing input topology file path')
    input_topology_filepath = property(get_input_topology_filepath, None, None, "Input topology file path (read only)")

    # Get the input topology file
    # If the file is not found try to download it
    def get_input_topology_file (self) -> File:
        # If we already have a value then return it
        if self._input_topology_file:
            return self._input_topology_file
        # Set the file
        self._input_topology_file = File(self.input_topology_filepath)
        # If the file already exists then we are done
        if self._input_topology_file.exists:
            return self._input_topology_file
        # Try to download it
        # If we do not have the required parameters to download it then we surrender here
        if not self.url:
            raise InputError('Missing input topology file "' + self._input_topology_file.filename + '"')
        # If the input topology name is the standard then download it from the topology endpoint
        if self._input_topology_file.filename == TOPOLOGY_FILENAME:
            # Check if the project has a topology and download it in json format if so
            topology_url = self.url + '/topology'
            try:
                sys.stdout.write('Downloading topology (' + self._input_topology_file.filename + ')\n')
                urllib.request.urlretrieve(topology_url, self._input_topology_file.path)
                # Rewrite the topology file in a pretty formatted way
                file_content = load_json(self._input_topology_file.path)
                save_json(file_content, self._input_topology_file.path, indent = 4)
            except:
                raise Exception('Something where wrong while downloading the topology')
            return self._input_topology_file
        # Otherwise, try to download it using the files endpoint
        # Note that this is not usually required
        topology_url = self.url + '/files/' + self._input_topology_file.filename
        try:
            sys.stdout.write('Downloading topology (' + self._input_topology_file.filename + ')\n')
            urllib.request.urlretrieve(topology_url, self._input_topology_file.path)
        except:
            raise Exception('Something where wrong while downloading the topology')
        # Before we finish, in case the topology is a '.top' file, we may need to download the itp files as well
        if self._input_topology_file.filename == 'topology.top':
            # Find available .itp files and download each of them
            itp_filenames = [filename for filename in self.available_files if filename[-4:] == '.itp']
            for itp_filename in itp_filenames:
                sys.stdout.write('Downloading itp file (' + itp_filename + ')\n')
                itp_url = self.url + '/files/' + itp_filename
                urllib.request.urlretrieve(itp_url, itp_filename)
        return self._input_topology_file
    input_topology_file = property(get_input_topology_file, None, None, "Input topology file (read only)")

    # Input structure filename ------------

    # Get the input structure filename
    # When calling this function make sure all MDs have the file or try to download it
    def get_input_structure_file (self) -> File:
        return self.reference_md._input_structure_file
    input_structure_file = property(get_input_structure_file, None, None, "Input structure filename for each MD (read only)")

    # Input trajectory filename ------------

    # Get the input trajectory filename(s) from the inputs
    # If file(s) are not found try to download it
    def get_input_trajectory_files (self) -> List[File]:
        return self.reference_md._input_trajectory_files
    input_trajectory_files = property(get_input_trajectory_files, None, None, "Input trajectory filenames for each MD (read only)")

    # Populations filename ------------

    def get_populations_file (self) -> File:
        if not self.get_file(self._populations_file):
            return None
        return self._populations_file
    populations_file = property(get_populations_file, None, None, "MSM equilibrium populations filename (read only)")

    # Transitions filename ------------

    def get_transitions_file (self) -> Optional[str]:
        if not self.get_file(self._transitions_file):
            return None
        return self._transitions_file
    transitions_file = property(get_transitions_file, None, None, "MSM transition probabilities filename (read only)")

    # ---------------------------------

    # DISCLAIMER:
    # Note that requesting files from any of the MDs already returns the project-specific files
    # For this reason these two functions are inherited from the MD class although 'get_file' is used for project files only

    # Remote available files
    # Note that requesting files from any of the MDs already return the project-specific files
    def get_available_files (self) -> List[str]:
        return self.reference_md.available_files
    available_files = property(get_available_files, None, None, "Remote available files (read only)")

    # Check if a file exists
    # If not, try to download it from the database
    # If the file is not found in the database it is fine, we do not even warn the user
    # Note that nowadays this function is used to get populations and transitions files, which are not common
    def get_file (self, target_file : File) -> bool:
        return self.reference_md.get_file(target_file)

    # Input file values -----------------------------------------

    # First of all set input themselves

    # Get inputs
    def get_inputs (self) -> dict:
        # If inputs are already loaded then return them
        if self._inputs:
            return self._inputs
        # Otherwise, load inputs from the inputs file
        inputs_data = None
        if self.inputs_file.format == 'json':
            inputs_data = load_json(self.inputs_file.path)
        elif self.inputs_file.format == 'yaml':
            inputs_data = load_yaml(self.inputs_file.path)
        else:
            raise InputError('Input file format is not supported. Please use json or yaml files.')
        if not inputs_data:
            raise InputError('Input file is empty')
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
                raise InputError('Missing input "' + name + '"')
            return value
        return getter

    # Assign the getters
    input_interactions = property(input_getter('interactions'), None, None, "Interactions to be analyzed (read only)")
    forced_references = property(input_getter('forced_references'), None, None, "Uniprot IDs to be used first when aligning protein sequences (read only)")
    pdb_ids = property(input_getter('pdbIds'), None, None, "Protein Data Bank IDs used for the setup of the system (read only)")
    input_type = property(input_getter('type'), None, None, "Set if its a trajectory or an ensemble (read only)")
    input_mds = property(input_getter('mds'), None, None, "Input MDs configuration (read only)")
    input_reference_md_index = property(input_getter('mdref'), None, None, "Input MD reference index (read only)")
    input_ligands = property(input_getter('ligands'), None, None, "Input ligand references (read only)")

    # PBC selection may come from the console or from the inputs file
    # Console has priority over the inputs file
    def get_pbc_selection (self) -> Optional[str]:
        # If we have an internal value then return it
        if self._pbc_selection:
            return self._pbc_selection
        # As an exception, we avoid asking for the inputs file if it is not available
        # This input is required for some early processing steps where we do not need the inputs file for anything else
        if not self.is_inputs_file_available():
            return None
        # Otherwise, find it in the inputs
        # Get the input value, whose key must exist
        self._pbc_selection = self.inputs.get('pbc_selection', missing_input_exception)
        if self._pbc_selection == missing_input_exception:
            raise InputError('Missing input "pbc_selection"')
        return self._pbc_selection
    pbc_selection = property(get_pbc_selection, None, None, "Selection of atoms which are still in periodic boundary conditions (read only)")

    # Set additional values infered from input values

    # Set if MDs are time dependent
    def check_is_time_dependent (self) -> bool:
        if self.input_type == 'trajectory':
            return True
        elif self.input_type == 'ensemble':
            return False
        raise InputError('Not supported input type value: ' + self.input_type)
    is_time_dependend = property(check_is_time_dependent, None, None, "Check if trajectory frames are time dependent (read only)")

    # Processed files ----------------------------------------------------

    # Set the expected output topology filename given the input topology filename
    # Note that topology formats are conserved
    def inherit_topology_filename (self) -> str:
        filename = self.input_topology_file.filename
        if not filename:
            return None
        if filename == TOPOLOGY_FILENAME:
            return filename
        if filename == RAW_CHARGES_FILENAME:
            return filename
        standard_format = self.input_topology_file.format
        return 'topology.' + standard_format

    # Get the processed topology file path
    def get_topology_filepath (self) -> str:
        # If we have a stored value then return it
        if self._topology_filepath:
            return self._topology_filepath
        # Otherwise we must find it
        self._topology_filepath = self.inherit_topology_filename()
        return self._topology_filepath
    topology_filepath = property(get_topology_filepath, None, None, "Topology file path (read only)")

    # Get the processed topology file
    def get_topology_file (self) -> str:
        # If we have a stored value then return it
        # This means we already found or generated this file
        if self._topology_file:
            return self._topology_file
        # If the file already exists then we are done
        self._topology_file = File(self.topology_filepath)
        # Run the processing logic
        self.reference_md.process_input_files()
        # Now that the file is sure to exist we return it
        return self._topology_file
    topology_file = property(get_topology_file, None, None, "Topology file (read only)")

    # Get the processed structure from the reference MD
    def get_structure_file (self) -> str:
        return self.reference_md.structure_file
    structure_file = property(get_structure_file, None, None, "Structure filename from the reference MD (read only)")

    # Get the processed trajectory from the reference MD
    def get_trajectory_file (self) -> str:
        return self.reference_md.trajectory_file
    trajectory_file = property(get_trajectory_file, None, None, "Trajectory filename from the reference MD (read only)")

    # ---------------------------------------------------------------------------------
    # Others values which may be found/calculated and files to be generated on demand
    # ---------------------------------------------------------------------------------

    # Parsed structure from reference MD
    def get_structure (self) -> 'Structure':
        return self.reference_md.structure
    structure = property(get_structure, None, None, "Parsed structure from the reference MD (read only)")

    # Indices of residues in periodic boundary conditions
    def get_pbc_residues (self) -> List[int]:
        # If we already have a stored value then return it
        if self._pbc_residues:
            return self._pbc_residues
        # Otherwise we must find the value
        self._pbc_residues = get_pbc_residues(self.structure, self.pbc_selection)
        return self._pbc_residues
    pbc_residues = property(get_pbc_residues, None, None, "Indices of residues in periodic boundary conditions (read only)")

     # Safe bonds
    def get_safe_bonds (self) -> List[List[int]]:
        # If we already have a stored value then return it
        # WARNING: Do not remove the self.topology_file checking
        # Note that checking if the topology file exists triggers all the processing logic
        # The processing logic is able to set the internal safe bonds value as well so this avoids repeating the process
        # Besides it generates the resorted files, if needed
        if self.topology_file and self._safe_bonds:
            return self._safe_bonds
        # If we have a resorted file then use it
        # Note that this is very excepcional
        if self.resorted_bonds_file.exists:
            print('Using resorted safe bonds')
            self._safe_bonds = load_json(self.resorted_bonds_file.path)
            return self._safe_bonds
        # Otherwise we must find safe bonds value
        # This should only happen if we are working with already processed files
        self._safe_bonds = get_safe_bonds(
            input_topology_file=self.topology_file,
            input_structure_file=self.structure_file,
            input_trajectory_file=self.trajectory_file,
            must_check_stable_bonds=(STABLE_BONDS_FLAG not in self.trust),
            snapshots=self.reference_md.snapshots,
            structure=self.structure,
        )
        return self._safe_bonds
    safe_bonds = property(get_safe_bonds, None, None, "Atom bonds to be trusted (read only)")

    # Atom charges
    def get_charges (self) -> List[float]:
        # If we already have a stored value then return it
        # WARNING: Do not remove the self.topology_file checking
        # Note that checking if the topology file exists triggers all the processing logic
        # The processing logic is able to set the internal atom charges value as well so this avoids repeating the process
        # Besides it generates the resorted files, if needed
        if self.topology_file and self._charges:
            return self._charges
        # If we have a resorted file then use it
        # Note that this is very excepcional
        if self.resorted_charges_file.exists:
            print('Using resorted atom charges')
            self._charges = load_json(self.resorted_charges_file.path)
            return self._charges
        # Otherwise we must find the value
        self._charges = get_charges(self.topology_file)
        return self._charges
    charges = property(get_charges, None, None, "Atom charges (read only)")

    # Equilibrium populations from a MSM
    def get_populations (self) -> Optional[List[float]]:
        # If we already have a stored value then return it
        if self._populations:
            return self._populations
        # Otherwise we must find the value
        if not self.populations_file:
            return None
        self._populations = read_file(self.populations_file)
        return self._populations
    populations = property(get_populations, None, None, "Equilibrium populations from a MSM (read only)")

    # Transition probabilities from a MSM
    def get_transitions (self) -> Optional[List[List[float]]]:
        # If we already have a stored value then return it
        if self._transitions:
            return self._transitions
        # Otherwise we must find the value
        if not self.transitions_file:
            return None
        self._transitions = read_file(self.transitions_file)
        return self._transitions
    transitions = property(get_transitions, None, None, "Transition probabilities from a MSM (read only)")

    # Protein residues mapping
    def get_protein_map (self) -> List[dict]:
        # If we already have a stored value then return it
        if self._protein_map:
            return self._protein_map
        print('-> Getting protein references')
        # Set the protein references file
        protein_references_filepath = self.project_pathify(PROTEIN_REFERENCES_FILENAME)
        protein_references_file = File(protein_references_filepath)
        # Otherwise we must find the value
        self._protein_map = generate_protein_mapping(
            structure = self.structure,
            protein_references_file = protein_references_file,
            register = self.register,
            mercy = self.mercy,
            forced_references = self.forced_references,
            pdb_ids = self.pdb_ids,
        )
        return self._protein_map
    protein_map = property(get_protein_map, None, None, "Residues mapping (read only)")

    # Define the output file of the protein mapping including protein references
    def get_protein_references_file (self, overwrite : bool = False) -> File:
        # Set the protein references file
        protein_references_filepath = self.project_pathify(PROTEIN_REFERENCES_FILENAME)
        protein_references_file = File(protein_references_filepath)
        # If the file already exists then return it
        # However if the overwrite argument is passed then delete it and proceed to produce it again
        if protein_references_file.exists:
            if not overwrite:
                return protein_references_file
            protein_references_file.remove()
        # Ask for the protein map thus producing the protein references file
        self.get_protein_map()
        return protein_references_file
    protein_references_file = property(get_protein_references_file, None, None, "File including protein refereces data mined from UniProt (read only)")

    # Ligand residues mapping
    def get_ligand_map (self) -> List[dict]:
        # If we already have a stored value then return it
        if self._ligand_map != None:
            return self._ligand_map
        print('-> Getting ligand references')
        # Set the ligand references file
        ligand_references_filepath = self.project_pathify(LIGAND_REFERENCES_FILENAME)
        ligand_references_file = File(ligand_references_filepath)
        # Otherwise we must find the value
        self._ligand_map, self.pubchem_name_list = generate_ligand_mapping(
            structure = self.structure,
            input_ligands = self.input_ligands,
            input_pdb_ids = self.pdb_ids,
            output_ligands_filepath = ligand_references_file.path, 
            mercy = self.mercy,
        )
        return self._ligand_map
    ligand_map = property(get_ligand_map, None, None, "Ligand references (read only)")

    # Define the output file of the ligand mapping including ligand references
    def get_ligand_references_file (self, overwrite : bool = False) -> File:
        # Set the ligand references file
        ligand_references_filepath = self.project_pathify(LIGAND_REFERENCES_FILENAME)
        ligand_references_file = File(ligand_references_filepath)
        # If the file already exists then return it
        # However if the overwrite argument is passed then delete it and proceed to produce it again
        if ligand_references_file.exists:
            if not overwrite:
                return ligand_references_file
            ligand_references_file.remove()
        # Ask for the ligand map thus producing the ligand references file
        self.get_ligand_map()
        return ligand_references_file
    ligand_references_file = property(get_ligand_references_file, None, None, "File including ligand refereces data mined from PubChem (read only)")

    # Build the residue map from both proteins and ligands maps
    # This is formatted as both the standard topology and metadata generators expect them
    def get_residue_map (self) -> List[dict]:
        # If we already have a stored value then return it
        if self._residue_map:
            return self._residue_map
        print('-> Getting residue references')
        # Otherwise we must find the value
        self._residue_map = generate_residue_mapping(
            protein_map = self.protein_map,
            ligand_map = self.ligand_map,
            structure = self.structure,
        )
        return self._residue_map
    residue_map = property(get_residue_map, None, None, "Residue map (read only)")

    # Metadata filename
    def get_metadata_file (self, overwrite : bool = False) -> File:
        # Set the metadata file
        metadata_filepath = self.project_pathify(OUTPUT_METADATA_FILENAME)
        metadata_file = File(metadata_filepath)
        # If the file already exists then send it
        if metadata_file.exists and not overwrite:
            return metadata_file
        print('-> Generating project metadata')
        # Set an input getter that gets the input as soon as called
        def get_input (name : str, optional : bool = False):
            if optional:
                return Project.input_getter(name, None)(self)
            return Project.input_getter(name)(self)
        # Otherwise, generate it
        generate_project_metadata(
            input_structure_filename = self.structure_file.path,
            input_trajectory_filename = self.trajectory_file.path,
            get_input = get_input,
            structure = self.structure,
            residue_map = self.residue_map,
            protein_references_file = self.protein_references_file,
            interactions = self.reference_md.interactions,
            register = self.register,
            output_metadata_filename = metadata_file.path,
            ligand_customized_names = self.pubchem_name_list,
        )
        return metadata_file
    metadata_file = property(get_metadata_file, None, None, "Project metadata filename (read only)")

    # Standard topology filename
    def get_standard_topology_file (self, overwrite : bool = False) -> File:
        # If we have a stored value then return it
        # This means we already found or generated this file
        if self._standard_topology_file:
            return self._standard_topology_file
        # Set the standard topology file
        standard_topology_filepath = self.project_pathify(TOPOLOGY_FILENAME)
        self._standard_topology_file = File(standard_topology_filepath)
        # If the file already exists and is not to be overwirtten then send it
        if self._standard_topology_file.exists and not overwrite:
            return self._standard_topology_file
        print('-> Generating topology')
        # Otherwise, generate it
        generate_topology(
            structure = self.structure,
            charges = self.charges,
            residue_map = self.residue_map,
            pbc_residues = self.pbc_residues,
            output_topology_filepath = self._standard_topology_file.path
        )
        # Register the last modification times of the new generated standard topology file
        # Note that this may be already the topology file, but it may be not
        self.register.update_mtime(self._standard_topology_file)
        return self._standard_topology_file
    standard_topology_file = property(get_standard_topology_file, None, None, "Standard topology filename (read only)")

    # Screenshot filename
    def get_screenshot_filename (self, overwrite : bool = False) -> str:
        # Set the screenshot file
        screenshot_filepath = self.project_pathify(OUTPUT_SCREENSHOT_FILENAME)
        screenshot_file = File(screenshot_filepath)
        # If the file already exists then send it
        if screenshot_file.exists and not overwrite:
            return screenshot_file
        print('-> Generating screenshot')
        # Otherwise, generate it
        get_screenshot(
            input_structure_filename = self.structure_file.path,
            output_screenshot_filename = screenshot_file.path,
        )
        return screenshot_file
    screenshot_filename = property(get_screenshot_filename, None, None, "Screenshot filename (read only)")


# AUXILIAR FUNCTIONS ---------------------------------------------------------------------------

# Set a function to read a file which may be in differen formats
# DANI: En cuanto se concrete el formato de los markov esta función no hará falta
def read_file (target_file : File) -> dict:
    # Get the file format
    file_format = target_file.filename.split('.')[-1]
    # Read numpy files
    if file_format == 'npy':
        return numpy.load(target_file.path)
    # Read JSON files
    if file_format == 'json':
        return load_json(target_file.path)

# Set a function to convert an MD name into an equivalent MD directory
def name_2_directory (name : str) -> str:
    # Replace white spaces with underscores
    directory = name.replace(' ', '_')
    # Remove problematic characters
    for character in FORBIDEN_DIRECTORY_CHARACTERS:
        directory = directory.replace(character, '')
    return directory

# Set a function to convert an MD directory into an equivalent MD name
def directory_2_name (directory : str) -> str:
    # Replace white spaces with underscores
    name = directory.replace('_', ' ')
    return name

# Remove the final slash if exists since it may cuse problems when recognizing input directories
def remove_final_slash (directory : str) -> str:
    if directory[-1] == '/':
        return directory[:-1]
    return directory

# Input files
input_files = {
    'istructure': MD.get_input_structure_file,
    'itrajectory': MD.get_input_trajectory_files,
    'itopology': Project.get_input_topology_file,
    'inputs': Project.get_inputs_file,
    'populations': Project.get_populations_file,
    'transitions': Project.get_transitions_file
}

# Processed files
processed_files = {
    'structure': MD.get_structure_file,
    'trajectory': MD.get_trajectory_file,
    'topology': Project.get_topology_file
}

# List of available analyses
analyses = {
    'clusters': MD.run_clusters_analysis,
    'dist': MD.run_dist_perres_analysis,
    'energies': MD.run_energies_analysis,
    'hbonds': MD.run_hbonds_analysis,
    #'helical': MD.run_helical_analysis,
    'markov': MD.run_markov_analysis,
    'pca': MD.run_pca_analysis,
    #'pcacons': MD.run_pca_contacts,
    'pockets': MD.run_pockets_analysis,
    'rgyr': MD.run_rgyr_analysis,
    'rmsds': MD.run_rmsds_analysis,
    'perres': MD.run_rmsd_perres_analysis,
    'pairwise': MD.run_rmsd_pairwise_analysis,
    'rmsf': MD.run_rmsf_analysis,
    'sas': MD.run_sas_analysis,
    'tmscore': MD.run_tmscores_analysis,
}

# List of requestables for the console
requestables = {
    **input_files,
    **processed_files,
    **analyses,
    'interactions': MD.get_interactions,
    'mapping': Project.get_protein_references_file,
    'ligands': Project.get_ligand_references_file,
    'screenshot': Project.get_screenshot_filename,
    'stopology': Project.get_standard_topology_file,
    'pmeta': Project.get_metadata_file,
    'mdmeta': MD.get_metadata_file,
}

# The actual main function
def workflow (
    # Project parameters
    project_parameters : dict = {},
    # The actual workflow parameters
    # The working directory
    working_directory : str = '.',
    # Download only
    download : bool = False,
    # Download and correct only
    setup : bool = False,
    # Run only specific analyses/processes
    include : Optional[List[str]] = None,
    # Run everything but specific analyses/processes
    exclude : Optional[List[str]] = None,
    # Overwrite already existing output files
    overwrite : Optional[ Union[ List[str], bool ] ] = None,
):

    # Check there are not input errors

    # Include and exclude are not compatible
    # This is to protect the user to do something which makes not sense
    if include and exclude:
        raise InputError('Include (-i) and exclude (-e) are not comaptible. Use one of these options.')

    # Move the current directory to the working directory
    chdir(working_directory)
    current_directory_name = getcwd().split('/')[-1]
    print('Running workflow for project at ' + current_directory_name)

    # Initiate the project project
    project = Project(**project_parameters)
    print(f'  {len(project.mds)} MDs are to be run')

    # Now iterate over the different MDs
    for md in project.mds:

        print('\n' + CYAN_HEADER + 'Running workflow for MD at ' + md.directory + COLOR_END)

        # Set a function to call getters with the proper instance
        def call_getter (getter : Callable, **args):
            instance = md if getter.__qualname__[0:3] == 'MD.' else project
            getter(instance, **args)

        # If download is passed as True then just download input files if they are missing and exit
        if download:
            for getter in input_files.values():
                call_getter(getter)
            continue

        # If setup is passed as True then process input files and exit here
        if setup:
            for getter in processed_files.values():
                call_getter(getter)
            continue

        # Set the list of processings and analyses to run
        tasks = None
        # If the user included specific tasks then add only these tasks to the list
        if include and len(include) > 0:
            sys.stdout.write(f"Executing specific dependencies: " + ', '.join(include) + '\n')
            tasks = include
        # Set the default tasks otherwise
        else:
            # Add all the analysis, the metadata, the standard topology and the screenshot
            tasks = [
                # Mapping is not included here
                # It does more things than generating the references file and it takes time
                'stopology',
                'screenshot',
                'pmeta',
                'ligands',
                'interactions',
                'mdmeta',
                *analyses.keys(),
            ]
            # If the exclude parameter was passed then remove excluded tasks from the default tasks
            if exclude and len(exclude) > 0:
                sys.stdout.write(f"Excluding specific dependencies: " + ', '.join(exclude) + '\n')
                tasks = [ name for name in tasks if name not in exclude ]

        # If the user requested to verwrite something, make sure it is in the tasks list
        if overwrite and type(overwrite) == list:
            for task in overwrite:
                if task not in tasks:
                    raise InputError(f'Task "{task}" is to be overwriten but it is not in the tasks list. Either include it or do not exclude it')
            
        # Run the tasks
        for task in tasks:
            # Get the function to be called
            getter = requestables[task]
            # Check if the current task is to run even if output files exists thus overriding them
            must_overwrite = overwrite == True or ( type(overwrite) == list and task in overwrite )
            # If we must overwrite then call the function with the overwrite parameter set as true
            if must_overwrite:
                # Make sure the function to be called has the overwrite argument
                function_arguments = getfullargspec(getter)[0]
                if 'overwrite' not in function_arguments:
                    raise InputError(f'Task "{task}" does not support overwrite')
                # Finally call the function
                call_getter(getter, overwrite=True)
                continue
            # Call the function with no additional arguments otherwise
            call_getter(getter)

        # Remove gromacs backups and other trash files
        remove_trash(md.directory)

    sys.stdout.write("Done!\n")