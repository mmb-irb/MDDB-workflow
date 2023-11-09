#!/usr/bin/env python

# This is the starter script

# Import python libraries
from os import chdir, remove, symlink, rename, walk, mkdir, getcwd
from os.path import exists, isabs
import sys
import io
import math
from pathlib import Path
import urllib.request
import json
import numpy
from glob import glob
from typing import Optional, Union, List, Callable

# Constants
from model_workflow.constants import *

# Import local tools
from model_workflow.tools.topology_manager import setup_structure
from model_workflow.tools.get_pytraj_trajectory import get_pytraj_trajectory
from model_workflow.tools.get_first_frame import get_first_frame
from model_workflow.tools.get_average import get_average
from model_workflow.tools.process_interactions import process_interactions
from model_workflow.tools.get_pbc_residues import get_pbc_residues
from model_workflow.tools.generate_metadata import generate_project_metadata, generate_md_metadata
from model_workflow.tools.generate_map import generate_map_online
from model_workflow.tools.generate_topology import generate_topology
from model_workflow.tools.get_summarized_trajectory import get_summarized_trajectory
from model_workflow.tools.get_frames_count import get_frames_count
from model_workflow.tools.get_charges import get_charges
from model_workflow.tools.remove_trash import remove_trash
from model_workflow.tools.get_screenshot import get_screenshot
from model_workflow.tools.filter_atoms import filter_atoms
from model_workflow.tools.image_and_fit import image_and_fit
from model_workflow.tools.structure_corrector import structure_corrector
#from model_workflow.tools.httpsf import mount
from model_workflow.tools.register import Register

# Import local analyses
from model_workflow.analyses.rmsds import rmsds
from model_workflow.analyses.tmscores import tmscores
from model_workflow.analyses.rmsf import rmsf
from model_workflow.analyses.rgyr import rgyr
from model_workflow.analyses.pca import pca
#from model_workflow.analyses.pca_contacts import pca_contacts
from model_workflow.analyses.rmsd_per_residue import rmsd_per_residue
from model_workflow.analyses.rmsd_pairwise import rmsd_pairwise
from model_workflow.analyses.distance_per_residue import distance_per_residue
#from model_workflow.analyses.hydrogen_bonds_2 import hydrogen_bonds
from model_workflow.analyses.hydrogen_bonds import hydrogen_bonds
from model_workflow.analyses.sasa import sasa
from model_workflow.analyses.energies import energies
from model_workflow.analyses.pockets import pockets
from model_workflow.analyses.rmsd_check import check_trajectory_integrity
#from model_workflow.analyses.helical_parameters import helical_parameters
from model_workflow.analyses.markov import markov

# Import mdtoolbelt tools
from mdtoolbelt.conversions import convert
from mdtoolbelt.structures import Structure
from mdtoolbelt.file import File
from mdtoolbelt.auxiliar import InputError

# Make the system output stream to not be buffered
# This is useful to make prints work on time in Slurm
# Otherwise, output logs are written after the script has fully run
# Note that this fix affects all modules and built-ins
unbuffered = io.TextIOWrapper(open(sys.stdout.fileno(), 'wb', 0), write_through=True)
sys.stdout = unbuffered

# Set an special exception for input errors
missing_input_exception = Exception('Missing input')

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
        self.directory = directory
        # If the directory does not exists then we may have 2 different scenarios
        if not exists(self.directory):
            # If we have an URL to donwload then it means we must download input files
            # Thus we simpy create the missing directroy and it be filled further
            if self.url:
                mkdir(self.directory)
            # Otherwise we are supposed to find input files locally
            # If the directory does not exist we have nothing to do so we raise an error
            else:
                raise InputError('MD directory ' + self.directory + ' does not exist')
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
            for parsed_path in glob(path):
                self._input_trajectory_files.append( File(parsed_path) )
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

        # Tests
        self._trajectory_integrity = None

        # Set a new MD specific register
        # In case the directory is the project directory itself, use the project register
        register_file = File(self.md_pathify(REGISTER_FILENAME))
        if register_file.path == self.project.register.file.path:
            self.register = self.project.register
        else:
            self.register = Register(file_path = register_file.path)

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
            raise InputError('Missing inputs file "' + self._input_structure_file.filename + '"')
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
            missing_filepaths = [ trajectory_file.relative_path for trajectory_file in missing_input_trajectory_files ]
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
        for md in self.project.input_mds:
            name = md['name']
            directory = name_2_directory(name)
            if directory == self.directory:
                self._md_inputs = md
                return self._md_inputs
        raise InputError('No MD input matches the current directory (' + self.directory + ')')

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
        output_topology_file = File(self.project._topology_filepath)

        print('-> Processing input files')

        # --- CONVERTING AND MERGING ------------------------------------------------------------

        # Set output filenames for the already converted structure
        input_structure_format = self.input_structure_file.format
        output_structure_format = output_structure_file.format
        converted_structure_filepath = self.md_pathify(CONVERTED_STRUCTURE)
        # If input structure already matches the output format then avoid the renaming
        if input_structure_format == output_structure_format:
            converted_structure_filepath = input_structure_file.path
        converted_structure_file = File(converted_structure_filepath)
        # Set output filenames for the already converted structure
        # Input trajectories should have all the same format
        input_trajectory_formats = set([ trajectory_file.format for trajectory_file in input_trajectory_files ])
        if len(input_trajectory_formats) > 1:
            raise InputError('All input trajectory files must have the same format')
        input_trajectories_format = list(input_trajectory_formats)[0]
        output_trajectory_format = output_trajectory_file.format
        converted_trajectory_filepath = self.md_pathify(CONVERTED_TRAJECTORY)
        # If input trajectory already matches the output format and is unique then avoid the renaming
        if input_trajectories_format == output_trajectory_format and len(input_trajectory_files) == 1:
            converted_trajectory_filepath = input_trajectory_files[0].path
        converted_trajectory_file = File(converted_trajectory_filepath)
        input_trajectory_absolute_paths = [ trajectory_file.path for trajectory_file in input_trajectory_files ]

        # Convert input structure and trajectories to output structure and trajectory
        if not converted_structure_file.exists or not converted_trajectory_file.exists:
            print(' * Converting and merging')
            convert(
                input_structure_filepath = input_structure_file.path,
                output_structure_filepath = converted_structure_file.path,
                input_trajectory_filepaths = input_trajectory_absolute_paths,
                output_trajectory_filepath = converted_trajectory_file.path,
            )

        # Topologies are never converted, but they are kept in their original format

        # --- FILTERING ATOMS ------------------------------------------------------------

        # Find out if we need to filter
        # i.e. check if there is a selection filter and it matches some atoms
        converted_structure = Structure.from_pdb_file(converted_structure_file.path)
        must_filter = bool(self.project.filter_selection and converted_structure.select(self.project.filter_selection))

        # Set output filenames for the already filtered structure and trajectory
        # Note that this is the only step affecting topology and thus here we output the definitive topology
        filtered_structure_file = File(self.md_pathify(FILTERED_STRUCTURE)) if must_filter else converted_structure_file
        filtered_trajectory_file = File(self.md_pathify(FILTERED_TRAJECTORY)) if must_filter else converted_trajectory_file
        filtered_topology_file = output_topology_file if must_filter else input_topology_file
        
        # Filter atoms in structure, trajectory and topology if required and not done yet
        if must_filter and (not filtered_structure_file.exists or not filtered_trajectory_file.exists):
            print(' * Filtering atoms')
            filter_atoms(
                input_structure_file = converted_structure_file,
                input_trajectory_file = converted_trajectory_file,
                input_topology_file = input_topology_file, # We use input topology
                output_structure_file = filtered_structure_file,
                output_trajectory_file = filtered_trajectory_file,
                output_topology_file = filtered_topology_file, # We genereate the definitive topology
                filter_selection = self.project.filter_selection
            )

        # --- IMAGING AND FITTING ------------------------------------------------------------

        # There is no logical way to know if the trajectory is already imaged or it must be imaged
        # We rely exclusively in input flags
        must_image = self.project.image or self.project.fit

        # Set output filenames for the already filtered structure and trajectory
        imaged_structure_file = File(self.md_pathify(IMAGED_STRUCTURE)) if must_image else filtered_structure_file
        imaged_trajectory_file = File(self.md_pathify(IMAGED_TRAJECTORY)) if must_image else filtered_trajectory_file

        # Image the trajectory if it is required
        # i.e. make the trajectory uniform avoiding atom jumps and making molecules to stay whole
        # Fit the trajectory by removing the translation and rotation if it is required
        if must_image and (not imaged_structure_file.exists or not imaged_trajectory_file.exists):
            print(' * Imaging and fitting')
            image_and_fit(
                input_structure_file = filtered_structure_file,
                input_trajectory_file = filtered_trajectory_file,
                input_topology_file = filtered_topology_file, # This is optional if there are no PBC residues
                output_structure_file = imaged_structure_file,
                output_trajectory_file = imaged_trajectory_file,
                image = self.project.image,
                fit = self.project.fit,
                translation = self.project.translation,
                pbc_selection = self.project.pbc_selection
            )

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
        self._snapshots = get_frames_count(imaged_structure_file.path, imaged_trajectory_file.path)

        print(' * Correcting structure')

        # Correct the structure
        structure_corrector(
            input_structure_file = imaged_structure_file,
            input_trajectory_file = imaged_trajectory_file,
            input_topology_file = filtered_topology_file,
            output_structure_file = output_structure_file,
            output_trajectory_file = output_trajectory_file,
            snapshots = self._snapshots,
            register = self.register,
            mercy = self.project.mercy,
            trust = self.project.trust
        )

        # Now we must rename files to match the output file in case there is any missmatch
        # Some processed files may remain with some intermediate filename
        input_and_output_files = [
            (input_structure_file, imaged_structure_file, output_structure_file),
            (input_trajectory_files[0], imaged_trajectory_file, output_trajectory_file),
            (input_topology_file, filtered_topology_file, output_topology_file)
        ]
        for input_file, processed_file, output_file in input_and_output_files:
            if not output_file.exists:
                # There is also a chance that the input files have not been modified
                # This means the input format has already the output format and it is not to be imaged, fitted or corrected
                # However we need the output files to exist and we dont want to rename the original ones to conserve them
                # In order to not duplicate data, we will setup a symbolic link to the input files with the output filepaths
                if processed_file == input_file:
                    symlink(input_file.relative_path, output_file.path)
                # Some processed files may remain with some intermediate filename
                else:
                    rename(processed_file.relative_path, output_file.path)

        # Save the internal variables
        self._structure_file = output_structure_file
        self._trajectory_file = output_trajectory_file
        self.project._topology_file = output_topology_file

        # --- Cleanup intermediate files

        intermediate_filenames = [
            CONVERTED_STRUCTURE, CONVERTED_TRAJECTORY,
            FILTERED_STRUCTURE, FILTERED_TRAJECTORY,
            IMAGED_STRUCTURE, IMAGED_TRAJECTORY
        ]
        for filename in intermediate_filenames:
            file_path = self.md_pathify(filename)
            if exists(file_path):
                remove(file_path)

        # --- RUNNING FINAL TESTS ------------------------------------------------------------

        # Note that some tests have been run already
        # e.g. stable bonds is run in the structure corrector function

        # Note that tests here do not modify any file

        # Check the trajectory has not sudden jumps
        self.is_trajectory_integral()

        # Make a final summary
        print('Tests summary:')
        for test_name in AVAILABLE_CHECKINGS:
            test_result = self.register.tests.get(test_name)
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

    # Check if any of the available tests is missing or failed
    def any_missing_processing_tests (self) -> bool:
        for checking in AVAILABLE_CHECKINGS:
            test_result = self.register.tests.get(checking, None)
            if not test_result:
                return True
        return False

    # Get the processed structure
    def get_structure_file (self) -> str:
        # If we have a stored value then return it
        # This means we already found or generated this file
        if self._structure_file:
            return self._structure_file
        # Set the file
        structure_filepath = self.md_pathify(STRUCTURE_FILENAME)
        self._structure_file = File(structure_filepath)
        # If file does not exist then run the processing logic to generate it
        # Also run it if any of the tests is missing or failed since tests are run in the processing logic
        if not self._structure_file.exists or self.any_missing_processing_tests():
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
        # If file does not exist then run the processing logic to generate it
        # Also run it if any of the tests is missing or failed since tests are run in the processing logic
        if not self._trajectory_file.exists or self.any_missing_processing_tests():
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
        # Note that checking if the trajectory filename exists triggers all the processing logic
        # The processing logic is able to set the internal snapshots value as well so this avoid repeating the process
        if self.trajectory_file and self._snapshots != None:
            return self._snapshots
        # If we already have a value in the register cache then use it
        field_name = 'snapshots'
        if field_name in self.register.cache:
            return self.register.cache[field_name]
        print('-> Counting snapshots')
        # Otherwise we must find the value
        # This happens when the input files are already porcessed and thus we did not yet count the frames
        self._snapshots = get_frames_count(self.structure_file.path, self.trajectory_file.path)
        # Save the snapshots value in the register cache as well
        self.register.cache[field_name] = self._snapshots
        return self._snapshots
    snapshots = property(get_snapshots, None, None, "Trajectory snapshots (read only)")

    # Parsed structure
    def get_structure (self) -> 'Structure':
        # If we already have a stored value then return it
        if self._structure:
            return self._structure
        # Otherwise we must set the structure
        # Note that this is not only the mdtoolbelt structure, but it also contains additional logic
        self._structure = setup_structure(self.structure_file.path)
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
    def get_metadata_filename (self) -> str:
        # If the file already exists then send it
        metadata_filepath = self.md_pathify(OUTPUT_METADATA_FILENAME)
        if exists(metadata_filepath):
            return metadata_filepath
        print('-> Generating MD metadata')
        # Otherwise, generate it
        generate_md_metadata(
            md_inputs = self.md_inputs, # DANI: No sería mejor pasarle los inputs?
            structure = self.structure,
            snapshots = self.snapshots,
            register = self.register,
            output_metadata_filename = metadata_filepath,
        )
        return metadata_filepath
    metadata_filename = property(get_metadata_filename, None, None, "Project metadata filename (read only)")

    # The interactions
    # This is a bit exceptional since it is a value to be used and an analysis file to be generated
    def get_interactions (self) -> List[dict]:
        # If we already have a stored value then return it
        if self._interactions != None:
            return self._interactions
        print('-> Processing interactions')
        # Otherwise, process interactions
        self._interactions = process_interactions(
            input_interactions = self.project.input_interactions,
            structure_filename = self.structure_file.path,
            trajectory_filename = self.trajectory_file.path,
            structure = self.structure,
            snapshots = self.snapshots,
            interactions_file = self.md_pathify(OUTPUT_INTERACTIONS_FILENAME),
            mercy = self.project.mercy,
            frames_limit = 1000,
            interaction_cutoff = self.project.interaction_cutoff
        )
        return self._interactions
    interactions = property(get_interactions, None, None, "Interactions (read only)")

    # Indices of residues in periodic boundary conditions
    # Inherited from project
    def get_pbc_residues (self) -> List[int]:
        return self.project.pbc_residues
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
    def get_residues_map (self) -> dict:
        return self.project.residues_map
    residues_map = property(get_residues_map, None, None, "Residues mapping (read only)")

    # ---------------------------------------------------------------------------------
    # Tests
    # ---------------------------------------------------------------------------------

    # Sudden jumps test
    def is_trajectory_integral (self) -> Optional[bool]:
        # If we already have a stored value then return it
        if self._trajectory_integrity != None:
            return self._trajectory_integrity
        # If we are missing the inputs file then we are missing the PBC residues se we skip this analysis by now
        if not self.project.is_inputs_file_available():
            return None
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
            check_selection = PROTEIN_AND_NUCLEIC,
            standard_deviations_cutoff = self.project.rmsd_cutoff,
        )
        return self._trajectory_integrity

    # ---------------------------------------------------------------------------------
    # Analyses
    # ---------------------------------------------------------------------------------

    # RMSDs
    def generate_rmsds_analysis (self):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_RMSDS_FILENAME)
        if exists(output_analysis_filepath):
            return
        print('-> Running RMSDs analysis')
        # WARNING: This analysis is fast enought to use the full trajectory
        # WARNING: However, the output file size depends on the trajectory size
        # WARNING: In very long trajectories the number of points may make the client go slow when loading data
        rmsds(
            trajectory_file = self.trajectory_file,
            first_frame_file = self.first_frame_file,
            average_structure_file = self.average_structure_file,
            output_analysis_filename = output_analysis_filepath,
            frames_limit = 5000,
            snapshots = self.snapshots,
            structure = self.structure,
            pbc_residues = self.pbc_residues,
        )

    # TM scores
    def generate_tmscores_analysis (self):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_TMSCORES_FILENAME)
        if exists(output_analysis_filepath):
            return
        print('-> Running TM scores analysis')
        # Here we set a small frames limit since this anlaysis is a bit slow
        tmscores(
            input_topology_filename = self.structure_file.path,
            input_trajectory_filename = self.trajectory_file.path,
            output_analysis_filename = output_analysis_filepath,
            first_frame_filename = self.first_frame_file.path,
            average_structure_filename = self.average_structure_file.path,
            structure = self.structure,
            pbc_residues = self.pbc_residues,
            snapshots = self.snapshots,
            frames_limit = 200,
        )

    # RMSF, atom fluctuation
    def generate_rmsf_analysis (self):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_RMSF_FILENAME)
        if exists(output_analysis_filepath):
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
    def generate_rgyr_analysis (self):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_RGYR_FILENAME)
        if exists(output_analysis_filepath):
            return
        print('-> Running RGYR analysis')
        # WARNING: This analysis is fast enought to use the full trajectory instead of the reduced one
        # WARNING: However, the output file size depends on the trajectory size
        # WARNING: In very long trajectories the number of points may make the client go slow when loading data
        rgyr(
            input_topology_filename = self.structure_file.path,
            input_trajectory_filename = self.trajectory_file.path,
            output_analysis_filename = output_analysis_filepath,
            snapshots = self.snapshots,
            frames_limit = 5000,
            structure = self.structure,
            pbc_residues = self.pbc_residues,
        )

    # PCA, principal component analysis
    def generate_pca_analysis (self):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_PCA_FILENAME)
        if exists(output_analysis_filepath):
            return
        print('-> Running PCA analysis')
        # WARNING: This analysis will generate several output files
        # File 'pca.average.pdb' is generated by the PCA and it was used by the client but not anymore
        # File 'covar.log' is generated by the PCA but never used
        pca(
            input_topology_filename = self.structure_file.path,
            input_trajectory_filename = self.trajectory_file.path,
            output_analysis_filename = output_analysis_filepath,
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
    # def generate_pca_contacts (self):
    #     # Do not run the analysis if the output file already exists
    #     output_analysis_filepath = self.md_pathify(OUTPUT_PCA_CONTACTS_FILENAME)
    #     if exists(output_analysis_filepath):
    #         return
    #     print('-> Running PCA contacts analysis')
    #     pca_contacts(
    #         trajectory = self.trajectory_file.path,
    #         topology = self.pdb_filename,
    #         interactions = self.interactions,
    #         output_analysis_filename = output_analysis_filepath
    #     )

    # RMSD per residue
    def genereate_rmsd_perres_analysis (self):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_RMSD_PERRES_FILENAME)
        if exists(output_analysis_filepath):
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
    def genereate_rmsd_pairwise_analysis (self):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_RMSD_PAIRWISE_FILENAME)
        if exists(output_analysis_filepath):
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

    # Distance per residue
    def generate_dist_perres_analysis (self):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_DIST_PERRES_FILENAME)
        if exists(output_analysis_filepath):
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
    def generate_hbonds_analysis (self):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_HBONDS_FILENAME)
        if exists(output_analysis_filepath):
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
    def generate_sas_analysis (self):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_SASA_FILENAME)
        if exists(output_analysis_filepath):
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
    def generate_energies_analysis (self):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_ENERGIES_FILENAME)
        if exists(output_analysis_filepath):
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
    def generate_pockets_analysis (self):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_POCKETS_FILENAME)
        if exists(output_analysis_filepath):
            return
        print('-> Running pockets analysis')
        # Run the analysis
        pockets(
            structure_file = self.structure_file,
            trajectory_file = self.trajectory_file,
            pockets_prefix = self.md_pathify(OUTPUT_POCKET_STRUCTURES_PREFIX),
            output_analysis_filename = output_analysis_filepath,
            mdpocket_folder = self.md_pathify(POCKETS_FOLDER),
            structure = self.structure,
            pbc_residues = self.pbc_residues,
            snapshots = self.snapshots,
            frames_limit = 100,
        )

    # Helical parameters
    # DANI: Al final lo reimplementará Subamoy (en python) osea que esto no lo hacemos de momento
    # def generate_helical_analysis (self):
    #     # Do not run the analysis if the output file already exists
    #     output_analysis_filepath = self.md_pathify(OUTPUT_HELICAL_PARAMETERS_FILENAME)
    #     if exists(output_analysis_filepath):
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
    def generate_markov_analysis (self):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_MARKOV_FILENAME)
        if exists(output_analysis_filepath):
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
        inputs_filepath : str = DEFAULT_INPUTS_FILENAME,
        # The input topology filename
        # Multiple formats are accepted but the default is our own parsed json topology
        input_topology_filepath : str = None,
        # Input structure filepath
        # It may be both relative to the project directory or to every MD directory
        input_structure_filepath : str = STRUCTURE_FILENAME,
        # Input trajectory filepaths
        # These files are searched in every MD directory so the path MUST be relative
        input_trajectory_filepaths : str = [ TRAJECTORY_FILENAME ],
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
        self.url = None
        if database_url and accession:
            self.url = database_url + '/rest/current/projects/' + accession

        # Save inputs for the register, even if they are not used in the class
        self.inputs_filepath = inputs_filepath
        self.input_topology_filepath = input_topology_filepath


        # Set internal variables for input filenames
        # Set the inputs file
        self._inputs_file = File(inputs_filepath)
        # Set the input topology file
        self._input_topology_file = File(input_topology_filepath)
        if not self._input_topology_file:
            self._input_topology_file = File(self.guess_input_topology_filename())
        # Input structure and trajectory filepaths
        # Do not parse them to files yet, let this to the MD class
        self.input_structure_filepath = input_structure_filepath
        # WARNING: Input trajectory filepaths may be both a list or a single string
        if type(input_trajectory_filepaths) == list:
            self.input_trajectory_filepaths = input_trajectory_filepaths 
        elif type(input_trajectory_filepaths) == str:
            self.input_trajectory_filepaths = [ input_trajectory_filepaths ]
        else:
            raise InputError('Input trajectory filepaths must be a list of strings or a string')
        # Check no trajectory path is absolute
        for path in self.input_trajectory_filepaths:
            if path[0] == '/':
                raise InputError('Trajectory paths MUST be relative, not absolute (' + path + ')')
        # Input populations and transitions for MSM
        self.populations_filepath = populations_filepath
        self._populations_file = File(self.populations_filepath)
        self.transitions_filepath = transitions_filepath
        self._transitions_file = File(self.transitions_filepath)

        # Set the processed topology filepath, which depends on the input topology filename
        # Note that this file is different from the standard topology, although it may be standard as well
        self._topology_filepath = self.inherit_topology_filename()
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
        self.filter_selection = filter_selection
        # Fix the filter selection input, if needed
        # If a boolean is passed instead of a string then we set its corresponding value
        if type(filter_selection) == bool:
            if filter_selection:
                self.filter_selection = PROTEIN_AND_NUCLEIC
            else:
                self.filter_selection = None
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
        self._charges = None
        self._populations = None
        self._transitions = None
        self._residues_map = None
        self._mds = None

        # Set a new entry for the register
        # This is useful to track previous workflow runs and problems
        register_inputs = {}
        for input_name in REGISTER_INPUTS:
            register_inputs[input_name] = getattr(self, input_name)
        self.register = Register(inputs=register_inputs)

    # Check MD directories to be right
    # If there is any problem then directly raise an input error
    def check_md_directories (self):
        # Check there is at least one MD
        if len(self._md_directories) < 1:
            raise InputError('There must be at least one MD')
        # Check there are not duplicated MD directories
        if len(set(self._md_directories)) != len(self._md_directories):
            raise InputError('There are duplicated MD directories')

    # Set a function to get MD directories
    def get_md_directories (self) -> list:
        # If MD directories are already declared then return them
        if self._md_directories:
            return self._md_directories
        # Otherwise use the default MDs
        self._md_directories = []
        # Use the MDs from the inputs file when available
        if self.is_inputs_file_available():
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
                    'You can either declare them using the "-mdir" option or by providing and inputs file')
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
            raise Exception('Missing inputs file "' + self._inputs_file.filename + '"')
        # Download the inputs json file if it does not exists
        sys.stdout.write('Downloading inputs (' + self._inputs_file.filename + ')\n')
        inputs_url = self.url + '/inputs/'
        urllib.request.urlretrieve(inputs_url, self._inputs_file.path)
        # Rewrite the inputs file in a pretty formatted way
        with open(self._inputs_file.path, 'r+') as file:
            file_content = json.load(file)
            file.seek(0)
            json.dump(file_content, file, indent=4)
            file.truncate()
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

    # Get the input charges filename from the inputs
    # If the file is not found try to download it
    def get_input_topology_file (self) -> File:
        # There must be a topology filename
        if not self._input_topology_file:
            return None
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
                with open(self._input_topology_file.path, 'r+') as file:
                    file_content = json.load(file)
                    file.seek(0)
                    json.dump(file_content, file, indent=4)
                    file.truncate()
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
    input_topology_file = property(get_input_topology_file, None, None, "Input topology filename (read only)")

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
        with open(self.inputs_file.path, 'r') as file:
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
                raise InputError('Missing input "' + name + '"')
            return value
        return getter

    # Assign the getters
    input_interactions = property(input_getter('interactions'), None, None, "Interactions to be analyzed (read only)")
    pbc_selection = property(input_getter('pbc_selection'), None, None, "Selection of atoms which are still in periodic boundary conditions (read only)")
    forced_references = property(input_getter('forced_references'), None, None, "Uniprot IDs to be used first when aligning protein sequences (read only)")
    pdb_ids = property(input_getter('type'), None, None, "Protein Data Bank IDs used for the setup of the system (read only)")
    input_type = property(input_getter('type'), None, None, "Set if its a trajectory or an ensemble (read only)")
    input_mds = property(input_getter('mds'), None, None, "Input MDs configuration (read only)")
    input_reference_md_index = property(input_getter('mdref'), None, None, "Input MD reference index (read only)")

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
        filename = self._input_topology_file.filename
        if not filename:
            return None
        if filename == TOPOLOGY_FILENAME:
            return filename
        if filename == RAW_CHARGES_FILENAME:
            return filename
        standard_format = self._input_topology_file.format
        return 'topology.' + standard_format

    # Get the processed topology
    def get_topology_file (self) -> str:
        # If we have a stored value then return it
        # This means we already found or generated this file
        if self._topology_file:
            return self._topology_file
        # If the file already exists then we are done
        self._topology_file = File(self._topology_filepath)
        # If file does not exist then run the processing logic to generate it
        # Also run it if any of the tests is missing or failed since tests are run in the processing logic
        if not self._topology_file.exists or self.reference_md.any_missing_processing_tests():
            self.reference_md.process_input_files()
        # Now that the file is sure to exist we return it
        return self._topology_file
    topology_file = property(get_topology_file, None, None, "Topology filename (read only)")

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

    # Atom charges
    def get_charges (self) -> List[float]:
        # If we already have a stored value then return it
        if self._charges:
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
        print('-> Generating project metadata')
        # Otherwise, generate it
        generate_project_metadata(
            input_structure_filename = self.structure_file.path,
            input_trajectory_filename = self.trajectory_file.path,
            inputs_filename = self.inputs_file.filename, # DANI: No sería mejor pasarle los inputs?
            structure = self.structure,
            residues_map = self.residues_map,
            interactions = self.reference_md.interactions,
            register = self.register,
            output_metadata_filename = OUTPUT_METADATA_FILENAME,
        )
        return OUTPUT_METADATA_FILENAME
    metadata_filename = property(get_metadata_filename, None, None, "Project metadata filename (read only)")

    # Standard topology filename
    def get_standard_topology_file (self) -> str:
        # If we have a stored value then return it
        # This means we already found or generated this file
        if self._standard_topology_file:
            return self._standard_topology_file
        # If the file already exists then send it
        self._standard_topology_file = File(TOPOLOGY_FILENAME)
        if self._standard_topology_file.exists:
            return self._standard_topology_file
        print('-> Generating topology')
        # Otherwise, generate it
        generate_topology(
            structure = self.structure,
            charges = self.charges,
            residues_map = self.residues_map,
            pbc_residues = self.pbc_residues,
            output_topology_filepath = self._standard_topology_file.path
        )
        return self._standard_topology_file
    standard_topology_file = property(get_standard_topology_file, None, None, "Standard topology filename (read only)")

    # Screenshot filename
    def get_screenshot_filename (self) -> str:
        # If the file already exists then send it
        if exists(OUTPUT_SCREENSHOT_FILENAME):
            return OUTPUT_SCREENSHOT_FILENAME
        print('-> Generating screenshot')
        # Otherwise, generate it
        get_screenshot(
            input_structure_filename = self.structure_file.path,
            output_screenshot_filename = OUTPUT_SCREENSHOT_FILENAME,
        )
        return OUTPUT_SCREENSHOT_FILENAME
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
        with open(target_file.path, 'r') as file:
            return json.load(file)

# Set a function to convert an MD name into an equivalent MD directory
def name_2_directory (name : str) -> str:
    # Make all letters lower and replace white spaces by underscores
    directory = name.lower().replace(' ', '_')
    # Remove problematic characters
    for character in FORBIDEN_DIRECTORY_CHARACTERS:
        directory = directory.replace(character, '')
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
    'dist': MD.generate_dist_perres_analysis,
    'energies': MD.generate_energies_analysis,
    'hbonds': MD.generate_hbonds_analysis,
    #'helical': MD.generate_helical_analysis,
    'markov': MD.generate_markov_analysis,
    'pca': MD.generate_pca_analysis,
    #'pcacons': MD.generate_pca_contacts,
    'pockets': MD.generate_pockets_analysis,
    'rgyr': MD.generate_rgyr_analysis,
    'rmsds': MD.generate_rmsds_analysis,
    'perres': MD.genereate_rmsd_pairwise_analysis,
    'pairwise': MD.genereate_rmsd_perres_analysis,
    'rmsf': MD.generate_rmsf_analysis,
    'sas': MD.generate_sas_analysis,
    'tmscore': MD.generate_tmscores_analysis,
}

# List of requestables for the console
requestables = {
    **input_files,
    **processed_files,
    **analyses,
    'interactions': MD.get_interactions,
    'snapshots': MD.get_snapshots,
    'charges': Project.get_charges,
    'mapping': Project.get_residues_map,
    'screenshot': Project.get_screenshot_filename,
    'stopology': Project.get_standard_topology_file,
    'pmeta': Project.get_metadata_filename,
    'mdmeta': MD.get_metadata_filename
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
):

    # Move the current directory to the working directory
    chdir(working_directory)
    current_absolute_path = getcwd()
    current_directory_name = current_absolute_path.split('/')[-1]
    print('Running workflow for project at ' + current_directory_name)

    # Initiate the project project
    project = Project(**project_parameters)
    print('  ' + str(len(project.mds)) + ' MDs are to be run')

    # Now iterate over the different MDs
    for md in project.mds:

        print('\n' + CYAN_HEADER + 'Running workflow for MD at ' + md.directory + COLOR_END)

        # Set a function to call getters with the proper instance
        def call_getter (getter : Callable):
            instance = md if getter.__qualname__[0:3] == 'MD.' else project
            getter(instance)

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

        # Run the requested analyses
        if include and len(include) > 0:
            sys.stdout.write(f"Executing specific dependencies: " + ', '.join(include) + '\n')
            # Include only the specified dependencies
            requested = [ getter for name, getter in requestables.items() if name in include ]
            for getter in requested:
                call_getter(getter)
            # Exit here
            continue

        # Set the default requests, when there are not specific requests
        # Request all the analysis, the metadata, the standard topology and the screenshot
        requests = [
            *analyses.keys(),
            'mdmeta',
            'pmeta',
            'stopology',
            'screenshot'
        ]

        # If the exclude parameter was passed then remove excluded requests from the default requests
        if exclude and len(exclude) > 0:
            sys.stdout.write(f"Excluding specific dependencies: " + ', '.join(exclude) + '\n')
            requests = [ name for name in requests if name not in exclude ]
            
        # Run the requests
        for request in requests:
            getter = requestables[request]
            call_getter(getter)

        # Remove gromacs backups
        # DANI: Esto iría mejor en otro sitio
        remove_trash()

    sys.stdout.write("Done!\n")