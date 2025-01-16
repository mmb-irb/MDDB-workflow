#!/usr/bin/env python

# This is the starter script

# Import python libraries
from os import chdir, rename, remove, walk, mkdir, getcwd
from os.path import exists, isdir, isabs
import sys
import io
import re
import numpy
from glob import glob
from inspect import getfullargspec

# Constants
from model_workflow.utils.constants import *

# Import local tools
from model_workflow.tools.topology_manager import setup_structure
from model_workflow.tools.get_pytraj_trajectory import get_pytraj_trajectory
from model_workflow.tools.get_first_frame import get_first_frame
from model_workflow.tools.get_average import get_average
from model_workflow.tools.get_bonds import find_safe_bonds, get_bonds_canonical_frame
from model_workflow.tools.process_interactions import process_interactions
from model_workflow.tools.generate_metadata import generate_project_metadata, generate_md_metadata
from model_workflow.tools.generate_ligands_desc import generate_ligand_mapping
from model_workflow.tools.chains import generate_chain_references
from model_workflow.tools.generate_pdb_references import generate_pdb_references
from model_workflow.tools.residue_mapping import generate_residue_mapping
from model_workflow.tools.generate_map import generate_protein_mapping
from model_workflow.tools.generate_topology import generate_topology
from model_workflow.tools.get_charges import get_charges
from model_workflow.tools.remove_trash import remove_trash
from model_workflow.tools.get_screenshot import get_screenshot
from model_workflow.tools.filter_atoms import filter_atoms
from model_workflow.tools.image_and_fit import image_and_fit
from model_workflow.tools.structure_corrector import structure_corrector
from model_workflow.tools.fix_gromacs_masses import fix_gromacs_masses
from model_workflow.tools.check_inputs import check_inputs

# Import local utils
#from model_workflow.utils.httpsf import mount
from model_workflow.utils.auxiliar import InputError, warn, load_json, load_yaml, list_files, is_glob, parse_glob
from model_workflow.utils.register import Register
from model_workflow.utils.conversions import convert
from model_workflow.utils.structures import Structure
from model_workflow.utils.file import File
from model_workflow.utils.remote import Remote
from model_workflow.utils.pyt_spells import get_frames_count
from model_workflow.utils.type_hints import *

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
from model_workflow.analyses.helical_parameters import helical_parameters
from model_workflow.analyses.markov import markov

# Make the system output stream to not be buffered
# This is useful to make prints work on time in Slurm
# Otherwise, output logs are written after the script has fully run
# Note that this fix affects all modules and built-ins
unbuffered = io.TextIOWrapper(open(sys.stdout.fileno(), 'wb', 0), write_through=True)
sys.stdout = unbuffered

# Set an special exception for input errors
MISSING_INPUT_EXCEPTION = Exception('Missing input')

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
        self.remote = None
        if self.project.database_url and self.project.accession:
            self.accession = f'{self.project.accession}.{self.number}'
            self.remote = Remote(f'{self.project.database_url}/rest/current/projects/{self.accession}')
        # Save the directory
        self.directory = remove_final_slash(directory)
        # If the directory does not exists then create it
        if not exists(self.directory):
            mkdir(self.directory)
        # Save the input structure filepath
        # They may be relative to the project directory (unique) or relative to the MD directory (one per MD)
        # If the path is absolute then it is considered unique
        # If the file does not exist and it is to be downloaded then it is downloaded for each MD
        # Priorize the MD directory over the project directory
        self.input_structure_filepath = input_structure_filepath
        # Set the internal variable for the input structure file, to be assigned later
        self._input_structure_file = None
        # Save the input trajectory filepaths
        self.input_trajectory_filepaths = input_trajectory_filepaths
        # Set the internal variable for the input trajectory files, to be assigned later
        self._input_trajectory_files = None

        # Processed structure and trajectory files
        self._structure_file = None
        self._trajectory_file = None

        # Other values which may be found/calculated on demand
        self._md_inputs = None
        self._snapshots = None
        self._structure = None
        self._pytraj_topology = None
        self._processed_interactions = None
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
        return f'<MD ({len(self.structure.atoms)} atoms)>'

    # Given a filename or relative path, add the MD directory path at the beginning
    def md_pathify (self, filename_or_relative_path : str) -> str:
        return self.directory + '/' + filename_or_relative_path

    # Input structure file ------------

    # Set a function to get input structure file path
    def get_input_structure_filepath (self) -> str:
        # Set a function to find out if a path is relative to MD directories or to the project directory
        # To do so just check if the file exists in any of those
        # In case it exists in both or none then assume it is relative to MD directory
        # Parse glob notation in the process
        def relativize_and_parse_paths (input_path : str, may_not_exist : bool = False) -> Optional[str]:
            # Check if it is an absolute path
            if isabs(input_path):
                abs_glob_parse = parse_glob(input_path)
                # If we had multiple results then we complain
                if len(abs_glob_parse) > 1:
                    raise InputError(f'Multiple structures found with "{input_path}": {", ".join(abs_glob_parse)}')
                # If we had no results then we complain
                if len(abs_glob_parse) == 0:
                    if self.remote:
                        warn('Spread syntax is not supported to download remote files')
                    raise InputError(f'No structure found with "{input_path}"')
                abs_parsed_filepath = abs_glob_parse[0]
                return abs_parsed_filepath
            # Check the MD directory
            md_relative_filepath = self.md_pathify(input_path)
            md_glob_parse = parse_glob(md_relative_filepath)
            if len(md_glob_parse) > 1:
                raise InputError(f'Multiple structures found with "{input_path}": {", ".join(md_glob_parse)}')
            md_parsed_filepath = md_glob_parse[0] if len(md_glob_parse) == 1 else None
            if md_parsed_filepath and File(md_parsed_filepath).exists:
                return md_parsed_filepath
            # Check the project directory
            project_relative_filepath = self.project.project_pathify(input_path)
            project_glob_parse = parse_glob(project_relative_filepath)
            if len(project_glob_parse) > 1:
                raise InputError(f'Multiple structures found with "{input_path}": {", ".join(project_glob_parse)}')
            project_parsed_filepath = project_glob_parse[0] if len(project_glob_parse) == 1 else None
            if project_parsed_filepath and File(project_parsed_filepath).exists:
                return project_parsed_filepath
            # At this point we can conclude the input structure file does not exist
            # If we have no paths at all then it means a glob pattern was passed and it didn't match
            # Note that if a glob pattern existed then it would mean the file actually existed
            if len(md_glob_parse) == 0 and len(project_glob_parse) == 0:
                # Warn the user in case it was trying to use glob syntax to donwload remote files
                if self.remote:
                    warn('Spread syntax is not supported to download remote files')
                raise InputError('No trajectory file was reached neither in the project directory or MD directories in path(s) ' + ', '.join(input_path))
            # If the path does not exist anywhere then we asume it will be downloaded and set it relative to the MD
            # However make sure we have a remote
            # As an exception, if the 'may not exist' flag is passed then we return the result even if there is no remote
            if not may_not_exist and not self.remote:
                raise InputError(f'Cannot find a structure file by "{input_path}" anywhere')
            return md_parsed_filepath
        # If we have a value passed through command line
        if self.input_structure_filepath:
            # Find out if it is relative to MD directories or to the project directory
            return relativize_and_parse_paths(self.input_structure_filepath)
        # If we have a value passed through the inputs file has the value
        if self.project.is_inputs_file_available():
            # Get the input value, whose key must exist
            inputs_value = self.project.get_input('input_structure_filepath')
            # If there is a valid input then use it
            if inputs_value: return relativize_and_parse_paths(inputs_value)
        # If there is not input structure then asume it is the default
        # Check the default structure file exists or it may be downloaded
        default_structure_filepath = relativize_and_parse_paths(STRUCTURE_FILENAME, may_not_exist=True)
        default_structure_file = File(default_structure_filepath)
        # AGUS: si default_structure_filepath es None, default_structure_file será un objeto File y no se puede evaluar como None
        # AGUS: de esta forma al evaluar directamente si default_structure_filepath es None, se evita el error
        if default_structure_filepath is not None:
            if default_structure_file.exists or self.remote:
                return default_structure_filepath
        # If there is not input structure anywhere then use the input topology
        # We will extract the structure from it using a sample frame from the trajectory
        # Note that topology input filepath must exist and an input error will raise otherwise
        # However if we are using the standard topology file we can not extract the PDB from it (yet)
        if self.project.input_topology_file.filename != STANDARD_TOPOLOGY_FILENAME:
            return self.project.input_topology_file.path
        # If we can not use the topology either then surrender
        raise InputError('There is not input structure at all')

    # Get the input pdb filename from the inputs
    # If the file is not found try to download it
    def get_input_structure_file (self) -> str:
        # If the input structure file is already defined then return it
        if self._input_structure_file:
            return self._input_structure_file
        # Otherwise we must set it
        # First set the input structure filepath
        input_structure_filepath = self.get_input_structure_filepath()
        # Now set the input structure file
        self._input_structure_file = File(input_structure_filepath)
        # If the file already exists then return it
        if self._input_structure_file.exists:
            return self._input_structure_file
        # Try to download it
        # If we do not have the required parameters to download it then we surrender here
        if not self.remote:
            raise InputError(f'Missing input structure file "{self._input_structure_file.path}"')
        # Download the structure
        # If the structure filename is the standard structure filename then use the structure endpoint instead
        if self._input_structure_file.filename == STRUCTURE_FILENAME:
            self.remote.download_standard_structure(self._input_structure_file)
        # Otherwise download the input strucutre file by its filename
        else:
            self.remote.download_file(self._input_structure_file)
        return self._input_structure_file
    input_structure_file = property(get_input_structure_file, None, None, "Input structure filename (read only)")

    # Input trajectory filename ------------

    # Set a function to get input trajectory file paths
    def get_input_trajectory_filepaths (self) -> str:
        # Set a function to check and fix input trajectory filepaths
        # Also relativize paths to the current MD directory and parse glob notation
        def relativize_and_parse_paths (input_paths : List[str]) -> List[str]:
            checked_paths = input_paths
            # Input trajectory filepaths may be both a list or a single string
            # However we must keep a list
            if type(checked_paths) == list:
                pass 
            elif type(checked_paths) == str:
                checked_paths = [ checked_paths ]
            else:
                raise InputError('Input trajectory filepaths must be a list of strings or a string')
            # Make sure all or none of the trajectory paths are absolute
            abs_count = sum([ isabs(path) for path in checked_paths ])
            if not (abs_count == 0 or abs_count == len(checked_paths)):
                raise InputError('All trajectory frames must be relative or absolute. Mixing is not supported')
            # Set a function to glob-parse and merge all paths
            def parse_all_glob (paths : List[str]) -> List[str]:
                parsed_paths = []
                for path in paths:
                    parsed_paths += parse_glob(path)
                return parsed_paths
            # In case trajectory paths are absolute
            if abs_count > 0:
                absolute_parsed_paths = parse_all_glob(checked_paths)
                # Check we successfully defined some trajectory file
                if len(absolute_parsed_paths) == 0:
                    # Warn the user in case it was trying to use glob syntax to donwload remote files
                    if self.remote:
                        warn('Spread syntax is not supported to download remote files')
                    raise InputError('No trajectory file was reached neither in the project directory or MD directories in path(s) ' + ', '.join(input_paths))
                return absolute_parsed_paths
            # If trajectory paths are not absolute then check if they are relative to the MD directory
            # Get paths relative to the current MD directory
            md_relative_paths = [ self.md_pathify(path) for path in checked_paths ]
            # In case there are glob characters we must parse the paths
            md_parsed_paths = parse_all_glob(md_relative_paths)
            # Check we successfully defined some trajectory file
            if len(md_parsed_paths) > 0:
                # If so, check at least one of the files do actually exist
                if any([ File(path).exists for path in md_parsed_paths ]):
                    return md_parsed_paths
            # If no trajectory files where found then asume they are relative to the project
            # Get paths relative to the project directory
            project_relative_paths = [ self.project.project_pathify(path) for path in checked_paths ]
            # In case there are glob characters we must parse the paths
            project_parsed_paths = parse_all_glob(project_relative_paths)
            # Check we successfully defined some trajectory file
            if len(project_parsed_paths) > 0:
                # If so, check at least one of the files do actually exist
                if any([ File(path).exists for path in project_parsed_paths ]):
                    return project_parsed_paths
            # At this point we can conclude the input trajectory file does not exist
            # If we have no paths at all then it means a glob pattern was passed and it didn't match
            # Note that if a glob pattern existed then it would mean the file actually existed
            if len(md_parsed_paths) == 0 and len(project_parsed_paths) == 0:
                # Warn the user in case it was trying to use glob syntax to donwload remote files
                if self.remote:
                    warn('Spread syntax is not supported to download remote files')
                raise InputError('No trajectory file was reached neither in the project directory or MD directories in path(s) ' + ', '.join(input_paths))
            # If we have a path however it may be downloaded from the database if we have a remote
            if not self.remote:
                raise InputError(f'Cannot find anywhere a trajectory file with path(s) "{", ".join(input_paths)}"')
            # In this case we set the path as MD relative
            # Note that if input path was not glob based it will be both as project relative and MD relative
            if len(md_parsed_paths) == 0: raise ValueError('This should never happen')
            return md_parsed_paths
        # If we have a value passed through command line
        if self.input_trajectory_filepaths:
            return relativize_and_parse_paths(self.input_trajectory_filepaths)
        # Check if the inputs file has the value
        if self.project.is_inputs_file_available():
            # Get the input value
            inputs_value = self.project.get_input('input_trajectory_filepaths')
            if inputs_value:
                return relativize_and_parse_paths(inputs_value)
        # If no input trajectory is passed then asume it is the default
        default_trajectory_filepath = self.md_pathify(TRAJECTORY_FILENAME)
        default_trajectory_file = File(default_trajectory_filepath)
        if default_trajectory_file.exists or self.remote:
            return relativize_and_parse_paths([ TRAJECTORY_FILENAME ])
        # If there is no trajectory available then we surrender
        raise InputError('There is not input trajectory at all')

    # Get the input trajectory filename(s) from the inputs
    # If file(s) are not found try to download it
    def get_input_trajectory_files (self) -> str:
        # If we already defined input trajectory files then return them
        if self._input_trajectory_files != None:
            return self._input_trajectory_files
        # Otherwise we must set the input trajectory files
        input_trajectory_filepaths = self.get_input_trajectory_filepaths()
        self._input_trajectory_files = [ File(path) for path in input_trajectory_filepaths ]
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
        if not self.remote:
            missing_filepaths = [ trajectory_file.path for trajectory_file in missing_input_trajectory_files ]
            raise InputError('Missing input trajectory files: ' + ', '.join(missing_filepaths))
        # Download each trajectory file (ususally it will be just one)
        for trajectory_file in self._input_trajectory_files:
            # If this is the main trajectory (the usual one) then use the dedicated endpoint
            if trajectory_file.filename == TRAJECTORY_FILENAME:
                frame_selection = '1:10:1' if self.project.sample_trajectory else None
                self.remote.download_trajectory(trajectory_file, frame_selection=frame_selection, format='xtc')
            # Otherwise, download it by its filename
            else:
                self.remote.download_file(trajectory_file)
        return self._input_trajectory_files
    input_trajectory_files = property(get_input_trajectory_files, None, None, "Input trajectory filenames (read only)")

    # MD specific inputs
    def get_md_inputs (self) -> dict:
        # If we already have a value stored then return it
        if self._md_inputs:
            return self._md_inputs
        # Otherwise we must find its value
        # If we have MD inputs in the inputs file then use them
        if self.project.input_mds:
            # Iterate over the different MD inputs to find out each directory
            # We must find the MD inputs whcih belong to this specific MD according to this directory
            for md in self.project.input_mds:
                # Get the directory according to the inputs
                directory = md.get(MD_DIRECTORY, None)
                if directory:
                    check_directory(directory)
                # If no directory is specified in the inputs then guess it from the MD name
                else:
                    name = md['name']
                    directory = name_2_directory(name)
                # If the directory matches then this is our MD inputs
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
        if not self.remote:
            return False
        # Check if the file is among the available remote files
        # If it is no then stop here
        if target_file.filename not in self.remote.available_files:
            return False
        # Download the file
        self.remote.download_file(target_file)
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
        # If the file exists but it is new then we must run the tests as well
        elif self.register.is_file_new(output_structure_file):
            required_tests.update(STRUCTURE_TESTS)
        # If the structure exists and it is not new but it was modified since the last time then we also run the tests
        elif self.register.is_file_modified(output_structure_file):
            message = 'Structure was modified since the last processing'
            warn(message)
            required_tests.update(STRUCTURE_TESTS)

        # If there is no trajectory then we must run some tests
        if not output_trajectory_file.exists:
            required_tests.update(TRAJECTORY_TESTS)
        # If the file exists but it is new then we must run the tests as well
        elif self.register.is_file_new(output_trajectory_file):
            required_tests.update(TRAJECTORY_TESTS)
            self.register.reset_cache()
        # If the trajectory was modified since the last time then we must run these tests as well
        elif self.register.is_file_modified(output_trajectory_file):
            message = 'Trajectory was modified since the last processing'
            warn(message)
            required_tests.update(TRAJECTORY_TESTS)
            self.register.reset_cache()

        # If there is no topology then we must run some tests
        if not output_topology_file.exists:
            required_tests.update(TOPOLOGY_TESTS)
        # If the file exists but it is new then we must run the tests as well
        elif self.project.register.is_file_new(output_topology_file):
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

        # Extend the required tests list with the base tests which are to be run by default and are not already passed
        for checking in AVAILABLE_CHECKINGS:
            test_result = self.register.tests.get(checking, None)
            # If the test result is true then it menas ir has already passed
            # If the test result is 'na' then it means it is not aplicable
            if test_result == True or test_result == 'na': continue
            # If the test result is None then it means it has never been run
            # If the test result is false then it means it failed
            required_tests.update(checking)

        # Check if the processing parameters (filter, image, etc.) have changed since the last time
        # If so, then we must reset all tests and rerun the processing
        previous_processed_parameters = self.register.cache.get(PROCESSED, None)
        current_processed_parameters = {
            'filter': self.project.filter_selection,
            'image': self.project.image,
            'fit': self.project.fit,
        }
        # Note that not passing any of these parameters is condiered as 'leave it as it is'
        # This means if we already filtered and now there is no filter parameter then we consider there is no change
        for key, value in current_processed_parameters.items():
            if not value and previous_processed_parameters != None:
                current_processed_parameters[key] = previous_processed_parameters[key]
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
        if outputs_exist and len(required_tests) == 0 and same_processed_paramaters:
            return

        print('-> Processing input files')

        # --- FIRST CHECK -----------------------------------------------------------------------

        check_inputs(input_structure_file, input_trajectory_files, input_topology_file)

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
            # Set a provisional PBC selection
            # Note that we can not rely in the standard PBC selection since it is parsed through the structure
            # However we still don't have the standard structure available
            provisional_pbc_selection = self.input_pbc_selection if self.project.is_inputs_file_available() else 'guess'
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
                input_pbc_selection = provisional_pbc_selection
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
        # Othwerise count the number of snapshots
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

        # Set output filenames for the already filtered structure and trajectory
        corrected_structure_file = File(self.md_pathify(CORRECTED_STRUCTURE))
        corrected_trajectory_file = File(self.md_pathify(CORRECTED_TRAJECTORY))

        # Correct the structure
        # This function reads and or modifies the following MD variables:
        #   snapshots, safe_bonds, register, mercy, trust
        structure_corrector(
            input_structure_file = imaged_structure_file,
            input_trajectory_file = imaged_trajectory_file,
            input_topology_file = filtered_topology_file,
            output_structure_file = corrected_structure_file,
            output_trajectory_file = corrected_trajectory_file,
            MD = self
        )

        # If the corrected output exists then use it
        # Otherwise use the previous step files
        # Corrected files are generated only when changes are made in these files
        corrected_structure_file = corrected_structure_file if corrected_structure_file.exists else imaged_structure_file
        corrected_trajectory_file = corrected_trajectory_file if corrected_trajectory_file.exists else imaged_trajectory_file

        # Set for every type of file (structure, trajectory and topology) tte input, the las processed step and the output files
        input_and_output_files = [
            (input_structure_file, corrected_structure_file, output_structure_file),
            (input_trajectory_files[0], corrected_trajectory_file, output_trajectory_file),
            (input_topology_file, filtered_topology_file, output_topology_file)
        ]
        # Set a list of intermediate files
        intermediate_files = set([
            converted_structure_file, converted_trajectory_file,
            filtered_structure_file, filtered_trajectory_file,
            imaged_structure_file, imaged_trajectory_file,
        ])
        # Now we must rename files to match the output file
        # Processed files remain with some intermediate filename
        for input_file, processed_file, output_file in input_and_output_files:
            # If the processed file is already the output file then there is nothing to do here
            # This means it was already the input file and no changes were made
            if processed_file == output_file:
                continue
            # There is a chance that the input files have not been modified
            # This means the input format has already the output format and it is not to be imaged, fitted or corrected
            # However we need the output files to exist and we dont want to rename the original ones to conserve them
            # In order to not duplicate data, we will setup a symbolic link to the input files with the output filepaths
            if processed_file == input_file:
                # If output file exists and its the same as the input file, we can not create a symlink from a file to the same file
                if output_file.exists:
                    output_file.remove()
                output_file.set_symlink_to(input_file)
            # Otherwise rename the last intermediate file as the output file
            else:
                
                # In case the processed file is a symlink we must make sure the symlink is not made to a intermediate step
                # Intermediate steps will be removed further and thus the symlink would break
                # If the symlinks points to the input file there is no problem though
                if processed_file.is_symlink():
                    target_file = processed_file.get_symlink()
                    if target_file in intermediate_files:
                        target_file.rename_to(output_file)
                    else:
                        processed_file.rename_to(output_file)
                # If the files is not a symlink then simply rename it
                else:
                    processed_file.rename_to(output_file)


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
            elif test_result == 'na':
                test_nice_result = BLUE_HEADER + 'Not applicable' + COLOR_END
            else:
                raise ValueError()
            
            print(f' - {test_nice_name} -> {test_nice_result}')

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
            
        # --- Cleanup intermediate files

        # Set also a list of input files
        inputs_files = set([ input_structure_file, *input_trajectory_files, input_topology_file ])
        # We must make sure an intermediate file is not actually an input file before deleting it
        removable_files = intermediate_files - inputs_files
        # Now delete every removable file
        for removable_file in removable_files:
            # Note that a broken symlink does not 'exists'
            if removable_file.exists or removable_file.is_symlink():
                removable_file.remove()

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

    # The processed interactions
    # This is a bit exceptional since it is a value to be used and an analysis file to be generated
    def get_processed_interactions (self, overwrite : bool = False) -> List[dict]:
        # If we already have a stored value then return it
        if self._processed_interactions != None:
            return self._processed_interactions
        # Set the processed interactions file
        processed_interactions_filepath = self.md_pathify(OUTPUT_PROCESSED_INTERACTIONS_FILENAME)
        processed_interactions_file = File(processed_interactions_filepath)
        # If the file already exists then processed interactions will be read from it
        # If the overwrite argument is passed we must delete it here
        if processed_interactions_file.exists and overwrite:
            processed_interactions_file.remove()
        # Otherwise, process interactions
        self._processed_interactions = process_interactions(
            input_interactions = self.input_interactions,
            structure_file = self.structure_file,
            trajectory_file = self.trajectory_file,
            structure = self.structure,
            snapshots = self.snapshots,
            processed_interactions_file = processed_interactions_file,
            mercy = self.project.mercy,
            register = self.register,
            frames_limit = 1000,
            interaction_cutoff = self.project.interaction_cutoff
        )
        return self._processed_interactions
    processed_interactions = property(get_processed_interactions, None, None, "Processed interactions (read only)")

    def count_valid_interactions (self) -> int:
        valid_interactions = [ interaction for interaction in self.processed_interactions if not interaction.get('failed', False) ]
        return len(valid_interactions)
    valid_interactions_count = property(count_valid_interactions, None, None, "Count of non-failed processed_interactions (read only)")

    # Set a function to get input values which may be MD specific
    # If the MD input is missing then we use the project input
    def input_getter (name : str):
        # Set the getter
        def getter (self):
            # Get the MD input
            value = self.md_inputs.get(name, None)
            if value != None:
                return value
            # If there is no MD input then return the project value
            return getattr(self.project, f'input_{name}')
        return getter

    # Assign the MD input getters
    input_interactions = property(input_getter('interactions'), None, None, "Interactions to be analyzed (read only)")
    input_pbc_selection = property(input_getter('pbc_selection'), None, None, "Selection of atoms which are still in periodic boundary conditions (read only)")

    # Periodic boundary conditions atom selection
    def get_pbc_selection (self) -> List[int]:
        # If there is no inputs file then asume there are no PBC atoms and warn the user
        if not self.project.is_inputs_file_available():
            warn('Since there is no inputs file we guess PBC atoms as solvent, counter ions and lipids')
            return self.structure.select_pbc_guess()
        # Otherwise use the input value
        return self.structure.select(self.input_pbc_selection)
    pbc_selection = property(get_pbc_selection, None, None, "Periodic boundary conditions atom selection (read only)")

    # Indices of residues in periodic boundary conditions
    # WARNING: Do not inherit project pbc residues
    # WARNING: It may trigger all the processing logic of the reference MD when there is no need
    def get_pbc_residues (self) -> List[int]:
        # If we already have a stored value then return it
        if self.project._pbc_residues:
            return self.project._pbc_residues
        # If there is no inputs file then asume there are no PBC residues and warn the user
        if not self.pbc_selection:
            self.project._pbc_residues = []
            return self.project._pbc_residues
        # Otherwise we parse the selection and return the list of residue indices     
        self.project._pbc_residues = self.structure.get_selection_residue_indices(self.pbc_selection)
        print(f'PBC residues "{self.input_pbc_selection}" -> {len(self.project._pbc_residues)} residues')
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
            pbc_selection = self.pbc_selection,
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
            pbc_selection = self.pbc_selection,
            ligand_map = self.project.ligand_map,
        )

    # TM scores
    def run_tmscores_analysis (self, overwrite : bool = False):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_TMSCORES_FILENAME)
        if exists(output_analysis_filepath) and not overwrite:
            return
        # Here we set a small frames limit since this anlaysis is a bit slow
        tmscores(
            input_trajectory_file = self.trajectory_file,
            output_analysis_filename = output_analysis_filepath,
            first_frame_file = self.first_frame_file,
            average_structure_file = self.average_structure_file,
            structure = self.structure,
            pbc_selection = self.pbc_selection,
            snapshots = self.snapshots,
            frames_limit = 200,
        )

    # RMSF, atom fluctuation
    def run_rmsf_analysis (self, overwrite : bool = False):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_RMSF_FILENAME)
        if exists(output_analysis_filepath) and not overwrite:
            return
        # This analysis is fast and the output size depends on the number of atoms only
        # For this reason here it is used the whole trajectory with no frames limit
        rmsf(
            input_topology_filename = self.structure_file.path,
            input_trajectory_filename = self.trajectory_file.path,
            output_analysis_filename = output_analysis_filepath,
            pbc_selection = self.pbc_selection,
        )

    # RGYR, radius of gyration
    def run_rgyr_analysis (self, overwrite : bool = False):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_RGYR_FILENAME)
        if exists(output_analysis_filepath) and not overwrite:
            return
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
            pbc_selection = self.pbc_selection,
        )

    # PCA, principal component analysis
    def run_pca_analysis (self, overwrite : bool = False):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_PCA_FILENAME)
        if exists(output_analysis_filepath) and not overwrite:
            return
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
            pbc_selection = self.pbc_selection,
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
    #     pca_contacts(
    #         trajectory = self.trajectory_file.path,
    #         topology = self.pdb_filename,
    #         interactions = self.processed_interactions,
    #         output_analysis_filename = output_analysis_filepath
    #     )

    # RMSD per residue
    def run_rmsd_perres_analysis (self, overwrite : bool = False):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_RMSD_PERRES_FILENAME)
        if exists(output_analysis_filepath) and not overwrite:
            return
        # WARNING: This analysis is fast enought to use the full trajectory instead of the reduced one
        # WARNING: However, the output file size depends on the trajectory size. It may be pretty big
        rmsd_per_residue(
            input_topology_filename = self.structure_file.path,
            input_trajectory_filename = self.trajectory_file.path,
            output_analysis_filename = output_analysis_filepath,
            structure = self.structure,
            pbc_selection = self.pbc_selection,
            snapshots = self.snapshots,
            frames_limit = 100,
        )

    # RMSD pairwise
    def run_rmsd_pairwise_analysis (self, overwrite : bool = False):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_RMSD_PAIRWISE_FILENAME)
        if exists(output_analysis_filepath) and not overwrite:
            return
        # WARNING: This analysis is fast enought to use the full trajectory instead of the reduced one
        # WARNING: However, the output file size depends on the trajectory size exponentially. It may be huge
        rmsd_pairwise(
            input_topology_filename = self.structure_file.path,
            input_trajectory_filename = self.trajectory_file.path,
            output_analysis_filename = output_analysis_filepath,
            interactions = self.processed_interactions,
            structure = self.structure,
            pbc_selection = self.pbc_selection,
            snapshots = self.snapshots,
            frames_limit = 200,
            overall_selection = "name CA or name C5"
        )

    # Clusters
    def run_clusters_analysis (self, overwrite : bool = False):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_CLUSTERS_FILENAME)
        # The number of output analyses should be the same number of valid processed interactions
        # Otherwise we are missing some analysis
        if exists(output_analysis_filepath) and not overwrite:
            return
        # Set the output filepaths for all runs
        output_run_filepath = self.md_pathify(OUTPUT_CLUSTERS_RUNS_FILENAME)
        # Set the output filepaths for additional images generated in this analysis
        output_screenshot_filepath = self.md_pathify(OUTPUT_CLUSTER_SCREENSHOT_FILENAMES)
        # In case the overwirte argument is passed delete all already existing outputs
        if overwrite:
            for outputs in [ output_run_filepath, output_screenshot_filepath ]:
                existing_outputs = glob(outputs)
                for existing_output in existing_outputs:
                    if exists(existing_output):
                        remove(existing_output)
        # Run the analysis
        clusters_analysis(
            input_structure_file = self.structure_file,
            input_trajectory_file = self.trajectory_file,
            interactions = self.processed_interactions,
            structure = self.structure,
            snapshots = self.snapshots,
            pbc_selection = self.pbc_selection,
            output_analysis_filename = output_analysis_filepath,
            output_run_filepath = output_run_filepath,
            output_screenshots_filename = output_screenshot_filepath,
        )

    # Distance per residue
    def run_dist_perres_analysis (self, overwrite : bool = False):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_DIST_PERRES_FILENAME)
        if exists(output_analysis_filepath) and not overwrite:
            return
        # WARNING: This analysis is not fast enought to use the full trajectory. It would take a while
        distance_per_residue(
            input_topology_filename = self.structure_file.path,
            input_trajectory_filename = self.trajectory_file.path,
            output_analysis_filename = output_analysis_filepath,
            interactions = self.processed_interactions,
            snapshots = self.snapshots,
            frames_limit = 200,
        )

    # Hydrogen bonds
    def run_hbonds_analysis (self, overwrite : bool = False):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_HBONDS_FILENAME)
        if exists(output_analysis_filepath) and not overwrite:
            return
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
            interactions = self.processed_interactions,
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
        # Run the analysis
        energies(
            input_trajectory_file = self.trajectory_file,
            output_analysis_filename = output_analysis_filepath,
            energies_folder = self.md_pathify(ENERGIES_FOLDER),
            structure = self.structure,
            interactions = self.processed_interactions,
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
        # Run the analysis
        pockets(
            structure_file = self.structure_file,
            trajectory_file = self.trajectory_file,
            pockets_prefix = self.md_pathify(OUTPUT_POCKET_STRUCTURES_PREFIX),
            output_analysis_filepath = output_analysis_filepath,
            mdpocket_folder = self.md_pathify(POCKETS_FOLDER),
            pbc_selection = self.pbc_selection,
            snapshots = self.snapshots,
            frames_limit = 100,
        )

    # Helical parameters
    def run_helical_analysis (self, overwrite : bool = False):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_HELICAL_PARAMETERS_FILENAME)
        if exists(output_analysis_filepath) and not overwrite:
            return
        # Run the analysis
        helical_parameters(
            input_topology_filename = self.structure_file.path,
            input_trajectory_filename = self.trajectory_file.path,
            output_analysis_filename = output_analysis_filepath,
            structure = self.structure,
            frames_limit = None,
        )
        
    # Markov
    def run_markov_analysis (self, overwrite : bool = False):
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.md_pathify(OUTPUT_MARKOV_FILENAME)
        if exists(output_analysis_filepath) and not overwrite:
            return
        # If there are no populations file then stop here to avoid the log and calculating dependencies
        if not self.populations:
            return
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
        # Set an alternative MD configuration input
        md_config : Optional[list] = None,
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
        #nassa_config: str = DEFAULT_NASSA_CONFIG_FILENAME,
        # Set it we must download just a few frames instead of the whole trajectory
        sample_trajectory : bool = False,
    ):
        # Save input parameters
        self.directory = remove_final_slash(directory)
        self.database_url = database_url
        self.accession = accession
        # Set the project URL in case we have the required data
        self.remote = None
        if self.database_url and self.accession:
            self.remote = Remote(f'{self.database_url}/rest/current/projects/{self.accession}')

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
        # Note that even if the input topology path is passed we do not check it exists
        # Never forget we can donwload some input files from the database on the fly
        self.input_topology_filepath = input_topology_filepath
        self._input_topology_file = None
        # Input structure and trajectory filepaths
        # Do not parse them to files yet, let this to the MD class
        self.input_structure_filepath = input_structure_filepath
        self.input_trajectory_filepaths = input_trajectory_filepaths

        # Make sure the new MD configuration (-md) was not passed as well as old MD inputs (-mdir, -stru, -traj)
        if md_config and (md_directories or input_structure_filepath or input_trajectory_filepaths):
            raise InputError('MD configurations (-md) is not compatible with old MD inputs (-mdir, -stru, -traj)')
        # Save the MD configurations
        self.md_config = md_config
        # Make sure MD configuration has the correct format
        if self.md_config:
            # Make sure all MD configurations have at least 3 values each
            for mdc in self.md_config:
                if len(mdc) < 3:
                    raise InputError('Wrong MD configuration: the patter is -md <directory> <structure> <trajectory> <trajectory 2> ...')
            # Make sure there are no duplictaed MD directories
            md_directories = [ mdc[0] for mdc in self.md_config ]
            if len(md_directories) > len(set(md_directories)):
                raise InputError('There are duplicated MD directories')

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
        # Check input MDs are correct to far
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
        self._input_pbc_selection = pbc_selection
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
        self._pbc_residues = None
        self._safe_bonds = None
        self._charges = None
        self._populations = None
        self._transitions = None
        self._pdb_ids = None
        self._pdb_references = None
        self._protein_map = None
        self._ligand_map = None
        self.pubchem_name_list = None
        self._residue_map = None
        self._mds = None

        # Force a couple of extraordinary files which is generated if atoms are resorted
        self.resorted_bonds_file = File(self.project_pathify(RESORTED_BONDS_FILENAME))
        self.resorted_charges_file = File(self.project_pathify(RESORTED_CHARGES_FILENAME))

        # Set a new entry for the register
        # This is useful to track previous workflow runs and problems
        register_file = File(self.project_pathify(REGISTER_FILENAME))
        self.register = Register(register_file)

    # Given a filename or relative path, add the project directory path at the beginning
    def project_pathify (self, filename_or_relative_path : str) -> str:
        return self.directory + '/' + filename_or_relative_path

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
        if self.is_inputs_file_available() and self.input_mds:
            for input_md in self.input_mds:
                # Get the directory according to the inputs
                directory = input_md.get(MD_DIRECTORY, None)
                if directory:
                    check_directory(directory)
                # If no directory is specified in the inputs then guess it from the MD name
                else:
                    name = input_md['name']
                    if not name:
                        name = 'unnamed'
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
            self._reference_md_index = self.get_input('mdref')
        # Otherwise we simply set the first MD as the reference and warn the user about this
        if self._reference_md_index == None:
            warn('No reference MD was specified. The first MD will be used as reference.')
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
        # New system with MD configurations (-md)
        if self.md_config:
            for n, config in enumerate(self.md_config, 1):
                md = MD(
                    project = self, number = n, directory = config[0],
                    input_structure_filepath = config[1],
                    input_trajectory_filepaths = config[2:],
                )
                self._mds.append(md)
        # Old system (-mdir, -stru -traj)
        else:
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
        if self.remote:
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
        if not self.remote:
            raise InputError(f'Missing inputs file "{self._inputs_file.filename}"')
        # Download the inputs json file if it does not exists
        self.remote.download_inputs_file(self._inputs_file)
        return self._inputs_file
    inputs_file = property(get_inputs_file, None, None, "Inputs filename (read only)")

    # Topology filename ------------

    # If there is not input topology filepath, we must try to guess it among the files in the project directory
    # Note that if we can download from the remote then we must check the remote available files as well
    def guess_input_topology_filepath (self) -> Optional[str]:
        # Find the first supported topology file according to its name and format
        def find_first_accepted_topology_filename (available_filenames : List[str]) -> Optional[str]:
            for filename in available_filenames:
                # Make sure it is a valid topology file candidate
                # i.e. topology.xxx
                filename_splits = filename.split('.')
                if len(filename_splits) != 2 or filename_splits[0] != 'topology':
                    continue
                # Then make sure its format is among the acceoted topology formats
                extension = filename_splits[1]
                format = EXTENSION_FORMATS[extension]
                if format in ACCEPTED_TOPOLOGY_FORMATS:
                    return filename
            return None
        # First check among the local available files
        local_files = list_files(self.directory)
        accepted_topology_filename = find_first_accepted_topology_filename(local_files)
        if accepted_topology_filename:
            return self.project_pathify(accepted_topology_filename)
        # In case we did not find a topology among the local files, repeat the process with the remote files
        if self.remote:
            remote_files = self.remote.available_files
            accepted_topology_filename = find_first_accepted_topology_filename(remote_files)
            if accepted_topology_filename:
                return self.project_pathify(accepted_topology_filename)
        # If no actual topology is to be found then try with the standard topology instead
        # Check if the standard topology file is available
        # Note that we do not use standard_topology_file property to avoid generating it at this point
        standard_topology_filepath = self.project_pathify(STANDARD_TOPOLOGY_FILENAME)
        standard_topology_file = File(standard_topology_filepath)
        if standard_topology_file.exists:
            return standard_topology_filepath
        # If not we may also try to download the standard topology
        if self.remote:
            self.remote.download_standard_topology(standard_topology_file)
            return standard_topology_filepath
        # DEPRECATED: Find if the raw charges file is present as a last resource
        if exists(RAW_CHARGES_FILENAME):
            return RAW_CHARGES_FILENAME
        # If we did not find any valid topology filepath at this point then return None
        return None

    # Get the input topology filepath from the inputs or try to guess it
    def get_input_topology_filepath (self) -> File:
        # Set a function to parse possible glob notation
        def parse (filepath : str) -> str:
            # If there is no glob pattern then just return the string as is
            if not is_glob(filepath):
                return filepath
            # If there is glob pattern then parse it
            parsed_filepaths = glob(filepath)
            if len(parsed_filepaths) == 0:
                # Warn the user in case it was trying to use glob syntax to donwload remote files
                if self.remote:
                    warn('Spread syntax is not supported to download remote files')
                raise InputError(f'No topologies found with "{filepath}"')
            if len(parsed_filepaths) > 1:
                raise InputError(f'Multiple topologies found with "{filepath}": {", ".join(parsed_filepaths)}')
            return parsed_filepaths[0]
        # If this value was passed through command line then it would be set as the internal value already
        if self.input_topology_filepath:
            return parse(self.input_topology_filepath)
        # Check if the inputs file has the value
        if self.is_inputs_file_available():
            # Get the input value, whose key must exist
            inputs_value = self.get_input('input_topology_filepath')
            # If there is a valid input then use it
            if inputs_value:
                return parse(inputs_value)
        # Otherwise we must guess which is the topology file
        guess = self.guess_input_topology_filepath()
        if guess:
            return guess
        # If nothing worked then surrender
        raise InputError('Missing input topology file path')

    # Get the input topology file
    # If the file is not found try to download it
    def get_input_topology_file (self) -> File:
        # If we already have a value then return it
        if self._input_topology_file:
            return self._input_topology_file
        # Set the input topology filepath
        input_topology_filepath = self.get_input_topology_filepath()
        # If no input is passed then we check the inputs file
        # Set the file
        self._input_topology_file = File(input_topology_filepath)
        # If the file already exists then we are done
        if self._input_topology_file.exists:
            return self._input_topology_file
        # Try to download it
        # If we do not have the required parameters to download it then we surrender here
        if not self.remote:
            raise InputError(f'Missing input topology file "{self._input_topology_file.filename}"')
        # Otherwise, try to download it using the files endpoint
        # Note that this is not usually required
        self.remote.download_file(self._input_topology_file)
        # In case the topology is a '.top' file we consider it is a Gromacs topology
        # It may come with additional itp files we must download as well
        if self._input_topology_file.format == 'top':
            # Find available .itp files and download each of them
            itp_filenames = [filename for filename in self.remote.available_files if filename[-4:] == '.itp']
            for itp_filename in itp_filenames:
                itp_filepath = self.project_pathify(itp_filename)
                itp_file = File(itp_filepath)
                self.remote.download_file(itp_file)
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
        # Legacy fixes
        old_pdb_ids = self._inputs.get('pdbIds', None)
        if old_pdb_ids:
            self._inputs['pdb_ids'] = old_pdb_ids
        # Finally return the updated inputs
        return self._inputs
    inputs = property(get_inputs, None, None, "Inputs from the inputs file (read only)")

    # Then set getters for every value in the inputs file

    # Get a specific 'input' value
    # Handle a possible missing keys
    def get_input (self, name: str):
        value = self.inputs.get(name, MISSING_INPUT_EXCEPTION)
        # If we had a value then return it
        if value != MISSING_INPUT_EXCEPTION:
            return value
        # If the field is not specified in the inputs file then set a defualt value
        default_value = DEFAULT_INPUT_VALUES.get(name, None)
        # Warn the user about this
        warn(f'Missing input "{name}" -> Using default value: {default_value}')
        return default_value

    # Set a function to get a specific 'input' value by its key/name
    # Note that we return the getter function but we do not call it just yet
    def input_getter (name : str):
        def getter (self):
            return self.get_input(name)
        return getter

    # Assign the getters
    input_interactions = property(input_getter('interactions'), None, None, "Interactions to be analyzed (read only)")
    forced_references = property(input_getter('forced_references'), None, None, "Uniprot IDs to be used first when aligning protein sequences (read only)")
    input_pdb_ids = property(input_getter('pdb_ids'), None, None, "Protein Data Bank IDs used for the setup of the system (read only)")
    input_type = property(input_getter('type'), None, None, "Set if its a trajectory or an ensemble (read only)")
    input_mds = property(input_getter('mds'), None, None, "Input MDs configuration (read only)")
    input_ligands = property(input_getter('ligands'), None, None, "Input ligand references (read only)")
    
    # PBC selection may come from the console or from the inputs file
    # Console has priority over the inputs file
    def get_input_pbc_selection (self) -> Optional[str]:
        # If we have an internal value then return it
        if self._input_pbc_selection:
            return self._input_pbc_selection
        # As an exception, we avoid asking for the inputs file if it is not available
        # This input is required for some early processing steps where we do not need the inputs file for anything else
        if not self.is_inputs_file_available():
            return None
        # Otherwise, find it in the inputs
        # Get the input value, whose key must exist
        self._input_pbc_selection = self.get_input('pbc_selection')
        return self._input_pbc_selection
    input_pbc_selection = property(get_input_pbc_selection, None, None, "Selection of atoms which are still in periodic boundary conditions (read only)")

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
        return self.reference_md.pbc_residues
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
        # Set if stable bonds have to be checked
        must_check_stable_bonds = STABLE_BONDS_FLAG not in self.trust
        # If this analysis has been already passed then we can trust structure bonds
        if self.register.tests.get(STABLE_BONDS_FLAG, None) == True:
            must_check_stable_bonds = False
        # Otherwise we must find safe bonds value
        # This should only happen if we are working with already processed files
        self._safe_bonds = find_safe_bonds(
            input_topology_file=self.topology_file,
            input_structure_file=self.structure_file,
            input_trajectory_file=self.trajectory_file,
            must_check_stable_bonds=must_check_stable_bonds,
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

    # Tested and standarized PDB ids
    def get_pdb_ids (self) -> List[str]:
        # If we already have a stored value then return it
        if self._pdb_ids != None:
            return self._pdb_ids
        # Otherwise test and standarize input PDB ids
        self._pdb_ids = []
        # If there is no input pdb ids (may be None) then stop here
        if not self.input_pdb_ids:
            return []
        # Iterate input PDB ids
        for input_pdb_id in self.input_pdb_ids:
            # First make sure this is a PDB id
            if not re.match(PDB_ID_FORMAT, input_pdb_id):
                raise InputError(f'Input PDB id "{input_pdb_id}" does not look like a PDB id')
            # Make letters upper
            pdb_id = input_pdb_id.upper()
            self._pdb_ids.append(pdb_id)
        return self._pdb_ids
    pdb_ids = property(get_pdb_ids, None, None, "Tested and standarized PDB ids (read only)")

    # PDB references
    def get_pdb_references (self) -> List[dict]:
        # If we already have a stored value then return it
        if self._pdb_references:
            return self._pdb_references
        # Set the PDB references file
        pdb_references_filepath = self.project_pathify(PDB_REFERENCES_FILENAME)
        pdb_references_file = File(pdb_references_filepath)
        # Otherwise we must find the value
        self._pdb_references = generate_pdb_references(
            pdb_ids = self.pdb_ids,
            pdb_references_file = pdb_references_file
        )
        return self._pdb_references
    pdb_references = property(get_pdb_references, None, None, "PDB references (read only)")

    # Define the PDB references output file
    def get_pdb_references_file (self, overwrite : bool = False) -> File:
        # Set the PDB references file
        pdb_references_filepath = self.project_pathify(PDB_REFERENCES_FILENAME)
        pdb_references_file = File(pdb_references_filepath)
        # If the file already exists then return it
        # However if the overwrite argument is passed then delete it and proceed to produce it again
        if pdb_references_file.exists:
            if not overwrite:
                return pdb_references_file
            pdb_references_file.remove()
        # Ask for the PDB references thus producing the PDB references file
        self.get_pdb_references()
        return pdb_references_file
    pdb_references_file = property(get_pdb_references_file, None, None, "File including PDB refereces data (read only)")

    # Protein residues mapping
    def get_protein_map (self) -> List[dict]:
        # If we already have a stored value then return it
        if self._protein_map:
            return self._protein_map
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

    # Get chain references
    def get_chain_references (self, overwrite : bool = False) -> List[str]:
        # Set the chains references file
        chains_references_filepath = self.project_pathify(OUTPUT_CHAINS_FILENAME)
        chains_references_file = File(chains_references_filepath)
        # If the file already exists and the overwrite option is passed then remove it
        if chains_references_file.exists and overwrite:
            chains_references_file.remove()
        # Call the function to generate the chain references
        chains = generate_chain_references(
            structure = self.structure,
            chains_references_file = chains_references_file,
            #chain_name=self.structure.chain_name,
        )
        return chains
    chains_data = property(get_chain_references, None, None, "Chain (read only)")

    # Ligand residues mapping
    def get_ligand_map (self) -> List[dict]:
        # If we already have a stored value then return it
        if self._ligand_map != None:
            return self._ligand_map
        # Set the ligand references file
        ligand_references_filepath = self.project_pathify(LIGAND_REFERENCES_FILENAME)
        ligand_references_file = File(ligand_references_filepath)
        # Otherwise we must find the value
        self._ligand_map, self.pubchem_name_list = generate_ligand_mapping(
            structure = self.structure,
            register = self.register,
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
        # Set an input getter that gets the input as soon as called
        def get_input (name : str):
            return Project.input_getter(name)(self)
        # Otherwise, generate it
        generate_project_metadata(
            input_structure_filename = self.structure_file.path,
            input_trajectory_filename = self.trajectory_file.path,
            get_input = get_input,
            structure = self.structure,
            residue_map = self.residue_map,
            protein_references_file = self.protein_references_file,
            pdb_ids = self.pdb_ids,
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
        standard_topology_filepath = self.project_pathify(STANDARD_TOPOLOGY_FILENAME)
        self._standard_topology_file = File(standard_topology_filepath)
        # If the file already exists and it is not to be overwirtten then send it
        if self._standard_topology_file.exists and not overwrite:
            return self._standard_topology_file
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
    for character in FORBIDDEN_DIRECTORY_CHARACTERS:
        directory = directory.replace(character, '')
    return directory

# Set a function to check for problematic characters in a directory path
def check_directory (directory : str) -> str:
    # Remove problematic characters
    for character in FORBIDDEN_DIRECTORY_CHARACTERS:
        if character in directory:
            raise InputError(f'Directory path "{directory}" includes the forbidden character "{character}"')

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

# Project input files
project_input_files = {
    'itopology': Project.get_input_topology_file,
    'inputs': Project.get_inputs_file,
    'populations': Project.get_populations_file,
    'transitions': Project.get_transitions_file
}
# MD input files
md_input_files = {
    'istructure': MD.get_input_structure_file,
    'itrajectory': MD.get_input_trajectory_files
}
# Both project and MD input files
input_files = { **project_input_files, **md_input_files }

# Project processed files
project_processed_files = {
    'topology': Project.get_topology_file
}
# MD processed files
md_processed_files = {
    'structure': MD.get_structure_file,
    'trajectory': MD.get_trajectory_file
}
# Both project and MD processed files
processed_files = { **project_processed_files, **md_processed_files }

# List of available analyses
analyses = {
    'clusters': MD.run_clusters_analysis,
    'dist': MD.run_dist_perres_analysis,
    'energies': MD.run_energies_analysis,
    'hbonds': MD.run_hbonds_analysis,
    'helical': MD.run_helical_analysis,
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

# Project requestable tasks
project_requestables = {
    **project_input_files,
    **project_processed_files,
    'pdbs': Project.get_pdb_references_file,
    'mapping': Project.get_protein_references_file,
    'ligands': Project.get_ligand_references_file,
    'screenshot': Project.get_screenshot_filename,
    'stopology': Project.get_standard_topology_file,
    'pmeta': Project.get_metadata_file,
    'chains': Project.get_chain_references,
}
# MD requestable tasks
md_requestables = {
    **md_input_files,
    **md_processed_files,
    **analyses,
    'interactions': MD.get_processed_interactions,
    'mdmeta': MD.get_metadata_file,
}
# List of requestables for the console
requestables = { **project_requestables, **md_requestables }

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
        raise InputError('Include (-i) and exclude (-e) are not compatible. Use one of these options.')

    # Make sure the working directory exists
    if not exists(working_directory):
        raise InputError(f'Working directory "{working_directory}" does not exist')

    # Make sure the working directory is actually a directory
    if not isdir(working_directory):
        raise InputError(f'Working directory "{working_directory}" is actually not a directory')

    # Move the current directory to the working directory
    chdir(working_directory)
    current_directory_name = getcwd().split('/')[-1]
    print(f'\n{CYAN_HEADER}Running workflow for project at {current_directory_name}{COLOR_END}')

    # Initiate the project project
    project = Project(**project_parameters)
    print(f'  {len(project.mds)} MDs are to be run')

    # Set the tasks to be run
    tasks = None
    # If the download argument is passed then just make sure input files are available
    if download:
        tasks = list(input_files.keys())
    # If the setup argument is passed then just process input files
    elif setup:
        tasks = list(processed_files.keys())
    # If the include argument then add only the specified tasks to the list
    elif include and len(include) > 0:
        tasks = include
    # Set the default tasks otherwise
    else:
        tasks = [
            # Project tasks
            'stopology',
            'screenshot',
            'pmeta',
            'pdbs',
            'chains',
            # MD tasks
            'mdmeta',
            'interactions',
            *analyses.keys(),
        ]
        # If the exclude parameter was passed then remove excluded tasks from the default tasks
        if exclude and len(exclude) > 0:
            tasks = [ name for name in tasks if name not in exclude ]

    # If the user requested to overwrite something, make sure it is in the tasks list
    if overwrite and type(overwrite) == list:
        for task in overwrite:
            if task not in tasks:
                raise InputError(f'Task "{task}" is to be overwriten but it is not in the tasks list. Either include it or do not exclude it')

    # Run the project tasks now
    project_tasks = [ task for task in tasks if task in project_requestables ]
    for task in project_tasks:
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
            getter(project, overwrite=True)
        # Call the function with no additional arguments otherwise
        else:
            getter(project)

    # Get the MD tasks
    md_tasks = [ task for task in tasks if task in md_requestables ]

    # If there are no MD tasks then we are done already
    if len(md_tasks) == 0:
        print("Finished!")
        return

    # Now iterate over the different MDs
    for md in project.mds:
        print(f'\n{CYAN_HEADER} Processing MD at {md.directory}{COLOR_END}')

        # Run the MD tasks
        for task in md_tasks:
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
                getter(md, overwrite=True)
            # Call the function with no additional arguments otherwise
            else:
                getter(md)

        # Remove gromacs backups and other trash files from this MD
        remove_trash(md.directory)

    # Remove gromacs backups and other trash files from the project
    remove_trash(project.directory)

    print("Done!")