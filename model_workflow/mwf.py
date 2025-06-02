#!/usr/bin/env python

# This is the starter script

# Import python libraries
from os import chdir, rename, remove, walk, mkdir, getcwd
from os.path import exists, isdir, isabs, relpath
from shutil import rmtree
import sys
import io
import re
import numpy
from glob import glob
from inspect import getfullargspec

# Constants
from model_workflow.utils.constants import *

# Import local utils
# Importing constants first is important
from model_workflow.utils.constants import *
#from model_workflow.utils.httpsf import mount
from model_workflow.utils.auxiliar import InputError, MISSING_TOPOLOGY
from model_workflow.utils.auxiliar import warn, load_json, load_yaml, list_files, is_directory_empty
from model_workflow.utils.auxiliar import is_glob, parse_glob, glob_filename, purge_glob
from model_workflow.utils.arg_cksum import get_cksum_id
from model_workflow.utils.register import Register
from model_workflow.utils.cache import Cache
from model_workflow.utils.conversions import convert
from model_workflow.utils.structures import Structure
from model_workflow.utils.topologies import Topology
from model_workflow.utils.file import File
from model_workflow.utils.remote import Remote
from model_workflow.utils.pyt_spells import get_frames_count, get_pytraj_trajectory
from model_workflow.utils.selections import Selection
from model_workflow.utils.type_hints import *

# Import local tools
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
from model_workflow.tools.generate_membrane_mapping import generate_membrane_mapping
from model_workflow.tools.generate_topology import generate_topology
from model_workflow.tools.get_charges import get_charges
from model_workflow.tools.remove_trash import remove_trash
from model_workflow.tools.get_screenshot import get_screenshot
from model_workflow.tools.filter_atoms import filter_atoms
from model_workflow.tools.image_and_fit import image_and_fit
from model_workflow.tools.structure_corrector import structure_corrector
from model_workflow.tools.fix_gromacs_masses import fix_gromacs_masses
from model_workflow.tools.check_inputs import check_inputs

# Import local analyses
from model_workflow.analyses.rmsds import rmsds
from model_workflow.analyses.tmscores import tmscores
from model_workflow.analyses.rmsf import rmsf
from model_workflow.analyses.rgyr import rgyr
from model_workflow.analyses.pca import pca
from model_workflow.analyses.density import density
from model_workflow.analyses.thickness import thickness
from model_workflow.analyses.area_per_lipid import area_per_lipid
from model_workflow.analyses.lipid_order import lipid_order
from model_workflow.analyses.lipid_interactions import lipid_interactions
#from model_workflow.analyses.pca_contacts import pca_contacts
from model_workflow.analyses.rmsd_per_residue import rmsd_per_residue
from model_workflow.analyses.rmsd_pairwise import rmsd_pairwise
from model_workflow.analyses.clusters import clusters_analysis
from model_workflow.analyses.distance_per_residue import distance_per_residue
from model_workflow.analyses.hydrogen_bonds import hydrogen_bonds
from model_workflow.analyses.sasa import sasa
from model_workflow.analyses.energies import energies
from model_workflow.analyses.dihedral_energies import compute_dihedral_energies
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

# Set a special exception for missing inputs
MISSING_INPUT_EXCEPTION = Exception('Missing input')

# Set a special exception for missing task function arguments
# This is used for easy debug when a new functions is added wrongly
MISSING_ARGUMENT_EXCEPTION = Exception('Missing argument')

# Set a special exception for when a value is missing
MISSING_VALUE_EXCEPTION = Exception('Missing value')

# Name of the argument used by all functions to know where to write output
OUTPUT_DIRECTORY_ARG = 'output_directory'

# Set the name of the section in the cache where argument value cksums are saved
CACHE_ARG_CKSUMS = 'arg_cksums'

# Run a fix in gromacs if not done before
# Note that this is run always at the moment the code is read, no matter the command or calling origin
fix_gromacs_masses()

# Set some variables which are filled at the end but are referred by previously defined functions
requestables = {}
inverted_requestables = {}

# Set a class to handle a generic task
class Task:
    def __init__ (self,
        # The task flag
        # This name is to be used by the include/exclude/overwrite arguments
        # It will also name the folder containing all analysis output
        flag : str,
        # The task nice name
        # This user-firendly name is to be used in the logs
        name : str,
        # The task function
        # Function argument names must correspond with Project/MD property names
        func : Callable,
        # The task function "additional" inputs
        # Project/MD properties are automatically sent to the function as arguments
        # However some analyses have additional arguments (e.g. frames limit, cutoffs, etc.)
        args : dict,
    ):
        # Save input arguments
        self.flag = flag
        self.name = name
        self.func = func
        self.args = args
        # Set internal values
        self._value = MISSING_VALUE_EXCEPTION

    # When a task is called
    def __call__(self, parent):
        # First of all check if this task has been already done in this very run
        # If so then return the stored vale
        if self._value != MISSING_VALUE_EXCEPTION: return self._value
        # Process the task function arguments
        processed_args = {}
        # Get the task function expected arguments
        specification = getfullargspec(self.func)
        expected_arguments = specification.args
        n_default_arguments = len(specification.defaults) if specification.defaults else 0
        # Find out which arguments are optional since they have default values
        default_arguments = set(expected_arguments[::-1][:n_default_arguments])
        # Find out if the function is to return output
        returns_output = bool(specification.annotations.get('return', None))
        # If one of the argument expected outputs is the output_directory then set it here
        # We will set a new directory with the flag name of the task, in the correspoding path
        # Note that while the task is beeing done the output directory has a different name
        # Thus the directory is hidden and marked as incomplete
        # The final output directory is the one without the incomplete prefix
        writes_output = OUTPUT_DIRECTORY_ARG in expected_arguments
        incomplete_output_directory = None
        final_output_directory = None
        if writes_output:
            # Set the output directory path
            incomplete_output_directory = parent.pathify(INCOMPLETE_PREFIX + self.flag)
            final_output_directory = incomplete_output_directory.replace(INCOMPLETE_PREFIX, '')
            # Add it to the processed args
            processed_args[OUTPUT_DIRECTORY_ARG] = incomplete_output_directory
            # Remove the expected argument from the list
            expected_arguments.remove(OUTPUT_DIRECTORY_ARG)
        # Iterate the reamining expected arguments
        for arg in expected_arguments:
            # First find the argument among the parent properties
            arg_value = self.find_arg_value(arg, parent, default_arguments)
            if arg_value == MISSING_ARGUMENT_EXCEPTION: continue
            # Add the processed argument
            processed_args[arg] = arg_value
        # Check if this dependency is to be overwriten
        forced_overwrite = self.flag in parent.overwritables
        # Get the list of inputs which have changed compared to a previous run
        # WARNING: Always get changed inputs, since this function updates the cache
        # If had_cache is false then it means this is the first time the task is ever done
        changed_inputs, had_cache = self.get_changed_inputs(parent, processed_args)
        any_input_changed = len(changed_inputs) > 0
        # We must overwrite outputs either if inputs changed or if it was forced by the user
        must_overwrite = forced_overwrite or any_input_changed
        # Check if output already exists
        # If the final directory already exists then it means the task was started in a previous run
        existing_incomplete_output = writes_output and exists(incomplete_output_directory)
        # If the final directory already exists then it means the task was done in a previous run
        existing_final_output = writes_output and exists(final_output_directory)
        # It should never happend that both directories exist
        if existing_incomplete_output and existing_final_output:
            raise RuntimeError('We have incomplete and finished output directories at the same time')
        # If we must overwrite then purge previous outputs
        if must_overwrite:
            if existing_incomplete_output: rmtree(incomplete_output_directory)
            if existing_final_output: rmtree(final_output_directory)
        # If the output already exists and it is not to be overwritten
        elif existing_final_output:
            # If we must proceed because the analysis returns output then use the final output directory
            if returns_output:
                processed_args[OUTPUT_DIRECTORY_ARG] = final_output_directory
            # If the task is not to return any output then we can skip the task
            else: 
                print(f'{GREY_HEADER}-> Task {self.flag} ({self.name}) already completed{COLOR_END}')
                return
        # Create the output directory, if necessary
        missing_incomplete_output = writes_output \
            and not exists(incomplete_output_directory) \
            and not exists(final_output_directory)
        if missing_incomplete_output: mkdir(incomplete_output_directory)
        # Finally call the function
        print(f'{GREEN_HEADER}-> Running task {self.flag} ({self.name}){COLOR_END}')
        # If the task is to be run again because an inputs changed then let the user know
        if any_input_changed and had_cache and not forced_overwrite:
            changes = ''.join([ '\n   - ' + inp for inp in changed_inputs ])
            print(f'{GREEN_HEADER}   The task is run again since the following inputs changed:{changes}{COLOR_END}')
        self._value = self.func(**processed_args)
        # Update the overwritables so this is not remade further in the same run
        parent.overwritables.discard(self.flag)
        # As a brief cleanup, if the output directory is empty at the end, then remove it
        # Otheriwse, change the incomplete directory name to its final name
        if writes_output and exists(incomplete_output_directory):
            if is_directory_empty(incomplete_output_directory): rmtree(incomplete_output_directory)
            else: rename(incomplete_output_directory, final_output_directory)
        # Now return the function result
        return self._value

    # Find argument values, thus running any dependency
    def find_arg_value (self, arg : str, parent : Union['Project', 'MD'], default_arguments : set):
        # First find the argument among the parent properties
        arg_value = getattr(parent, arg, MISSING_ARGUMENT_EXCEPTION)
        if arg_value != MISSING_ARGUMENT_EXCEPTION: return arg_value
        # If the parent is an MD then it may happen the property is from the Project
        if isinstance(parent, MD):
            arg_value = getattr(parent.project, arg, MISSING_ARGUMENT_EXCEPTION)
            if arg_value != MISSING_ARGUMENT_EXCEPTION: return arg_value
        # If the property is missing then search among the additional arguments
        arg_value = self.args.get(arg, MISSING_ARGUMENT_EXCEPTION)
        if arg_value != MISSING_ARGUMENT_EXCEPTION: return arg_value
        # It may also happen that the argument has a default value
        # If this is the case then we can skip it
        if arg in default_arguments: return MISSING_ARGUMENT_EXCEPTION
        # NEVER FORGET: Function arguments must have the same name that the Project/MD property
        # If the argument is still missing then you programmed the function wrongly or...
        # You may have forgotten the additional argument in the task args
        raise RuntimeError(f'Function "{self.func.__name__}" expects argument "{arg}" but it is missing')
    
    # Find out if inputs changed regarding the last run
    def get_changed_inputs (self,
        parent : Union['Project', 'MD'],
        processed_args : dict) -> Tuple[ List[str], bool ]:
        # Get cache argument references
        all_cksums = parent.cache.retrieve(CACHE_ARG_CKSUMS, {})
        task_cksums = all_cksums.get(self.func.__name__, None)
        had_cache = False if task_cksums == None else True
        if task_cksums == None:
            task_cksums = {}
            all_cksums[self.func.__name__] = task_cksums
        # Check argument by argument
        # Keep a list with arguments which have changed
        unmatched_arguments = []
        for arg_name, arg_value in processed_args.items():
            # Skip the output directory argument
            # Changes in this argument are not actual changes
            if arg_name == OUTPUT_DIRECTORY_ARG: continue
            # Get the cksum from the new argument value
            new_cksum = get_cksum_id(arg_value)
            # Retrieve the cksum from the old argument value
            old_cksum = task_cksums.get(arg_name, None)
            # Compare new and old cksums
            if new_cksum != old_cksum:
                # If we found a missmatch then add it to the list
                unmatched_arguments.append(arg_name)
                # Update the references
                task_cksums[arg_name] = new_cksum
        # If we found no missmatch then stop here
        if len(unmatched_arguments) > 0:
            # If there were differences then update the cache
            parent.cache.update(CACHE_ARG_CKSUMS, all_cksums)
        return unmatched_arguments, had_cache

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
            self.remote = Remote(self.project.database_url, self.accession)
        # Save the directory
        # If it is an absolute then make it relative to the project
        if isabs(directory):
            # This function already removes the final slash
            self.directory = relpath(directory, self.project.directory)
        # Otherwise save it as is but just removing the final slash (if any)
        else:
            self.directory = remove_final_slash(directory)
        self.directory = self.project.pathify(self.directory)
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
        register_filepath = self.pathify(REGISTER_FILENAME)
        register_file = File(register_filepath)
        if register_file.path == self.project.register.file.path:
            self.register = self.project.register
        else:
            self.register = Register(register_file)

        # Set a new MD specific cache
        # In case the directory is the project directory itself, use the project cache
        cache_filepath = self.pathify(CACHE_FILENAME)
        cache_file = File(cache_filepath)
        if cache_file.path == self.project.cache.file.path:
            self.cache = self.project.cache
        else:
            self.cache = Cache(cache_file)

        # Set tasks whose output is to be overwritten
        self.overwritables = set()

    def __repr__ (self):
        return 'MD'
    
    # This function is able to find its caller "self" function
    # Then it finds its associated label in the requestables
    def _get_task (self) -> str:
        caller_data = sys._getframe().f_back
        caller_name = caller_data.f_code.co_name
        caller_func = getattr(self, caller_name).__func__
        return inverted_requestables[caller_func]

    # Given a filename or relative path, add the MD directory path at the beginning
    def pathify (self, filename_or_relative_path : str) -> str:
        return self.directory + '/' + filename_or_relative_path

    # Input structure file ------------

    def get_input_structure_filepath (self) -> str:
        """Set a function to get input structure file path"""
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
            md_relative_filepath = self.pathify(input_path)
            md_glob_parse = parse_glob(md_relative_filepath)
            if len(md_glob_parse) > 1:
                raise InputError(f'Multiple structures found with "{input_path}": {", ".join(md_glob_parse)}')
            md_parsed_filepath = md_glob_parse[0] if len(md_glob_parse) == 1 else None
            if md_parsed_filepath and File(md_parsed_filepath).exists:
                return md_parsed_filepath
            # Check the project directory
            project_relative_filepath = self.project.pathify(input_path)
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

    def get_input_structure_file (self) -> str:
        """Get the input pdb filename from the inputs.
        If the file is not found try to download it."""
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

    def get_input_trajectory_filepaths (self) -> str:
        """Set a function to get input trajectory file paths."""
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
            md_relative_paths = [ self.pathify(path) for path in checked_paths ]
            # In case there are glob characters we must parse the paths
            md_parsed_paths = parse_all_glob(md_relative_paths)
            # Check we successfully defined some trajectory file
            if len(md_parsed_paths) > 0:
                # If so, check at least one of the files do actually exist
                if any([ File(path).exists for path in md_parsed_paths ]):
                    return md_parsed_paths
            # If no trajectory files where found then asume they are relative to the project
            # Get paths relative to the project directory
            project_relative_paths = [ self.project.pathify(path) for path in checked_paths ]
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
        default_trajectory_filepath = self.pathify(TRAJECTORY_FILENAME)
        default_trajectory_file = File(default_trajectory_filepath)
        if default_trajectory_file.exists or self.remote:
            return relativize_and_parse_paths([ TRAJECTORY_FILENAME ])
        # If there is no trajectory available then we surrender
        raise InputError('There is not input trajectory at all')

    def get_input_trajectory_files (self) -> str:
        """Get the input trajectory filename(s) from the inputs.
        If file(s) are not found try to download it."""
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
                frame_selection = f'1:{self.project.sample_trajectory}:1' if self.project.sample_trajectory else None
                self.remote.download_trajectory(trajectory_file, frame_selection=frame_selection, format='xtc')
            # Otherwise, download it by its filename
            else:
                self.remote.download_file(trajectory_file)
        return self._input_trajectory_files
    input_trajectory_files = property(get_input_trajectory_files, None, None, "Input trajectory filenames (read only)")

    def get_md_inputs (self) -> dict:
        """MD specific inputs."""
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

    def get_file (self, target_file : File) -> bool:
        """Check if a file exists. If not, try to download it from the database.
        If the file is not found in the database it is fine, we do not even warn the user.
        Note that this function is used to get populations and transitions files, which are not common."""
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

    def process_input_files (self):
        """Process input files to generate the processed files.
        This process corrects and standarizes the topology, the trajectory and the structure."""
        # Set the input filepaths
        input_structure_file = self.input_structure_file
        input_trajectory_files = self.input_trajectory_files
        input_topology_file = self.project.input_topology_file

        # Set the output filepaths
        output_structure_filepath = self.pathify(STRUCTURE_FILENAME)
        output_structure_file = File(output_structure_filepath)
        output_trajectory_filepath = self.pathify(TRAJECTORY_FILENAME)
        output_trajectory_file = File(output_trajectory_filepath)
        output_topology_filepath = self.project.topology_filepath
        output_topology_file = File(self.pathify(output_topology_filepath)) if output_topology_filepath else MISSING_TOPOLOGY

        # If all output files already exist we may skip the processing
        topology_already_processed = output_topology_file == MISSING_TOPOLOGY or output_topology_file.exists
        outputs_exist = output_structure_file.exists and output_trajectory_file.exists and topology_already_processed

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
            self.cache.reset()
        # If the trajectory was modified since the last time then we must run these tests as well
        elif self.register.is_file_modified(output_trajectory_file):
            message = 'Trajectory was modified since the last processing'
            warn(message)
            required_tests.update(TRAJECTORY_TESTS)
            self.cache.reset()

        # In case there is a topology to be processed...
        if output_topology_file != MISSING_TOPOLOGY:
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
            required_tests.update([checking])

        # Check if the processing parameters (filter, image, etc.) have changed since the last time
        # If so, then we must reset all tests and rerun the processing
        previous_processed_parameters = self.cache.retrieve(PROCESSED)
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

        # Make sure we do not enter in a loop
        # This may happen when we read/call an output value/file by mistake
        if hasattr(self, '_processed'): raise RuntimeError('Looped processing')
        self._processed = True

        print('-> Processing input files')

        # --- FIRST CHECK -----------------------------------------------------------------------

        check_inputs(input_structure_file, input_trajectory_files, input_topology_file)

        # --- CONVERTING AND MERGING ------------------------------------------------------------

        # Set the output format for the already converted structure
        input_structure_format = self.input_structure_file.format
        output_structure_format = output_structure_file.format
        converted_structure_filepath = self.pathify(CONVERTED_STRUCTURE)
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
        converted_trajectory_filepath = self.pathify(CONVERTED_TRAJECTORY)
        # If input trajectory already matches the output format and is unique then avoid the renaming
        if input_trajectories_format == output_trajectory_format and len(input_trajectory_files) == 1:
            converted_trajectory_filepath = input_trajectory_files[0].path
        converted_trajectory_file = File(converted_trajectory_filepath)
        # Join all input trajectory paths
        input_trajectory_paths = [ trajectory_file.path for trajectory_file in input_trajectory_files ]

        # Set an intermeidate file for the trajectory while it is being converted
        # This prevents using an incomplete trajectory in case the workflow is suddenly interrupted while converting
        incompleted_converted_trajectory_filepath = self.pathify(INCOMPLETE_PREFIX + CONVERTED_TRAJECTORY)
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

        # --- provisional reference structure ---

        # Now that we MUST have a PDB file we can set a provisional structure instance
        # Note that this structure is not yet corrected so it must be used with care
        # Otherwise we could have silent errors
        provisional_structure = Structure.from_pdb_file(converted_structure_file.path)
        # Now we can set a provisional coarse grain selection
        # This selection is useful to avoid problems with CG atom elements
        # Since this is proviosonal we will make it silent
        provisional_cg_selection = self._set_cg_selection(provisional_structure, verbose=False)
        for atom_index in provisional_cg_selection.atom_indices:
            provisional_structure.atoms[atom_index].element = CG_ATOM_ELEMENT

        # --- FILTERING ATOMS ------------------------------------------------------------

        # Find out if we need to filter
        # i.e. check if there is a selection filter and it matches some atoms
        must_filter = bool(self.project.filter_selection)

        # Set output filenames for the already filtered structure and trajectory
        # Note that this is the only step affecting topology and thus here we output the definitive topology
        filtered_structure_file = File(self.pathify(FILTERED_STRUCTURE)) if must_filter else converted_structure_file
        filtered_trajectory_file = File(self.pathify(FILTERED_TRAJECTORY)) if must_filter else converted_trajectory_file
        filtered_topology_file = output_topology_file if must_filter else input_topology_file

        # Set an intermeidate file for the trajectory while it is being filtered
        # This prevents using an incomplete trajectory in case the workflow is suddenly interrupted while filtering
        incompleted_filtered_trajectory_filepath = self.pathify(INCOMPLETE_PREFIX + FILTERED_TRAJECTORY)
        incompleted_filtered_trajectory_file = File(incompleted_filtered_trajectory_filepath)
        # If there is an incomplete trajectory then remove it
        if incompleted_filtered_trajectory_file.exists:
            incompleted_filtered_trajectory_file.remove()

        # Check if any output file is missing
        missing_filter_output = not filtered_structure_file.exists or not filtered_trajectory_file.exists

        # Check if parameters have changed
        # Note that for this specific step only filtering is important
        previous_filtered_parameters = self.cache.retrieve(FILTERED)
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
                reference_structure = provisional_structure,
                filter_selection = self.project.filter_selection,
            )
            # Once filetered, rename the trajectory file as completed
            rename(incompleted_filtered_trajectory_file.path, filtered_trajectory_file.path)
            # Update the cache
            self.cache.update(FILTERED, current_filtered_parameters)

        # --- provisional reference structure ---

        # Now that we have a filtered PDB file we have to update provisional structure instance
        # Note that this structure is not yet corrected so it must be used with care
        # Otherwise we could have silent errors
        provisional_structure = Structure.from_pdb_file(filtered_structure_file.path)
        # Again, set the coarse grain atoms
        # Since elements may be needed to guess PBC selection we must solve them right before
        # Since this is proviosonal we will make it silent
        provisional_cg_selection = self._set_cg_selection(provisional_structure, verbose=False)
        for atom_index in provisional_cg_selection.atom_indices:
            provisional_structure.atoms[atom_index].element = CG_ATOM_ELEMENT
        # Also we can set a provisional PBC selection
        # This selection is useful both for imaging/fitting and for the correction
        # We will make sure that the provisonal and the final PBC selections match
        # Since this is proviosonal we will make it silent
        provisional_pbc_selection = self._set_pbc_selection(provisional_structure, verbose=False)

        # --- IMAGING AND FITTING ------------------------------------------------------------

        # There is no logical way to know if the trajectory is already imaged or it must be imaged
        # We rely exclusively in input flags
        must_image = self.project.image or self.project.fit

        # Set output filenames for the already filtered structure and trajectory
        imaged_structure_file = File(self.pathify(IMAGED_STRUCTURE)) if must_image else filtered_structure_file
        imaged_trajectory_file = File(self.pathify(IMAGED_TRAJECTORY)) if must_image else filtered_trajectory_file

        # Set an intermeidate file for the trajectory while it is being imaged
        # This prevents using an incomplete trajectory in case the workflow is suddenly interrupted while imaging
        incompleted_imaged_trajectory_filepath = self.pathify(INCOMPLETE_PREFIX + IMAGED_TRAJECTORY)
        incompleted_imaged_trajectory_file = File(incompleted_imaged_trajectory_filepath)
        # If there is an incomplete trajectory then remove it
        if incompleted_imaged_trajectory_file.exists:
            incompleted_imaged_trajectory_file.remove()

        # Check if any output file is missing
        missing_imaged_output = not imaged_structure_file.exists or not imaged_trajectory_file.exists

        # Check if parameters have changed
        # Note that for this step the filter parameters is also important
        previous_imaged_parameters = self.cache.retrieve(IMAGED)
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
                structure = provisional_structure,
                pbc_selection = provisional_pbc_selection
            )
            # Once imaged, rename the trajectory file as completed
            rename(incompleted_imaged_trajectory_file.path, imaged_trajectory_file.path)
            # Update the cache
            self.cache.update(IMAGED, current_imaged_parameters)
            # Update the provisional strucutre coordinates
            imaged_structure = Structure.from_pdb_file(imaged_structure_file.path)
            imaged_structure_coords = [ atom.coords for atom in imaged_structure.atoms ]
            provisional_structure.set_new_coordinates(imaged_structure_coords)

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
        # If we already have a value in the cache then use it
        cached_snapshots = self.cache.retrieve(SNAPSHOTS_FLAG)
        if cached_snapshots != None:
            self._snapshots = cached_snapshots
        # Othwerise count the number of snapshots
        else:
            self._snapshots = get_frames_count(imaged_structure_file, imaged_trajectory_file)
            # Save the snapshots value in the cache as well
            self.cache.update(SNAPSHOTS_FLAG, self._snapshots)

        # WARNING:
        # We may need to resort atoms in the structure corrector function
        # In such case, bonds and charges must be resorted as well and saved apart to keep values coherent
        # Bonds are calculated during the structure corrector but atom charges must be extracted no
        self.project._charges = get_charges(filtered_topology_file)

        print(' * Correcting structure')

        # Set output filenames for the already filtered structure and trajectory
        corrected_structure_file = File(self.pathify(CORRECTED_STRUCTURE))
        corrected_trajectory_file = File(self.pathify(CORRECTED_TRAJECTORY))

        # Correct the structure
        # This function reads and or modifies the following MD variables:
        #   snapshots, safe_bonds, register, cache, mercy, trust
        structure_corrector(
            structure = provisional_structure,
            input_trajectory_file = imaged_trajectory_file,
            input_topology_file = filtered_topology_file,
            output_structure_file = corrected_structure_file,
            output_trajectory_file = corrected_trajectory_file,
            MD = self,
            pbc_selection = provisional_pbc_selection
        )

        # If the corrected output exists then use it
        # Otherwise use the previous step files
        # Corrected files are generated only when changes are made in these files
        corrected_structure_file = corrected_structure_file if corrected_structure_file.exists else imaged_structure_file
        corrected_trajectory_file = corrected_trajectory_file if corrected_trajectory_file.exists else imaged_trajectory_file

        # Set for every type of file (structure, trajectory and topology) the input, the last processed step and the output files
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
        if output_topology_file != MISSING_TOPOLOGY:
            self.project.register.update_mtime(output_topology_file)

        # Update the parameters used to get the last processed structure and trajectory files
        self.cache.update(PROCESSED, current_processed_parameters)

        # --- Definitive PBC selection ---

        # Now that we have the corrected structure we can set the definitive PBC atoms
        # Make sure the selection is identical to the provisional selection
        if self.pbc_selection != provisional_pbc_selection:
            raise InputError('PBC selection is not consistent after correcting the structure. '
                'Please consider using a different PBC selection. '
                'Avoid relying in atom distances or elements to avoid this problem.')

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
                # Remove previous warnings
                self.register.remove_warnings(test_name)
                # Get test pretty name
                test_nice_name = NICE_NAMES[test_name]
                # Issue the corresponding warning            
                self.register.add_warning(test_name, test_nice_name + ' was skipped and never run before')
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
        structure_filepath = self.pathify(STRUCTURE_FILENAME)
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
        trajectory_filepath = self.pathify(TRAJECTORY_FILENAME)
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
        # If we already have a value in the cache then use it
        cached_value = self.cache.retrieve(SNAPSHOTS_FLAG)
        if cached_value != None:
            return cached_value
        # Otherwise we must find the value
        # This happens when the input files are already porcessed and thus we did not yet count the frames
        self._snapshots = get_frames_count(self.structure_file, self.trajectory_file)
        # Save the snapshots value in the cache as well
        self.cache.update(SNAPSHOTS_FLAG, self._snapshots)
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
        # Make sure the structure file exists at this point
        if not self.structure_file.exists:
            raise ValueError('Trying to set standard structure but file '
                f'{self.structure_file.path} does not exist yet. Are you trying '
                'to access the standard structure before processing input files?')
        # Note that this is not only the structure class, but it also contains additional logic
        self._structure = Structure.from_pdb_file(self.structure_file.path)
        # If the stable bonds test failed and we had mercy then it is sure our structure will have wrong bonds
        # In order to make it coherent with the topology we will mine topology bonds from here and force them in the structure
        # If we fail to get bonds from topology then just go along with the default structure bonds
        if not self.register.tests.get(STABLE_BONDS_FLAG, None):
            self._structure.bonds = self.safe_bonds
        # Same procedure if we have coarse grain atoms
        elif self.cg_selection:
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
        first_frame_filepath = self.pathify(FIRST_FRAME_FILENAME)
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
        average_structure_filepath = self.pathify(AVERAGE_STRUCTURE_FILENAME)
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

    def get_metadata_file (self) -> File:
        """Generate the MD metadata file."""
        # Get the task name
        task = self._get_task()
        # Check if this dependency is to be overwriten
        must_overwrite = task in self.overwritables
        # Update the overwritables so this is not remade further in the same run
        self.overwritables.discard(task)
        # Set the metadata file
        metadata_filepath = self.pathify(OUTPUT_METADATA_FILENAME)
        metadata_file = File(metadata_filepath)
        # If the file already exists
        if metadata_file.exists and not must_overwrite:
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
    get_processed_interactions = Task('inter', 'Interaccions processing',
        process_interactions, { 'frames_limit': 1000 })
    interactions = property(get_processed_interactions, None, None, "Processed interactions (read only)")

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
    input_cg_selection = property(input_getter('cg_selection'), None, None, "Selection of atoms which are not actual atoms but coarse grain beads (read only)")

    # Internal function to set PBC selection
    # It may parse the inputs file selection string if it is available or guess it otherwise
    def _set_pbc_selection (self, reference_structure : 'Structure', verbose : bool = True) -> 'Selection':
        # Otherwise we must set the PBC selection
        if verbose: print('Setting Periodic Boundary Conditions (PBC) atoms selection')
        selection_string = None
        # If there is inputs file then get the input pbc selection
        if self.project.is_inputs_file_available():
            if verbose: print(' Using selection string in the inputs file')
            selection_string = self.input_pbc_selection
        # If there is no inputs file we guess PBC atoms automatically
        else:
            if verbose: print(' No inputs file -> Selection string will be set automatically')
            selection_string = 'auto'
        # Parse the selection string using the reference structure
        parsed_selection = None
        # If the input PBC selection is 'auto' then guess it automatically
        if selection_string == 'auto':
            # To guess PBC atoms (with the current implementation) we must make sure ther eis no CG
            if reference_structure.has_cg():
                raise InputError('We can not guess PBC atoms in CG systems. Please set PBC atoms manually.\n'
                    ' Use the "-pbc" argument or set the inputs file "pbc_selection" field.')
            if verbose: print(' Guessing PBC atoms as solvent, counter ions and lipids')
            parsed_selection = reference_structure.select_pbc_guess()
        # If we have a valid input value then use it
        elif selection_string:
            if verbose: print(f' Selecting PBC atoms "{selection_string}"')
            parsed_selection = reference_structure.select(selection_string)
            if not parsed_selection:
                raise InputError(f'PBC selection "{selection_string}" selected no atoms')
        # If we have an input value but it is empty then we set an empty selection
        else:
            if verbose: print(' No PBC atoms selected')
            parsed_selection = Selection()
        # Log a few of the selected residue names
        if verbose and parsed_selection:
            print(f' Parsed PBC selection has {len(parsed_selection)} atoms')
            selected_residues = reference_structure.get_selection_residues(parsed_selection)
            selected_residue_names = list(set([ residue.name for residue in selected_residues ]))
            limit = 3 # Show a maximum of 3 residue names
            example_residue_names = ', '.join(selected_residue_names[0:limit])
            if len(selected_residue_names) > limit: example_residue_names += ', etc.'
            print('  e.g. ' + example_residue_names)
        return parsed_selection

    # Periodic boundary conditions atom selection
    def get_pbc_selection (self) -> 'Selection':
        # If we already have a stored value then return it
        if self.project._pbc_selection != None:
            return self.project._pbc_selection
        # Otherwise we must set the PBC selection
        self.project._pbc_selection = self._set_pbc_selection(self.structure)
        return self.project._pbc_selection
    pbc_selection = property(get_pbc_selection, None, None, "Periodic boundary conditions atom selection (read only)")

    # Indices of residues in periodic boundary conditions
    # WARNING: Do not inherit project pbc residues
    # WARNING: It may trigger all the processing logic of the reference MD when there is no need
    def get_pbc_residues (self) -> List[int]:
        # If we already have a stored value then return it
        if self.project._pbc_residues:
            return self.project._pbc_residues
        # If there is no inputs file then asume there are no PBC residues
        if not self.pbc_selection:
            self.project._pbc_residues = []
            return self.project._pbc_residues
        # Otherwise we parse the selection and return the list of residue indices     
        self.project._pbc_residues = self.structure.get_selection_residue_indices(self.pbc_selection)
        print(f'PBC residues "{self.input_pbc_selection}" -> {len(self.project._pbc_residues)} residues')
        return self.project._pbc_residues
    pbc_residues = property(get_pbc_residues, None, None, "Indices of residues in periodic boundary conditions (read only)")

    # Set the coare grain selection
    # DANI: Esto algún día habría que tratar de automatizarlo
    def _set_cg_selection (self, reference_structure : 'Structure', verbose : bool = True) -> 'Selection':
        if verbose: print('Setting Coarse Grained (CG) atoms selection')
        # If there is no inputs file then asum there is no CG selection
        if not self.project.is_inputs_file_available():
            if verbose: print(' No inputs file -> Asuming there is no CG at all')
            return Selection()
        # Otherwise we use the selection string from the inputs
        if verbose: print(' Using selection string in the inputs file')
        selection_string = self.input_cg_selection
        # If the selection is empty, again, assume there is no CG selection
        if not selection_string:
            print(' Empty selection -> There is no CG at all')
            return Selection()
        # Otherwise, process it
        # If we have a valid input value then use it
        elif selection_string:
            if verbose: print(f' Selecting CG atoms "{selection_string}"')
            parsed_selection = reference_structure.select(selection_string)
        # If we have an input value but it is empty then we set an empty selection
        else:
            if verbose: print(' No CG atoms selected')
            parsed_selection = Selection()
        # Lof the parsed selection size
        if verbose: print(f' Parsed CG selection has {len(parsed_selection)} atoms')
        # Log a few of the selected residue names
        if verbose and parsed_selection:
            selected_residues = reference_structure.get_selection_residues(parsed_selection)
            selected_residue_names = list(set([ residue.name for residue in selected_residues ]))
            limit = 3 # Show a maximum of 3 residue names
            example_residue_names = ', '.join(selected_residue_names[0:limit])
            if len(selected_residue_names) > limit: example_residue_names += ', etc.'
            print('  e.g. ' + example_residue_names)
        return parsed_selection

    # Coarse grain atom selection
    def get_cg_selection (self) -> 'Selection':
        # If we already have a stored value then return it
        if self.project._cg_selection:
            return self.project._cg_selection
        # Otherwise we must set the PBC selection
        self.project._cg_selection = self._set_cg_selection(self.structure)
        return self.project._cg_selection
    cg_selection = property(get_cg_selection, None, None, "Periodic boundary conditions atom selection (read only)")

    # Indices of residues in coarse grain
    # WARNING: Do not inherit project cg residues
    # WARNING: It may trigger all the processing logic of the reference MD when there is no need
    def get_cg_residues (self) -> List[int]:
        # If we already have a stored value then return it
        if self.project._cg_residues:
            return self.project._cg_residues
        # If there is no inputs file then asume there are no cg residues
        if not self.cg_selection:
            self.project._cg_residues = []
            return self.project._cg_residues
        # Otherwise we parse the selection and return the list of residue indices     
        self.project._cg_residues = self.structure.get_selection_residue_indices(self.cg_selection)
        print(f'CG residues "{self.input_cg_selection}" -> {len(self.project._cg_residues)} residues')
        return self.project._cg_residues
    cg_residues = property(get_cg_residues, None, None, "Indices of residues in coarse grain (read only)")

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
        non_pbc_ions_selection = self.structure.select_ions() - self.pbc_selection
        excluded_atoms_selection = non_pbc_ions_selection + self.structure.select_cg()
        # If all atoms are to be excluded then set the first frame as the reference frame and stop here
        if len(excluded_atoms_selection) == len(self.structure.atoms):
            self._reference_frame = 0
            return self._reference_frame
        # Find the first frame in the whole trajectory where safe bonds are respected
        self._reference_frame = get_bonds_canonical_frame(
            structure_filepath = structure_filepath,
            trajectory_filepath = trajectory_filepath,
            snapshots = self.snapshots,
            reference_bonds = self.safe_bonds,
            excluded_atoms_selection = excluded_atoms_selection
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
            snapshots = self.snapshots,
        )
        return self._trajectory_integrity

    # ---------------------------------------------------------------------------------
    # Analyses
    # ---------------------------------------------------------------------------------

    def run_rmsds_analysis (self):
        """RMSDs analysis."""
        # Get the task name
        task = self._get_task()
        # Check if this dependency is to be overwriten
        must_overwrite = task in self.overwritables
        # Update the overwritables so this is not remade further in the same run
        self.overwritables.discard(task)
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.pathify(OUTPUT_RMSDS_FILENAME)
        if exists(output_analysis_filepath) and not must_overwrite:
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

    def run_tmscores_analysis (self):
        """TM scores analysis."""
        # Get the task name
        task = self._get_task()
        # Check if this dependency is to be overwriten
        must_overwrite = task in self.overwritables
        # Update the overwritables so this is not remade further in the same run
        self.overwritables.discard(task)
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.pathify(OUTPUT_TMSCORES_FILENAME)
        if exists(output_analysis_filepath) and not must_overwrite:
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

    def run_rmsf_analysis (self):
        """RMSF, atom fluctuation analysis."""
        # Get the task name
        task = self._get_task()
        # Check if this dependency is to be overwriten
        must_overwrite = task in self.overwritables
        # Update the overwritables so this is not remade further in the same run
        self.overwritables.discard(task)
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.pathify(OUTPUT_RMSF_FILENAME)
        if exists(output_analysis_filepath) and not must_overwrite:
            return
        # This analysis is fast and the output size depends on the number of atoms only
        # For this reason here it is used the whole trajectory with no frames limit
        rmsf(
            input_topology_filename = self.structure_file.path,
            input_trajectory_filename = self.trajectory_file.path,
            output_analysis_filename = output_analysis_filepath,
            pbc_selection = self.pbc_selection,
        )

    def run_rgyr_analysis (self):
        """Radius of gyration analysis."""
        # Get the task name
        task = self._get_task()
        # Check if this dependency is to be overwriten
        must_overwrite = task in self.overwritables
        # Update the overwritables so this is not remade further in the same run
        self.overwritables.discard(task)
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.pathify(OUTPUT_RGYR_FILENAME)
        if exists(output_analysis_filepath) and not must_overwrite:
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

    def run_pca_analysis (self):
        """PCA, principal component analysis."""
        # Get the task name
        task = self._get_task()
        # Check if this dependency is to be overwriten
        must_overwrite = task in self.overwritables
        # Update the overwritables so this is not remade further in the same run
        self.overwritables.discard(task)
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.pathify(OUTPUT_PCA_FILENAME)
        if exists(output_analysis_filepath) and not must_overwrite:
            return
        # WARNING: This analysis will generate several output files
        # File 'pca.average.pdb' is generated by the PCA and it was used by the client but not anymore
        # File 'covar.log' is generated by the PCA but never used
        pca(
            input_topology_file = self.structure_file,
            input_trajectory_file = self.trajectory_file,
            output_analysis_filepath = output_analysis_filepath,
            output_trajectory_projections_prefix = self.pathify(OUTPUT_PCA_PROJECTION_PREFIX),
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
    # def run_pca_contacts (self):
    #     # Get the task name
    #     task = self._get_task()
    #     # Check if this dependency is to be overwriten
    #     must_overwrite = task in self.overwritables
    #     # Update the overwritables so this is not remade further in the same run
    #     self.overwritables.discard(task)
    #     # Do not run the analysis if the output file already exists
    #     output_analysis_filepath = self.pathify(OUTPUT_PCA_CONTACTS_FILENAME)
    #     if exists(output_analysis_filepath) and not must_overwrite:
    #         return
    #     pca_contacts(
    #         trajectory = self.trajectory_file.path,
    #         topology = self.pdb_filename,
    #         interactions = self.interactions,
    #         output_analysis_filename = output_analysis_filepath
    #     )

    def run_rmsd_perres_analysis (self):
        """RMSD per residue analysis."""
        # Get the task name
        task = self._get_task()
        # Check if this dependency is to be overwriten
        must_overwrite = task in self.overwritables
        # Update the overwritables so this is not remade further in the same run
        self.overwritables.discard(task)
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.pathify(OUTPUT_RMSD_PERRES_FILENAME)
        if exists(output_analysis_filepath) and not must_overwrite:
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
    def run_rmsd_pairwise_analysis (self):
        """Perform an analysis for the overall structure and then one more analysis for each interaction."""
        # Get the task name
        task = self._get_task()
        # Check if this dependency is to be overwriten
        must_overwrite = task in self.overwritables
        # Update the overwritables so this is not remade further in the same run
        self.overwritables.discard(task)
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.pathify(OUTPUT_RMSD_PAIRWISE_FILENAME)
        if must_overwrite: purge_glob(output_analysis_filepath)
        # WARNING: This analysis is fast enought to use the full trajectory instead of the reduced one
        # WARNING: However, the output file size depends on the trajectory size exponentially. It may be huge
        rmsd_pairwise(
            input_topology_filepath = self.structure_file.path,
            input_trajectory_filepath = self.trajectory_file.path,
            output_analysis_filepath = output_analysis_filepath,
            interactions = self.interactions,
            structure = self.structure,
            pbc_selection = self.pbc_selection,
            snapshots = self.snapshots,
            frames_limit = 200,
            overall_selection = "name CA or name C5"
        )

    def run_clusters_analysis (self):
        """Run the cluster analysis."""
        # Get the task name
        task = self._get_task()
        # Check if this dependency is to be overwriten
        must_overwrite = task in self.overwritables
        # Update the overwritables so this is not remade further in the same run
        self.overwritables.discard(task)
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.pathify(OUTPUT_CLUSTERS_FILENAME)
        if exists(output_analysis_filepath) and not must_overwrite:
            return
        # Set the output filepaths for additional images generated in this analysis
        output_screenshot_filepath = self.pathify(OUTPUT_CLUSTER_SCREENSHOT_FILENAMES)
        # In case the overwirte argument is passed delete all already existing outputs
        if must_overwrite:
            # Delete the summary if it exists
            if exists(output_analysis_filepath):
                remove(output_analysis_filepath)
            # Get the glob-patterned path of the output analyses
            output_analysis_glob_pattern = glob_filename(output_analysis_filepath)
            # Iterate glob matches
            for outputs in [ output_analysis_glob_pattern, output_screenshot_filepath ]:
                existing_outputs = glob(outputs)
                for existing_output in existing_outputs:
                    if exists(existing_output):
                        remove(existing_output)
        # Run the analysis
        clusters_analysis(
            input_structure_file = self.structure_file,
            input_trajectory_file = self.trajectory_file,
            interactions = self.interactions,
            structure = self.structure,
            snapshots = self.snapshots,
            pbc_selection = self.pbc_selection,
            output_analysis_filepath = output_analysis_filepath,
            output_screenshots_filename = output_screenshot_filepath,
        )

    def run_dist_perres_analysis (self):
        """Calculate the distance mean and standard deviation of each pair of residues*."""
        # Get the task name
        task = self._get_task()
        # Check if this dependency is to be overwriten
        must_overwrite = task in self.overwritables
        # Update the overwritables so this is not remade further in the same run
        self.overwritables.discard(task)
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.pathify(OUTPUT_DIST_PERRES_FILENAME)
        if must_overwrite: purge_glob(output_analysis_filepath)
        # WARNING: This analysis is not fast enought to use the full trajectory. It would take a while
        distance_per_residue(
            input_topology_filepath = self.structure_file.path,
            input_trajectory_filepath = self.trajectory_file.path,
            output_analysis_filepath = output_analysis_filepath,
            structure = self.structure,
            interactions = self.interactions,
            snapshots = self.snapshots,
            frames_limit = 200,
        )

    # Hydrogen bonds
    # WARNING: the output file size depends on the number of hydrogen bonds
    # WARNING: analyses must be no heavier than 16Mb in BSON format
    # WARNING: In case of large surface interaction the output analysis may be larger than the limit
    run_hbonds_analysis = Task('hbonds', 'Hydrogen bonds analysis',
        hydrogen_bonds, { 'time_splits': 100 })

    # SASA, solvent accessible surface analysis
    run_sas_analysis = Task('sas', 'Solvent accessible surface analysis',
        sasa, { 'frames_limit': 100 })
    
    # Perform the electrostatic and vdw energies analysis for each pair of interaction agents
    run_energies_analysis = Task('energies', 'Energies analysis',
        energies, { 'frames_limit': 100 })

    def run_dihedral_energies (self):
        """Calculate torsions and then dihedral energies for every dihedral along the trajectory."""
        # Get the task name
        task = self._get_task()
        # Check if this dependency is to be overwriten
        must_overwrite = task in self.overwritables
        # Update the overwritables so this is not remade further in the same run
        self.overwritables.discard(task)
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.pathify(OUTPUT_DIHEDRAL_ENERGIES_FILENAME)
        if exists(output_analysis_filepath) and not must_overwrite:
            return
        # Run the analysis
        compute_dihedral_energies(
            input_structure_file = self.structure_file,
            input_trajectory_file = self.trajectory_file,
            output_analysis_filepath = output_analysis_filepath,
            dihedrals_data = self.project.dihedrals,
            snapshots = self.snapshots,
            frames_limit = 100,
        )

    def run_pockets_analysis (self):
        """Perform the pockets analysis."""
        # Get the task name
        task = self._get_task()
        # Check if this dependency is to be overwriten
        must_overwrite = task in self.overwritables
        # Update the overwritables so this is not remade further in the same run
        self.overwritables.discard(task)
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.pathify(OUTPUT_POCKETS_FILENAME)
        if exists(output_analysis_filepath) and not must_overwrite:
            return
        # If we must overwritte pockets then delete the pockets folder and everything inside
        mdpocket_folder = self.pathify(POCKETS_FOLDER)
        if must_overwrite and exists(mdpocket_folder):
            rmtree(mdpocket_folder)
        # Run the analysis
        pockets(
            structure_file = self.structure_file,
            trajectory_file = self.trajectory_file,
            pockets_prefix = self.pathify(OUTPUT_POCKET_STRUCTURES_PREFIX),
            output_analysis_filepath = output_analysis_filepath,
            mdpocket_folder = mdpocket_folder,
            pbc_selection = self.pbc_selection,
            snapshots = self.snapshots,
            frames_limit = 100,
        )

    # Helical parameters
    def run_helical_analysis (self):
        # Get the task name
        task = self._get_task()
        # Check if this dependency is to be overwriten
        must_overwrite = task in self.overwritables
        # Update the overwritables so this is not remade further in the same run
        self.overwritables.discard(task)
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.pathify(OUTPUT_HELICAL_PARAMETERS_FILENAME)
        if exists(output_analysis_filepath) and not must_overwrite:
            return
        # Run the analysis
        helical_parameters(
            input_topology_filename = self.topology_file.path,
            input_trajectory_filename = self.trajectory_file.path,
            output_analysis_filename = output_analysis_filepath,
            structure = self.structure,
            structure_filename= self.structure_file.path,
            frames_limit = None,
        )
        
    # Markov
    def run_markov_analysis (self):
        # Get the task name
        task = self._get_task()
        # Check if this dependency is to be overwriten
        must_overwrite = task in self.overwritables
        # Update the overwritables so this is not remade further in the same run
        self.overwritables.discard(task)
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.pathify(OUTPUT_MARKOV_FILENAME)
        if exists(output_analysis_filepath) and not must_overwrite:
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

    # MEMBRANE ANALYSES    
    def run_density_analysis (self):
        """Membrane density analysis."""
        # Get the task name
        task = self._get_task()
        # Check if this dependency is to be overwriten
        must_overwrite = task in self.overwritables
        # Update the overwritables so this is not remade further in the same run
        self.overwritables.discard(task)
        # Do not run the analysis if the output file already exists
        output_analysis_filepath = self.pathify(OUTPUT_DENSITY_FILENAME)
        if exists(output_analysis_filepath) and not must_overwrite:
            return
        # Run the analysis
        density(
            input_structure_filepath = self.structure_file.path,
            input_trajectory_filepath = self.trajectory_file.path,
            output_analysis_filepath = output_analysis_filepath,
            membrane_map = self.project.membrane_map,
            structure = self.structure,
            snapshots = self.snapshots,
        )

    def run_thickness_analysis (self):
        """Membrane thickness analysis."""
        # Get the task name
        task = self._get_task()
        # Check if this dependency is to be overwriten
        must_overwrite = task in self.overwritables
        # Update the overwritables so this is not remade further in the same run
        self.overwritables.discard(task)
        # Do not run the analysis if the output file already exists
        output_thickness_filepath = self.pathify(OUTPUT_THICKNESS_FILENAME)
        if exists(output_thickness_filepath) and not must_overwrite:
            return
        # Run the analysis
        thickness(
            input_structure_filepath = self.structure_file.path,
            input_trajectory_filepath = self.trajectory_file.path,
            output_analysis_filepath = output_thickness_filepath,
            membrane_map = self.project.membrane_map,
            snapshots = self.snapshots,
        )

    def run_apl_analysis (self):
        """Area per lipid analysis."""
        # Get the task name
        task = self._get_task()
        # Check if this dependency is to be overwriten
        must_overwrite = task in self.overwritables
        # Update the overwritables so this is not remade further in the same run
        self.overwritables.discard(task)
        # Do not run the analysis if the output file already exists
        output_apl_filepath = self.pathify(OUTPUT_APL_FILENAME)
        if exists(output_apl_filepath) and not must_overwrite:
            return
        # Run the analysis
        area_per_lipid(
            input_structure_filepath = self.structure_file.path,
            input_trajectory_filepath = self.trajectory_file.path,
            output_analysis_filepath = output_apl_filepath,
            membrane_map = self.project.membrane_map,
        )

    def run_lipid_order_analysis (self):
        """Calculate lipid order parameters for membranes."""
        # Get the task name
        task = self._get_task()
        # Check if this dependency is to be overwriten
        must_overwrite = task in self.overwritables
        # Update the overwritables so this is not remade further in the same run
        self.overwritables.discard(task)
        # Do not run the analysis if the output file already exists
        output_lipid_order_filepath = self.pathify(OUTPUT_LIPID_ORDER_FILENAME)
        if exists(output_lipid_order_filepath) and not must_overwrite:
            return
        # Run the analysis
        lipid_order(
            input_trajectory_filepath = self.trajectory_file.path,
            topology_file=self.project.standard_topology_file,
            output_analysis_filepath = output_lipid_order_filepath,
            membrane_map = self.project.membrane_map,
            snapshots = self.snapshots,
        )

    def run_lipid_interactions_analysis (self):
        """Lipid-protein interactions analysis."""
        # Get the task name
        task = self._get_task()
        # Check if this dependency is to be overwriten
        must_overwrite = task in self.overwritables
        # Update the overwritables so this is not remade further in the same run
        self.overwritables.discard(task)
        # Do not run the analysis if the output file already exists
        output_lipid_interactions_filepath = self.pathify(OUTPUT_LIPID_INTERACTIONS_FILENAME)
        if exists(output_lipid_interactions_filepath) and not must_overwrite:
            return
        # Run the analysis
        lipid_interactions(
            input_trajectory_filepath = self.trajectory_file.path,
            topology_file=self.project.standard_topology_file,
            output_analysis_filepath = output_lipid_interactions_filepath,
            membrane_map = self.project.membrane_map,
            snapshots = self.snapshots,
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
        cg_selection : Optional[str] = None,
        image : bool = False,
        fit : bool = False,
        translation : List[float] = [0, 0, 0],
        mercy : Union[ List[str], bool ] = [],
        trust : Union[ List[str], bool ] = [],
        pca_selection : str = PROTEIN_AND_NUCLEIC_BACKBONE,
        pca_fit_selection : str = PROTEIN_AND_NUCLEIC_BACKBONE,
        rmsd_cutoff : float = DEFAULT_RMSD_CUTOFF,
        interaction_cutoff : float = DEFAULT_INTERACTION_CUTOFF,
        interactions_auto : Optional[str] = None,
        # Set it we must download just a few frames instead of the whole trajectory
        sample_trajectory : Optional[int] = None,
    ):
        # Save input parameters
        self.directory = remove_final_slash(directory)
        self.database_url = database_url
        self.accession = accession
        # Set the project URL in case we have the required data
        self.remote = None
        if self.database_url and self.accession:
            self.remote = Remote(self.database_url, self.accession)

        # Set the inputs file
        # Set the expected default name in case there is no inputs file since it may be downloaded
        self._inputs_file = File(self.pathify(DEFAULT_INPUTS_FILENAME))
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
        self._input_cg_selection = cg_selection
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
        self.interactions_auto = interactions_auto
        # Set the inputs, where values from the inputs file will be stored
        self._inputs = None

        # Other values which may be found/calculated on demand
        self._pbc_selection = None
        self._pbc_residues = None
        self._cg_selection = None
        self._cg_residues = None
        self._safe_bonds = None
        self._charges = None
        self._topology_reader = None
        self._dihedrals = None
        self._populations = None
        self._transitions = None
        self._pdb_ids = None
        self._pdb_references = None
        self._protein_map = None
        self._ligand_map = None
        self._membrane_map = None
        self.pubchem_name_list = None
        self._residue_map = None
        self._mds = None

        # Force a couple of extraordinary files which is generated if atoms are resorted
        self.resorted_bonds_file = File(self.pathify(RESORTED_BONDS_FILENAME))
        self.resorted_charges_file = File(self.pathify(RESORTED_CHARGES_FILENAME))

        # Set a new entry for the register
        # This is useful to track previous workflow runs and problems
        register_filepath = self.pathify(REGISTER_FILENAME)
        register_file = File(register_filepath)
        self.register = Register(register_file)

        # Set the cache
        cache_filepath = self.pathify(CACHE_FILENAME)
        cache_file = File(cache_filepath)
        self.cache = Cache(cache_file)

        # Set tasks whose output is to be overwritten
        self.overwritables = set()

    def __repr__ (self):
        return 'Project'

    # This function is able to find its caller "self" function
    # Then it finds its associated label in the requestables
    def _get_task (self) -> str:
        caller_data = sys._getframe().f_back
        caller_name = caller_data.f_code.co_name
        caller_func = getattr(self, caller_name).__func__
        return inverted_requestables[caller_func]

    # Given a filename or relative path, add the project directory path at the beginning
    def pathify (self, filename_or_relative_path : str) -> str:
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

    def is_inputs_file_available (self) -> bool:
        """Set a function to check if inputs file is available.
        Note that asking for it when it is not available will lead to raising an input error."""
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

    def get_inputs_file (self) -> File:
        """Set a function to load the inputs file"""
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

    
    def guess_input_topology_filepath (self) -> Optional[str]:
        """If there is not input topology filepath, we try to guess it among the files in the project directory.
        Note that if we can download from the remote then we must check the remote available files as well."""
        # Find the first supported topology file according to its name and format
        def find_first_accepted_topology_filename (available_filenames : List[str]) -> Optional[str]:
            for filename in available_filenames:
                # Make sure it is a valid topology file candidate
                # i.e. topology.xxx
                filename_splits = filename.split('.')
                if len(filename_splits) != 2 or filename_splits[0] != 'topology':
                    continue
                # Then make sure its format is among the accepted topology formats
                extension = filename_splits[1]
                format = EXTENSION_FORMATS[extension]
                if format in ACCEPTED_TOPOLOGY_FORMATS:
                    return filename
            return None
        # First check among the local available files
        local_files = list_files(self.directory)
        accepted_topology_filename = find_first_accepted_topology_filename(local_files)
        if accepted_topology_filename:
            return self.pathify(accepted_topology_filename)
        # In case we did not find a topology among the local files, repeat the process with the remote files
        if self.remote:
            remote_files = self.remote.available_files
            accepted_topology_filename = find_first_accepted_topology_filename(remote_files)
            if accepted_topology_filename:
                return self.pathify(accepted_topology_filename)
        # If no actual topology is to be found then try with the standard topology instead
        # Check if the standard topology file is available
        # Note that we do not use standard_topology_file property to avoid generating it at this point
        standard_topology_filepath = self.pathify(STANDARD_TOPOLOGY_FILENAME)
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

    
    def get_input_topology_filepath (self) -> Optional[str]:
        """Get the input topology filepath from the inputs or try to guess it.
        If the input topology filepath is a 'no' flag then we consider there is no topology at all
        So far we extract atom charges and atom bonds from the topology file
        In this scenario we can keep working but there are some consecuences:
        1 - Analysis using atom charges such as 'energies' will be skipped
        2 - The standard topology file will not include atom charges
        3 - Bonds will be guessed"""
        if type(self.input_topology_filepath) == str and self.input_topology_filepath.lower() in { 'no', 'not', 'na' }:
            return MISSING_TOPOLOGY
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
        raise InputError('Missing input topology file path. Please provide a topology file using the "-top" argument.\n' +
            '  Note that you may run the workflow without a topology file. To do so, use the "-top no" argument.\n' +
            '  However this has implications since we usually mine atom charges and bonds from the topology file.\n' +
            '  Some analyses such us the interaction energies will be skiped')

    def get_input_topology_file (self) -> Optional[File]:
        """Get the input topology file.
        If the file is not found try to download it."""
        # If we already have a value then return it
        if self._input_topology_file != None:
            return self._input_topology_file
        # Set the input topology filepath
        input_topology_filepath = self.get_input_topology_filepath()
        # If the input filepath is None then it menas we must proceed without a topology
        if input_topology_filepath == MISSING_TOPOLOGY:
            self._input_topology_file = MISSING_TOPOLOGY
            return self._input_topology_file
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
                itp_filepath = self.pathify(itp_filename)
                itp_file = File(itp_filepath)
                self.remote.download_file(itp_file)
        return self._input_topology_file
    input_topology_file = property(get_input_topology_file, None, None, "Input topology file (read only)")

    # Input structure filename ------------
    def get_input_structure_file (self) -> File:
        """Get the input structure filename."""
        # When calling this function make sure all MDs have the file or try to download it
        return self.reference_md._input_structure_file
    input_structure_file = property(get_input_structure_file, None, None, "Input structure filename for each MD (read only)")

    # Input trajectory filename ------------

    def get_input_trajectory_files (self) -> List[File]:
        """Get the input trajectory filename(s) from the inputs.
        If file(s) are not found try to download it."""
        return self.reference_md._input_trajectory_files
    input_trajectory_files = property(get_input_trajectory_files, None, None, "Input trajectory filenames for each MD (read only)")

    # Populations filename ------------

    def get_populations_file (self) -> File:
        """Get the MSM equilibrium populations filename."""
        if not self.get_file(self._populations_file):
            return None
        return self._populations_file
    populations_file = property(get_populations_file, None, None, "MSM equilibrium populations filename (read only)")

    # Transitions filename ------------

    def get_transitions_file (self) -> Optional[str]:
        """Get the MSM transition probabilities filename."""
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

    # CG selection may come from the console or from the inputs file
    # Console has priority over the inputs file
    def get_input_cg_selection (self) -> Optional[str]:
        # If we have an internal value then return it
        if self._input_cg_selection:
            return self._input_cg_selection
        # As an exception, we avoid asking for the inputs file if it is not available
        # This input is required for some early processing steps where we do not need the inputs file for anything else
        if not self.is_inputs_file_available():
            return None
        # Otherwise, find it in the inputs
        # Get the input value, whose key must exist
        self._input_cg_selection = self.get_input('cg_selection')
        return self._input_cg_selection
    input_cg_selection = property(get_input_cg_selection, None, None, "Selection of atoms which are not acutal atoms but Coarse Grained beads (read only)")

    # Set additional values infered from input values

    def check_is_time_dependent (self) -> bool:
        """Set if MDs are time dependent."""
        if self.input_type == 'trajectory':
            return True
        elif self.input_type == 'ensemble':
            return False
        raise InputError('Not supported input type value: ' + self.input_type)
    is_time_dependent = property(check_is_time_dependent, None, None, "Check if trajectory frames are time dependent (read only)")

    # Processed files ----------------------------------------------------

    # Set the expected output topology filename given the input topology filename
    # Note that topology formats are conserved
    def inherit_topology_filename (self) -> Optional[str]:
        if self.input_topology_file == MISSING_TOPOLOGY:
            return None
        filename = self.input_topology_file.filename
        if not filename:
            return None
        if filename == RAW_CHARGES_FILENAME:
            return filename
        standard_format = self.input_topology_file.format
        return 'topology.' + standard_format

    def get_topology_filepath (self) -> str:
        """Get the processed topology file path."""
        # If we have a stored value then return it
        if self._topology_filepath:
            return self._topology_filepath
        # Otherwise we must find it
        self._topology_filepath = self.inherit_topology_filename()
        return self._topology_filepath
    topology_filepath = property(get_topology_filepath, None, None, "Topology file path (read only)")

    def get_topology_file (self) -> str:
        """Get the processed topology file."""
        # If we have a stored value then return it
        # This means we already found or generated this file
        if self._topology_file != None:
            return self._topology_file
        # If the file already exists then we are done
        self._topology_file = File(self.topology_filepath) if self.topology_filepath != None else MISSING_TOPOLOGY
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

    # Indices of residues in coarse grain
    def get_cg_residues (self) -> List[int]:
        return self.reference_md.cg_residues
    cg_residues = property(get_cg_residues, None, None, "Indices of residues in coarse grain (read only)")


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

    # Topolody data reader
    def get_topology_reader (self) -> 'Topology':
        # If we already have a stored value then return it
        if self._topology_reader: return self._topology_reader
        # Instantiate the topology reader
        self._topology_reader = Topology(self.topology_file)
        return self._topology_reader
    topology_reader = property(get_topology_reader, None, None, "Topology reader (read only)")

    # Dihedrals data
    def get_dihedrals (self) -> List[dict]:
        # If we already have a stored value then return it
        if self._dihedrals: return self._dihedrals
        # Calculate the dihedrals otherwise
        self._dihedrals = self.topology_reader.get_dihedrals_data()
        return self._dihedrals
    dihedrals = property(get_dihedrals, None, None, "Topology dihedrals (read only)")

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

    def get_pdb_references (self) -> List[dict]:
        """Get PDB references."""
        # If we already have a stored value then return it
        if self._pdb_references:
            return self._pdb_references
        # Get the task name
        task = self._get_task()
        # Check if this dependency is to be overwriten
        must_overwrite = task in self.overwritables
        # Update the overwritables so this is not remade further in the same run
        self.overwritables.discard(task)
        # Set the PDB references file
        pdb_references_filepath = self.pathify(PDB_REFERENCES_FILENAME)
        pdb_references_file = File(pdb_references_filepath)
        # If the file already exists and we must overwrite it then delete it here
        if pdb_references_file.exists and must_overwrite:
            pdb_references_file.remove()
        # Otherwise we must find the value
        self._pdb_references = generate_pdb_references(
            pdb_ids = self.pdb_ids,
            pdb_references_file = pdb_references_file
        )
        return self._pdb_references
    pdb_references = property(get_pdb_references, None, None, "PDB references (read only)")

    def get_pdb_references_file (self) -> File:
        """Define the PDB references output file."""
        # Set the PDB references file
        pdb_references_filepath = self.pathify(PDB_REFERENCES_FILENAME)
        pdb_references_file = File(pdb_references_filepath)
        # Ask for the references thus producing the requested output file
        # WARNING: We must always run the references, no matter if the file already exists
        # WARNING: There could be a pending overwrite to run
        # WARNING: Also note that if it was run already it won't run again
        self.get_pdb_references()
        return pdb_references_file
    pdb_references_file = property(get_pdb_references_file, None, None, "File including PDB refereces data (read only)")

    def get_protein_map (self) -> List[dict]:
        """Map the structure aminoacids sequences against the Uniprot reference sequences."""
        # If we already have a stored value then return it
        if self._protein_map != None:
            return self._protein_map
        # Get the task name
        task = self._get_task()
        # Check if this dependency is to be overwriten
        must_overwrite = task in self.overwritables
        # Update the overwritables so this is not remade further in the same run
        self.overwritables.discard(task)
        # Set the protein references file
        protein_references_filepath = self.pathify(PROTEIN_REFERENCES_FILENAME)
        protein_references_file = File(protein_references_filepath)
        # If the file already exists and we must overwrite it then delete it here
        if protein_references_file.exists and must_overwrite:
            protein_references_file.remove()
        # Otherwise we must find the value
        self._protein_map = generate_protein_mapping(
            structure = self.structure,
            protein_references_file = protein_references_file,
            database_url = self.database_url,
            register = self.register,
            mercy = self.mercy,
            forced_references = self.forced_references,
            pdb_ids = self.pdb_ids,
        )
        return self._protein_map
    protein_map = property(get_protein_map, None, None, "Residues mapping (read only)")

    # Define the output file of the protein mapping including protein references
    def get_protein_references_file (self) -> File:
        # Set the protein references file
        protein_references_filepath = self.pathify(PROTEIN_REFERENCES_FILENAME)
        protein_references_file = File(protein_references_filepath)
        # Ask for the map thus producing the requested output file
        # WARNING: We must always run the map, no matter if the file already exists
        # WARNING: There could be a pending overwrite to run
        # WARNING: Also note that if it was run already it won't run again
        self.get_protein_map()
        # Return the file
        return protein_references_file
    protein_references_file = property(get_protein_references_file, None, None, "File including protein refereces data mined from UniProt (read only)")

    def get_chain_references (self) -> List[str]:
        """Get chain references."""
        # Get the task name
        task = self._get_task()
        # Check if this dependency is to be overwriten
        must_overwrite = task in self.overwritables
        # Update the overwritables so this is not remade further in the same run
        self.overwritables.discard(task)
        # Set the chains references file
        chains_references_filepath = self.pathify(OUTPUT_CHAINS_FILENAME)
        chains_references_file = File(chains_references_filepath)
        # If the file already exists and we must overwrite then remove it
        if chains_references_file.exists and must_overwrite:
            chains_references_file.remove()
        # Call the function to generate the chain references
        chains = generate_chain_references(
            structure = self.structure,
            chains_references_file = chains_references_file,
            database_url=self.database_url
            #chain_name=self.structure.chain_name,
        )
        return chains
    chains_data = property(get_chain_references, None, None, "Chain (read only)")

    def get_ligand_map (self) -> List[dict]:
        """Get the ligand residues mapping."""
        # If we already have a stored value then return it
        if self._ligand_map != None:
            return self._ligand_map
        # Get the task name
        task = self._get_task()
        # Check if this dependency is to be overwriten
        must_overwrite = task in self.overwritables
        # Update the overwritables so this is not remade further in the same run
        self.overwritables.discard(task)
        # Set the ligand references file
        ligand_references_filepath = self.pathify(LIGAND_REFERENCES_FILENAME)
        ligand_references_file = File(ligand_references_filepath)
        # If the file already exists and we must overwrite it then delete it here
        if ligand_references_file.exists and must_overwrite:
            ligand_references_file.remove()
        # Otherwise we must find the value
        self._ligand_map, self.pubchem_name_list = generate_ligand_mapping(
            structure = self.structure,
            cache = self.cache,
            input_ligands = self.input_ligands,
            input_pdb_ids = self.pdb_ids,
            output_ligands_filepath = ligand_references_file.path, 
            mercy = self.mercy,
        )
        return self._ligand_map
    ligand_map = property(get_ligand_map, None, None, "Ligand references (read only)")

    # Define the output file of the ligand mapping including ligand references
    def get_ligand_references_file (self) -> File:
        # Set the ligand references file
        ligand_references_filepath = self.pathify(LIGAND_REFERENCES_FILENAME)
        ligand_references_file = File(ligand_references_filepath)
        # Ask for the map thus producing the requested output file
        # WARNING: We must always run the map, no matter if the file already exists
        # WARNING: There could be a pending overwrite to run
        # WARNING: Also note that if it was run already it won't run again
        self.get_ligand_map()
        return ligand_references_file
    ligand_references_file = property(get_ligand_references_file, None, None, "File including ligand refereces data mined from PubChem (read only)")

    def get_membrane_map (self) -> List[dict]:
        """Get mapping of residues in the membrane."""
        # If we already have a stored value then return it
        if self._membrane_map:
            return self._membrane_map
        # Get the task name
        task = self._get_task()
        # Check if this dependency is to be overwriten
        must_overwrite = task in self.overwritables
        # Update the overwritables so this is not remade further in the same run
        self.overwritables.discard(task)
        # Set the membrane mapping file
        mem_map_filepath = self.pathify(MEMBRANE_MAPPING_FILENAME)
        mem_map_file = File(mem_map_filepath)
        # If the file already exists then send it
        if mem_map_file.exists and not must_overwrite:
            self._membrane_map =  load_json(mem_map_file.path)
        else:
            self._membrane_map = generate_membrane_mapping(
                structure = self.structure,
                topology_file=self.standard_topology_file,
                structure_file=self.structure_file,
            )   
        if self._membrane_map is None or self._membrane_map['n_mems'] == 0:
            print('No membrane available. Related analyses will be skipped.')
        return self._membrane_map
    membrane_map = property(get_membrane_map, None, None, "Membrane mapping (read only)")

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

    
    def get_metadata_file (self) -> File:
        """Generate the project metadata file to be upload to the database."""
        # Get the task name
        task = self._get_task()
        # Check if this dependency is to be overwriten
        must_overwrite = task in self.overwritables
        # Update the overwritables so this is not remade further in the same run
        self.overwritables.discard(task)
        # Set the metadata file
        metadata_filepath = self.pathify(OUTPUT_METADATA_FILENAME)
        metadata_file = File(metadata_filepath)
        # If the file already exists then send it
        if metadata_file.exists and not must_overwrite:
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
            interactions = self.reference_md.interactions
        )
        return metadata_file
    metadata_file = property(get_metadata_file, None, None, "Project metadata filename (read only)")

    def get_standard_topology_file (self) -> File:
        """Generate the standardized topology file."""
        # Get the task name
        task = self._get_task()
        # Check if this dependency is to be overwriten
        must_overwrite = task in self.overwritables
        # Update the overwritables so this is not remade further in the same run
        self.overwritables.discard(task)
        # If we have a stored value then return it
        # This means we already found or generated this file
        if self._standard_topology_file:
            return self._standard_topology_file
        # Set the standard topology file
        standard_topology_filepath = self.pathify(STANDARD_TOPOLOGY_FILENAME)
        self._standard_topology_file = File(standard_topology_filepath)
        # If the file already exists and it is not to be overwirtten then send it
        if self._standard_topology_file.exists and not must_overwrite:
            return self._standard_topology_file
        # Otherwise, generate it
        generate_topology(
            structure = self.structure,
            charges = self.charges,
            residue_map = self.residue_map,
            pbc_residues = self.pbc_residues,
            cg_residues = self.cg_residues,
            output_topology_filepath = self._standard_topology_file.path
        )
        # Register the last modification times of the new generated standard topology file
        # Note that this may be already the topology file, but it may be not
        self.register.update_mtime(self._standard_topology_file)
        return self._standard_topology_file
    standard_topology_file = property(get_standard_topology_file, None, None, "Standard topology filename (read only)")

    def get_screenshot_filename (self) -> str:
        """Generate a screenshot of the system."""
        # Get the task name
        task = self._get_task()
        # Check if this dependency is to be overwriten
        must_overwrite = task in self.overwritables
        # Update the overwritables so this is not remade further in the same run
        self.overwritables.discard(task)
        # Set the screenshot file
        screenshot_filepath = self.pathify(OUTPUT_SCREENSHOT_FILENAME)
        screenshot_file = File(screenshot_filepath)
        # If the file already exists then send it
        if screenshot_file.exists and not must_overwrite:
            return screenshot_file
        # Otherwise, generate it
        get_screenshot(
            structure = self.structure,
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

def check_directory (directory : str) -> str:
    """Check for problematic characters in a directory path."""
    # Remove problematic characters
    for character in FORBIDDEN_DIRECTORY_CHARACTERS:
        if character in directory:
            raise InputError(f'Directory path "{directory}" includes the forbidden character "{character}"')

def directory_2_name (directory : str) -> str:
    """Convert an MD directory into an equivalent MD name."""
    # Replace white spaces with underscores
    name = directory.replace('_', ' ')
    return name

def remove_final_slash (directory : str) -> str:
    """Remove the final slash if exists since it may cause 
    problems when recognizing input directories."""
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
default_analyses = {
    'clusters': MD.run_clusters_analysis,
    'dist': MD.run_dist_perres_analysis,
    #'energies': MD.run_energies_analysis,
    #'hbonds': MD.run_hbonds_analysis,
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
    #'sas': MD.run_sas_analysis,
    'tmscore': MD.run_tmscores_analysis,
    'density': MD.run_density_analysis,
    'thickness': MD.run_thickness_analysis,
    'apl': MD.run_apl_analysis,
    'lorder': MD.run_lipid_order_analysis,
    'linter': MD.run_lipid_interactions_analysis,
}

extra_analyses = {
    'dihedrals': MD.run_dihedral_energies,
}

analyses = { **default_analyses, **extra_analyses }

# Project requestable tasks
project_requestables = {
    **project_input_files,
    **project_processed_files,
    'pdbs': Project.get_pdb_references,
    'mapping': Project.get_protein_map,
    'ligands': Project.get_ligand_map,
    'screenshot': Project.get_screenshot_filename,
    'stopology': Project.get_standard_topology_file,
    'pmeta': Project.get_metadata_file,
    'chains': Project.get_chain_references,
    'membrane': Project.get_membrane_map, 
    'membranes': Project.get_membrane_map, 
}
# Add available tasks to project requestables
for callable in vars(Project).values():
    if isinstance(callable, Task): project_requestables[callable.flag] = callable
# MD requestable tasks
md_requestables = {
    **md_input_files,
    **md_processed_files,
    **analyses,
    #'inter': MD.get_processed_interactions,
    'mdmeta': MD.get_metadata_file,
}
# Add available tasks to project requestables
for callable in vars(MD).values():
    if isinstance(callable, Task): md_requestables[callable.flag] = callable
# Requestables for the console
# Note that this constant is global
requestables.update({ **project_requestables, **md_requestables })
# Inverted requestables for every function to know which is its 'label'
inverted_requestables.update({ v: k for k, v in requestables.items() })

# Set groups of dependencies to be requested together using only one flag
DEPENDENCY_FLAGS = {
    'download': list(input_files.keys()),
    'setup': list(processed_files.keys()),
    'network': [ 'mapping', 'ligands', 'chains', 'pdbs', 'membrane' ],
    'minimal': [ 'pmeta', 'mdmeta', 'stopology' ],
    'interdeps': [ 'inter', 'pairwise', 'hbonds', 'energies', 'perres', 'clusters' ]
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
        warn('The "-d" or "--download" argument is deprecated. Please use "-i download" instead.')
        tasks = list(input_files.keys())
    # If the setup argument is passed then just process input files
    elif setup:
        warn('The "-s" or "--setup" argument is deprecated. Please use "-i setup" instead.')
        tasks = list(processed_files.keys())
    # If the include argument then add only the specified tasks to the list
    elif include and len(include) > 0:
        tasks = [ *include ]
        # Find for special flags among included
        for flag, dependencies, in DEPENDENCY_FLAGS.items():
            if flag not in tasks: continue
            # If the flag is found then remove it and write the corresponding dependencie instead
            # Make sure not to duplicate a dependency if it was already included
            tasks.remove(flag)
            for dep in dependencies:
                if dep in tasks: continue
                tasks.append(dep)
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
            'inter',
            *default_analyses.keys(),
        ]
        # WARNING: Do not run helical by default, it will fail in the default environment
        # There is a separated enviornment to run this analysis
        tasks.remove('helical')
        # If the exclude parameter was passed then remove excluded tasks from the default tasks
        if exclude and len(exclude) > 0:
            tasks = [ name for name in tasks if name not in exclude ]

    # If the user requested to overwrite something, make sure it is in the tasks list

    # Update the overwritable variable with the requested overwrites
    overwritables = set()
    if overwrite:
        # If the overwrite argument is simply true then add all requestables to the overwritable
        if type(overwrite) == bool:
            for task in tasks:
                overwritables.add(task)
        # If the overwrite argument is a list of tasks then iterate them
        elif type(overwrite) == list:
            for task in overwrite:
                # Make sure the task to be overwriten is among the tasks to be run
                if task not in tasks:
                    raise InputError(f'Task "{task}" is to be overwriten but it is not in the tasks list. Either include it or do not exclude it')
                # Add it to the global variable
                overwritables.add(task)
        else: raise ValueError('Not supported overwrite type')

    # Get project tasks
    project_tasks = [ task for task in tasks if task in project_requestables ]
    # Get the MD tasks
    md_tasks = [ task for task in tasks if task in md_requestables ]

    # Set project overwritables
    project.overwritables = set([ task for task in project_tasks if task in overwritables ])
    # Set MD overwrittables
    # Note that this must be done before running project tasks
    # Some project tasks rely in MD tasks
    for md in project.mds:
        md.overwritables = set([ task for task in md_tasks if task in overwritables ])

    # Run the project tasks now
    for task in project_tasks:
        # Get the function to be called and call it
        getter = requestables[task]
        getter(project)

    # If there are no MD tasks then we are done already
    if len(md_tasks) == 0:
        print("Finished!")
        return

    # Now iterate over the different MDs
    for md in project.mds:
        print(f'\n{CYAN_HEADER} Processing MD at {md.directory}{COLOR_END}')
        # Run the MD tasks
        for task in md_tasks:
            # Get the function to be called and call it
            getter = requestables[task]
            getter(md)

        # Remove gromacs backups and other trash files from this MD
        remove_trash(md.directory)

    # Remove gromacs backups and other trash files from the project
    remove_trash(project.directory)

    print("Done!")
