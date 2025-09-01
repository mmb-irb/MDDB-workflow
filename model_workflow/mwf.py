#!/usr/bin/env python

# This is the starter script

# Import python libraries
from os import chdir, rename, remove, walk, mkdir, getcwd
from os.path import exists, isdir, isabs, relpath, normpath
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
from model_workflow.utils.auxiliar import is_glob, parse_glob, safe_getattr
from model_workflow.utils.arg_cksum import get_cksum_id
from model_workflow.utils.register import Register
from model_workflow.utils.cache import Cache
from model_workflow.utils.structures import Structure
from model_workflow.utils.topologies import Topology
from model_workflow.utils.file import File
from model_workflow.utils.remote import Remote
from model_workflow.utils.pyt_spells import get_frames_count, get_average_structure
from model_workflow.utils.selections import Selection
from model_workflow.utils.mda_spells import get_mda_universe
from model_workflow.utils.type_hints import *

# Import local tools
from model_workflow.tools.get_first_frame import get_first_frame
from model_workflow.tools.get_bonds import find_safe_bonds, get_bonds_canonical_frame
from model_workflow.tools.process_interactions import process_interactions
from model_workflow.tools.find_interaction_types import find_interaction_types
from model_workflow.tools.generate_metadata import prepare_project_metadata, generate_md_metadata
from model_workflow.tools.generate_ligands_desc import generate_ligand_mapping
from model_workflow.tools.chains import prepare_chain_references
from model_workflow.tools.generate_pdb_references import prepare_pdb_references
from model_workflow.tools.residue_mapping import generate_residue_mapping
from model_workflow.tools.generate_map import generate_protein_mapping
from model_workflow.tools.generate_lipid_references import generate_lipid_references
from model_workflow.tools.generate_membrane_mapping import generate_membrane_mapping
from model_workflow.tools.generate_topology import generate_topology
from model_workflow.tools.get_charges import get_charges
from model_workflow.tools.remove_trash import remove_trash
from model_workflow.tools.get_screenshot import get_screenshot
from model_workflow.tools.process_input_files import process_input_files

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
# Only when not in a Jupyter notebook or using pytest
# Check if we're in an interactive Python shell like Jupyter
if not hasattr(sys, 'ps1') and not 'pytest' in sys.modules:  
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
OUTPUT_FILEPATH_ARG = 'output_filepath'
OUTPUT_DIRECTORY_ARG = 'output_directory'

# Set some variables which are filled at the end but are referred by previously defined functions
requestables = {}
inverted_requestables = {}

# Set a class to handle a generic task
# WARNING: Note that a task is static in relation to its Project/MD
# This means that all MDs share the same task and thus it can not store internal values
# Instead all its internal vallues are stored in the parent using a key name
# This is not easy to change with the curent implementation
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
        args : dict = {},
        # In case this task is to produce an output file, set here its name
        # The actual path relative to its project/MD will be set automatically
        # For those tasks which generate a directory with multiple outputs this is not necessary
        # However this may come in handy by tasks with a single file output
        # Specillay when this output file is used later in this workflow
        output_filename : Optional[str] = None,
        # Set if the returned output is to be cached
        # Note that argument values are always cached, this is not optional
        use_cache : bool = True
    ):
        # Save input arguments
        self.flag = flag
        self.name = name
        self.func = func
        self.args = args
        self.output_filename = output_filename
        self.use_cache = use_cache
        # Set the key used to store and retireve data in the parent and cache
        self.parent_output_key = f'_{self.flag}_task_output'
        self.parent_output_filepath_key = f'{self.flag}_task_output_filepath'
        self.cache_output_key = f'{self.flag}_task_output'
        self.cache_arg_cksums = f'{self.flag}_task_arg_cksums'
        # Para el arg_cksum
        self.__name__ = self.flag

    # Set internal functions to handle parent saved output
    # This output is not saved in the task itself, but in the parent, because the task is static
    def _get_parent_output (self, parent):
        return safe_getattr(parent, self.parent_output_key, MISSING_VALUE_EXCEPTION)
    def _set_parent_output (self, parent, new_output):
        return setattr(parent, self.parent_output_key, new_output)
    def _get_parent_output_filepath (self, parent):
        return safe_getattr(parent, self.parent_output_filepath_key, MISSING_VALUE_EXCEPTION)
    def _set_parent_output_filepath (self, parent, new_filepath):
        return setattr(parent, self.parent_output_filepath_key, new_filepath)
    # Get the task output, running the task if necessary
    def get_output (self, parent):
        # If we already have a value stored from this run then return it
        output = self._get_parent_output(parent)
        if output != MISSING_VALUE_EXCEPTION: return output
        # Otherwise run the task and return the output
        return self(parent)
    output = property(get_output, None, None, "Task output (read only)")
    # Asking for the output file or filepath implies running the Task, then returning the file/filepath
    def get_output_filepath (self, parent) -> str:
        # If we already have a filepath stored from this run then return it
        filepath = self._get_parent_output_filepath(parent)
        if filepath != MISSING_VALUE_EXCEPTION: return filepath
        # Otherwise run the task and return the filepath
        self(parent)
        filepath = self._get_parent_output_filepath(parent)
        if filepath != MISSING_VALUE_EXCEPTION: return filepath
        raise ValueError(f'Task {self.flag} has no output filepath after running')
    output_filepath = property(get_output_filepath, None, None, "Task output filepath (read only)")
    def get_output_file (self, parent) -> str:
        filepath = self.get_output_filepath(parent)
        return File(filepath)
    output_file = property(get_output_file, None, None, "Task output file (read only)")

    # When the task is printed, show the flag
    def __repr__ (self): return f'<Task ({self.flag})>'

    # When a task is called
    def __call__(self, parent):
        # First of all check if this task has been already done in this very run
        # If so then return the stored vale
        output = self._get_parent_output(parent)
        if output != MISSING_VALUE_EXCEPTION: return output
        # Process the task function arguments
        processed_args = {}
        # Get the task function expected arguments
        specification = getfullargspec(self.func)
        expected_arguments = specification.args
        n_default_arguments = len(specification.defaults) if specification.defaults else 0
        # Find out which arguments are optional since they have default values
        default_arguments = set(expected_arguments[::-1][:n_default_arguments])
        # If one of the expected arguments is the output_filename then set it here
        output_filepath = None
        writes_output_file = OUTPUT_FILEPATH_ARG in expected_arguments
        if writes_output_file:
            # The task should have a defined output file
            if not self.output_filename:
                raise RuntimeError(f'Task {self.flag} must have an "output_filename"')
            # Set the output file path
            output_filepath = parent.pathify(self.output_filename)
            self._set_parent_output_filepath(parent, output_filepath)
            # Add it to the processed args
            processed_args[OUTPUT_FILEPATH_ARG] = output_filepath
            # Remove the expected argument from the list
            expected_arguments.remove(OUTPUT_FILEPATH_ARG)
        # If one of the expected arguments is the output_directory then set it here
        # We will set a new directory with the flag name of the task, in the correspoding path
        # Note that while the task is beeing done the output directory has a different name
        # Thus the directory is hidden and marked as incomplete
        # The final output directory is the one without the incomplete prefix
        writes_output_dir = OUTPUT_DIRECTORY_ARG in expected_arguments
        incomplete_output_directory = None
        final_output_directory = None
        if writes_output_dir:
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
        # Check again if the task has output already
        # It may happend that some dependencies assign output on their own
        # e.g. charges, bonds
        # If so then return the stored vale
        output = self._get_parent_output(parent)
        if output != MISSING_VALUE_EXCEPTION: return output
        # Find if we have cached output
        if self.use_cache:
            output = parent.cache.retrieve(self.cache_output_key, MISSING_VALUE_EXCEPTION)
            self._set_parent_output(parent, output)
        # Check if this dependency is to be overwriten
        forced_overwrite = self.flag in parent.overwritables
        # Get the list of inputs which have changed compared to a previous run
        # WARNING: Always get changed inputs, since this function updates the cache
        # If had_cache is false then it means this is the first time the task is ever done
        changed_inputs, had_cache, cache_cksums = self.get_changed_inputs(parent, processed_args)
        any_input_changed = len(changed_inputs) > 0
        # Update the cache inputs
        parent.cache.update(self.cache_arg_cksums, cache_cksums)
        # We must overwrite outputs either if inputs changed or if it was forced by the user
        must_overwrite = forced_overwrite or any_input_changed
        # Check if output already exists
        # If the final directory already exists then it means the task was started in a previous run
        existing_incomplete_output = writes_output_dir and exists(incomplete_output_directory)
        # If the final directory already exists then it means the task was done in a previous run
        existing_final_output = writes_output_dir and exists(final_output_directory)
        # If the output file already exists then it also means the task was done in a previous run
        existing_output_file = writes_output_file and exists(output_filepath)
        # If we already have a cached output result
        existing_output_data = output != MISSING_VALUE_EXCEPTION
        # If we must overwrite then purge previous outputs
        if must_overwrite:
            if existing_incomplete_output: rmtree(incomplete_output_directory)
            if existing_final_output: rmtree(final_output_directory)
            if existing_output_file: remove(output_filepath)
            if existing_output_data: parent.cache.delete(self.cache_output_key)
        # If already existing output is not to be overwritten then check if it is already what we need
        else:
            # If output files/directories are expected then they must exist
            # If output data is expected then it must be cached
            satisfied_output = (not writes_output_dir or exists(final_output_directory)) \
                and (not writes_output_file or exists(output_filepath)) \
                and (output != MISSING_VALUE_EXCEPTION)
            # If we already have the expected output then we can skip the task at all
            if satisfied_output:
                print(f'{GREY_HEADER}-> Task {self.flag} ({self.name}) already completed{COLOR_END}')
                return output
        # If we are at this point then we are missing some output so we must proceed to run the task
        # Use the final output directory instead of the incomplete one if exists
        # Note that we must check if it exists again since it may have been deleted since the last check
        if writes_output_dir and exists(final_output_directory):
            processed_args[OUTPUT_DIRECTORY_ARG] = final_output_directory
        # Create the incomplete output directory, if necessary
        missing_incomplete_output = writes_output_dir \
            and not exists(incomplete_output_directory) \
            and not exists(final_output_directory)
        if missing_incomplete_output: mkdir(incomplete_output_directory)
        # Finally call the function
        print(f'{GREEN_HEADER}-> Running task {self.flag} ({self.name}){COLOR_END}')
        # If the task is to be run again because an inputs changed then let the user know
        if any_input_changed and had_cache and not forced_overwrite:
            changes = ''.join([ '\n   - ' + inp for inp in changed_inputs ])
            print(f'{GREEN_HEADER}   The task is run again since the following inputs changed:{changes}{COLOR_END}')
        # Save a few internal values the task although the task is static
        # We save it right before calling the function in case the function uses this task as input
        self.changed_inputs = changed_inputs
        self.cache_cksums = cache_cksums
        # Run the actual task
        output = self.func(**processed_args)
        self._set_parent_output(parent, output)
        # Update cache output unless it is marked to not save it
        if self.use_cache: parent.cache.update(self.cache_output_key, output)
        # Update the overwritables so this is not remade further in the same run
        parent.overwritables.discard(self.flag)
        # As a brief cleanup, if the output directory is empty at the end, then remove it
        # Otheriwse, change the incomplete directory name to its final name
        if writes_output_dir and exists(incomplete_output_directory):
            if is_directory_empty(incomplete_output_directory): rmtree(incomplete_output_directory)
            else: rename(incomplete_output_directory, final_output_directory)
        # Now return the function result
        return output

    # Find argument values, thus running any dependency
    def find_arg_value (self, arg : str, parent : Union['Project', 'MD'], default_arguments : set):
        # Word 'task' is reserved for getting the task itself
        if arg == 'task': return self
        # Word 'self' is reserved for getting the caller Project/MD
        if arg == 'self': return parent
        # Check if the argument is an MD property
        arg_value = safe_getattr(parent, arg, MISSING_ARGUMENT_EXCEPTION)
        if arg_value != MISSING_ARGUMENT_EXCEPTION: return arg_value
        # If the parent is an MD then it may happen the property is from the Project
        if isinstance(parent, MD):
            arg_value = safe_getattr(parent.project, arg, MISSING_ARGUMENT_EXCEPTION)
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
        raise RuntimeError(f'Function "{self.func.__name__}" from task "{self.flag}" expects argument "{arg}" but it is missing')
    
    # Find out if inputs changed regarding the last run
    def get_changed_inputs (self,
        parent : Union['Project', 'MD'],
        processed_args : dict) -> Tuple[ List[str], bool ]:
        # Get cache argument references
        cache_cksums = parent.cache.retrieve(self.cache_arg_cksums, MISSING_VALUE_EXCEPTION)
        had_cache = False if cache_cksums == MISSING_VALUE_EXCEPTION else True
        if cache_cksums == MISSING_VALUE_EXCEPTION: cache_cksums = {}
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
            old_cksum = cache_cksums.get(arg_name, None)
            # Compare new and old cksums
            if new_cksum != old_cksum:
                # If we found a missmatch then add it to the list
                unmatched_arguments.append(arg_name)
                # Update the references
                cache_cksums[arg_name] = new_cksum
        return unmatched_arguments, had_cache, cache_cksums
    
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
        # Make sure the MD directory does not equal the project directory
        if normpath(self.directory) == normpath(self.project.directory):
            raise InputError(f'MD {self.number} has the same directory as the project: {self.directory}')
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
        self._structure = None

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
        # Save also warnings apart since they are to be used as an input for metadata tasks
        self.warnings = self.register.warnings

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
                    name = md.get('name', None)
                    if not name: raise InputError('There is a MD with no name and no directory. Please define at least one of them.')
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
    
    # Make a summary of tests and their status
    def print_tests_summary (self):
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
    # This is run after processing input files
    def _issue_required_test_warnings (self):
        for test_name in AVAILABLE_CHECKINGS:
            # If test was not skipped then proceed
            if test_name not in self.project.trust: continue
            # If test passed in a previous run the proceed
            test_result = self.register.tests.get(test_name)
            if test_result == True: continue
            # If test failed in a previous run we can also proceed
            # The failing warning must be among the inherited warnings, so there is no need to add more warnings here
            if test_result == False: continue
            # If the test has been always skipped then issue a warning
            if test_result == None:
                # Remove previous warnings
                self.register.remove_warnings(test_name)
                # Get test pretty name
                test_nice_name = NICE_NAMES[test_name]
                # Issue the corresponding warning            
                self.register.add_warning(test_name, test_nice_name + ' was skipped and never run before')
                continue
            raise ValueError('Test value is not supported')

    # Processed files ----------------------------------------------------

    # Run the actual processing to generate output processed files out of input raw files
    # And by "files" I mean structure, trajectory and topology
    process_input_files = Task('inpro', 'Process input files', process_input_files)

    def get_structure_file (self) -> str:
        """Get the processed structure."""
        # If we have a stored value then return it
        # This means we already found or generated this file
        if self._structure_file:
            return self._structure_file
        # Set the file
        structure_filepath = self.pathify(STRUCTURE_FILENAME)
        self._structure_file = File(structure_filepath)
        # If the faith flag was passed then simply make sure the input file makes sense
        if self.project.faith:
            if self.input_structure_file != self._structure_file:
                raise InputError('Input structure file is not equal to output structure file but the "faith" flag was used.\n'
                    '  Please refrain from using the faith argument (-f) if you ignore its effect.')
            if not self.input_structure_file.exists:
                raise InputError('Input structure file does not exist but the "faith" flag was used.\n'
                    '  Please refrain from using the faith argument (-f) if you ignore its effect.')
            return self._structure_file
        # Run the processing logic
        self.process_input_files(self)
        # Now that the file is sure to exist we return it
        return self._structure_file
    structure_file = property(get_structure_file, None, None, "Structure file (read only)")

    def get_trajectory_file (self) -> str:
        """Get the processed trajectory."""
        # If we have a stored value then return it
        # This means we already found or generated this file
        if self._trajectory_file:
            return self._trajectory_file
        # If the file already exists then we are done
        trajectory_filepath = self.pathify(TRAJECTORY_FILENAME)
        self._trajectory_file = File(trajectory_filepath)
        # If the faith flag was passed then simply make sure the input file makes sense
        if self.project.faith:
            if len(self.input_trajectory_files) > 1:
                raise InputError('Several input trajectory files but the "faith" flag was used.\n'
                    '  Please refrain from using the faith argument (-f) if you ignore its effect.')
            sample = self.input_trajectory_files[0]
            if sample != self._trajectory_file:
                raise InputError('Input trajectory file is not equal to output trajectory file but the "faith" flag was used.\n'
                    '  Please refrain from using the faith argument (-f) if you ignore its effect.')
            if not self._trajectory_file.exists:
                raise InputError('Input trajectory file does not exist but the "faith" flag was used.\n'
                    '  Please refrain from using the faith argument (-f) if you ignore its effect.')
            return self._trajectory_file
        # Run the processing logic
        self.process_input_files(self)
        # Now that the file is sure to exist we return it
        return self._trajectory_file
    trajectory_file = property(get_trajectory_file, None, None, "Trajectory file (read only)")

    # Get the processed topology from the project
    def get_topology_file (self) -> str:
        return self.project.get_topology_file()
    topology_file = property(get_topology_file, None, None, "Topology filename from the project (read only)")

    # ---------------------------------------------------------------------------------
    # Others values which may be found/calculated and files to be generated on demand
    # ---------------------------------------------------------------------------------

    # Trajectory snapshots
    get_snapshots = Task('frames', 'Count trajectory frames', get_frames_count)
    snapshots = property(get_snapshots, None, None, "Trajectory snapshots (read only)")

    # Safe bonds
    def get_reference_bonds (self) -> List[List[int]]:
        return self.project.reference_bonds
    reference_bonds = property(get_reference_bonds, None, None, "Atom bonds to be trusted (read only)")

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
            self._structure.bonds = self.reference_bonds
        # Same procedure if we have coarse grain atoms
        elif self.cg_selection:
            self._structure.bonds = self.reference_bonds
        return self._structure
    structure = property(get_structure, None, None, "Parsed structure (read only)")

    # First frame PDB file
    get_first_frame = Task('firstframe', 'Get first frame structure',
        get_first_frame, output_filename = FIRST_FRAME_FILENAME)
    get_first_frame_file = get_first_frame.get_output_file
    first_frame_file = property(get_first_frame_file, None, None, "First frame (read only)")

    # Average structure filename
    get_average_structure = Task('average', 'Get average structure',
        get_average_structure, output_filename = AVERAGE_STRUCTURE_FILENAME)
    get_average_structure_file = get_average_structure.get_output_file
    average_structure_file = property(get_average_structure_file, None, None, "Average structure filename (read only)")

    # Produce the MD metadata file to be uploaded to the database
    prepare_metadata = Task('mdmeta', 'Prepare MD metadata',
        generate_md_metadata, output_filename=OUTPUT_METADATA_FILENAME)

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
    def _set_pbc_selection (self, reference_structure : 'Structure', verbose : bool = False) -> 'Selection':
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
    def _set_cg_selection (self, reference_structure : 'Structure', verbose : bool = False) -> 'Selection':
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
            if verbose: print(' Empty selection -> There is no CG at all')
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
    get_reference_frame = Task('reframe', 'Reference frame', get_bonds_canonical_frame)
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

    # RMSDs analysis
    run_rmsds_analysis = Task('rmsds', 'RMSDs analysis',
        rmsds, { 'frames_limit': 5000 })

    # TM scores analysis
    run_tmscores_analysis = Task('tmscore', 'TM scores analysis',
        tmscores, { 'frames_limit': 200 })

    # RMSF, atom fluctuation analysis
    run_rmsf_analysis = Task('rmsf', 'Fluctuation analysis', rmsf)

    # Radius of gyration analysis
    run_rgyr_analysis = Task('rgyr', 'Radius of gyration analysis',
        rgyr, { 'frames_limit': 5000 })

    # PCA, principal component analysis
    run_pca_analysis = Task('pca', 'Principal component analysis',
        pca, { 'frames_limit': 2000, 'projection_frames': 20 })

    # PCA contacts
    # DANI: Intenta usar mucha memoria, hay que revisar
    # DANI: Puede saltar un error de imposible alojar tanta memoria
    # DANI: Puede comerse toda la ram y que al final salte un error de 'Terminado (killed)'
    # DANI: De momento me lo salto
    # DANI: Lleva mucho tiempo sin mantenerse, habrá que cambiar varias cosas al recuperarlo
    # run_pca_contacts('pcacons', 'PCA contacts', pca_contacts)

    # RMSD per residue analysis
    run_rmsd_perres_analysis = Task('perres', 'RMSD per residue analysis',
        rmsd_per_residue, { 'frames_limit': 100 })

    # RMSD pairwise
    # Perform an analysis for the overall structure and then one more analysis for each interaction
    run_rmsd_pairwise_analysis = Task('pairwise', 'RMSD pairwise',
        rmsd_pairwise, { 'frames_limit': 200, 'overall_selection': "name CA or name C5'" })

    # Run the cluster analysis
    run_clusters_analysis = Task('clusters', 'Clusters analysis',
        clusters_analysis, { 'frames_limit': 1000, 'desired_n_clusters': 20 })

    # Calculate the distance mean and standard deviation of each pair of residues
    run_dist_perres_analysis = Task('dist', 'Distance per residue',
        distance_per_residue, { 'frames_limit': 200 })

    # Hydrogen bonds
    run_hbonds_analysis = Task('hbonds', 'Hydrogen bonds analysis',
        hydrogen_bonds, { 'time_splits': 100 })

    # SASA, solvent accessible surface analysis
    run_sas_analysis = Task('sas', 'Solvent accessible surface analysis',
        sasa, { 'frames_limit': 100 })
    
    # Perform the electrostatic and vdw energies analysis for each pair of interaction agents
    run_energies_analysis = Task('energies', 'Energies analysis',
        energies, { 'frames_limit': 100 })

    # Calculate torsions and then dihedral energies for every dihedral along the trajectory
    run_dihedral_energies = Task('dihedrals', 'Dihedral energies analysis',
        compute_dihedral_energies, { 'frames_limit': 100 })

    # Perform the pockets analysis
    run_pockets_analysis = Task('pockets', 'Pockets analysis',
        pockets, { 'frames_limit': 100, 'maximum_pockets_number': 10 })

    # Helical parameters
    run_helical_analysis = Task('helical', 'Helical parameters', helical_parameters)
        
    # Markov
    run_markov_analysis = Task('markov', 'Markov', markov, { 'rmsd_selection': PROTEIN_AND_NUCLEIC })

    # Membrane density analysis
    run_density_analysis = Task('density', 'Membrane density analysis',
        density, { 'frames_limit': 1000 })

    # Membrane thickness analysis
    run_thickness_analysis = Task('thickness', 'Membrane thickness analysis',
        thickness, { 'frames_limit': 100 })

    # Area per lipid analysis
    run_apl_analysis = Task('apl', 'Membrane area per lipid analysis', area_per_lipid)

    # Calculate lipid order parameters for membranes
    run_lipid_order_analysis = Task('lorder', 'Membrane lipid order analysis',
        lipid_order, { 'frames_limit': 100 })

    # Lipid-protein interactions analysis
    run_lipid_interactions_analysis = Task('linter', 'Membrane lipid-protein interactions analysis',
        lipid_interactions, { 'frames_limit': 100 })
        
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
        faith : bool = False,
        pca_analysis_selection : str = PROTEIN_AND_NUCLEIC_BACKBONE,
        pca_fit_selection : str = PROTEIN_AND_NUCLEIC_BACKBONE,
        rmsd_cutoff : float = DEFAULT_RMSD_CUTOFF,
        interaction_cutoff : float = DEFAULT_INTERACTION_CUTOFF,
        interactions_auto : Optional[str] = None,
        guess_bonds : bool = False,
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
        if md_config and (md_directories or input_trajectory_filepaths):
            raise InputError('MD configurations (-md) is not compatible with old MD inputs (-mdir, -traj)')
        # Save the MD configurations
        self.md_config = md_config
        # Make sure MD configuration has the correct format
        if self.md_config:
            # Make sure all MD configurations have at least 3 values each
            for mdc in self.md_config:
                if len(mdc) < 2:
                    raise InputError('Wrong MD configuration: the patter is -md <directory> <trajectory> <trajectory 2> ...')
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
        self.faith = faith
        self.pca_analysis_selection = pca_analysis_selection
        self.pca_fit_selection = pca_fit_selection
        self.rmsd_cutoff = rmsd_cutoff
        self.interaction_cutoff = interaction_cutoff
        self.sample_trajectory = sample_trajectory
        self.interactions_auto = interactions_auto
        self.guess_bonds = guess_bonds
        # Set the inputs, where values from the inputs file will be stored
        self._inputs = None

        # Other values which may be found/calculated on demand
        self._pbc_selection = None
        self._pbc_residues = None
        self._cg_selection = None
        self._cg_residues = None
        self._reference_bonds = None
        self._topology_reader = None
        self._dihedrals = None
        self._populations = None
        self._transitions = None
        self._pdb_ids = None
        self._mds = None

        # Force a couple of extraordinary files which is generated if atoms are resorted
        self.resorted_bonds_file = File(self.pathify(RESORTED_BONDS_FILENAME))
        self.resorted_charges_file = File(self.pathify(RESORTED_CHARGES_FILENAME))

        # Set a new entry for the register
        # This is useful to track previous workflow runs and problems
        register_filepath = self.pathify(REGISTER_FILENAME)
        register_file = File(register_filepath)
        self.register = Register(register_file)
        # Save also warnings apart since they are to be used as an input for metadata tasks
        self.warnings = self.register.warnings

        # Set the cache
        cache_filepath = self.pathify(CACHE_FILENAME)
        cache_file = File(cache_filepath)
        self.cache = Cache(cache_file)

        # Set tasks whose output is to be overwritten
        self.overwritables = set()

    def __repr__ (self):
        return 'Project'

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
                directory = config[0]
                # LEGACY 
                # In a previous version, the md config argument also holded the structure
                # This was the second argument, so we check if we have more than 2 arguments
                # If this is the case, then check if the second argument has different format
                # Note that PDB format is also a trajectory supported format
                has_structure = False
                if len(config) > 2:
                    first_sample = File(config[1])
                    second_sample = File(config[2])
                    if first_sample.format != second_sample.format:
                        has_structure = True
                # Finally set the input structure and trajectories
                input_structure_filepath = config[1] if has_structure else self.input_structure_filepath
                input_trajectory_filepaths = config[2:] if has_structure else config[1:]
                # Define the MD
                md = MD(
                    project = self, number = n, directory = directory,
                    input_structure_filepath = input_structure_filepath,
                    input_trajectory_filepaths = input_trajectory_filepaths,
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
        """Set a function to load the inputs file."""
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
    input_protein_references = property(input_getter('forced_references'), None, None, "Uniprot IDs to be used first when aligning protein sequences (read only)")
    input_pdb_ids = property(input_getter('pdb_ids'), None, None, "Protein Data Bank IDs used for the setup of the system (read only)")
    input_type = property(input_getter('type'), None, None, "Set if its a trajectory or an ensemble (read only)")
    input_mds = property(input_getter('mds'), None, None, "Input MDs configuration (read only)")
    input_ligands = property(input_getter('ligands'), None, None, "Input ligand references (read only)")
    input_force_fields = property(input_getter('ff'), None, None, "Input force fields (read only)")
    input_collections = property(input_getter('collections'), None, None, "Input collections (read only)")
    input_chain_names = property(input_getter('chainnames'), None, None, "Input chain names (read only)")
    input_framestep = property(input_getter('framestep'), None, None, "Input framestep (read only)")
    input_name = property(input_getter('name'), None, None, "Input name (read only)")
    input_description = property(input_getter('description'), None, None, "Input description (read only)")
    input_authors = property(input_getter('authors'), None, None, "Input authors (read only)")
    input_groups = property(input_getter('groups'), None, None, "Input groups (read only)")
    input_contact = property(input_getter('contact'), None, None, "Input contact (read only)")
    input_program = property(input_getter('program'), None, None, "Input program (read only)")
    input_version = property(input_getter('version'), None, None, "Input version (read only)")
    input_method = property(input_getter('method'), None, None, "Input method (read only)")
    input_license = property(input_getter('license'), None, None, "Input license (read only)")
    input_linkcense = property(input_getter('linkcense'), None, None, "Input license link (read only)")
    input_citation = property(input_getter('citation'), None, None, "Input citation (read only)")
    input_thanks = property(input_getter('thanks'), None, None, "Input acknowledgements (read only)")
    input_links = property(input_getter('links'), None, None, "Input links (read only)")
    input_timestep = property(input_getter('timestep'), None, None, "Input timestep (read only)")
    input_temperature = property(input_getter('temp'), None, None, "Input temperature (read only)")
    input_ensemble = property(input_getter('ensemble'), None, None, "Input ensemble (read only)")
    input_water = property(input_getter('wat'), None, None, "Input water force field (read only)")
    input_boxtype = property(input_getter('boxtype'), None, None, "Input boxtype (read only)")
    input_pbc_selection = property(input_getter('pbc_selection'), None, None, "Input Periodic Boundary Conditions (PBC) selection (read only)")
    input_cg_selection = property(input_getter('cg_selection'), None, None, "Input Coarse Grained (CG) selection (read only)")
    input_customs = property(input_getter('customs'), None, None, "Input custom representations (read only)")
    input_orientation = property(input_getter('orientation'), None, None, "Input orientation (read only)")
    input_multimeric = property(input_getter('multimeric'), None, None, "Input multimeric labels (read only)")
    # Additional topic-specific inputs
    input_cv19_unit = property(input_getter('cv19_unit'), None, None, "Input Covid-19 Unit (read only)")
    input_cv19_startconf = property(input_getter('cv19_startconf'), None, None, "Input Covid-19 starting conformation (read only)")
    input_cv19_abs = property(input_getter('cv19_abs'), None, None, "Input Covid-19 antibodies (read only)")
    input_cv19_nanobs = property(input_getter('cv19_nanobs'), None, None, "Input Covid-19 nanobodies (read only)")

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
        raise InputError(f'Not supported input "type" value: {self.input_type}. It must be "trajectory" or "ensemble"')
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
        inherited_filename = self.inherit_topology_filename()
        self._topology_filepath = self.pathify(inherited_filename) if inherited_filename else None
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
        # If the faith flag was passed then simply make sure the input file makes sense
        if self.faith:
            if self.input_topology_file != self._topology_file:
                raise InputError('Input topology file is not equal to output topology file but the "faith" flag was used.\n'
                    '  Please refrain from using the faith argument (-f) if you ignore its effect.')
            if not self.input_topology_file.exists:
                raise InputError('Input topology file does not exist but the "faith" flag was used.\n'
                    '  Please refrain from using the faith argument (-f) if you ignore its effect.')
            return self._topology_file
        # Run the processing logic
        self.reference_md.process_input_files(self.reference_md)
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

    # Reference MD spanshots
    # Used next to the reference MD trajectory data
    def get_snaphsots (self) -> int:
        return self.reference_md.snapshots
    snapshots = property(get_snaphsots, None, None, "Reference MD snapshots (read only)")
    
    # Check if we must check stable bonds
    def get_check_stable_bonds (self) -> bool:
        # Set if stable bonds have to be checked
        must_check = STABLE_BONDS_FLAG not in self.trust
        # If this analysis has been already passed then we can trust structure bonds
        if self.register.tests.get(STABLE_BONDS_FLAG, None) == True:
            must_check = False
        return must_check
    must_check_stable_bonds = property(get_cg_residues, None, None, "Indices of residues in coarse grain (read only)")

    # Reference bonds
    get_reference_bonds = Task('refbonds', 'Reference bonds', find_safe_bonds)
    reference_bonds = property(get_reference_bonds, None, None, "Atom bonds to be trusted (read only)")

    # Atom charges
    get_charges = Task('charges', 'Getting atom charges', get_charges)
    charges = property(get_charges, None, None, "Atom charges (read only)")

    # Topolody data reader
    def get_topology_reader (self) -> 'Topology':
        # If we already have a stored value then return it
        if self._topology_reader: return self._topology_reader
        # Instantiate the topology reader
        self._topology_reader = Topology(self.topology_file)
        return self._topology_reader
    topology_reader = property(get_topology_reader, None, None, "Topology reader (read only)")

    def get_dihedrals (self) -> List[dict]:
        """Get the topology dihedrals."""
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
        # If input PDB ids is a string instead of a list then fix it
        input_pdb_ids = self.input_pdb_ids
        if type(input_pdb_ids) == str:
            input_pdb_ids = [ input_pdb_ids ]
        # Iterate input PDB ids
        for input_pdb_id in input_pdb_ids:
            # First make sure this is a PDB id
            if not re.match(PDB_ID_FORMAT, input_pdb_id):
                raise InputError(f'Input PDB id "{input_pdb_id}" does not look like a PDB id')
            # Make letters upper
            pdb_id = input_pdb_id.upper()
            self._pdb_ids.append(pdb_id)
        return self._pdb_ids
    pdb_ids = property(get_pdb_ids, None, None, "Tested and standarized PDB ids (read only)")

    # Prepare the PDB references file to be uploaded to the database
    get_pdb_references = Task('pdbs', 'Prepare PDB references',
        prepare_pdb_references, output_filename = PDB_REFERENCES_FILENAME)
    pdb_references_file = get_pdb_references.get_output_file

    # Map the structure aminoacids sequences against the Uniprot reference sequences
    get_protein_map = Task('protmap', 'Protein residues mapping',
        generate_protein_mapping, output_filename=PROTEIN_REFERENCES_FILENAME)
    protein_map = property(get_protein_map, None, None, "Protein residues mapping (read only)")

    # Define the output file of the protein mapping including protein references
    get_protein_references_file = get_protein_map.get_output_file
    protein_references_file = property(get_protein_references_file, None, None, "File including protein refereces data mined from UniProt (read only)")

    # Get chain references
    get_chain_references = Task('chains', 'Chain references',
        prepare_chain_references, output_filename = OUTPUT_CHAINS_FILENAME)

    # Get the ligand residues mapping
    get_ligand_map = Task('ligmap', 'Ligand residues mapping',
        generate_ligand_mapping, output_filename = LIGAND_REFERENCES_FILENAME)
    ligand_map = property(get_ligand_map, None, None, "Ligand references (read only)")

    # Define the output file of the ligand mapping including ligand references
    get_ligand_references_file = get_ligand_map.get_output_file
    ligand_references_file = property(get_ligand_references_file, None, None, "File including ligand refereces data mined from PubChem (read only)")

    # MDAnalysis Universe object
    get_MDAnalysis_Universe = Task('mda_univ', 'MDAnalysis Universe object',
        get_mda_universe, use_cache = False)
    universe = property(get_MDAnalysis_Universe, None, None, "MDAnalysis Universe object (read only)")

    # Get the lipid references
    get_lipid_map = Task('lipmap', 'Lipid mapping',
        generate_lipid_references, output_filename = LIPID_REFERENCES_FILENAME)
    lipid_map = property(get_lipid_map, None, None, "Lipid mapping (read only)")

    # Define the output file of the lipid mapping including lipid references
    get_lipid_references_file = get_lipid_map.get_output_file
    lipid_references_file = property(get_lipid_references_file, None, None, "File including lipid references data mined from PubChem (read only)")

    # Get mapping of residues in the membrane
    get_membrane_map = Task('membranes', 'Membrane mapping',
        generate_membrane_mapping, output_filename = MEMBRANE_MAPPING_FILENAME) 
    membrane_map = property(get_membrane_map, None, None, "Membrane mapping (read only)")

    # Build the residue map from both proteins and ligands maps
    # This is formatted as both the standard topology and metadata producers expect them
    get_residue_map = Task('resmap', 'Residue mapping', generate_residue_mapping)
    residue_map = property(get_residue_map, None, None, "Residue map (read only)")

    # Get interaction types
    get_interaction_types = Task('intertypes', 'Finding interaction types', find_interaction_types)
    interaction_types = property(get_interaction_types, None, None, "Interaction types (read only)")

    # Prepare the project metadata file to be upload to the database
    prepare_metadata = Task('pmeta', 'Prepare project metadata',
        prepare_project_metadata, output_filename=OUTPUT_METADATA_FILENAME)

    # Prepare the standard topology file to be uploaded to the database
    prepare_standard_topology = Task('stopology', 'Standard topology file',
        generate_topology, output_filename = STANDARD_TOPOLOGY_FILENAME)
    get_standard_topology_file = prepare_standard_topology.get_output_file
    standard_topology_file = property(get_standard_topology_file, None, None, "Standard topology filename (read only)")

    # Get a screenshot of the system
    get_screenshot_filename = Task('screenshot', 'Screenshot file',
        get_screenshot, output_filename = OUTPUT_SCREENSHOT_FILENAME)
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
    # Remove a possible starting './'
    # Replace white spaces with underscores
    name = directory.split('/')[-1].replace('_', ' ')
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

# Project requestable tasks
project_requestables = {
    **project_input_files,
    **project_processed_files,
}
# Add available tasks to project requestables
for callable in vars(Project).values():
    if isinstance(callable, Task): project_requestables[callable.flag] = callable
# MD requestable tasks
md_requestables = {
    **md_input_files,
    **md_processed_files,
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
    'meta': ['pmeta', 'mdmeta'],
    'network': [ 'mapping', 'ligands', 'chains', 'pdbs', 'membrane' ],
    'minimal': [ 'pmeta', 'mdmeta', 'stopology' ],
    'interdeps': [ 'interactions', 'pairwise', 'hbonds', 'energies', 'perres', 'clusters', 'dist' ],
    'membs': ['membranes', 'density',  'thickness', 'apl', 'lorder', 'linter']
}

# Set the default analyses to be run when no task is specified
DEFAULT_ANALYSES = ['clusters', 'dist', 'energies', 'hbonds', 'pca', 'pockets',
    'rgyr', 'rmsds', 'perres', 'pairwise', 'rmsf', 'sas', 'tmscore', 'density',
    'thickness', 'apl', 'lorder', 'linter']

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
            *DEFAULT_ANALYSES,
        ]
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
