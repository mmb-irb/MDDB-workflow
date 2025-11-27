# This is the starter script

from os import chdir, walk, mkdir, getcwd
from os.path import exists, isdir, isabs, relpath, normpath, split, basename
import sys
import io
import re
import numpy
from glob import glob

# Constants
# Importing constants first is important
from mddb_workflow.utils.constants import *

# Import local utils
from mddb_workflow.utils.auxiliar import InputError, MISSING_TOPOLOGY
from mddb_workflow.utils.auxiliar import warn, load_json, save_json, load_yaml, save_yaml
from mddb_workflow.utils.auxiliar import is_glob, parse_glob, is_url, url_to_source_filename
from mddb_workflow.utils.auxiliar import read_ndict, write_ndict, get_git_version, download_file
from mddb_workflow.utils.auxiliar import is_standard_topology
from mddb_workflow.utils.register import Register
from mddb_workflow.utils.cache import Cache
from mddb_workflow.utils.structures import Structure
from mddb_workflow.utils.topologies import Topology
from mddb_workflow.utils.file import File
from mddb_workflow.utils.remote import Remote
from mddb_workflow.utils.pyt_spells import get_frames_count, get_average_structure
from mddb_workflow.utils.selections import Selection
from mddb_workflow.utils.mda_spells import get_mda_universe
from mddb_workflow.utils.tasks import Task
from mddb_workflow.utils.type_hints import *

# Import local tools
from mddb_workflow.tools.get_first_frame import get_first_frame
from mddb_workflow.tools.get_bonds import find_safe_bonds, get_bonds_reference_frame
from mddb_workflow.tools.process_interactions import process_interactions
from mddb_workflow.tools.generate_metadata import prepare_project_metadata, generate_md_metadata
from mddb_workflow.tools.chains import prepare_chain_references
from mddb_workflow.tools.generate_pdb_references import prepare_pdb_references
from mddb_workflow.tools.generate_map import generate_protein_mapping
from mddb_workflow.tools.get_inchi_keys import generate_inchikeys, generate_inchi_references
from mddb_workflow.tools.get_lipids import generate_lipid_references
from mddb_workflow.tools.membrane_mapping import generate_membrane_mapping
from mddb_workflow.tools.generate_ligands_desc import generate_ligand_references
from mddb_workflow.tools.residue_mapping import generate_residue_mapping
from mddb_workflow.tools.generate_topology import generate_topology
from mddb_workflow.tools.get_charges import get_charges
from mddb_workflow.tools.remove_trash import remove_trash
from mddb_workflow.tools.get_screenshot import get_screenshot
from mddb_workflow.tools.process_input_files import process_input_files
from mddb_workflow.tools.provenance import produce_provenance
from mddb_workflow.tools.get_reduced_trajectory import calculate_frame_step

# Import local analyses
from mddb_workflow.analyses.rmsds import rmsds
from mddb_workflow.analyses.tmscores import tmscores
from mddb_workflow.analyses.rmsf import rmsf
from mddb_workflow.analyses.rgyr import rgyr
from mddb_workflow.analyses.pca import pca
from mddb_workflow.analyses.density import density
from mddb_workflow.analyses.thickness import thickness
from mddb_workflow.analyses.area_per_lipid import area_per_lipid
from mddb_workflow.analyses.lipid_order import lipid_order
from mddb_workflow.analyses.lipid_interactions import lipid_interactions
from mddb_workflow.analyses.channels import channels
# from mddb_workflow.analyses.pca_contacts import pca_contacts
from mddb_workflow.analyses.rmsd_per_residue import rmsd_per_residue
from mddb_workflow.analyses.rmsd_pairwise import rmsd_pairwise
from mddb_workflow.analyses.clusters import clusters_analysis
from mddb_workflow.analyses.distance_per_residue import distance_per_residue
from mddb_workflow.analyses.hydrogen_bonds import hydrogen_bonds
from mddb_workflow.analyses.sasa import sasa
from mddb_workflow.analyses.energies import energies
from mddb_workflow.analyses.dihedral_energies import compute_dihedral_energies
from mddb_workflow.analyses.pockets import pockets
from mddb_workflow.analyses.rmsd_check import check_trajectory_integrity
from mddb_workflow.analyses.helical_parameters import helical_parameters
from mddb_workflow.analyses.markov import markov

# Make the system output stream to not be buffered
# Only when not in a Jupyter notebook or using pytest
# Check if we're in an interactive Python shell like Jupyter
if not hasattr(sys, 'ps1') and 'pytest' not in sys.modules:
    # This is useful to make prints work on time in Slurm
    # Otherwise, output logs are written after the script has fully run
    # Note that this fix affects all modules and built-ins
    unbuffered = io.TextIOWrapper(open(sys.stdout.fileno(), 'wb', 0), write_through=True)
    sys.stdout = unbuffered

# Set a special exception for missing inputs
MISSING_INPUT_EXCEPTION = Exception('Missing input')


class MD:
    """A Molecular Dynamics (MD) is the union of a structure and a trajectory.
    Having this data several analyses are possible.
    Note that an MD is always defined inside of a Project and thus it has additional topology and metadata.
    """

    def __init__(self,
        project: 'Project',
        number: int,
        directory: str,
        input_structure_filepath: str,
        input_trajectory_filepaths: list[str],
    ):
        """Initialize the MD object.

        Args:
            project (Project): The parent project this MD belongs to.
            number (int): The number of the MD according to its accession.
            directory (str): The local directory where the MD takes place.
            input_structure_filepath (str): The input structure file path.
            input_trajectory_filepaths (list[str]): The input trajectory file paths.

        """
        # Save the inputs
        self.project = project
        if not project:
            raise Exception('Project is mandatory to instantiate a new MD')
        # Save the MD number and index
        self.number = number
        self.index = number - 1
        # Set the MD accession and request URL
        self.accession = None
        self.remote = None
        if self.project.database_url and self.project.accession:
            self.accession = f'{self.project.accession}.{self.number}'
            self.remote = Remote(self.project.database_url, self.accession)
        # Save the directory
        self.directory = normpath(directory)
        # Now set the director relative to the project
        self.directory = self.project.pathify(self.directory)
        if normpath(self.directory) == normpath(self.project.directory):
            raise InputError(f'MD {self.number} has the same directory as the project: {self.directory}')
        # Save the directory name alone apart
        self.directory_location, self.directory_name = split(self.directory)
        # If the directory does not exists then create it
        if not exists(self.directory):
            mkdir(self.directory)
        # Save the input structure filepath
        # They may be relative to the project directory (unique) or relative to the MD directory (one per MD)
        # If the path is absolute then it is considered unique
        # If the file does not exist and it is to be downloaded then it is downloaded for each MD
        # Priorize the MD directory over the project directory
        self.arg_input_structure_filepath = input_structure_filepath
        self._input_structure_filepath = None
        self._input_structure_url = None
        # Set the internal variable for the input structure file, to be assigned later
        self._input_structure_file = None
        # Save the input trajectory filepaths
        self.arg_input_trajectory_filepaths = input_trajectory_filepaths
        self._input_trajectory_filepaths = None
        self._input_trajectory_urls = None
        # Set the internal variable for the input trajectory files, to be assigned later
        self._input_trajectory_files = None

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

        # Get MD inputs just to fill the inputs' "mds" value
        # Some functions may fail otherwise when its value is missing
        if self.project.is_inputs_file_available():
            self.get_md_inputs()

    def __repr__(self):
        """Return a string representation of the MD object."""
        return 'MD'

    def pathify(self, filename_or_relative_path: str) -> str:
        """Given a filename or relative path, add the MD directory path at the beginning."""
        return normpath(self.directory + '/' + filename_or_relative_path)

    # Input structure file ------------

    def get_input_structure_filepath(self) -> str:
        """Set a function to get input structure file path."""
        # Return the internal value if it is already assigned
        if self._input_structure_filepath is not None:
            return self._input_structure_filepath

        def relativize_and_parse_paths(input_path: str, may_not_exist: bool = False) -> Optional[str]:
            """Find out if a path is relative to MD directories or to the project directory.

            To do so just check if the file exists in any of those.
            In case it exists in both or none then assume it is relative to MD directory.
            Parse glob notation in the process.
            """
            # Check if it is a URL
            if is_url(input_path):
                self._input_structure_url = input_path
                # Set the paths for the further download
                source_filename = url_to_source_filename(input_path)
                md_relative_filepath = self.pathify(source_filename)
                return md_relative_filepath
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
            md_parsed_filepath = self.project.pathify(input_path) if f'{self.directory_name}/' in md_parsed_filepath else self.pathify(input_path)
            return md_parsed_filepath
        # If we have a value passed through command line
        if self.arg_input_structure_filepath:
            # Find out if it is relative to MD directories or to the project directory
            self._input_structure_filepath = relativize_and_parse_paths(self.arg_input_structure_filepath)
            # Save the parsed value in the inputs file
            self.project.update_inputs(
                f'mds.{self.index}.input_structure_filepath',
                self._input_structure_filepath)
            return self._input_structure_filepath
        # If we have a value passed through the inputs file has the value
        if self.project.is_inputs_file_available():
            # Get the input value, whose key must exist
            inputs_value = self.get_input('input_structure_filepath')
            # If there is a valid input then use it
            if inputs_value:
                self._input_structure_filepath = relativize_and_parse_paths(inputs_value)
                return self._input_structure_filepath
        # If there is not input structure anywhere then use the input topology
        # We will extract the structure from it using a sample frame from the trajectory
        # Note that topology input filepath must exist and an input error will raise otherwise
        # However if we are using the standard topology file we can not extract the PDB from it (yet)
        if self.project.input_topology_file != MISSING_TOPOLOGY and \
            not is_standard_topology(self.project.input_topology_file):
            return self.project.input_topology_file.path
        # If we can not use the topology either then surrender
        raise InputError('There is not input structure at all')

    def get_input_structure_file(self) -> str:
        """Get the input pdb filename from the inputs.
        If the file is not found try to download it.
        """
        # If the input structure file is already defined then return it
        if self._input_structure_file:
            return self._input_structure_file
        # Otherwise we must set it
        # First set the input structure filepath
        input_structure_filepath = self.get_input_structure_filepath()
        # Now set the input structure file
        input_structure_file = File(input_structure_filepath)
        # If there is a structure URL then we must donwload the structure first
        input_structure_url = self._input_structure_url
        if input_structure_url and not input_structure_file.exists:
            original_filename = input_structure_url.split('/')[-1]
            # If there is a remote then use it
            if self.remote:
                # If the structure filename is the standard structure filename then use the structure endpoint instead
                if original_filename == STRUCTURE_FILENAME:
                    self.remote.download_standard_structure(input_structure_file)
                # Otherwise download the input strucutre file by its filename
                else:
                    self.remote.download_file(original_filename, input_structure_file)
            # Otherwise use the URL as is
            else:
                download_file(input_structure_url, input_structure_file)
        # If the file already exists then return it
        if input_structure_file.exists:
            self._input_structure_file = input_structure_file
            return self._input_structure_file
        raise InputError(f'Missing input structure file "{input_structure_file.path}"')
    input_structure_file = property(get_input_structure_file, None, None, "Input structure filename (read only)")

    # Input trajectory filename ------------

    def get_input_trajectory_filepaths(self) -> str:
        """Get the input trajectory file paths."""
        # Return the internal value if it is already assigned
        if self._input_trajectory_filepaths is not None:
            return self._input_trajectory_filepaths

        def relativize_and_parse_paths(input_paths: list[str]) -> list[str]:
            """Check and fix input trajectory filepaths.
            Also relativize paths to the current MD directory and parse glob notation.
            """
            checked_paths = input_paths
            # Input trajectory filepaths may be both a list or a single string
            # However we must keep a list
            if type(checked_paths) is list:
                pass
            elif type(checked_paths) is str:
                checked_paths = [checked_paths]
            else:
                raise InputError('Input trajectory filepaths must be a list of strings or a string')
            # Make sure all or none of the trajectory paths are URLs
            url_count = sum([is_url(path) for path in checked_paths])
            if not (url_count == 0 or url_count == len(checked_paths)):
                raise InputError('All trajectory paths must be paths or URLs. Mixing is not supported')
            # In case trajectory paths are URLs
            if url_count > 0:
                self._input_trajectory_urls = checked_paths
                # Set the paths for the further download
                parsed_paths = []
                for path in checked_paths:
                    source_filename = url_to_source_filename(path)
                    md_relative_filepath = self.pathify(source_filename)
                    parsed_paths.append(md_relative_filepath)
                return parsed_paths
            # Make sure all or none of the trajectory paths are absolute
            abs_count = sum([isabs(path) for path in checked_paths])
            if not (abs_count == 0 or abs_count == len(checked_paths)):
                raise InputError('All trajectory paths must be relative or absolute. Mixing is not supported')

            def parse_all_glob(paths: list[str]) -> list[str]:
                """Glob-parse and merge all paths."""
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
            md_relative_paths = [self.pathify(path) for path in checked_paths]
            # In case there are glob characters we must parse the paths
            md_parsed_paths = parse_all_glob(md_relative_paths)
            # Check we successfully defined some trajectory file
            if len(md_parsed_paths) > 0:
                # If so, check at least one of the files do actually exist
                if any([File(path).exists for path in md_parsed_paths]):
                    return md_parsed_paths
            # If no trajectory files where found then asume they are relative to the project
            # Get paths relative to the project directory
            project_relative_paths = [self.project.pathify(path) for path in checked_paths]
            # In case there are glob characters we must parse the paths
            project_parsed_paths = parse_all_glob(project_relative_paths)
            # Check we successfully defined some trajectory file
            if len(project_parsed_paths) > 0:
                # If so, check at least one of the files do actually exist
                if any([File(path).exists for path in project_parsed_paths]):
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
            # Note that if input path was not glob based it will be both as project relative and MD relative
            if len(md_parsed_paths) == 0: raise ValueError('This should never happen')
            # If file is to be downloaded then we must make sure the path is relative to the project
            project_relative_paths = [
                self.project.pathify(path) if f'{self.directory_name}/' in path else self.pathify(path) for path in checked_paths
            ]
            return project_relative_paths
        # If we have a value passed through command line
        if self.arg_input_trajectory_filepaths:
            self._input_trajectory_filepaths = relativize_and_parse_paths(self.arg_input_trajectory_filepaths)
            # Save the parsed value in the inputs file
            self.project.update_inputs(
                f'mds.{self.index}.input_trajectory_filepaths',
                self._input_trajectory_filepaths)
            return self._input_trajectory_filepaths
        # Check if the inputs file has the value
        if self.project.is_inputs_file_available():
            # Get the input value
            inputs_value = self.get_input('input_trajectory_filepaths')
            if inputs_value:
                self._input_trajectory_filepaths = relativize_and_parse_paths(inputs_value)
                return self._input_trajectory_filepaths
        # If there is no trajectory available then we surrender
        raise InputError('There is not input trajectory at all')

    def get_input_trajectory_files(self) -> str:
        """Get the input trajectory filename(s) from the inputs.
        If file(s) are not found try to download it.
        """
        # If we already defined input trajectory files then return them
        if self._input_trajectory_files is not None:
            return self._input_trajectory_files
        # Otherwise we must set the input trajectory files
        input_trajectory_filepaths = self.get_input_trajectory_filepaths()
        input_trajectory_files = [File(path) for path in input_trajectory_filepaths]
        # If there are input trajectory URLs then download the trajectories first
        input_trajectory_urls = self._input_trajectory_urls
        if input_trajectory_urls:
            for trajectory_url, trajectory_file in zip(input_trajectory_urls, input_trajectory_files):
                # If the trajectory file already exists then skip it
                if trajectory_file.exists: continue
                original_filename = trajectory_url.split('/')[-1]
                # If there is a remote then use it
                if self.remote:
                    # If the trajectory filename is the standard trajectory filename then use the trajectory endpoint instead
                    if original_filename == TRAJECTORY_FILENAME:
                        frame_selection = None
                        if self.project.sample_trajectory:
                            remote_frames = self.remote.snapshots
                            maximum_desired_frames = self.project.sample_trajectory
                            step, final_frames = calculate_frame_step(remote_frames, maximum_desired_frames)
                            frame_selection = f'1:{remote_frames}:{step}'
                        self.remote.download_trajectory(trajectory_file, frame_selection=frame_selection, format='xtc')
                    # Otherwise download the input trajectory file by its filename
                    else:
                        if self.project.sample_trajectory:
                            raise InputError('The "-smp" argument is supported only when asking for the standard trajectory')
                        self.remote.download_file(original_filename, trajectory_file)
                # Otherwise use the URL as is
                else:
                    if self.project.sample_trajectory:
                        raise InputError('The "-smp" argument is supported only when using the "-proj" argument')
                    download_file(trajectory_url, trajectory_file)
        # Find missing trajectory files
        missing_input_trajectory_files = []
        for trajectory_file in input_trajectory_files:
            if not trajectory_file.exists:
                missing_input_trajectory_files.append(trajectory_file)
        # If all files already exists then we are done
        if len(missing_input_trajectory_files) == 0:
            self._input_trajectory_files = input_trajectory_files
            return self._input_trajectory_files
        missing_filepaths = [trajectory_file.path for trajectory_file in missing_input_trajectory_files]
        raise InputError('Missing input trajectory files: ' + ', '.join(missing_filepaths))
    input_trajectory_files = property(get_input_trajectory_files, None, None, "Input trajectory filenames (read only)")

    def get_md_inputs(self) -> dict:
        """Get MD specific inputs."""
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
                if self.project.pathify(directory) == self.directory:
                    self._md_inputs = md
                    return self._md_inputs
        # If this MD directory has not associated inputs then it means it was passed through command line
        # We set a new MD inputs for it
        new_md_name = directory_2_name(self.directory)
        self._md_inputs = {'name': new_md_name, 'mdir': self.directory}
        # Update the inputs file with the new MD inputs
        mds = self.project.inputs.get('mds', [])
        new_mds_inputs = [*mds, self._md_inputs]
        self.project.update_inputs('mds', new_mds_inputs)
        return self._md_inputs

    md_inputs = property(get_md_inputs, None, None, "MD specific inputs (read only)")

    def get_input(self, name: str):
        """Get a specific 'input' value from MD inputs."""
        value = self.md_inputs.get(name, MISSING_INPUT_EXCEPTION)
        # If we had a value then return it
        if value != MISSING_INPUT_EXCEPTION:
            return value
        return self.project.get_input(name)

    # ---------------------------------

    def get_file(self, target_file: File) -> bool:
        """Check if a file exists. If not, try to download it from the database.
        If the file is not found in the database it is fine, we do not even warn the user.
        Note that this function is used to get populations and transitions files, which are not common.
        """
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
        self.remote.download_file(target_file.filename, target_file)
        return True

    def print_tests_summary(self):
        """Make a summary of tests and their status."""
        print('Tests summary:')
        for test_name in AVAILABLE_CHECKINGS:
            test_result = self.register.tests.get(test_name, None)
            # Print things pretty
            test_nice_name = NICE_NAMES[test_name]
            test_nice_result = None
            if test_result is None:
                test_nice_result = YELLOW_HEADER + 'Not run' + COLOR_END
            elif test_result is False:
                test_nice_result = RED_HEADER + 'Failed' + COLOR_END
            elif test_result is True:
                test_nice_result = GREEN_HEADER + 'Passed' + COLOR_END
            elif test_result == 'na':
                test_nice_result = BLUE_HEADER + 'Not applicable' + COLOR_END
            else:
                raise ValueError()

            print(f' - {test_nice_name} -> {test_nice_result}')

    # Issue some warnings if failed or never run tests are skipped
    # This is run after processing input files
    # RUBEN: mover a inpro?
    def _issue_required_test_warnings(self):
        for test_name in AVAILABLE_CHECKINGS:
            # If test was not skipped then proceed
            if test_name not in self.project.trust: continue
            # If test passed in a previous run the proceed
            test_result = self.register.tests.get(test_name)
            if test_result is True: continue
            # If test failed in a previous run we can also proceed
            # The failing warning must be among the inherited warnings, so there is no need to add more warnings here
            if test_result is False: continue
            # If the test has been always skipped then issue a warning
            if test_result is None:
                # Remove previous warnings
                self.register.remove_warnings(test_name)
                # Get test pretty name
                test_nice_name = NICE_NAMES[test_name]
                # Issue the corresponding warning
                self.register.add_warning(test_name, test_nice_name + ' was skipped and never run before')
                continue
            raise ValueError('Test value is not supported')

    # Processed files ----------------------------------------------------

    def get_topology_filepath(self) -> str:
        """Get the processed topology file path."""
        return self.project.get_topology_filepath()
    topology_filepath = property(get_topology_filepath, None, None, "Topology file path (read only)")

    # Run the actual processing to generate output processed files out of input raw files
    # And by "files" I mean structure, trajectory and topology
    input_files_processing = Task('inpro', 'Input files processing', process_input_files,
        output_filenames={
            'output_topology_filepath': get_topology_filepath,
            'output_structure_filepath': STRUCTURE_FILENAME,
            'output_trajectory_filepath': TRAJECTORY_FILENAME,
        })

    # Output main files
    get_structure_file = input_files_processing.get_output_file_getter('output_structure_filepath')
    structure_file = property(get_structure_file, None, None, "Structure file (read only)")
    get_trajectory_file = input_files_processing.get_output_file_getter('output_trajectory_filepath')
    trajectory_file = property(get_trajectory_file, None, None, "Trajectory file (read only)")

    def get_topology_file(self) -> 'File':
        """Get the processed topology from the project."""
        return self.project.get_topology_file()
    topology_file = property(get_topology_file, None, None,
                             "Topology filename from the project (read only)")

    # ---------------------------------------------------------------------------------
    # Others values which may be found/calculated and files to be generated on demand
    # ---------------------------------------------------------------------------------

    # Trajectory snapshots
    get_snapshots = Task('frames', 'Count trajectory frames', get_frames_count)
    snapshots = property(get_snapshots, None, None, "Trajectory snapshots (read only)")

    def get_reference_bonds(self) -> list[list[int]]:
        """Get the reference bonds."""
        return self.project.reference_bonds
    reference_bonds = property(get_reference_bonds, None, None, "Atom bonds to be trusted (read only)")

    def get_structure(self) -> 'Structure':
        """Get the parsed structure."""
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
        get_first_frame, output_filenames={'output_filepath': FIRST_FRAME_FILENAME})
    get_first_frame_file = get_first_frame.get_output_file_getter('output_filepath')
    first_frame_file = property(get_first_frame_file, None, None, "First frame (read only)")

    # Average structure filename
    get_average_structure = Task('average', 'Get average structure',
        get_average_structure, output_filenames={'output_filepath': AVERAGE_STRUCTURE_FILENAME})
    get_average_structure_file = get_average_structure.get_output_file_getter('output_filepath')
    average_structure_file = property(get_average_structure_file, None, None, "Average structure filename (read only)")

    # Produce the MD metadata file to be uploaded to the database
    prepare_metadata = Task('mdmeta', 'Prepare MD metadata',
        generate_md_metadata, output_filenames={'output_filepath': OUTPUT_METADATA_FILENAME})

    # The processed interactions
    get_processed_interactions = Task('inter', 'Interactions processing',
        process_interactions, {'frames_limit': 1000})
    interactions = property(get_processed_interactions, None, None, "Processed interactions (read only)")

    # MDAnalysis Universe object
    get_MDAnalysis_Universe = Task('mda_univ', 'MDAnalysis Universe object',
        get_mda_universe, use_cache=False)
    universe = property(get_MDAnalysis_Universe, None, None, "MDAnalysis Universe object (read only)")

    def input_getter(name: str):
        """Get input values which may be MD specific.
        If the MD input is missing then we use the project input value.
        """
        # Set the getter
        def getter(self):
            # Get the MD input
            value = self.md_inputs.get(name, None)
            if value is not None:
                return value
            # If there is no MD input then return the project value
            return getattr(self.project, f'input_{name}')
        return getter

    # Assign the MD input getters
    input_interactions = property(input_getter('interactions'), None, None, "Interactions to be analyzed (read only)")
    input_pbc_selection = property(input_getter('pbc_selection'), None, None, "Selection of atoms which are still in periodic boundary conditions (read only)")
    input_cg_selection = property(input_getter('cg_selection'), None, None, "Selection of atoms which are not actual atoms but coarse grain beads (read only)")

    def _set_pbc_selection(self, reference_structure: 'Structure', verbose: bool = False) -> 'Selection':
        """Set PBC selection.
        It may parse the inputs file selection string if it is available or guess it otherwise.
        """
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
            selected_residue_names = list(set([residue.name for residue in selected_residues]))
            limit = 3  # Show a maximum of 3 residue names
            example_residue_names = ', '.join(selected_residue_names[0:limit])
            if len(selected_residue_names) > limit: example_residue_names += ', etc.'
            print('  e.g. ' + example_residue_names)
        return parsed_selection

    def get_pbc_selection(self) -> 'Selection':
        """Get the periodic boundary conditions atom selection."""
        # If we already have a stored value then return it
        if self.project._pbc_selection is not None:
            return self.project._pbc_selection
        # Otherwise we must set the PBC selection
        self.project._pbc_selection = self._set_pbc_selection(self.structure)
        return self.project._pbc_selection
    pbc_selection = property(get_pbc_selection, None, None, "Periodic boundary conditions atom selection (read only)")

    # WARNING: Do not inherit project pbc residues
    # WARNING: It may trigger all the processing logic of the reference MD when there is no need
    def get_pbc_residues(self) -> list[int]:
        """Get indices of residues in periodic boundary conditions."""
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

    # DANI: Esto algún día habría que tratar de automatizarlo
    def _set_cg_selection(self, reference_structure: 'Structure', verbose: bool = False) -> 'Selection':
        """Set the coarse grain selection."""
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
            selected_residue_names = list(set([residue.name for residue in selected_residues]))
            limit = 3  # Show a maximum of 3 residue names
            example_residue_names = ', '.join(selected_residue_names[0:limit])
            if len(selected_residue_names) > limit: example_residue_names += ', etc.'
            print('  e.g. ' + example_residue_names)
        return parsed_selection

    def get_cg_selection(self) -> 'Selection':
        """Get the coarse grain atom selection."""
        # If we already have a stored value then return it
        if self.project._cg_selection:
            return self.project._cg_selection
        # Otherwise we must set the PBC selection
        self.project._cg_selection = self._set_cg_selection(self.structure)
        return self.project._cg_selection
    cg_selection = property(get_cg_selection, None, None, "Periodic boundary conditions atom selection (read only)")

    # WARNING: Do not inherit project cg residues
    # WARNING: It may trigger all the processing logic of the reference MD when there is no need
    def get_cg_residues(self) -> list[int]:
        """Get indices of residues in coarse grain."""
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

    def get_populations(self) -> list[float]:
        """Get equilibrium populations from a MSM from the project."""
        return self.project.populations
    populations = property(get_populations, None, None, "Equilibrium populations from a MSM (read only)")

    def get_transitions(self) -> list[list[float]]:
        """Get transition probabilities from a MSM from the project."""
        return self.project.transitions
    transitions = property(get_transitions, None, None, "Transition probabilities from a MSM (read only)")

    def get_protein_map(self) -> dict:
        """Get the residues mapping from the project."""
        return self.project.protein_map
    protein_map = property(get_protein_map, None, None, "Residues mapping (read only)")

    def get_charges(self) -> dict:
        """Get the residues mapping from the project."""
        return self.project.charges
    charges = property(get_charges, None, None, "Residues charges (read only)")

    # Reference frame
    get_reference_frame = Task('reframe', 'Reference frame', get_bonds_reference_frame)
    reference_frame = property(get_reference_frame, None, None, "Reference frame to be used to represent the MD (read only)")

    # ---------------------------------------------------------------------------------
    # Tests
    # ---------------------------------------------------------------------------------

    def is_trajectory_integral(self) -> Optional[bool]:
        """Sudden jumps test."""
        # If we already have a stored value then return it
        if self._trajectory_integrity is not None:
            return self._trajectory_integrity
        # Otherwise we must find the value
        self._trajectory_integrity = check_trajectory_integrity(
            input_structure_filename=self.structure_file.path,
            input_trajectory_filename=self.trajectory_file.path,
            structure=self.structure,
            pbc_selection=self.pbc_selection,
            mercy=self.project.mercy,
            trust=self.project.trust,
            register=self.register,
            # time_length = self.time_length,
            check_selection=ALL_ATOMS,
            standard_deviations_cutoff=self.project.rmsd_cutoff,
            snapshots=self.snapshots,
        )
        return self._trajectory_integrity

    # ---------------------------------------------------------------------------------
    # Analyses
    # ---------------------------------------------------------------------------------

    # RMSDs analysis
    run_rmsds_analysis = Task('rmsds', 'RMSDs analysis',
        rmsds, {'frames_limit': 5000})

    # TM scores analysis
    run_tmscores_analysis = Task('tmscore', 'TM scores analysis',
        tmscores, {'frames_limit': 200})

    # RMSF, atom fluctuation analysis
    run_rmsf_analysis = Task('rmsf', 'Fluctuation analysis', rmsf)

    # Radius of gyration analysis
    run_rgyr_analysis = Task('rgyr', 'Radius of gyration analysis',
        rgyr, {'frames_limit': 5000})

    # PCA, principal component analysis
    run_pca_analysis = Task('pca', 'Principal component analysis',
        pca, {'frames_limit': 2000, 'projection_frames': 20})

    # PCA contacts
    # DANI: Intenta usar mucha memoria, hay que revisar
    # DANI: Puede saltar un error de imposible alojar tanta memoria
    # DANI: Puede comerse toda la ram y que al final salte un error de 'Terminado (killed)'
    # DANI: De momento me lo salto
    # DANI: Lleva mucho tiempo sin mantenerse, habrá que cambiar varias cosas al recuperarlo
    # run_pca_contacts('pcacons', 'PCA contacts', pca_contacts)

    # RMSD per residue analysis
    run_rmsd_perres_analysis = Task('perres', 'RMSD per residue analysis',
        rmsd_per_residue, {'frames_limit': 100})

    # RMSD pairwise
    # Perform an analysis for the overall structure and then one more analysis for each interaction
    run_rmsd_pairwise_analysis = Task('pairwise', 'RMSD pairwise',
        rmsd_pairwise, {'frames_limit': 200, 'overall_selection': "name CA or name C5'"})

    # Run the cluster analysis
    run_clusters_analysis = Task('clusters', 'Clusters analysis',
        clusters_analysis, {'frames_limit': 1000, 'desired_n_clusters': 20})

    # Calculate the distance mean and standard deviation of each pair of residues
    run_dist_perres_analysis = Task('dist', 'Distance per residue',
        distance_per_residue, {'frames_limit': 200})

    # Hydrogen bonds
    run_hbonds_analysis = Task('hbonds', 'Hydrogen bonds analysis',
        hydrogen_bonds, {'time_splits': 100})

    # SASA, solvent accessible surface analysis
    run_sas_analysis = Task('sas', 'Solvent accessible surface analysis',
        sasa, {'frames_limit': 100})

    # Perform the electrostatic and vdw energies analysis for each pair of interaction agents
    run_energies_analysis = Task('energies', 'Energies analysis',
        energies, {'frames_limit': 100})

    # Calculate torsions and then dihedral energies for every dihedral along the trajectory
    run_dihedral_energies = Task('dihedrals', 'Dihedral energies analysis',
        compute_dihedral_energies, {'frames_limit': 100})

    # Perform the pockets analysis
    run_pockets_analysis = Task('pockets', 'Pockets analysis',
        pockets, {'frames_limit': 100, 'maximum_pockets_number': 10})

    # Helical parameters
    run_helical_analysis = Task('helical', 'Helical parameters', helical_parameters)

    # Markov
    run_markov_analysis = Task('markov', 'Markov', markov, {'rmsd_selection': PROTEIN_AND_NUCLEIC})

    # Membrane density analysis
    run_density_analysis = Task('density', 'Membrane density analysis',
        density, {'frames_limit': 1000})

    # Membrane thickness analysis
    run_thickness_analysis = Task('thickness', 'Membrane thickness analysis',
        thickness, {'frames_limit': 100})

    # Area per lipid analysis
    run_apl_analysis = Task('apl', 'Membrane area per lipid analysis', area_per_lipid)

    # Calculate lipid order parameters for membranes
    run_lipid_order_analysis = Task('lorder', 'Membrane lipid order analysis',
        lipid_order, {'frames_limit': 100})

    # Lipid-protein interactions analysis
    run_lipid_interactions_analysis = Task('linter', 'Membrane lipid-protein interactions analysis',
        lipid_interactions, {'frames_limit': 100})

    run_channels_analysis = Task('channels', 'Membrane channels analysis',
        channels, {'frames_limit': 10})

    def get_warnings(self) -> list:
        """Get the warnings.

        The warnings list should not be reasigned, but it was back in the day.
        To avoid silent bugs, we read it directly from the register every time.
        """
        return self.register.warnings
    warnings = property(get_warnings, None, None, "MD warnings to be written in metadata")


class Project:
    """Class for the main project of an MDDB accession.
    A project is a set of related MDs.
    These MDs share all or most topology and metadata.
    """

    def __init__(self,
        directory: str = '.',
        accession: Optional[str] = None,
        database_url: str = DEFAULT_API_URL,
        inputs_filepath: str = None,
        input_topology_filepath: Optional[str] = None,
        input_structure_filepath: Optional[str] = None,
        input_trajectory_filepaths: Optional[list[str]] = None,
        md_directories: Optional[list[str]] = None,
        md_config: Optional[list[list[str]]] = None,
        reference_md_index: Optional[int] = None,
        forced_inputs: Optional[list[list[str]]] = None,
        populations_filepath: str = DEFAULT_POPULATIONS_FILENAME,
        transitions_filepath: str = DEFAULT_TRANSITIONS_FILENAME,
        aiida_data_filepath: Optional[str] = None,
        filter_selection: bool | str = False,
        pbc_selection: Optional[str] = None,
        cg_selection: Optional[str] = None,
        image: bool = False,
        fit: bool = False,
        translation: list[float] = [0, 0, 0],
        mercy: list[str] | bool = [],
        trust: list[str] | bool = [],
        faith: bool = False,
        pca_analysis_selection: str = PROTEIN_AND_NUCLEIC_BACKBONE,
        pca_fit_selection: str = PROTEIN_AND_NUCLEIC_BACKBONE,
        rmsd_cutoff: float = DEFAULT_RMSD_CUTOFF,
        interaction_cutoff: float = DEFAULT_INTERACTION_CUTOFF,
        interactions_auto: Optional[str] = None,
        guess_bonds: bool = False,
        ignore_bonds: bool = False,
        sample_trajectory: Optional[int] = None,
    ):
        """Initialize a Project.

        Args:
            directory (str):
                Local directory where the project takes place.
            accession (Optional[str]):
                Project accession to download missing input files from the database (if already uploaded).
            database_url (str):
                API URL to download missing data. when an accession is provided.
            inputs_filepath (str):
                Path to a file with inputs for metadata, simulation parameters and analysis config.
            input_topology_filepath (Optional[str]):
                Path to input topology file relative to the project directory.
                Multiple formats accepted; default is our parsed JSON topology.
            input_structure_filepath (Optional[str]):
                Path to input structure file. It may be relative to the project or to each MD directory.
                If this value is not passed then the standard structure file is used as input by default.
            input_trajectory_filepaths (Optional[list[str]]):
                Paths to input trajectory files relative to each MD directory.
                If this value is not passed then the standard trajectory file path is used as input by default.
            md_directories (Optional[list[str]]):
                Path to the different MD directories.
                Each directory is to contain an independent trajectory and structure.
                Several output files will be generated in every MD directory.
            md_config (Optional[list]):
                Configuration of a specific MD. You may declare as many as you want.
                Every MD requires a directory name and at least one trajectory path.
                The structure is -md <directory> <trajectory_1> <trajectory_2> ...
                Note that all trajectories from the same MD will be merged.
                For legacy reasons, you may also provide a specific structure for an MD.
                e.g. -md <directory> <structure> <trajectory_1> <trajectory_2> ...
            reference_md_index (Optional[int]):
                Index of the reference MD (used by project-level functions; defaults to first MD).
            forced_inputs (Optional[list]):
                Force a specific input through the command line.
                Inputs passed through command line have priority over the ones from the inputs file.
                In fact, these values will overwritten or be appended in the inputs file.
                Every forced input requires an input name and a value.
                The structure is -fi <input name> <new input value>
            populations_filepath (str):
                Path to equilibrium populations file (Markov State Model only)
            transitions_filepath (str):
                Path to transition probabilities file (Markov State Model only).
            aiida_data_filepath (Optional[str]):
                Path to the AiiDA data file.
                This file may be generated by the aiida-gromacs plugin and contains provenance data.
            filter_selection (bool|str):
                Atoms selection to be filtered in VMD format.
                If the argument is passed alone (i.e. with no selection) then water and counter ions are filtered.
            pbc_selection (Optional[str]):
                Selection of atoms which stay in Periodic Boundary Conditions even after imaging the trajectory.
                e.g. remaining solvent, ions, membrane lipids, etc.
                Selection passed through console overrides the one in inputs file.
            cg_selection (Optional[str]):
                Selection of atoms which are not actual atoms but Coarse Grained beads.
                Selection passed through console overrides the one in inputs file.
            image (bool):
                Set if the trajectory is to be imaged so atoms stay in the PBC box. See -pbc for more information.
            fit (bool):
                Set if the trajectory is to be fitted (both rotation and translation) to minimize the RMSD to PROTEIN_AND_NUCLEIC_BACKBONE selection.
            translation (list[float]):
                Set the x y z translation for the imaging process.
                e.g. -trans 0.5 -1 0
            mercy (list[str]|bool):
                Failures to be tolerated (or boolean to set all/none).
            trust (list[str]|bool):
                Tests to skip/trust (or boolean to set all/none).
            faith (bool):
                If True, require input files to match expected output files and skip processing.
            pca_analysis_selection (str):
                Atom selection for PCA analysis in VMD syntax.
            pca_fit_selection (str):
                Atom selection for the PCA fitting in VMD syntax.
            rmsd_cutoff (float):
                Set the cutoff for the RMSD sudden jumps analysis to fail.
                This cutoff stands for the number of standard deviations away from the mean an RMSD value is to be.
            interaction_cutoff (float):
                Set the cutoff for the interactions analysis to fail.
                This cutoff stands for percent of the trajectory where the interaction happens (from 0 to 1).
            interactions_auto (Optional[str]):
                Guess input interactions automatically. A VMD selection may be passed to limit guessed interactions to a specific subset of atoms.
            guess_bonds (bool):
                Force the workflow to guess atom bonds based on distance and atom radii in different frames along the trajectory instead of mining topology bonds.
            ignore_bonds (bool):
                Force the workflow to ignore atom bonds. This will result in many check-ins being skipped
            sample_trajectory (Optional[int]):
                If passed, download the first 10 (by default) frames from the trajectory.
                You can specify a different number by providing an integer value.

        """
        # Save input parameters
        # RUBEN: directory nunca se usa, eliminar?
        self.directory = normpath(directory)
        # If it is an absolute path then make it relative to the project
        if isabs(self.directory):
            self.directory = relpath(self.directory)
        # Save the directory name alone apart
        if self.directory == '.':
            self.directory_name = basename(getcwd())
        else:
            self.directory_name = basename(self.directory)

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
        self.arg_input_topology_filepath = input_topology_filepath
        self._input_topology_filepath = None
        self._input_topology_file = None
        self._input_topology_url = None
        # Input structure and trajectory filepaths
        # Do not parse them to files yet, let this to the MD class
        self.input_structure_filepath = input_structure_filepath
        self.input_trajectory_filepaths = input_trajectory_filepaths

        # Make sure the new MD configuration (-md) was not passed as well as old MD inputs (-mdir, -traj)
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
            md_directories = [mdc[0] for mdc in self.md_config]
            if len(md_directories) > len(set(md_directories)):
                raise InputError('There are duplicated MD directories')

        # Input populations and transitions for MSM
        self.populations_filepath = populations_filepath
        self._populations_file = File(self.populations_filepath)
        self.transitions_filepath = transitions_filepath
        self._transitions_file = File(self.transitions_filepath)
        # Input AiiDA data
        self.aiida_data_filepath = aiida_data_filepath
        self._aiida_data_file = File(self.aiida_data_filepath) if aiida_data_filepath else None

        # Set the processed topology filepath, which depends on the input topology filename
        # Note that this file is different from the standard topology, although it may be standard as well
        self._topology_filepath = None

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
        if type(mercy) is bool:
            if mercy:
                self.mercy = AVAILABLE_FAILURES
            else:
                self.mercy = []
        self.trust = trust
        # Fix the trust input, if needed
        # If a boolean is passed instead of a list then we set its corresponding value
        if type(trust) is bool:
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
        self.ignore_bonds = ignore_bonds
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

        # Save forced inputs and use them when necessary
        self.forced_inputs = {}
        if forced_inputs:
            # Make sure the format is respected
            for forced_input in forced_inputs:
                n_values = len(forced_input)
                if n_values == 0:
                    raise InputError('There is an empty "-fin". Please remove it from the command line.')
                if n_values == 1:
                    only_value = forced_input[0]
                    raise InputError(f'There is a "-fin {only_value}" which is missing the new input value.')
                if n_values > 2:
                    suggested_fix = f'{forced_input[0]} "{" ".join(forced_input[1:])}"'
                    raise InputError(f'Too many values in "-fin {" ".join(forced_input)}".\n' +
                        ' Note that only two values are expected: -fin <input name> <new input value>\n' +
                        f' Did you forget the quotes maybe? Try this: -fin {suggested_fix}')

            # Save forced inputs as a dict
            self.forced_inputs = {name: value for name, value in forced_inputs}

            # Check that forced inputs exist
            # This is to prevent the user from loosing a lot of time for a silly typo
            template_inputs = load_yaml(INPUTS_TEMPLATE_FILEPATH)
            for input_name in self.forced_inputs.keys():
                if input_name in template_inputs: continue
                available_inputs = ', '.join(template_inputs.keys())
                raise InputError(f'Unrecognized forced input "{input_name}". Available inputs: {available_inputs}')

            # Overwrite input file values
            for input_name, input_value in self.forced_inputs.items():
                self.update_inputs(input_name, input_value)

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

    def __repr__(self):
        """Return a string representation of the Project object."""
        return 'Project'

    def pathify(self, filename_or_relative_path: str) -> str:
        """Given a filename or relative path, add the project directory path at the beginning."""
        return normpath(self.directory + '/' + filename_or_relative_path)

    def check_md_directories(self):
        """Check MD directories to be right.
        If there is any problem then directly raise an input error.
        """
        # Check there is at least one MD
        if len(self._md_directories) < 1:
            raise InputError('There must be at least one MD')
        # Check there are not duplicated MD directories
        if len(set(self._md_directories)) != len(self._md_directories):
            raise InputError('There are duplicated MD directories')

    def get_md_directories(self) -> list:
        """Get MD directories."""
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

    def get_reference_md_index(self) -> int:
        """Get the reference MD index."""
        # If we are already have a value then return it
        if self._reference_md_index:
            return self._reference_md_index
        # Otherwise we must find the reference MD index
        # If the inputs file is available then it must declare the reference MD index
        if self.is_inputs_file_available():
            self._reference_md_index = self.get_input('mdref')
        # Otherwise we simply set the first MD as the reference and warn the user about this
        if self._reference_md_index is None:
            warn('No reference MD was specified. The first MD will be used as reference.')
            self._reference_md_index = 0
        return self._reference_md_index
    reference_md_index = property(get_reference_md_index, None, None, "Reference MD index (read only)")

    def get_reference_md(self) -> MD:
        """Get the reference MD."""
        # If we are already have a value then return it
        if self._reference_md:
            return self._reference_md
        # Otherwise we must find the reference MD
        self._reference_md = self.mds[self.reference_md_index]
        return self._reference_md
    reference_md: MD = property(get_reference_md, None, None, "Reference MD (read only)")

    def get_mds(self) -> list:
        """Get the available MDs (read only)."""
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
                    project=self, number=n, directory=directory,
                    input_structure_filepath=input_structure_filepath,
                    input_trajectory_filepaths=input_trajectory_filepaths,
                )
                self._mds.append(md)
        # Old system (-mdir, -stru -traj)
        else:
            for n, md_directory in enumerate(self.md_directories, 1):
                md = MD(
                    project=self, number=n, directory=md_directory,
                    input_structure_filepath=self.input_structure_filepath,
                    input_trajectory_filepaths=self.input_trajectory_filepaths,
                )
                self._mds.append(md)
        return self._mds
    mds: list[MD] = property(get_mds, None, None, "Available MDs (read only)")

    # Check input files exist when their filenames are read
    # If they do not exist then try to download them
    # If the download is not possible then raise an error

    # Inputs filename ------------

    def is_inputs_file_available(self) -> bool:
        """Set a function to check if inputs file is available.
        Note that asking for it when it is not available will lead to raising an input error.
        """
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

    def get_inputs_file(self) -> File:
        """Set a function to load the inputs yaml file."""
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

    def get_input_topology_filepath(self) -> Optional[str]:
        """Get the input topology filepath from the inputs or try to guess it.

        If the input topology filepath is a 'no' flag then we consider there is no topology at all
        So far we extract atom charges and atom bonds from the topology file
        In this scenario we can keep working but there are some consecuences:
          1. Analysis using atom charges such as 'energies' will be skipped
          2. The standard topology file will not include atom charges
          3. Bonds will be guessed
        """
        # If we already have an internal value calculated then return it
        if self._input_topology_filepath is not None:
            return self._input_topology_filepath

        def parse(filepath: str) -> str:
            """Parse possible glob notation."""
            # If it is a URL then set the paths for the further download
            if is_url(filepath):
                self._input_topology_url = filepath
                source_filename = url_to_source_filename(filepath)
                return source_filename
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
        if self.arg_input_topology_filepath:
            if self.arg_input_topology_filepath.lower() in {'no', 'not', 'na'}:
                self._input_topology_filepath = MISSING_TOPOLOGY
                return self._input_topology_filepath
            self._input_topology_filepath = parse(self.arg_input_topology_filepath)
            # Update the input topology fielpath in the inputs file, in case it is not matching
            self.update_inputs('input_topology_filepath', relpath(self._input_topology_filepath, self.directory))
            return self._input_topology_filepath
        # Check if the inputs file has the value
        if self.is_inputs_file_available():
            # Get the input value, whose key must exist
            inputs_value = self.get_input('input_topology_filepath')
            # If there is a valid input then use it
            if inputs_value is not None:
                # WARNING: the yaml parser automatically converts 'no' to False
                if inputs_value is False or inputs_value.lower() in {'no', 'not', 'na'}:
                    self._input_topology_filepath = MISSING_TOPOLOGY
                    return self._input_topology_filepath
                parsed_input_value = parse(inputs_value)
                self._input_topology_filepath = self.pathify(parsed_input_value)
                return self._input_topology_filepath
        # If nothing worked then surrender
        raise InputError('Missing input topology file path. Please provide a topology file using the "-top" argument.\n' +
            '  Note that you may run the workflow without a topology file. To do so, use the "-top no" argument.\n' +
            '  However this has implications since we usually mine atom charges and bonds from the topology file.\n' +
            '  Some analyses such us the interaction energies will be skiped')

    def get_input_topology_file(self) -> Optional[File]:
        """Get the input topology file.
        If the file is not found try to download it.
        """
        # If we already have a value then return it
        if self._input_topology_file is not None:
            return self._input_topology_file
        # Set the input topology filepath
        input_topology_filepath = self.get_input_topology_filepath()
        # If the input filepath is None then it menas we must proceed without a topology
        if input_topology_filepath == MISSING_TOPOLOGY:
            self._input_topology_file = MISSING_TOPOLOGY
            return self._input_topology_file
        # If no input is passed then we check the inputs file
        # Set the file
        input_topology_file = File(input_topology_filepath)
        # If there is a topology URL then we must donwload the topology first
        input_topology_url = self._input_topology_url
        if input_topology_url and not input_topology_file.exists:
            original_filename = input_topology_url.split('/')[-1]
            # If there is a remote then use it
            if self.remote:
                # If the topology filename is the standard topology filename then use the topology endpoint instead
                if original_filename == STANDARD_TOPOLOGY_FILENAME:
                    self.remote.download_standard_topology(input_topology_file)
                # Otherwise download the input strucutre file by its filename
                else:
                    self.remote.download_file(original_filename, input_topology_file)
                # In case the topology is a '.top' file we consider it is a Gromacs topology
                # It may come with additional itp files we must download as well
                if input_topology_file.format == 'top':
                    # Find available .itp files and download each of them
                    itp_filenames = [filename for filename in self.remote.available_files if filename[-4:] == '.itp']
                    for itp_filename in itp_filenames:
                        itp_filepath = self.pathify(itp_filename)
                        itp_file = File(itp_filepath)
                        self.remote.download_file(itp_file.filename, itp_file)
            # Otherwise use the URL as is
            else:
                if input_topology_file.format == 'top':
                    warn('Automatic download of itp files is not supported without the "-proj" argument.' +
                         ' Thus if the topology has associated itp files they will not be downloaded.')
                download_file(input_topology_url, input_topology_file)
        # If the file already exists then we are done
        if input_topology_file.exists:
            self._input_topology_file = input_topology_file
            return self._input_topology_file
        raise InputError(f'Missing input topology file "{input_topology_file.filename}"')
    input_topology_file = property(get_input_topology_file, None, None, "Input topology file (read only)")

    def get_input_structure_file(self) -> File:
        """Get the input structure filename."""
        # When calling this function make sure all MDs have the file or try to download it
        return self.reference_md._input_structure_file
    input_structure_file = property(get_input_structure_file, None, None, "Input structure filename for each MD (read only)")

    def get_input_trajectory_files(self) -> list[File]:
        """Get the input trajectory filename(s) from the inputs.
        If file(s) are not found try to download it.
        """
        return self.reference_md._input_trajectory_files
    input_trajectory_files = property(get_input_trajectory_files, None, None, "Input trajectory filenames for each MD (read only)")

    def get_populations_file(self) -> Optional[File]:
        """Get the MSM equilibrium populations file."""
        if not self.get_file(self._populations_file):
            return None
        return self._populations_file
    populations_file = property(get_populations_file, None, None, "MSM equilibrium populations file (read only)")

    def get_transitions_file(self) -> Optional[File]:
        """Get the MSM transition probabilities file."""
        if not self.get_file(self._transitions_file):
            return None
        return self._transitions_file
    transitions_file = property(get_transitions_file, None, None, "MSM transition probabilities file (read only)")

    def get_aiida_data_file(self) -> Optional[File]:
        """Get the AiiDA data file."""
        if not self._aiida_data_file: return None
        if not self.get_file(self._aiida_data_file): return None
        return self._aiida_data_file
    aiida_data_file = property(get_aiida_data_file, None, None, "AiiDA data file (read only)")

    # ---------------------------------

    def get_file(self, target_file: File) -> bool:
        """Check if a file exists.
        If not, try to download it from the database.
        If the file is not found in the database it is fine, we do not even warn the user.
        Note that nowadays this function is used to get populations and transitions files, which are not common.
        """
        return self.reference_md.get_file(target_file)

    # Input file values -----------------------------------------

    # First of all set input themselves

    def get_inputs(self) -> dict:
        """Get inputs."""
        # If inputs are already loaded then return them
        if self._inputs:
            return self._inputs
        # When loading the inuts file, replace some values automatically
        replaces = [
            ('$DIR', self.directory_name)
        ]
        # Otherwise, load inputs from the inputs file
        inputs_data = None
        if self.inputs_file.format == 'json':
            inputs_data = load_json(self.inputs_file.path, replaces)
        elif self.inputs_file.format == 'yaml':
            inputs_data = load_yaml(self.inputs_file.path, replaces)
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

    def update_inputs(self, nested_key: str, new_value):
        """Permanently update the inputs file.
        This may be done when command line inputs do not match file inputs.
        """
        # If there is no inputs file then rhere is nothing to update
        if not self.is_inputs_file_available(): return
        # If the input already matches then do nothing
        current_value = read_ndict(self.inputs, nested_key, MISSING_INPUT_EXCEPTION)
        if current_value == new_value: return
        # Set the new value
        write_ndict(self.inputs, nested_key, new_value)
        # If there is no inputs file then do not try to save anything
        if not self.is_inputs_file_available(): return
        print(f'* Field "{nested_key}" in the inputs file will be permanently modified')
        # Write the new inputs to disk
        if self.inputs_file.format == 'json':
            save_json(self.inputs, self.inputs_file.path)
        elif self.inputs_file.format == 'yaml':
            # Note that comments in the original YAML file will be not kept
            save_yaml(self.inputs, self.inputs_file.path)
        else:
            raise InputError('Input file format is not supported. Please use json or yaml files.')

    # Then set getters for every value in the inputs file
    def get_input(self, name: str) -> any:
        """Get a specific 'input' value."""
        # Check if the value of this input was forced from command line
        if name in self.forced_inputs:
            return self.forced_inputs[name]
        # Get the input value from the inputs file
        value = self.inputs.get(name, MISSING_INPUT_EXCEPTION)
        # If we had a value then return it
        if value != MISSING_INPUT_EXCEPTION:
            return value
        # If the field is not specified in the inputs file then set a defualt value
        default_value = DEFAULT_INPUT_VALUES.get(name, None)
        # Warn the user about this
        warn(f'Missing input "{name}" -> Using default value: {default_value}')
        return default_value

    def input_property(name: str, doc: str = ""):
        """Set a function to get a specific 'input' value by its key/name.
        Note that we return the property without calling the getter.
        """
        def getter(self: 'Project'):
            return self.get_input(name)
        return property(getter, doc=doc)

    # Assign the getters
    input_interactions = input_property('interactions', "Interactions to be analyzed (read only)")
    input_protein_references = input_property('forced_references', "Uniprot IDs to be used first when aligning protein sequences (read only)")
    input_pdb_ids = input_property('pdb_ids', "Protein Data Bank IDs used for the setup of the system (read only)")
    input_type = input_property('type', "Set if its a trajectory or an ensemble (read only)")
    input_mds = input_property('mds', "Input MDs configuration (read only)")
    input_ligands = input_property('ligands', "Input ligand references (read only)")
    input_force_fields = input_property('ff', "Input force fields (read only)")
    input_collections = input_property('collections', "Input collections (read only)")
    input_chain_names = input_property('chainnames', "Input chain names (read only)")
    input_framestep = input_property('framestep', "Input framestep (read only)")
    input_name = input_property('name', "Input name (read only)")
    input_description = input_property('description', "Input description (read only)")
    input_authors = input_property('authors', "Input authors (read only)")
    input_groups = input_property('groups', "Input groups (read only)")
    input_contact = input_property('contact', "Input contact (read only)")
    input_program = input_property('program', "Input program (read only)")
    input_version = input_property('version', "Input version (read only)")
    input_method = input_property('method', "Input method (read only)")
    input_license = input_property('license', "Input license (read only)")
    input_linkcense = input_property('linkcense', "Input license link (read only)")
    input_citation = input_property('citation', "Input citation (read only)")
    input_thanks = input_property('thanks', "Input acknowledgements (read only)")
    input_links = input_property('links', "Input links (read only)")
    input_timestep = input_property('timestep', "Input timestep (read only)")
    input_temperature = input_property('temp', "Input temperature (read only)")
    input_ensemble = input_property('ensemble', "Input ensemble (read only)")
    input_water = input_property('wat', "Input water force field (read only)")
    input_boxtype = input_property('boxtype', "Input boxtype (read only)")
    input_pbc_selection = input_property('pbc_selection', "Input Periodic Boundary Conditions (PBC) selection (read only)")
    input_cg_selection = input_property('cg_selection', "Input Coarse Grained (CG) selection (read only)")
    input_customs = input_property('customs', "Input custom representations (read only)")
    input_orientation = input_property('orientation', "Input orientation (read only)")
    input_multimeric = input_property('multimeric', "Input multimeric labels (read only)")
    # Additional topic-specific inputs
    input_cv19_unit = input_property('cv19_unit', "Input Covid-19 Unit (read only)")
    input_cv19_startconf = input_property('cv19_startconf', "Input Covid-19 starting conformation (read only)")
    input_cv19_abs = input_property('cv19_abs', "Input Covid-19 antibodies (read only)")
    input_cv19_nanobs = input_property('cv19_nanobs', "Input Covid-19 nanobodies (read only)")

    def get_input_pbc_selection(self) -> Optional[str]:
        """PBC selection may come from the console or from the inputs file.
        Console has priority over the inputs file.
        """
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

    def get_input_cg_selection(self) -> Optional[str]:
        """CG selection may come from the console or from the inputs file.
        Console has priority over the inputs file.
        """
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

    def check_is_time_dependent(self) -> bool:
        """Set if MDs are time dependent."""
        if self.input_type == 'trajectory':
            return True
        elif self.input_type == 'ensemble':
            return False
        raise InputError(f'Not supported input "type" value: {self.input_type}. It must be "trajectory" or "ensemble"')
    is_time_dependent = property(check_is_time_dependent, None, None, "Check if trajectory frames are time dependent (read only)")

    # Processed files ----------------------------------------------------

    def inherit_topology_filename(self) -> Optional[str]:
        """Set the expected output topology filename given the input topology filename.
        Note that topology formats are conserved.
        """
        if self.input_topology_file == MISSING_TOPOLOGY:
            return None
        filename = self.input_topology_file.filename
        if not filename:
            return None
        if filename == RAW_CHARGES_FILENAME:
            return filename
        standard_format = self.input_topology_file.format
        return 'topology.' + standard_format

    def get_topology_filepath(self) -> str:
        """Get the processed topology file path."""
        # If we have a stored value then return it
        if self._topology_filepath:
            return self._topology_filepath
        # Otherwise we must find it
        inherited_filename = self.inherit_topology_filename()
        self._topology_filepath = self.pathify(inherited_filename) if inherited_filename else None
        return self._topology_filepath
    topology_filepath = property(get_topology_filepath, None, None, "Topology file path (read only)")

    def get_topology_file(self) -> 'File':
        """Get the processed topology from the reference MD."""
        getter = self.reference_md.input_files_processing.get_output_file_getter('output_topology_filepath')
        return getter(self.reference_md)
    topology_file = property(get_topology_file, None, None, "Topology file (read only)")

    def get_structure_file(self) -> 'File':
        """Get the processed structure from the reference MD."""
        return self.reference_md.structure_file
    structure_file = property(get_structure_file, None, None, "Structure filename from the reference MD (read only)")

    def get_trajectory_file(self) -> 'File':
        """Get the processed trajectory from the reference MD."""
        return self.reference_md.trajectory_file
    trajectory_file = property(get_trajectory_file, None, None, "Trajectory filename from the reference MD (read only)")

    # ---------------------------------------------------------------------------------
    # Others values which may be found/calculated and files to be generated on demand
    # ---------------------------------------------------------------------------------

    def get_structure(self) -> 'Structure':
        """Get the parsed structure from the reference MD."""
        return self.reference_md.structure
    structure = property(get_structure, None, None, "Parsed structure from the reference MD (read only)")

    def get_pbc_residues(self) -> list[int]:
        """Get the indices of residues in periodic boundary conditions."""
        return self.reference_md.pbc_residues
    pbc_residues = property(get_pbc_residues, None, None, "Indices of residues in periodic boundary conditions (read only)")

    def get_cg_residues(self) -> list[int]:
        """Get the indices of residues in coarse grain."""
        return self.reference_md.cg_residues
    cg_residues = property(get_cg_residues, None, None, "Indices of residues in coarse grain (read only)")

    def get_snapshots(self) -> int:
        """Get the reference MD snapshots."""
        return self.reference_md.snapshots
    snapshots = property(get_snapshots, None, None, "Reference MD snapshots (read only)")

    def get_universe(self) -> int:
        """Get the MDAnalysis Universe from the reference MD."""
        return self.reference_md.universe
    universe = property(get_universe, None, None, "MDAnalysis Universe object (read only)")

    def get_processed_interactions(self) -> dict:
        """Get the processed interactions from the reference replica, which are the same for all replicas."""
        return self.reference_md.interactions
    interactions = property(get_processed_interactions, None, None, "Processed interactions (read only)")

    def get_check_stable_bonds(self) -> bool:
        """Check if we must check stable bonds."""
        # Set if stable bonds have to be checked
        must_check = STABLE_BONDS_FLAG not in self.trust
        # If this analysis has been already passed then we can trust structure bonds
        if self.register.tests.get(STABLE_BONDS_FLAG, None) is True:
            must_check = False
        return must_check
    must_check_stable_bonds = property(get_check_stable_bonds, None, None, "Check if we must check stable bonds (read only)")

    # Reference bonds
    get_reference_bonds = Task('refbonds', 'Reference bonds', find_safe_bonds)
    reference_bonds = property(get_reference_bonds, None, None, "Atom bonds to be trusted (read only)")

    # Atom charges
    get_charges = Task('charges', 'Getting atom charges', get_charges)
    charges = property(get_charges, None, None, "Atom charges (read only)")

    def get_topology_reader(self) -> 'Topology':
        """Get the topology data reader."""
        # If we already have a stored value then return it
        if self._topology_reader: return self._topology_reader
        # Instantiate the topology reader
        self._topology_reader = Topology(self.topology_file)
        return self._topology_reader
    topology_reader = property(get_topology_reader, None, None, "Topology reader (read only)")

    def get_dihedrals(self) -> list[dict]:
        """Get the topology dihedrals."""
        # If we already have a stored value then return it
        if self._dihedrals: return self._dihedrals
        # Calculate the dihedrals otherwise
        self._dihedrals = self.topology_reader.get_dihedrals_data()
        return self._dihedrals
    dihedrals = property(get_dihedrals, None, None, "Topology dihedrals (read only)")

    def get_populations(self) -> Optional[list[float]]:
        """Get the equilibrium populations from a MSM."""
        # If we already have a stored value then return it
        if self._populations:
            return self._populations
        # Otherwise we must find the value
        if not self.populations_file:
            return None
        self._populations = read_file(self.populations_file)
        return self._populations
    populations = property(get_populations, None, None, "Equilibrium populations from a MSM (read only)")

    def get_transitions(self) -> Optional[list[list[float]]]:
        """Get the transition probabilities from a MSM."""
        # If we already have a stored value then return it
        if self._transitions:
            return self._transitions
        # Otherwise we must find the value
        if not self.transitions_file:
            return None
        self._transitions = read_file(self.transitions_file)
        return self._transitions
    transitions = property(get_transitions, None, None, "Transition probabilities from a MSM (read only)")

    def get_pdb_ids(self) -> list[str]:
        """Get the tested and standardized PDB ids."""
        # If we already have a stored value then return it
        if self._pdb_ids is not None:
            return self._pdb_ids
        # Otherwise test and standarize input PDB ids
        self._pdb_ids = []
        # If there is no input pdb ids (may be None) then stop here
        if not self.input_pdb_ids:
            return []
        # If input PDB ids is a string instead of a list then fix it
        input_pdb_ids = self.input_pdb_ids
        if type(input_pdb_ids) is str:
            input_pdb_ids = [input_pdb_ids]
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
        prepare_pdb_references, output_filenames={'output_filepath': PDB_REFERENCES_FILENAME})
    pdb_references_file = get_pdb_references.get_output_file_getter('output_filepath')

    # Map the structure aminoacids sequences against the Uniprot reference sequences
    get_protein_map = Task('protmap', 'Protein residues mapping',
        generate_protein_mapping, output_filenames={'output_filepath': PROTEIN_REFERENCES_FILENAME})
    protein_map = property(get_protein_map, None, None, "Protein residues mapping (read only)")

    # Define the output file of the protein mapping including protein references
    get_protein_references_file = get_protein_map.get_output_file_getter('output_filepath')
    protein_references_file = property(get_protein_references_file, None, None, "File including protein refereces data mined from UniProt (read only)")

    get_chain_references = Task('chains', 'Chain references',
        prepare_chain_references, output_filenames={'output_filepath': OUTPUT_CHAINS_FILENAME})

    get_inchikeys = Task('inchikeys', 'InChI keys residues mapping', generate_inchikeys)
    inchikeys = property(get_inchikeys, None, None, "InChI keys (read only)")

    get_lipid_references = Task('lipmap', 'Lipid references', generate_lipid_references)
    lipid_references = property(get_lipid_references, None, None, "Lipid references (read only)")

    get_membrane_map = Task('memmap', 'Membrane residues mapping', generate_membrane_mapping,
        output_filenames={'output_filepath': MEMBRANE_MAPPING_FILENAME})
    membrane_map = property(get_membrane_map, None, None, "Membrane mapping (read only)")

    get_ligand_references = Task('ligmap', 'Ligand residues mapping', generate_ligand_references)
    ligand_references = property(get_ligand_references, None, None, "Ligand references (read only)")

    get_inchi_references = Task('inchimap', 'InChI references', generate_inchi_references,
        output_filenames={'output_filepath': INCHIKEY_REFERENCES_FILENAME})
    inchikey_map = property(get_inchi_references, None, None, "InChI references (read only)")

    # Build the residue map from both proteins and ligands maps
    # This is formatted as both the standard topology and metadata producers expect them
    get_residue_map = Task('resmap', 'Residue mapping', generate_residue_mapping)
    residue_map = property(get_residue_map, None, None, "Residue map (read only)")

    # Prepare the project metadata file to be upload to the database
    prepare_metadata = Task('pmeta', 'Prepare project metadata',
        prepare_project_metadata, output_filenames={'output_filepath': OUTPUT_METADATA_FILENAME})

    # Prepare the standard topology file to be uploaded to the database
    prepare_standard_topology = Task('stopology', 'Standard topology file',
        generate_topology, output_filenames={'output_filepath': STANDARD_TOPOLOGY_FILENAME})
    get_standard_topology_file = prepare_standard_topology.get_output_file_getter('output_filepath')
    standard_topology_file = property(get_standard_topology_file, None, None, "Standard topology filename (read only)")

    # Get a screenshot of the system
    get_screenshot_filename = Task('screenshot', 'Screenshot file',
        get_screenshot, output_filenames={'output_filepath': OUTPUT_SCREENSHOT_FILENAME})
    screenshot_filename = property(get_screenshot_filename, None, None, "Screenshot filename (read only)")

    # Provenance data
    produce_provenance = Task('aiidata', 'Produce provenance', produce_provenance)

    def get_warnings(self) -> list:
        """Get the warnings.

        The warnings list should not be reasigned, but it was back in the day.
        To avoid silent bugs, we read it directly from the register every time.
        """
        return self.register.warnings
    warnings = property(get_warnings, None, None, "Project warnings to be written in metadata")


# AUXILIAR FUNCTIONS ---------------------------------------------------------------------------

# DANI: En cuanto se concrete el formato de los markov esta función no hará falta
def read_file(target_file: File) -> dict:
    """Set a function to read a file which may be in different formats."""
    # Get the file format
    file_format = target_file.filename.split('.')[-1]
    # Read numpy files
    if file_format == 'npy':
        return numpy.load(target_file.path)
    # Read JSON files
    if file_format == 'json':
        return load_json(target_file.path)


def name_2_directory(name: str) -> str:
    """Set a function to convert an MD name into an equivalent MD directory."""
    # Replace white spaces with underscores
    directory = name.replace(' ', '_')
    # Remove problematic characters
    for character in FORBIDDEN_DIRECTORY_CHARACTERS:
        directory = directory.replace(character, '')
    return directory


def check_directory(directory: str) -> str:
    """Check for problematic characters in a directory path."""
    # Remove problematic characters
    directory_characters = set(directory)
    for character in FORBIDDEN_DIRECTORY_CHARACTERS:
        if character in directory_characters:
            raise InputError(f'Directory path "{directory}" includes the forbidden character "{character}"')


def directory_2_name(directory: str) -> str:
    """Convert an MD directory into an equivalent MD name."""
    # Remove a possible starting './'
    # Replace white spaces with underscores
    name = directory.split('/')[-1].replace('_', ' ')
    return name


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
input_files = {**project_input_files, **md_input_files}

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
processed_files = {**project_processed_files, **md_processed_files}

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
requestables = {**project_requestables, **md_requestables}

# Set groups of dependencies to be requested together using only one flag
DEPENDENCY_FLAGS = {
    'download': list(input_files.keys()),
    'setup': list(processed_files.keys()),
    'meta': ['pmeta', 'mdmeta'],
    'network': ['resmap', 'ligands', 'chains', 'pdbs', 'memmap'],
    'minimal': ['pmeta', 'mdmeta', 'stopology'],
    'interdeps': ['interactions', 'pairwise', 'hbonds', 'energies', 'perres', 'clusters', 'dist'],
    'membs': ['memmap', 'density', 'thickness', 'apl', 'lorder', 'linter', 'channels']
}

# Set the default analyses to be run when no task is specified
DEFAULT_ANALYSES = ['clusters', 'dist', 'energies', 'hbonds', 'pca', 'pockets',
    'rgyr', 'rmsds', 'perres', 'pairwise', 'rmsf', 'sas', 'tmscore', 'density',
    'thickness', 'apl', 'lorder', 'linter']


def workflow(
    project_parameters: dict = {},
    working_directory: str = '.',
    download: bool = False,
    setup: bool = False,
    include: Optional[list[str]] = None,
    exclude: Optional[list[str]] = None,
    overwrite: Optional[list[str] | bool] = None,
):
    """Run the MDDB workflow.

    Args:
        project_parameters (dict): Parameters to initiate the project.
        working_directory (str): Working directory where the project is located.
        download (bool): (Deprecated: use -i download) If passed, only download required files. Then exits.
        setup (bool): (Deprecated: use -i setup) If passed, only download required files and run mandatory dependencies. Then exits.
        include (list[str] | None): Set the unique analyses or tools to be run. All other steps will be skipped.
        exclude (list[str] | None): List of analyses or tools to be skipped. All other steps will be run. Note that if we exclude a dependency of something else then it will be run anyway.
        overwrite (list[str] | bool | None): List of tasks that will be re-run, overwriting previous output files. Use this flag alone to overwrite everything.

    """
    # Check there are not input errors

    # Include and exclude are not compatible
    # This is to protect the user to do something which makes not sense
    if include and exclude:
        raise InputError('Include (-i) and exclude (-e) are not compatible. Use one of these options.')

    # Save the directory where the workflow is called from so we can come back at the very end
    workflow_call_directory = getcwd()

    # Make sure the working directory exists
    if not exists(working_directory):
        raise InputError(f'Working directory "{working_directory}" does not exist')

    # Make sure the working directory is actually a directory
    if not isdir(working_directory):
        raise InputError(f'Working directory "{working_directory}" is actually not a directory')

    # Move the current directory to the working directory
    chdir(working_directory)
    current_directory_name = getcwd().split('/')[-1]
    git_version = get_git_version()
    print(f'\n{CYAN_HEADER}Running MDDB workflow ({git_version}) for project at {current_directory_name}{COLOR_END}')

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
        tasks = [*include]
        # Search for special flags among included
        for flag, dependencies in DEPENDENCY_FLAGS.items():
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
            excluded_dependencies = [*exclude]
            # Search for special flags among excluded
            for flag, dependencies in DEPENDENCY_FLAGS.items():
                if flag not in exclude: continue
                # If the flag is found then exclude the dependencies instead
                # Make sure not to duplicate a dependency if it was already included
                excluded_dependencies.remove(flag)
                for dep in dependencies:
                    if dep in exclude: continue
                    excluded_dependencies.append(dep)
            tasks = [name for name in tasks if name not in excluded_dependencies]

    # If the user requested to overwrite something, make sure it is in the tasks list

    # Update the overwritable variable with the requested overwrites
    overwritables = set()
    if overwrite:
        # If the overwrite argument is simply true then add all requestables to the overwritable
        if type(overwrite) is bool:
            for task in tasks:
                overwritables.add(task)
        # If the overwrite argument is a list of tasks then iterate them
        elif type(overwrite) is list:
            for task in overwrite:
                # Make sure the task to be overwriten is among the tasks to be run
                if task not in tasks:
                    raise InputError(f'Task "{task}" is to be overwriten but it is not in the tasks list. Either include it or do not exclude it')
                # Add it to the global variable
                overwritables.add(task)
        else: raise ValueError('Not supported overwrite type')

    # Get project tasks
    project_tasks = [task for task in tasks if task in project_requestables]
    # Get the MD tasks
    md_tasks = [task for task in tasks if task in md_requestables]

    # Set project overwritables
    project.overwritables = set([task for task in project_tasks if task in overwritables])
    # Set MD overwrittables
    # Note that this must be done before running project tasks
    # Some project tasks rely in MD tasks
    for md in project.mds:
        md.overwritables = set([task for task in md_tasks if task in overwritables])

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

    # Return back to the place where the workflow was called originally
    # This is not important for many applications
    # But if you call the workflow function from a python script then this is important
    chdir(workflow_call_directory)

    print("Done!")
