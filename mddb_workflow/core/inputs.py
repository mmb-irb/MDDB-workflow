"""Input document loading and value resolution.

This module contains the orchestration layer around the input schema.  The
schema itself remains in :mod:`inputs_schema`; this module is responsible for
loading a document and resolving values supplied by the caller.
"""

from collections.abc import Callable
from copy import deepcopy
from os.path import normpath
from typing import Any

from mddb_workflow.core.inputs_schema import validate_inputs
import mddb_workflow.utils.auxiliar as aux
from mddb_workflow.utils.auxiliar import InputError, warn

from mddb_workflow.utils.constants import (
    DEFAULT_INPUT_VALUE_WARNINGS,
    DEFAULT_INPUT_VALUES,
    FORBIDDEN_DIRECTORY_CHARACTERS,
    MD_DIRECTORY,
    MD_NAME,
    MD_REMOVED_FLAG,
    MD_INPUT_STRUCTURE_FILEPATH,
    MD_INPUT_TOPOLOGY_FILEPATH,
    MD_INPUT_TRAJECTORY_FILEPATHS,
    TOPOLOGY_SUPPORTED_FORMATS,
)
from mddb_workflow.utils.file import File


_UNSET = object()
_MISSING_INPUT = object()


class ProjectInputs:
    """Load and resolve the inputs document for one project.

    This is the inputs abstraction for one project.  It provides the shared
    behavior needed by the workflow: loading, schema validation, forced-value
    precedence, defaults, input-file persistence, MD configuration merging,
    and autofilling.

    Args:
        inputs_file: The input document to load.
        directory_name: Value used when replacing ``$DIR`` in the document.
        forced_inputs: Explicit values that take precedence over file values,
            supplied as raw ``-fin`` pairs.
        input_loader: Optional callback accepting the input ``File`` and
            downloading or otherwise making it available in place.
    """

    def __init__(self,
        inputs_file: File,
        directory_name: str,
        forced_inputs: list[list[str]] | None = None,
        input_loader: Callable[[File], None] | None = None,
    ):
        self.inputs_file = inputs_file
        self.directory_name = directory_name
        self.forced_inputs = self._normalize_forced_inputs(forced_inputs)
        self.defaults = DEFAULT_INPUT_VALUES
        self.default_warnings = DEFAULT_INPUT_VALUE_WARNINGS
        self.input_loader = input_loader

        # Set the file inputs, where values from the inputs file will be stored
        self._file_inputs: dict[str, Any] | None = None

        # Keep track of warnings already issued in order to not repeat them
        self._already_warned: set[str] = set()

        # Validate explicit overrides immediately.  This keeps malformed
        # command-line field names from surviving until the first lookup.
        if self.forced_inputs:
            validate_inputs(self.forced_inputs, strict_unknown=True)

        # Check if the inputs file is available at this point
        if not self.is_file_available:
            warn('Missing inputs file. Allowed tasks will be very limited.')

        # Persist forced inputs after ProjectInputs has validated them.  This
        # prevents the user from losing a lot of time for a silly typo.
        self.apply_forced_inputs()

    @staticmethod
    def _normalize_forced_inputs(forced_inputs: list[list[str]] | None) -> dict[str, Any]:
        """Normalize raw ``-fin`` pairs into one input mapping."""
        if not forced_inputs:
            return {}
        # Save forced inputs and use them when necessary
        for forced_input in forced_inputs:
            # Make sure the format is respected
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
        return {name: value for name, value in forced_inputs}

    def apply_forced_inputs(self) -> None:
        """Persist forced values in the input document when it is available."""
        # Overwrite input file values
        for input_name, input_value in self.forced_inputs.items():
            self.update_file_inputs(input_name, input_value)

    def _ensure_file_available(self) -> bool:
        """Make the input file available, downloading it when configured."""
        if self.file_exists:
            return True
        if self.input_loader is None:
            return False
        self.input_loader(self.inputs_file)
        return self.inputs_file.exists

    @property
    def file_exists(self) -> bool:
        """Whether the input document currently exists locally."""
        return self.inputs_file.exists

    @property
    def is_file_available(self) -> bool:
        """Whether the input document exists or can be loaded on demand."""
        # A loader means the file can be made available on demand, matching
        # the project's existing remote-input behavior.
        return self.file_exists or self.input_loader is not None

    def get_inputs_file(self) -> File:
        """Return the input file, loading it when a loader is configured."""
        if not self._ensure_file_available():
            raise InputError(f'Missing inputs file "{self.inputs_file.filename}"')
        return self.inputs_file

    @property
    def file_inputs(self) -> dict[str, Any]:
        """Return the validated values loaded from the input document."""
        # If inputs are already loaded then return them
        if self._file_inputs is not None:
            return self._file_inputs
        # If we have no inputs file at this point, try to make it available
        # through the configured loader before reporting it as missing.
        self.get_inputs_file()
        # When loading the inputs file, replace some values automatically
        replaces = [('$DIR', self.directory_name)]
        # If we have an inputs file stored locally, load it according to its format
        if self.inputs_file.format == 'json':
            inputs_data = aux.load_json(self.inputs_file.path, replaces)
        elif self.inputs_file.format == 'yaml':
            inputs_data = aux.load_yaml(self.inputs_file.path, replaces)
        else:
            raise InputError('Input file format is not supported. Please use json or yaml files.')
        if not inputs_data:
            raise InputError('Input file is empty')
        # Legacy fixes (applied before validation so legacy values are validated too)
        old_pdb_ids = inputs_data.pop('pdbIds', None)
        if old_pdb_ids:
            warn('Inputs "pdbIds" field is deprecated. Replacing it by "pdb_ids".')
            inputs_data['pdb_ids'] = old_pdb_ids
        # Validate the inputs against the workflow schema before using them
        # This raises an InputError with a clear per-field message on any problem
        validate_inputs(inputs_data)
        self._file_inputs = inputs_data
        # Finally return the updated inputs
        return self._file_inputs

    def update_file_inputs(self, nested_key: str, new_value: Any) -> bool:
        """Permanently update the inputs file.
        This may be done when command line inputs do not match file inputs.
        Return True if the inputs is updated correctly. Return False if there is no update.
        """
        # If there is no inputs file then there is nothing to update
        if not self.is_file_available: return False
        # If the input already matches then do nothing
        current_value = aux.read_ndict(self.file_inputs, nested_key, _MISSING_INPUT)
        if current_value is not _MISSING_INPUT and current_value == new_value: return False
        # Set the new value
        aux.write_ndict(self.file_inputs, nested_key, new_value)
        print(f'* Field "{nested_key}" in the inputs file will be permanently modified')
        # Write the new inputs to disk
        if self.inputs_file.format == 'json':
            aux.save_json(self.file_inputs, self.inputs_file.path)
        elif self.inputs_file.format in ['yaml', 'yml']:
            # Note that comments in the original YAML file will be not kept
            aux.save_yaml(self.file_inputs, self.inputs_file.path)
        else:
            raise InputError('Input file format is not supported. Please use json or yaml files.')
        return True

    def merge_md_config(
        self,
        input_md_config: list[list[str]] | None,
        input_md_directories: list[str] | None = None,
        input_topology_filepath: str | None = None,
        input_structure_filepath: str | None = None,
        input_trajectory_filepaths: list[str] | None = None,
    ) -> list[dict]:
        """Load, merge, and persist the project's CLI and file MD configuration."""
        # Get MD configuration from the inputs file, if any
        # An absent or null optional value is represented internally as an empty list
        # Detached copy so changes remain distinguishable for update_file_inputs equality checks
        md_config: list[dict] = deepcopy(self.get('mds') or [])
        # Purge removed MDs from the list
        active_md_config = [c for c in md_config if not c.get(MD_REMOVED_FLAG)]
        # Merge MD configs from inputs file and the CLI arguments
        # First scenario argument: the new way
        if input_md_config:
            # Make sure all MD configurations have at least 2 values each
            for mdc in input_md_config:
                if len(mdc) >= 2: continue
                raise InputError('Wrong MD configuration: the patter is -md <directory> <trajectory> <trajectory 2> ...')
            # Make sure there are no duplicated MD directories
            arg_md_directories = [mdc[0] for mdc in input_md_config]
            if len(arg_md_directories) > len(set(arg_md_directories)):
                raise InputError('There are duplicated MD directories. Every directory behind every "-md" must be unique.')
            # Iterate MD config entries and add / merge them with the already existing MD config
            for arg_config in input_md_config:
                md_dir = arg_config[0]
                # Find if there is already a MD configuration for this directory
                config = next((c for c in active_md_config if c[MD_DIRECTORY] == md_dir), None)
                # Otherwise create a new one and add it to the list
                if config is None:
                    config = {MD_DIRECTORY: md_dir}
                    md_config.append(config)
                # An input topology/structure for a specific MD may be passed before the trajectory
                # In order to tell if the topology/structure was passed we check input file formats
                # Note that PDB format is both a structure and a trajectory supported format
                has_structure = False
                has_topology = False
                if len(arg_config) > 2:
                    first_sample = File(arg_config[1])
                    second_sample = File(arg_config[2])
                    if first_sample.format != second_sample.format:
                        if first_sample.format in TOPOLOGY_SUPPORTED_FORMATS:
                            has_topology = True
                        else:
                            has_structure = True
                # Finally set the input topology, structure and trajectories files
                md_input_topology_filepath = arg_config[1] if has_topology else input_topology_filepath
                md_input_structure_filepath = arg_config[1] if has_structure else input_structure_filepath
                md_input_trajectory_filepaths = arg_config[2:] if has_structure or has_topology else arg_config[1:]
                # Add the input files to the MD configuration
                config.update({
                    MD_INPUT_TOPOLOGY_FILEPATH: md_input_topology_filepath,
                    MD_INPUT_STRUCTURE_FILEPATH: md_input_structure_filepath,
                    MD_INPUT_TRAJECTORY_FILEPATHS: md_input_trajectory_filepaths,
                })
        # Second scenario argument: the old way
        elif input_md_directories:
            warn('The "-mdir" argument is deprecated. Please consider using the "-md" argument. Use the "--help" argument to see how it works.')
            parsed_md_directories = []
            # Parse any glob notation
            for md_dir in input_md_directories:
                if aux.is_glob(md_dir):
                    parsed_md_directories.extend(aux.parse_glob(md_dir))
                else:
                    parsed_md_directories.append(md_dir)
            # Match the modern -md behavior and reject repeated directories
            if len(parsed_md_directories) > len(set(parsed_md_directories)):
                raise InputError('There are duplicated MD directories. Every directory behind "-mdir" must be unique.')
            # For each parsed MD directory, check if it is already defined in the config file
            # If not then add a new config for this
            for md_dir in parsed_md_directories:
                # Find if there is already a MD configuration for this directory
                config = next((c for c in active_md_config if c[MD_DIRECTORY] == md_dir), None)
                # Otherwise create a new one and add it to the list
                if config is None:
                    config = {MD_DIRECTORY: md_dir}
                    md_config.append(config)
        # Add the generic topology, structure or trajectory arguments, if any, to the MD config
        for config in md_config:
            if config.get(MD_REMOVED_FLAG): continue
            if input_topology_filepath:
                config[MD_INPUT_TOPOLOGY_FILEPATH] = input_topology_filepath
            if input_structure_filepath:
                config[MD_INPUT_STRUCTURE_FILEPATH] = input_structure_filepath
            if input_trajectory_filepaths:
                config[MD_INPUT_TRAJECTORY_FILEPATHS] = input_trajectory_filepaths
        # Persist the complete merged MD configuration before the Project checks
        # and fills generated names and directories.
        self.update_file_inputs('mds', md_config)
        if self.is_file_available:
            # Return the file_inputs dict so Project and MDs stay in sync
            return self.file_inputs['mds']
        return md_config

    @staticmethod
    def _name_to_directory(name: str) -> str:
        """Convert an MD name into an equivalent MD directory."""
        # Replace white spaces with underscores
        directory = name.replace(' ', '_')
        # Remove problematic characters
        for character in FORBIDDEN_DIRECTORY_CHARACTERS:
            directory = directory.replace(character, '')
        return directory

    @staticmethod
    def _directory_to_name(directory: str) -> str:
        """Convert an MD directory into an equivalent MD name."""
        # The normpath prevents a possible ending '/' which would make this not work
        # Replace underscores with white spaces
        return normpath(directory).split('/')[-1].replace('_', ' ')

    def validate_md_config(self, md_config: list[dict]) -> list[dict]:
        """Validate MD presence, fill identifiers, and enforce uniqueness."""
        # There must be at least one MD configuration
        if len(md_config) == 0:
            raise InputError('There must be at least one MD')
        # Run a few checks to make sure all inputs are coherent
        # Make sure all input MDs have unique name and directory
        names = {}
        directories = {}
        for md_index, md_inputs in enumerate(md_config):
            # Skip removed MDs
            if md_inputs.get(MD_REMOVED_FLAG):
                continue
            # Make sure the MD has at least a name or a directory
            directory = md_inputs.get(MD_DIRECTORY, None)
            name = md_inputs.get(MD_NAME, None)
            if not directory and not name:
                raise InputError(f'There is a MD (index {md_index}) with no name and no directory.' +
                    ' Please define at least one of them.')
            # Now fill the gaps
            # If there is a name and not a directory then issue a directory using the name
            if directory is None and name is not None:
                directory = self._name_to_directory(name)
                # Update before md_inputs
                self.update_file_inputs(f'mds.{md_index}.{MD_DIRECTORY}', directory)
                md_inputs[MD_DIRECTORY] = directory
            # If there is a directory and not a name then issue a name using the directory
            if name is None and directory is not None:
                name = self._directory_to_name(directory)
                # Update before md_inputs
                self.update_file_inputs(f'mds.{md_index}.{MD_NAME}', name)
                md_inputs[MD_NAME] = name
            # Add current names and directories to the counts
            current_name_count = names.get(name, 0)
            names[name] = current_name_count + 1
            current_directory_count = directories.get(directory, 0)
            directories[directory] = current_directory_count + 1
        # If there were any duplicates then report it
        repeats = False
        for name, name_count in names.items():
            if name_count == 1:
                continue
            warn(f'There are {name_count} MDs with the same name: {name}')
            repeats = True
        for directory, directory_count in directories.items():
            if directory_count == 1:
                continue
            warn(f'There are {directory_count} MDs with the same directory: {directory}')
            repeats = True
        if repeats:
            raise InputError('Duplicated values in MD inputs (see warnings above).' +
                ' All MD names and directories must be unique.')
        return md_config

    def validate(self, inputs_data: dict[str, Any] | None = None) -> dict[str, Any]:
        """Validate a document or the current effective values.

        The returned dictionary is the same object when ``inputs_data`` is
        supplied, matching the behavior of ``validate_inputs``.
        """
        if inputs_data is None:
            inputs_data = self.effective_inputs
        if not isinstance(inputs_data, dict):
            inputs_data = dict(inputs_data)
        validate_inputs(inputs_data)
        return inputs_data

    def ensure_file_from_template(self, template_file: File) -> bool:
        """Create the inputs file from a template when it is not available.

        Return True when a file was created and False when the file was
        already available, including when a loader can provide it lazily.
        """
        # Make a copy of the template in the local directory if there is not an inputs file yet
        if self.is_file_available:
            print(" Inputs file already exists")
            return False
        template_file.copy_to(self.inputs_file)
        print(f" File {self.inputs_file.path} has been created from a template")
        return True

    def autofill(self, values: dict[str, Any]) -> dict[str, bool]:
        """Persist generated input values and report which fields changed."""
        # Now fill the inputs file with the values generated by the workflow
        return {
            nested_key: self.update_file_inputs(nested_key, new_value)
            for nested_key, new_value in values.items()
        }

    @property
    def effective_inputs(self) -> dict[str, Any]:
        """Return explicit file values with forced values applied.

        Defaults are deliberately not materialized here.  They are resolved
        by :meth:`get`, so merely inspecting effective inputs does not turn
        implicit defaults into persistent file values.
        """
        values = deepcopy(self.file_inputs) if self.is_file_available else {}
        values.update(self.forced_inputs)
        return values

    def get(self, name: str, default: Any = _UNSET) -> Any:
        """Resolve one value using forced, file, and default precedence."""
        # Check if the value was forced from the command line
        if name in self.forced_inputs:
            return self.forced_inputs[name]
        # If the inputs file is available, get the input value from it
        if self.is_file_available:
            value = self.file_inputs.get(name, _UNSET)
            # If we had a value then return it
            if value is not _UNSET:
                return value

        # If the field is not specified in the inputs file then return the
        # explicitly provided fallback, when one exists
        if default is not _UNSET:
            return default
        # Otherwise return the workflow default value
        return self.get_default(name)

    def get_default(self, name: str) -> Any:
        """Return a default value and issue its warning at most once."""
        # If the field is not specified in the inputs file then set a default value
        default_value = self.defaults.get(name, None)
        # Warn the user about this value we are about to assign automatically
        # If the value is None then it means there is no default value
        if default_value is not None and name not in self._already_warned:
            warn(f'Missing input "{name}" -> Using default value: {default_value}')

        # If there is an additional warning then display it here
        additional_warning = self.default_warnings.get(name)
        if additional_warning is not None and name not in self._already_warned:
            warn(additional_warning)

        self._already_warned.add(name)
        return default_value
