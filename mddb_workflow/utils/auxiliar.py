# Auxiliar generic functions and classes used along the workflow

from mddb_workflow import __path__, __version__
from mddb_workflow.utils.constants import RESIDUE_NAME_LETTERS, PROTEIN_RESIDUE_NAME_LETTERS
from mddb_workflow.utils.constants import YELLOW_HEADER, COLOR_END
from mddb_workflow.utils.constants import STANDARD_TOPOLOGY_FILENAME
from mddb_workflow.utils.type_hints import *

import os
from os.path import isfile, exists
import re
import sys
import json
import yaml
from glob import glob
from struct import pack
# NEVER FORGET: GraphQL has a problem with urllib.parse -> It will always return error 400 (Bad request)
# We must use requests instead
import requests
import urllib.request
from subprocess import run, PIPE
from dataclasses import asdict, is_dataclass


def is_imported(module_name: str) -> bool:
    """Check if a module has been imported."""
    return module_name in sys.modules


class QuietException (Exception):
    """Exception which is not to print traceback.
    They are used when the problem is not in our code.
    """
    pass


class InputError (QuietException):
    """Quite exception for when user input is wrong."""
    pass


class TestFailure (QuietException):
    """Quite exception for when MD data has not passed a quality check test."""
    pass


class EnvironmentError (QuietException):
    """Quite exception for when the problem is not in the code but in the environment."""
    pass


class ToolError (QuietException):
    """Quite exception for when the problem comes from a third party dependency."""
    pass


class RemoteServiceError (QuietException):
    """Quite exception for when the problem comes from a remote service."""
    pass


class ForcedStop (QuietException):
    """Quite exception for when we stop the workflow in purpose."""
    pass


class NoReferableException (Exception):
    """No referable exception for PDB synthetic constructs or chimeric entities."""
    def __str__(self): return f'No referable sequence {self.sequence}'
    def __repr__(self): return self.__str__()
    def get_sequence(self) -> str:
        return self.args[0]
    sequence = property(get_sequence, None, None, 'Aminoacids sequence')


def custom_excepthook(exception_class, message, traceback):
    """Handle a custom exception where our input error exception has a quiet behaviour."""
    # Quite behaviour if it is our input error exception
    if QuietException in exception_class.__bases__:
        print('{0}: {1}'.format(exception_class.__name__, message))  # Only print Error Type and Message
        return
    # Default behaviour otherwise
    sys.__excepthook__(exception_class, message, traceback)


sys.excepthook = custom_excepthook

# Set a special exceptions for when the topology is missing
MISSING_TOPOLOGY = Exception('Missing topology')
MISSING_CHARGES = Exception('Missing atom charges')
MISSING_BONDS = Exception('Missing atom bonds')
JSON_SERIALIZABLE_MISSING_BONDS = 'MB'

# Keep all exceptions in a set
STANDARD_EXCEPTIONS = {MISSING_TOPOLOGY, MISSING_CHARGES, MISSING_BONDS}

# Set a function to get the next letter from an input letter in alphabetic order
# Return None if we run out of letters
letters = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
    'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'}


def get_new_letter(current_letters: set) -> Optional[str]:
    return next((letter for letter in letters if letter not in current_letters), None)


def all_cases(names: list[str] | set[str]) -> set[str]:
    """Given a list or set of names.

    Return a set with all case-posibilites:
    - All upper case
    - All lower case
    - First upper case and the rest lower case

    """
    all_names = []
    for name in names:
        all_upper = name.upper()
        all_lower = name.lower()
        one_upper = name[0].upper() + name[1:].lower()
        all_names += [all_upper, all_lower, one_upper]
    return set(all_names)


def residue_name_to_letter(residue_name: str) -> str:
    """Given a residue name, return its single letter."""
    return RESIDUE_NAME_LETTERS.get(residue_name, 'X')


def protein_residue_name_to_letter(residue_name: str) -> str:
    """Given a protein residue name, return its single letter."""
    return PROTEIN_RESIDUE_NAME_LETTERS.get(residue_name, 'X')


OBJECT_TYPES = {dict, list, tuple}


def recursive_transformer(target_object: dict | list | tuple, transformer: Optional[Callable] = None) -> dict | list | tuple:
    """Recursive transformer for nested dicts and lists.
    Note that a new object is created to prevent mutating the original.
    If no transformer function is passed then it becomes a recursive cloner.
    WARNING: the transformer function could mutate some types of the original object if not done properly.
    """
    object_type = type(target_object)
    # Get a starting object and entries iterator depending on the object type
    if object_type is dict:
        clone = {}
        entries = target_object.items()
    elif object_type is list or object_type is tuple:
        # Note that if it is a tuple we make a list anyway and we convert it to tuple at the end
        # WARNING: Note that tuples do not support reassigning items
        clone = [None for i in range(len(target_object))]
        entries = enumerate(target_object)
    else: ValueError(f'The recursive cloner should only be applied to object and lists, not to {object_type}')
    # Iterate the different entries in the object
    for index_or_key, value in entries:
        # Get he value type
        value_type = type(value)
        # If it is a dict or list then call the transformer recursively
        if value_type in OBJECT_TYPES: clone[index_or_key] = recursive_transformer(value, transformer)
        # If it is not an object type then apply the transformer to it
        else: clone[index_or_key] = transformer(value) if transformer else value
    # If it was a tuple then make the conversion now
    if object_type == tuple: clone = tuple(clone)
    return clone


# Set some headers for the serializer
EXCEPTION_HEADER = 'Exception: '


# LORE: This was originally intended to support exceptions in the cache
def json_serializer(object: dict | list | tuple) -> dict | list | tuple:
    """Serialize a standard JSON with support for additional types."""
    def serializer(value):
        # If we have exceptions then convert them to text with an appropiate header
        if type(value) is Exception:
            return f'{EXCEPTION_HEADER}{value}'
        # This must be done before the set check because asdict() will create dicts with sets inside
        if is_dataclass(value) and not isinstance(value, type):
            dict_value = asdict(value)
            return recursive_transformer(dict_value, serializer)
        if isinstance(value, set):
            return list(value)
        # If the type is not among the ones we check then assume it is already serializable
        return value
    object_clone = recursive_transformer(object, serializer)
    return object_clone


def json_deserializer(object: dict | list | tuple) -> dict | list | tuple:
    """Deserialize a standard JSON with support for additional types."""
    def deserializer(value):
        # Check if there is any value which was adapted to become JSON serialized and restore it
        if type(value) is str and value[0:11] == EXCEPTION_HEADER:
            # WARNING: Do not declare new exceptions here but use the constant ones
            # WARNING: Otherwise further equality comparisions will fail
            exception_message = value[11:]
            standard_exception = next((exception for exception in STANDARD_EXCEPTIONS if str(exception) == exception_message), None)
            if standard_exception is None:
                raise ValueError(f'Exception "{exception_message}" is not among standard exceptions')
            return standard_exception
        # If the type is not among the ones we check then assume it is already deserialized
        return value
    object_clone = recursive_transformer(object, deserializer)
    return object_clone


def load_json(filepath: str, replaces: Optional[list[tuple]] = []) -> dict:
    """Load a JSON with additional logic to better handle problems."""
    try:
        with open(filepath, 'r') as file:
            content = file.read()
            # Make pure text replacements
            for replace in replaces:
                target, replacement = replace
                content = content.replace(target, replacement)
            # Parse the content to JSON
            parsed_content = json.loads(content)
            # Deserialize some types like Exceptions, which are stored in JSONs in this context
            deserialized_content = json_deserializer(parsed_content)
        return deserialized_content
    except Exception as error:
        raise Exception(f'Something went wrong when loading JSON file {filepath}: {str(error)}')


def save_json(content, filepath: str, indent: Optional[int] = None):
    """Save a JSON with additional logic to better handle problems."""
    try:
        with open(filepath, 'w') as file:
            serialized_content = json_serializer(content)
            json.dump(serialized_content, file, indent=indent)
    except Exception as error:
        # Rename the JSON file since it will be half written thus giving problems when loaded
        os.rename(filepath, filepath + '.wrong')
        raise Exception(f'Something went wrong when saving JSON file {filepath}: {str(error)}')


# DANI: Por algún motivo yaml.load también funciona con archivos en formato JSON
def load_yaml(filepath: str, replaces: Optional[list[tuple]] = []) -> dict:
    """Load a YAML with additional logic to better handle problems.
    The argument replaces allows to replace file content before beeing processed.
    Every replace is a tuple whith two values: the target and the replacement.
    """
    try:
        with open(filepath, 'r') as file:
            content = file.read()
            # Pre-process fields that commonly contain colons (DOIs, URLs, etc.)
            # Match field: value where value isn't already quoted and may span multiple lines
            content = re.sub(
                r'^(citation:\s*)(?!")([^\n]+(?:\n(?!\w+:)[^\n]+)*)$',
                r'\1"\2"',
                content,
                flags=re.MULTILINE
            )
            for replace in replaces:
                target, replacement = replace
                content = content.replace(target, replacement)
            parsed_content = yaml.load(content, Loader=yaml.CLoader)
        return parsed_content
    except Exception as error:
        warn(str(error).replace('\n', ' '))
        raise InputError('Something went wrong when loading YAML file ' + filepath)


def save_yaml(content, filepath: str):
    """Save a YAML with additional logic to better handle problems."""
    with open(filepath, 'w') as file:
        yaml.dump(content, file)


# Set a few constants to erase previou logs in the terminal
CURSOR_UP_ONE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'
ERASE_PREVIOUS_LINES = CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE


def delete_previous_log():
    """Remove previous line in the terminal."""
    print(ERASE_PREVIOUS_LINES)


def reprint(text: str):
    """Reprint text in the same line in the terminal."""
    delete_previous_log()
    print(text)


def warn(message: str):
    """Print a message with a colored warning header."""
    print(YELLOW_HEADER + '⚠  WARNING: ' + COLOR_END + message)


def mean(values: list[float]) -> float:
    """Get the mean/average of a list of values."""
    return sum(values) / len(values)


def round_to_hundredths(number: float) -> float:
    """Round a number to hundredths."""
    return round(number * 100) / 100


def round_to_thousandths(number: float) -> float:
    """Round a number to thousandths."""
    return round(number * 1000) / 1000


def ranger(numbers: list[int]) -> str:
    """Given a list with numbers, create a string where number in a row are represented rangedly.

    Example:
        [1, 3, 5, 6, 7, 8] => "1, 3, 5-8"

    """
    # Remove duplicates and sort numbers
    sorted_numbers = sorted(list(set(numbers)))
    # Get the number of numbers in the list
    number_count = len(sorted_numbers)
    # If there is only one number then finish here
    if number_count == 1:
        return str(sorted_numbers[0])
    # Start the parsing otherwise
    ranged = ''
    last_number = -1
    # Iterate numbers
    for i, number in enumerate(sorted_numbers):
        # Skip this number if it was already included in a previous serie
        if i <= last_number: continue
        # Add current number to the ranged string
        ranged += ',' + str(number)
        # Now iterate numbers after the current number
        next_index = i+1
        for j, next_number in enumerate(sorted_numbers[next_index:], next_index):
            # Set the length of the serie
            length = j - i
            # End of the serie
            if next_number - number != length:
                # The length here is the length which broke the serie
                # i.e. if the length here is 2 the actual serie length is 1
                serie_length = length - 1
                if serie_length > 1:
                    last_serie_number = j - 1
                    previous_number = sorted_numbers[last_serie_number]
                    ranged += '-' + str(previous_number)
                    last_number = last_serie_number
                break
            # End of the selection
            if j == number_count - 1:
                if length > 1:
                    ranged += '-' + str(next_number)
                    last_number = j
    # Remove the first coma before returning the ranged string
    return ranged[1:]


def otherwise(values: list) -> Generator[tuple, None, None]:
    """Set a special iteration system.
    Return one value of the array and a new array with all other values for each value.
    """
    for v, value in enumerate(values):
        others = values[0:v] + values[v+1:]
        yield value, others


# List files in a directory
def list_files(directory: str) -> list[str]:
    """List files in a directory."""
    return [f for f in os.listdir(directory) if isfile(f'{directory}/{f}')]


# Check if a directory is empty
def is_directory_empty(directory: str) -> bool:
    """Check if a directory is empty."""
    return len(os.listdir(directory)) == 0


GLOB_CHARACTERS = ['*', '?', '[']


def is_glob(path: str) -> bool:
    """Check if a string has patterns to be parsed by a glob function.

    Note that this is not trivial, but this function should be good enough for our case.
    https://stackoverflow.com/questions/42283009/check-if-string-is-a-glob-pattern
    """
    # Find unescaped glob characters
    for c, character in enumerate(path):
        if character not in GLOB_CHARACTERS:
            continue
        if c == 0:
            return True
        previous_characters = path[c-1]
        if previous_characters != '\\':
            return True
    return False


def parse_glob(path: str) -> list[str]:
    """Parse a glob path into one or several results.

    If the path has no glob characters then return it as it is.
    Otherwise parse the glob pattern.
    """
    # If there is no glob pattern then just return the string as is
    if not is_glob(path):
        return [path]
    # If there is glob pattern then parse it
    parsed_filepaths = glob(path)
    return parsed_filepaths


def is_url(path: str) -> bool:
    """Return whether the passed string is a URL or not."""
    return path[0:4] == 'http'


def url_to_source_filename(url: str) -> str:
    """Set the filename of an input file downloaded from an input URL.

    In this scenario we are free to set our own paths or filenames.
    Note that the original name will usually be the very same output filename.
    In order to avoid both filenames being the same we will add a header here.
    """
    original_filename = url.split('/')[-1]
    return 'source_' + original_filename


def download_file(request_url: str, output_file: 'File'):
    """Download files from a specific URL."""
    print(f'Downloading file "{output_file.path}" from {request_url}\n')
    try:
        urllib.request.urlretrieve(request_url, output_file.path)
    except urllib.error.HTTPError as error:
        if error.code == 404:
            raise Exception(f'Missing remote file "{output_file.filename}"')
        # If we don't know the error then simply say something went wrong
        raise Exception(f'Something went wrong when downloading file "{output_file.filename}" from {request_url}')


def is_standard_topology(file: 'File') -> bool:
    """Check if a file is a standard topology.
    Note that the filename may include the source header.
    """
    return file.filename.endswith(STANDARD_TOPOLOGY_FILENAME)


# Supported byte sizes
SUPPORTED_BYTE_SIZES = {
    2: 'e',
    4: 'f',
    8: 'd'
}


def store_binary_data(data: list[float], byte_size: int, filepath: str):
    """Store binary data to a file.

    Args:
        data: A list of numeric values
        byte_size: The number of bytes for each value in data to be occupied
        filepath: The output file path

    """
    # Check bit size to make sense
    letter = SUPPORTED_BYTE_SIZES.get(byte_size, None)
    if not letter:
        raise ValueError(f'Not supported byte size {byte_size}, please select one of these: {", ".join(SUPPORTED_BYTE_SIZES.keys())}')
    # Set the binary format
    # '<' stands for little endian
    byte_flag = f'<{letter}'
    # Start writting the output file
    with open(filepath, 'wb') as file:
        # Iterate over data list values
        for value in data:
            value = float(value)
            file.write(pack(byte_flag, value))


class CaptureOutput (object):
    """Capture all stdout or stderr within a code region even if it comes from another non-python threaded process.

    https://stackoverflow.com/questions/24277488/in-python-how-to-capture-the-stdout-from-a-c-shared-library-to-a-variable
    """
    escape_char = "\b"
    def __init__(self, stream: str = 'stdout'):
        # Get sys stdout or stderr
        if not hasattr(sys, stream):
            raise ValueError(f'Unknown stream "{stream}". Expected stream value is "stdout" or "stderr"')
        self.original_stream = getattr(sys, stream)
        self.original_streamfd = self.original_stream.fileno()
        self.captured_text = ""
        # Create a pipe so the stream can be captured:
        self.pipe_out, self.pipe_in = os.pipe()
    def __enter__(self):
        self.captured_text = ""
        # Save a copy of the stream:
        self.streamfd = os.dup(self.original_streamfd)
        # Replace the original stream with our write pipe:
        os.dup2(self.pipe_in, self.original_streamfd)
        return self
    def __exit__(self, type, value, traceback):
        # Print the escape character to make the readOutput method stop:
        self.original_stream.write(self.escape_char)
        # Flush the stream to make sure all our data goes in before
        # the escape character:
        self.original_stream.flush()
        self.readOutput()
        # Close the pipe:
        os.close(self.pipe_in)
        os.close(self.pipe_out)
        # Restore the original stream:
        os.dup2(self.streamfd, self.original_streamfd)
        # Close the duplicate stream:
        os.close(self.streamfd)
    def readOutput(self):
        while True:
            char = os.read(self.pipe_out, 1).decode(self.original_stream.encoding)
            if not char or self.escape_char in char:
                break
            self.captured_text += char


def request_pdb_data(pdb_id: str, query: str) -> dict:
    """Request data to the PDB GraphQL API.

    Note that this function may be used for either PDB ids or PDB molecule ids, depending on the query.
    The query parameter may be constructed using the following page:
    https://data.rcsb.org/graphql/index.html
    """
    # Make sure the PDB id is valid as we set the correct key to mine the response data
    if len(pdb_id) == 4: data_key = 'entry'
    elif len(pdb_id) < 4 or len(pdb_id) > 4: data_key = 'chem_comp'
    else: raise ValueError(f'Wrong PDB id "{pdb_id}". It must be 4 (entries) or less (ligands) characters long')
    # Set the request URL
    request_url = 'https://data.rcsb.org/graphql'
    # Set the POST data
    post_data = {
        "query": query,
        "variables": {"id": pdb_id}
    }
    # Send the request
    try:
        response = requests.post(request_url, json=post_data)
    except requests.exceptions.ConnectionError as error:
        raise ConnectionError('No internet connection :(') from None
    # Get the response
    parsed_response = json.loads(response.text)['data'][data_key]
    if parsed_response is None:
        new_pdb_id = request_replaced_pdb(pdb_id)
        if new_pdb_id:
            parsed_response = request_pdb_data(new_pdb_id, query)
        else:
            print(f'PDB id {pdb_id} not found')
    return parsed_response


def request_replaced_pdb(pdb_id):
    """Use the RCSB REST API to get the replaced PDB id.

    This is useful when the PDB is obsolete and has been replaced.
    """
    query_url = 'https://data.rcsb.org/rest/v1/holdings/removed/' + pdb_id
    response = requests.get(query_url, headers={'Content-Type': 'application/json'})
    # Check if the response is OK
    if response.status_code == 200:
        try:
            return response.json()['rcsb_repository_holdings_removed']['id_codes_replaced_by'][0]
        except:
            print(f'Error when mining replaced PDB id for {pdb_id}')
            return None
    else:
        return None


def numerate_filename(filename: str, number: int, zeros: int = 2, separator: str = '_') -> str:
    """Given a filename, set a suffix number on it, right before the extension.

    Args:
        filename: The original filename
        number: The number to add as suffix
        zeros: The number of zeros to fill the name
        separator: The separator between filename and number

    """
    splits = filename.split('.')
    sufix = separator + str(number).zfill(zeros)
    return '.'.join(splits[0:-1]) + sufix + '.' + splits[-1]


def glob_filename(filename: str, separator: str = '_') -> str:
    """Given a filename, set a suffix including '*', right before the extension.
    This should match all filenames obtained through the numerate_filename function when used in bash.
    """
    splits = filename.split('.')
    sufix = separator + '*'
    return '.'.join(splits[0:-1]) + sufix + '.' + splits[-1]


def purge_glob(filename: str):
    """Delete all files matched by the glob_filename function."""
    glob_pattern = glob_filename(filename)
    existing_outputs = glob(glob_pattern)
    for existing_output in existing_outputs:
        if exists(existing_output): os.remove(existing_output)


def get_analysis_name(filename: str) -> str:
    """Given a filename with the pattern 'mda.xxxx.json', get the 'xxxx' out of it."""
    name_search = re.search(r'/mda.([A-Za-z0-9_-]*).json$', filename)
    if not name_search:
        raise ValueError(f'Wrong expected format in filename {filename}')
    # To make it coherent with the rest of analyses, the analysis name become parsed when loaded in the database
    # Every '_' is replaced by '-' so we must keep the analysis name coherent or the web client will not find it
    return name_search[1].replace('_', '-')


# DANI: Do not use getattr with a default argument or hasattr
# DANI: If you do, you will loose any further AtributeError(s)
# DANI: Thus you will have very silent errors every time you have a silly typo
# DANI: This is a python itself unresolved error https://bugs.python.org/issue39865
def safe_hasattr(instance, attribute_name: str) -> bool:
    """Use a safe alternative to hasattr."""
    return attribute_name in set(dir(instance))


def safe_getattr(instance, attribute_name: str, default):
    """Use a safe alternative to getattr."""
    if not safe_hasattr(instance, attribute_name): return default
    return getattr(instance, attribute_name)


def read_ndict(nested_dict: dict, nested_key: str, placeholder=KeyError('Missing nested key')):
    """Read a value in a nested dictionary using a single combined key.
    Return the placeholder if any key in the path does not exist.
    """
    keys = nested_key.split('.')
    value = nested_dict
    for key in keys:
        # support list indices
        if key.isdigit():
            if type(value) is not list: return placeholder
            index = int(key)
            value = value[index]
        # support dict keys
        else:
            if type(value) is not dict: return placeholder
            value = value.get(key, placeholder)
            if value == placeholder: return placeholder
    return value


def write_ndict(nested_dict: dict, nested_key: str, value):
    """Write a value in a nested dictionary using a single combined key.
    Raise an error if any key in the path is missing.
    """
    keys = nested_key.split('.')
    nested_keys = keys[0:-1]
    next_target = nested_dict
    for k, key in enumerate(nested_keys):
        # support list indices
        if key.isdigit():
            if type(next_target) is not list:
                raise ValueError(f'{".".join(nested_keys[0:k])} should be a list, but it is {next_target}')
            index = int(key)
            next_target = next_target[index]
        # support dict keys
        else:
            if type(next_target) is not dict:
                raise ValueError(f'{".".join(nested_keys[0:k])} should be a dict, but it is {next_target}')
            missing_key_error = KeyError(f'Missing nested key {key}')
            next_target = next_target.get(key, missing_key_error)
            if next_target == missing_key_error: raise missing_key_error
    field = keys[-1]
    next_target[field] = value


def get_git_version() -> str:
    """Get the current git version."""
    git_command = f"git -C {__path__[0]} describe"
    process = run(git_command, shell=True, stdout=PIPE)
    return process.stdout.decode().replace('\n', '') or __version__
