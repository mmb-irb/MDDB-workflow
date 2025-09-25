# Auxiliar generic functions and classes used along the workflow

from model_workflow import __path__
from model_workflow.utils.constants import RESIDUE_NAME_LETTERS, PROTEIN_RESIDUE_NAME_LETTERS
from model_workflow.utils.constants import YELLOW_HEADER, COLOR_END

import os
from os import rename, listdir, remove
from os.path import isfile, exists
import re
import sys
import json
import yaml
from glob import glob
from typing import Optional, List, Set, Generator, Union
from struct import pack
# NEVER FORGET: GraphQL has a problem with urllib.parse -> It will always return error 400 (Bad request)
# We must use requests instead
import requests
from subprocess import run, PIPE


# Check if a module has been imported
def is_imported (module_name : str) -> bool:
    return module_name in sys.modules

# Set custom exception which is not to print traceback
# They are used when the problem is not in our code
class QuietException (Exception):
    pass

# Set a custom quite exception for when user input is wrong
class InputError (QuietException):
    pass

# Set a custom quite exception for when MD data has not passed a quality check test
class TestFailure (QuietException):
    pass

# Set a custom quite exception for when the problem comes from a third party dependency
class ToolError (QuietException):
    pass

# Set a custom quite exception for when the problem comes from a remote service
class RemoteServiceError (QuietException):
    pass

# Set a no referable exception for PDB synthetic constructs or chimeric entities
class NoReferableException (Exception):
    def __str__ (self): return f'No referable sequence {self.sequence}'
    def __repr__ (self): return self.__str__()
    def get_sequence (self) -> str:
        return self.args[0]
    sequence = property(get_sequence, None, None, 'Aminoacids sequence')

# Set a custom exception handler where our input error exception has a quiet behaviour
def custom_excepthook (exception_class, message, traceback):
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

# Set a function to get the next letter from an input letter in alphabetic order
# Return None if we run out of letters
letters = { 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
    'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z' }
def get_new_letter(current_letters : set) -> Optional[str]:
    return next((letter for letter in letters if letter not in current_letters), None)

# Given a list or set of names, return a set with all case-posibilites:
# All upper case
# All lower case
# First upper case and the rest lower case
def all_cases (names : Union[List[str], Set[str]]) -> Set[str]:
    all_names = []
    for name in names:
        all_upper = name.upper()
        all_lower = name.lower()
        one_upper = name[0].upper() + name[1:].lower()
        all_names += [ all_upper, all_lower, one_upper ]
    return set(all_names)

# Given a residue name, return its single letter
def residue_name_to_letter (residue_name : str) -> str:
    return RESIDUE_NAME_LETTERS.get(residue_name, 'X')

# Given a protein residue name, return its single letter
def protein_residue_name_to_letter (residue_name : str) -> str:
    return PROTEIN_RESIDUE_NAME_LETTERS.get(residue_name, 'X')

# Set a JSON loader with additional logic to better handle problems
def load_json (filepath : str) -> dict: 
    try:
        with open(filepath, 'r') as file:
            content = json.load(file)
        return content
    except Exception as error:
        raise Exception(f'Something went wrong when loading JSON file {filepath}: {str(error)}')

# Set a JSON saver with additional logic to better handle problems
def save_json (content, filepath : str, indent : Optional[int] = None):
    try:
        with open(filepath, 'w') as file:
            json.dump(content, file, indent=indent)
    except Exception as error:
        # Rename the JSON file since it will be half written thus giving problems when loaded
        rename(filepath, filepath + '.wrong')
        raise Exception(f'Something went wrong when saving JSON file {filepath}: {str(error)}')

# Set a YAML loader with additional logic to better handle problems
# DANI: Por algún motivo yaml.load también funciona con archivos en formato JSON
def load_yaml (filepath : str):
    try:
        with open(filepath, 'r') as file:
            content = yaml.load(file, Loader=yaml.CLoader)
        return content
    except Exception as error:
        warn(str(error).replace('\n', ' '))
        raise InputError('Something went wrong when loading YAML file ' + filepath)

# Set a YAML saver with additional logic to better handle problems
def save_yaml (content, filepath : str):
    with open(filepath, 'w') as file:
        yaml.dump(content, file)

# Set a few constants to erase previou logs in the terminal
CURSOR_UP_ONE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'
ERASE_PREVIOUS_LINES = CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE

# Set a function to remove previous line
def delete_previous_log ():
    print(ERASE_PREVIOUS_LINES)

# Set a function to reprint in the same line
def reprint (text : str):
    delete_previous_log()
    print(text)

# Set a function to print a messahe with a colored warning header
def warn (message : str):
    print(YELLOW_HEADER + '⚠  WARNING: ' + COLOR_END + message)

# Get the mean/average of a list of values
def mean(values : List[float]) -> float:
    return sum(values) / len(values)

# Round a number to hundredths
def round_to_hundredths (number : float) -> float:
    return round(number * 100) / 100

# Round a number to hundredths
def round_to_thousandths (number : float) -> float:
    return round(number * 1000) / 1000

# Given a list with numbers,  create a string where number in a row are represented rangedly
# e.g. [1, 3, 5, 6, 7, 8] => "1, 3, 5-8"
def ranger (numbers : List[int]) -> str:
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

# Set a special iteration system
# Return one value of the array and a new array with all other values for each value
def otherwise (values : list) -> Generator[tuple, None, None]:
    for v, value in enumerate(values):
        others = values[0:v] + values[v+1:]
        yield value, others

# List files in a directory
def list_files (directory : str) -> List[str]:
    return [f for f in listdir(directory) if isfile(f'{directory}/{f}')]

# Check if a directory is empty
def is_directory_empty (directory : str) -> bool:
    return len(listdir(directory)) == 0

# Set a function to check if a string has patterns to be parsed by a glob function
# Note that this is not trivial, but this function should be good enough for our case
# https://stackoverflow.com/questions/42283009/check-if-string-is-a-glob-pattern
GLOB_CHARACTERS = ['*', '?', '[']
def is_glob (path : str) -> bool:
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

# Parse a glob path into one or several results
# If the path has no glob characters then return it as it is
# Otherwise make sure
def parse_glob (path : str) -> List[str]:
    # If there is no glob pattern then just return the string as is
    if not is_glob(path):
        return [ path ]
    # If there is glob pattern then parse it
    parsed_filepaths = glob(path)
    return parsed_filepaths

# Supported byte sizes
SUPPORTED_BYTE_SIZES = {
    2: 'e',
    4: 'f',
    8: 'd'
}

# Data is a list of numeric values
# Bit size is the number of bits for each value in data to be occupied
def store_binary_data (data : List[float], byte_size : int, filepath : str):
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

# Capture all stdout or stderr within a code region even if it comes from another non-python threaded process
# https://stackoverflow.com/questions/24277488/in-python-how-to-capture-the-stdout-from-a-c-shared-library-to-a-variable
class CaptureOutput (object):
    escape_char = "\b"
    def __init__(self, stream : str = 'stdout'):
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
            char = os.read(self.pipe_out,1).decode(self.original_stream.encoding)
            if not char or self.escape_char in char:
                break
            self.captured_text += char

# Set a function to request data to the PDB GraphQL API
# Note that this function may be used for either PDB ids or PDB molecule ids, depending on the query
# The query parameter may be constructed using the following page:
# https://data.rcsb.org/graphql/index.html
def request_pdb_data (pdb_id : str, query : str) -> dict:
    # Make sure the PDB id is valid as we set the correct key to mine the response data
    if len(pdb_id) == 4: data_key = 'entry'
    elif len(pdb_id) < 4 or len(pdb_id) > 4: data_key = 'chem_comp'
    else: raise ValueError(f'Wrong PDB id "{pdb_id}". It must be 4 (entries) or less (ligands) characters long')
    # Set the request URL
    request_url = 'https://data.rcsb.org/graphql'
    # Set the POST data
    post_data = {
        "query": query,
        "variables": { "id": pdb_id }
    }
    # Send the request
    try:
        response = requests.post(request_url, json=post_data)
    except requests.exceptions.ConnectionError as error:
        raise ConnectionError('No internet connection :(') from None
    # Get the response
    parsed_response = json.loads(response.text)['data'][data_key]
    if parsed_response == None:
        new_pdb_id = request_replaced_pdb(pdb_id)
        if new_pdb_id:
            parsed_response = request_pdb_data(new_pdb_id, query)
        else:
            print(f'PDB id {pdb_id} not found')
    return parsed_response

# Use the RCSB REST API to get the replaced PDB id
# This is useful when the PDB is obsolete and has been replaced
def request_replaced_pdb(pdb_id):
    query_url = 'https://data.rcsb.org/rest/v1/holdings/removed/'+pdb_id
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

# Given a filename, set a sufix number on it, right before the extension
# Set also the number of zeros to fill the name
def numerate_filename (filename : str, number : int, zeros : int = 2, separator : str = '_') -> str:
    splits = filename.split('.')
    sufix = separator + str(number).zfill(zeros)
    return '.'.join(splits[0:-1]) + sufix + '.' + splits[-1]

# Given a filename, set a sufix including '*', right before the extension
# This should match all filenames obtanied through the numerate_filename function when used in bash
def glob_filename (filename : str, separator : str = '_') -> str:
    splits = filename.split('.')
    sufix = separator + '*'
    return '.'.join(splits[0:-1]) + sufix + '.' + splits[-1]

# Delete all files matched by the glob_filename function
def purge_glob (filename : str):
    glob_pattern = glob_filename(filename)
    existing_outputs = glob(glob_pattern)
    for existing_output in existing_outputs:
        if exists(existing_output): remove(existing_output)

# Given a filename with the the pattern 'mda.xxxx.json', get the 'xxxx' out of it
def get_analysis_name (filename : str) -> str:
    name_search = re.search(r'/mda.([A-Za-z0-9_-]*).json$', filename)
    if not name_search:
        raise ValueError(f'Wrong expected format in filename {filename}')
    # To make it coherent with the rest of analyses, the analysis name become parsed when loaded in the database
    # Every '_' is replaced by '-' so we must keep the analysis name coherent or the web client will not find it
    return name_search[1].replace('_', '-')

# Use a safe alternative to hasattr/getattr
# DANI: Do not use getattr with a default argument or hasattr
# DANI: If you do, you will loose any further AtributeError(s)
# DANI: Thus you will have very silent errors every time you have a silly typo
# DANI: This is a python itself unresolved error https://bugs.python.org/issue39865
def safe_hasattr (instance, attribute_name : str) -> bool:
    return attribute_name in set(dir(instance))
def safe_getattr (instance, attribute_name : str, defualt):
    if not safe_hasattr(instance, attribute_name): return defualt
    return getattr(instance, attribute_name)
    
# Function to read and write a dict nested value with using a single combined key

# Read a value in a nested dictionary and return the placeholder if any key in the path does not exist
def read_ndict (nested_dict : dict, nested_key : str, placeholder = KeyError('Missing nested key')):
    keys = nested_key.split('.')
    value = nested_dict
    for key in keys:
        # support list indices
        if key.isdigit():
            index = int(key)
            value = value[index]
        # support dict keys
        else:
            value = value.get(key, placeholder)
            if value == placeholder: return placeholder
    return value

# Write a value in a nested dictionary and raise an error if any key in the path s missing
def write_ndict (nested_dict : dict, nested_key : str, value):
    keys = nested_key.split('.')
    nested_keys = keys[0:-1]
    next_target = nested_dict
    for key in nested_keys:
        # support list indices
        if key.isdigit():
            index = int(key)
            next_target = next_target[index]
        # support dict keys
        else:
            missing_key_error = KeyError(f'Missing nested key {key}')
            next_target = next_target.get(key, missing_key_error)
            if next_target == missing_key_error: return missing_key_error
    field = keys[-1]
    next_target[field] = value
    
# Get the current git version
def get_git_version () -> str:
    git_command = f"git -C {__path__[0]} describe"
    process = run(git_command, shell=True, stdout=PIPE)
    return process.stdout.decode().replace('\n','')