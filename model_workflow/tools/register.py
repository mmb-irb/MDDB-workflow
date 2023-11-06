from os.path import exists
from datetime import datetime
from typing import Optional
import json
import atexit

from model_workflow.constants import REGISTER_FILENAME
from model_workflow.tools.file import File

# The register tracks activity along multiple runs and thus avoids repeating some already succeeded tests
# It is also responsible for storing test failure warnings to be written in metadata
class Register:
    def __init__ (self, inputs : Optional[dict] = None, file_path : str = REGISTER_FILENAME):
        # Save the previous register
        self.file = File(file_path)
        # Save the current date
        self.date = datetime.today().strftime('%d-%m-%Y %H:%M:%S')
        # Save the inputs, if passed
        self.inputs = inputs
        # Set a cache for some already calculated values
        self.cache = {}
        # Set the tests tracker
        self.tests = {}
        # Set the warnings list, which will be filled by failing tests
        self.warnings = []
        # Inherit cache and test results from the register last entry
        # Note that warnings are not inherited from the last entry
        # Falied tests are to be repeated and they raise the warnings
        self.entries = []
        self.last_entry = None
        if exists(self.file.path):
            # Read the register in disk
            with open(self.file.path, 'r') as file:
                self.entries = json.load(file)
                self.last_entry = self.entries[-1]
            # Inherit the cache
            for field_name, field_value in self.last_entry['cache'].items():
                self.cache[field_name] = field_value
            # Inherit test results
            for test_name, test_result in self.last_entry['tests'].items():
                self.tests[test_name] = test_result
        # Set the save function to bre unned when the process is over
        atexit.register(self.save)

    def __repr__ (self):
        return str(self.to_dict())

    def to_dict (self) -> dict:
        # Set a dictionary with the current values
        dictionary = {
            'inputs': self.inputs,
            'date': self.date,
            'cache': self.cache,
            'tests': self.tests,
            'warnings': self.warnings,
        }
        # Remove inputs if it is empty
        if not dictionary['inputs']:
            del dictionary['inputs']
        return dictionary

    # Save the register to a json file
    def save (self):
        # Set a new entry for the current run
        current_entry = self.to_dict()
        # Add the new entry to the list
        self.entries.append(current_entry)
        # Write entries to disk
        with open(self.file.path, 'w') as file:
            json.dump(self.entries, file, indent=4)