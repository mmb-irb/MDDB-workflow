from sys import argv
from os.path import exists
from datetime import datetime
from typing import Optional, List
import atexit

from model_workflow.utils.auxiliar import load_json, save_json
from model_workflow.utils.constants import REGISTER_FILENAME, YELLOW_HEADER, COLOR_END, AVAILABLE_CHECKINGS
from model_workflow.utils.file import File

# The register tracks activity along multiple runs and thus avoids repeating some already succeeded tests
# It is also responsible for storing test failure warnings to be written in metadata
class Register:
    def __init__ (self, file_path : str = REGISTER_FILENAME):
        # Save the previous register
        self.file = File(file_path)
        # Save the current workflow call
        self.call = ' '.join(argv)
        # Save the current date
        self.date = datetime.today().strftime('%d-%m-%Y %H:%M:%S')
        # Set a cache for some already calculated values
        self.cache = {}
        # Set the tests tracker
        self.tests = { checking: None for checking in AVAILABLE_CHECKINGS }
        # Set the warnings list, which will be filled by failing tests
        self.warnings = []
        # Inherit cache, test results and warning from the register last entry
        self.entries = []
        self.last_entry = None
        if exists(self.file.path):
            # Read the register in disk
            self.entries = load_json(self.file.path)
            self.last_entry = self.entries[-1]
            # Inherit the cache
            for field_name, field_value in self.last_entry['cache'].items():
                self.cache[field_name] = field_value
            # Inherit test results
            for test_name, test_result in self.last_entry['tests'].items():
                self.tests[test_name] = test_result
            # Inherit warnings
            for warning in self.last_entry['warnings']:
                # DANI: Para quitarnos de encima warnings con el formato antiguo
                if not warning.get('tag', None):
                    continue
                self.warnings.append(warning)
        # Set the save function to be runned when the process is over
        atexit.register(self.save)

    def __repr__ (self):
        return str(self.to_dict())

    def to_dict (self) -> dict:
        # Set a dictionary with the current values
        dictionary = {
            'call': self.call,
            'date': self.date,
            'cache': self.cache,
            'tests': self.tests,
            'warnings': self.warnings,
        }
        return dictionary

    # Get current warnings filtered by tag
    def get_warnings (self, tag : str) -> List[dict]:
        return [ warning for warning in self.warnings if warning['tag'] == tag ]

    # Add warnings with the right format
    # A flag is to be passed to handle further removal of warnings
    def add_warning (self, tag : str, message : str):
        print(YELLOW_HEADER + 'WARNING: ' + COLOR_END + message)
        warning = { 'tag': tag, 'message': message }
        self.warnings.append(warning)

    # Remove warnings filtered by tag
    def remove_warnings (self, tag : str):
        self.warnings = [ warning for warning in self.warnings if warning['tag'] != tag ]

    # Save the register to a json file
    def save (self):
        # If path does not exist then do nothing
        # WARNING: I know this looks a bit silent 
        # WARNING: Otherwise it is a constant spam when something goes wrong close to beginning
        if not exists(self.file.basepath):
            return
        # Set a new entry for the current run
        current_entry = self.to_dict()
        # Add the new entry to the list
        self.entries.append(current_entry)
        # Write entries to disk
        save_json(self.entries, self.file.path, indent = 4)