from sys import argv
from os.path import exists, getmtime
from datetime import datetime
from time import strftime, gmtime

from model_workflow.utils.auxiliar import load_json, save_json, warn
from model_workflow.utils.type_hints import *

# Set dates format
date_style = '%d-%m-%Y %H:%M:%S'

# The register tracks activity along multiple runs and thus avoids repeating some already succeeded tests
# It is also responsible for storing test failure warnings to be written in metadata
class Register:
    def __init__ (self, register_file : 'File'):
        # Save the previous register
        self.file = register_file
        # Save the current workflow call
        # Quote those arguments including space, since it means they were quoted when inputed
        quoted_argv = [ f"'{arg}'" if ' ' in arg else arg for arg in argv ]
        self.call = ' '.join(quoted_argv)
        # Save the current date
        self.date = datetime.today().strftime(date_style)
        # Set record for the modification times of processed input files
        # This allows to know if any of those files have been modified and thus we must reset some register fields
        self.mtimes = {}
        # Set a cache for some already calculated values
        self.cache = {}
        # Set the tests tracker
        self.tests = {}
        # Set the warnings list, which will be filled by failing tests
        self.warnings = []
        # Inherit cache, test results and warning from the register last entry
        self.entries = []
        if self.file.exists:
            # Read the register in disk
            self.entries = load_json(self.file.path)
            last_entry = self.entries[-1]
            # Inherit modification times
            self.mtimes = last_entry.get('mtimes', {})
            # Inherit the cache
            for field_name, field_value in last_entry['cache'].items():
                self.cache[field_name] = field_value
            # Inherit test results
            for test_name, test_result in last_entry['tests'].items():
                self.tests[test_name] = test_result
            # Inherit warnings
            for warning in last_entry['warnings']:
                # DANI: Para quitarnos de encima warnings con el formato antiguo
                if not warning.get('tag', None):
                    continue
                self.warnings.append(warning)
        # Save the entry for the first time
        self.save()

    def __repr__ (self):
        return str(self.to_dict())

    def to_dict (self) -> dict:
        # Set a dictionary with the current values
        dictionary = {
            'call': self.call,
            'date': self.date,
            'mtimes': self.mtimes,
            'cache': self.cache,
            'tests': self.tests,
            'warnings': self.warnings,
        }
        return dictionary

    # Get new and previous mtimes of a target file
    def get_mtime (self, target_file : 'File') -> tuple:
        new_raw_mtime = getmtime(target_file.path)
        new_mtime = strftime(date_style, gmtime(new_raw_mtime))
        previous_mtime = self.mtimes.get(target_file.filename, None)
        return new_mtime, previous_mtime

    # Update a modification time
    def update_mtime (self, target_file : 'File'):
        # Get the new and the previous value
        new_mtime, previous_mtime = self.get_mtime(target_file)
        # If the new value is already the previous value then do nothing
        if new_mtime == previous_mtime:
            return
        # Overwrite previous value and save the register
        self.mtimes[target_file.filename] = new_mtime
        self.save()

    # Check if a file is new
    def is_file_new (self, target_file : 'File') -> bool:
        # Get the new and the previous value
        new_mtime, previous_mtime = self.get_mtime(target_file)
        # If the previous value is None then it is new
        return previous_mtime == None

    # Check if a file does not match the already registered modification time or it is new
    def is_file_modified (self, target_file : 'File') -> bool:
        # Get the new and the previous value
        new_mtime, previous_mtime = self.get_mtime(target_file)
        # If the new value is already the previous value then it has not been modified
        if new_mtime == previous_mtime:
            return False
        return True

    # Update the cache and save the register
    def update_cache (self, key : str, value):
        self.cache[key] = value
        self.save()

    # Reset the cache
    # This is called when some input files are modified
    def reset_cache (self):
        self.cache = {}
        self.save()

    # Update a test result and save the register
    def update_test (self, key : str, value : Optional[bool]):
        self.tests[key] = value
        self.save()

    # Get current warnings filtered by tag
    def get_warnings (self, tag : str) -> List[dict]:
        return [ warning for warning in self.warnings if warning['tag'] == tag ]

    # Add warnings with the right format and save the register
    # A flag is to be passed to handle further removal of warnings
    def add_warning (self, tag : str, message : str):
        warn(message)
        # If we already had this exact warning then do not repeat it
        for warning in self.warnings:
            if warning['tag'] == tag and warning['message'] == message : return
        # Add a new warning
        warning = { 'tag': tag, 'message': message }
        self.warnings.append(warning)
        self.save()

    # Remove warnings filtered by tag and save the register
    def remove_warnings (self, tag : str):
        self.warnings = [ warning for warning in self.warnings if warning['tag'] != tag ]
        self.save()

    # Save the register to a json file
    def save (self):
        # If path does not exist then do nothing
        # WARNING: I know this looks a bit silent 
        # WARNING: Otherwise it is a constant spam when something goes wrong close to beginning
        if not exists(self.file.basepath):
            return
        # Set a new entry for the current run
        current_entry = self.to_dict()
        # Write entries to disk
        save_json(self.entries + [ current_entry ], self.file.path, indent = 4)