from sys import argv
from os.path import exists
from datetime import datetime

from model_workflow.utils.auxiliar import load_json, save_json, warn, get_git_version
from model_workflow.utils.constants import DATE_STYLE
from model_workflow.utils.type_hints import *


class Register:
    """ The register tracks activity along multiple runs and thus avoids repeating some already succeeded tests
    It is also responsible for storing test failure warnings to be written in metadata. """
    def __init__ (self, register_file : 'File'):
        # Save the previous register
        self.file = register_file
        # Save the current workflow call
        # Quote those arguments including space, since it means they were quoted when inputed
        quoted_argv = [ f"'{arg}'" if ' ' in arg else arg for arg in argv ]
        self.call = ' '.join(quoted_argv)
        # Save the current date
        self.date = datetime.today().strftime(DATE_STYLE)
        # Save the current version
        self.version = get_git_version()
        # Set the tests tracker
        self.tests = {}
        # Set the warnings list, which will be filled by failing tests
        self.warnings = []
        # Inherit test results and warning from the register last entry
        self.entries = []
        if self.file.exists:
            # Read the register in disk
            self.entries = load_json(self.file.path)
            last_entry = self.entries[-1]
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
            'version': self.version,
            'tests': self.tests,
            'warnings': self.warnings,
        }
        return dictionary

    def update_test (self, key : str, value : Optional[bool]):
        """ Update a test result and save the register. """
        self.tests[key] = value
        self.save()

    def get_warnings (self, tag : str) -> list[dict]:
        """ Get current warnings filtered by tag. """
        return [ warning for warning in self.warnings if warning['tag'] == tag ]

    def add_warning (self, tag : str, message : str):
        """ Add warnings with the right format and save the register. 
        A flag is to be passed to handle further removal of warnings. """
        warn(message)
        # If we already had this exact warning then do not repeat it
        for warning in self.warnings:
            if warning['tag'] == tag and warning['message'] == message : return
        # Add a new warning
        warning = { 'tag': tag, 'message': message }
        self.warnings.append(warning)
        self.save()

    def remove_warnings (self, tag : str):
        """ Remove warnings filtered by tag and save the register. """
        self.warnings = [ warning for warning in self.warnings if warning['tag'] != tag ]
        self.save()

    def save (self):
        """ Save the register to a json file. """
        # If path does not exist then do nothing
        # WARNING: I know this looks a bit silent 
        # WARNING: Otherwise it is a constant spam when something goes wrong close to beginning
        if not exists(self.file.basepath):
            return
        # Set a new entry for the current run
        current_entry = self.to_dict()
        # Write entries to disk
        save_json(self.entries + [ current_entry ], self.file.path, indent = 4)