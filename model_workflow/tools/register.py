from os.path import exists
from datetime import datetime
import json

from model_workflow.constants import REGISTER_FILENAME

REGISTER_INPUTS = [
    'directory',
    'accession',
    'database_url',
    'inputs_filepath',
    'input_topology_filepath',
    'input_structure_filepath',
    'input_trajectory_filepaths',
    'populations_filepath',
    'transitions_filepath',
    'md_directories',
    'reference_md_directory',
    'filter_selection',
    'image',
    'fit',
    'translation',
    'mercy',
    'trust',
    'pca_selection',
    'pca_fit_selection',
    'rmsd_cutoff',
    'interaction_cutoff',
    'sample_trajectory',
]

# The register tracks activity along multiple runs and thus avoids repeating some already succeeded tests
# It is also responsible for storing test failure warnings to be written in metadata
class Register:
    def __init__ (self, current_project : 'Project', file_path : str = REGISTER_FILENAME):
        self.current_project = current_project
        self.file_path = file_path
        # Load previous entries if any
        self.entries = []
        if exists(self.file_path):
            with open(self.file_path, 'r') as file:
                self.entries = json.load(file)
        # Set the current run tracked values
        self.date = datetime.today().strftime('%d-%m-%Y %H:%M:%S')
        self.inputs = {}
        for register_input in REGISTER_INPUTS:
            self.inputs[register_input] = getattr(current_project, register_input)
        # Set the tests tracker
        self.tests = {}
        # Inherit test results from the last entry
        self.last_entry = None
        if len(self.entries) > 0:
            self.last_entry = self.entries[-1]
            for test_name, test_result in self.last_entry['tests'].items():
                self.tests[test_name] = test_result
        # Set the warnings list, which will be filled by failing tests
        # Note that warnings are not inherited from the previous entry since falied tests are to be repeated
        self.warnings = []
        # Set subsections for individual MDs to have their own register
        # MD registers have tests and warnings
        self.mds = {}
        for md_directory in self.current_project.md_directories:
            # Set the basic MD register
            md_register = {
                'tests': [],
                'warnings': []
            }
            # Inherit test results from the last entry
            if self.last_entry:
                last_md_register = self.last_entry['mds'][md_directory]
                for test_name, test_result in last_md_register['tests'].items():
                    md_register['tests'][test_name] = test_result
            # Save the current MD register using the MD directory as key
            self.mds[md_directory] = md_register


    def save (self):
        # Set a new entry for the current run
        current_entry = {
            'date': self.date,
            'inputs': self.inputs,
            'tests': self.tests,
            'warnings': self.warnings,
            'mds': self.mds
        }
        # Add the new entry to the list
        self.entries.append(current_entry)
        # Write entries to disk
        with open(self.file_path, 'w') as file:
            json.dump(self.entries, file, indent=4)