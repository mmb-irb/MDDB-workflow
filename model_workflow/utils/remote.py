import urllib.request
import json

from model_workflow.utils.auxiliar import load_json, save_json, InputError
from model_workflow.utils.type_hints import *

class Remote:
    def __init__ (self, database_url : str, accession : str):
        # Save input arguments
        self.database_url = database_url
        self.accession = accession
        # Set the URL
        self.url = f'{database_url}/rest/current/projects/{accession}'
        # Set internal variables
        self._project_data = None
        self._available_files = None
        # Download project data to make sure we have database access and the project exists
        self.get_project_data()

    # Get project data
    # This is only used to make sure the project exists by now
    def get_project_data (self):
        # Return the internal value if we already have it
        if self._project_data != None:
            return self._project_data
        # Otherwise request the project data to the API
        try:
            response = urllib.request.urlopen(self.url)
            self._project_data = json.loads(response.read())
        except urllib.error.HTTPError as error:
            # Try to provide comprehensive error logs depending on the error
            # If project was not found
            if error.code == 404:
                raise InputError(f'Remote project "{self.accession}" not found')
            # If we don't know the error then simply say something went wrong
            raise Exception('Error when downloading project data: ' + self.url)
        except:
            raise Exception('Something went wrong when requesting project data: ' + self.url)

    # Get available files in the remove project
    def get_available_files (self):
        # Return the internal value if we already have it
        if self._available_files != None:
            return self._available_files
        # Otherwise request the available files to the API
        request_url = self.url + '/files'
        print('GETTIG REMOTE FILES')
        try:
            response = urllib.request.urlopen(request_url)
            self._available_files = json.loads(response.read())
        except:
            raise Exception('Something went wrong when requesting available files: ' + request_url)
        return self._available_files
    available_files = property(get_available_files, None, None, "Remote available files (read only)")

    # Download a specific file
    def download_file (self, output_file : 'File'):
        request_url = f'{self.url}/files/{output_file.filename}'
        print(f'Downloading file "{output_file.filename}" ({output_file.path})\n')
        try:
            urllib.request.urlretrieve(request_url, output_file.path)
        except urllib.error.HTTPError as error:
            # Try to provide comprehensive error logs depending on the error
            # If file was not found
            if error.code == 404:
                raise Exception(f'Missing remote file "{output_file.filename}"')
            # If we don't know the error then simply say something went wrong
            raise Exception(f'Something went wrong when downloading file "{output_file.filename}": ' + request_url)

    # Download the project standard topology
    def download_standard_topology (self, output_file : 'File'):
        request_url = self.url + '/topology'
        print(f'Downloading standard topology ({output_file.path})\n')
        try:
            urllib.request.urlretrieve(request_url, output_file.path)
        except:
            raise Exception('Something went wrong when downloading the standard topology: ' + request_url)
        
    # Download the standard structure
    def download_standard_structure (self, output_file : 'File'):
        request_url = self.url + '/structure'
        print(f'Downloading standard structure ({output_file.path})\n')
        try:
            urllib.request.urlretrieve(request_url, output_file.path)
        except:
            raise Exception('Something went wrong when downloading the standard structure: ' + request_url)
        
    # Download the main trajectory
    def download_trajectory (self,
        output_file : 'File',
        frame_selection : Optional[str] = None,
        atom_selection : Optional[str] = None,
        format : Optional[str] = None
    ):
        # Set the base URL
        request_url = self.url + '/trajectory'
        # Additional arguments to be included in the URL
        arguments = []
        if frame_selection:
            arguments.append(f'frames={frame_selection}')
        if atom_selection:
            arguments.append(f'atoms={atom_selection}')
        if format:
            arguments.append(f'format={format}')
        if len(arguments) > 0:
            request_url += '?' + '&'.join(arguments)
        # Send the request
        print(f'Downloading main trajectory ({output_file.path})\n')
        try:
            urllib.request.urlretrieve(request_url, output_file.path)
        except:
            raise Exception('Something went wrong when downloading the main trajectory: ' + request_url)
        
    # Download the inputs file
    def download_inputs_file (self, output_file : 'File'):
        request_url = self.url + '/inputs'
        # In case this is a json file we must specify the format in the query
        is_json = output_file.format == 'json'
        if is_json:
            request_url += '?format=json'
        # Send the request
        print(f'Downloading input files ({output_file.path})\n')
        try:
            urllib.request.urlretrieve(request_url, output_file.path)
        except:
            raise Exception('Something went wrong when downloading the inputs file: ' + request_url)
        # If this is a json file then rewrite the inputs file in a pretty formatted way (with indentation)
        if is_json:
            file_content = load_json(output_file.path)
            save_json(file_content, output_file.path, indent = 4)