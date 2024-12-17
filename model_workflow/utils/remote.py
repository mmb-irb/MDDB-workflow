import urllib.request
import json
from typing import Optional

from model_workflow.utils.auxiliar import load_json, save_json
from model_workflow.utils.type_hints import *

class Remote:
    def __init__ (self, url : str):
        # Save arguments
        self.url = url
        # Set internal variables
        self._available_files = None

    # Get available files in the remove project
    def get_available_files (self):
        # Return the internal value if we already have it
        if self._available_files != None:
            return self._available_files
        # Otherwise request the available files to the API
        request_url = self.url + '/files'
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
            raise Exception('Something went wrong when downloading the input files: ' + request_url)
        # If this is a json file then rewrite the inputs file in a pretty formatted way (with indentation)
        if is_json:
            file_content = load_json(output_file.path)
            save_json(file_content, output_file.path, indent = 4)