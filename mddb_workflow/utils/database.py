import urllib.request
import urllib.error
import ssl
import json
from tqdm import tqdm
from mddb_workflow.utils.auxiliar import load_json, save_json, InputError, RemoteServiceError, warn
from mddb_workflow.utils.constants import INCOMPLETE_PREFIX, DEFAULT_INPUTS_FILENAME
from mddb_workflow.utils.type_hints import *

# When downloading files, set the chunk size in bytes
CHUNK_SIZE = 1024 * 1024  # 1 MB

# Create a system to skip SSL certificates authentication
NO_SSL_CONTEXT = ssl.create_default_context()
NO_SSL_CONTEXT.check_hostname = False
NO_SSL_CONTEXT.verify_mode = ssl.CERT_NONE

# Header for workflow request source
WORKFLOW_REQUEST_SOURCE_HEADER = 'x-mddb-request-source'
WORKFLOW_REQUEST_SOURCE_VALUE = 'workflow'


def workflow_urlopen(url, *args, **kwargs):
    """Set a workflow header for metrics urllib.request.urlopen."""
    req = urllib.request.Request(url)
    req.add_header(WORKFLOW_REQUEST_SOURCE_HEADER, WORKFLOW_REQUEST_SOURCE_VALUE)
    return urllib.request.urlopen(req, *args, **kwargs)


class Remote:
    """Class to handle remote projects in the database."""
    def __init__(self, database: 'Database', accession: str, context=None):
        """Initialize the remote project handler.

        Args:
            database: Database handler to access the database.
            accession: Accession of the remote project to be handled.
            context: SSL context to be used in the requests.

        """
        # Set the URL
        self.database = database
        self.accession = accession
        self.project_url = f'{self.database.url}rest/current/projects/{accession}'
        # Set the context
        self.context = context
        # Set internal variables
        self._project_data = None
        self._available_files = None
        # Download project data to make sure we have database access and the project exists
        self.get_project_data()

    def __str__(self) -> str:
        return f'< Remote {self.project_url} >'

    def get_project_data(self) -> dict:
        """Get project data.
        This is only used to make sure the project exists by now.
        """
        # Return the internal value if we already have it
        if self._project_data != None:
            return self._project_data
        # Make sure the database is alive (and thus the provided database URL is valid)
        if not self.database.is_alive():
            raise RemoteServiceError('Database not available')
        # Otherwise request the project data to the API
        try:
            response = workflow_urlopen(self.project_url, context=self.context)
            self._project_data = json.loads(response.read())
            return self._project_data
        except urllib.error.HTTPError as error:
            # Try to provide comprehensive error logs depending on the error
            # If project was not found
            if error.code == 404:
                raise InputError(f'Remote project "{self.accession}" not found in {self.project_url}')
            # If we don't know the error then simply say something went wrong
            raise Exception(f'Error when downloading project data: {self.project_url} with error: {error}')
        except urllib.error.URLError as error:
            # If we don't know the error then simply say something went wrong
            raise Exception(f'Error when downloading project data: {self.project_url} with error: {error}')
        except Exception as error:
            raise Exception(f'Something went wrong when requesting project data: {self.project_url} with error: {error}')
    project_data = property(get_project_data, None, None, "Project data (read only)")

    def get_snaphsots(self):
        """Get the number of snapshots in the remote trajectory."""
        return self._project_data['metadata']['mdFrames']
    snapshots = property(get_snaphsots, None, None, "Number of snapshots in the remote trajectory (read only)")

    def _get_file(self, request_url: str, description: str, output_file: 'File' = None):
        """Get the content of a request and optionally save it to a file."""
        file_str = '' if output_file is None else f' ({output_file.path})'
        print(f'Downloading {description}{file_str}\n')
        try:
            response = workflow_urlopen(request_url, context=self.context)
            if output_file:
                with open(output_file.path, 'wb') as file:
                    file.write(response.read())
            else:
                content = json.loads(response.read())
                return content
        except Exception as error:
            raise Exception(f'Something went wrong when downloading {description}: {request_url} with error: {error}')

    def get_available_files(self):
        """Get the available files in the remote project."""
        # Return the internal value if we already have it
        if self._available_files != None:
            return self._available_files
        # Otherwise request the available files to the API
        request_url = self.project_url + '/files'
        self._available_files = self._get_file(request_url, 'available files')
        return self._available_files
    available_files = property(get_available_files, None, None, "Remote available files (read only)")

    def download_file(self, target_filename: str, output_file: 'File'):
        """Download a specific file from the project/files endpoint."""
        request_url = f'{self.project_url}/files/{target_filename}'
        print(f'Downloading file "{target_filename}" in {output_file.path}\n')
        try:
            response = workflow_urlopen(request_url, context=self.context)
            with open(output_file.path, 'wb') as file:
                while True:
                    chunk = response.read(CHUNK_SIZE)
                    if not chunk: break
                    file.write(chunk)
        except urllib.error.HTTPError as error:
            if error.code == 404:
                raise Exception(f'Missing remote file "{target_filename}"')
            # If we don't know the error then simply say something went wrong
            raise Exception(f'Something went wrong when downloading file "{target_filename}" in {request_url} with error: {error}')

    def get_standard_topology(self):
        """Get the standard topology of the project."""
        # Download the project standard topology
        request_url = self.project_url + '/topology'
        return self._get_file(request_url, 'standard topology')

    def download_standard_topology(self, output_file: 'File'):
        """Download the standard topology of the project."""
        request_url = self.project_url + '/topology'
        self._get_file(request_url, 'standard topology', output_file)

    def download_standard_structure(self, output_file: 'File'):
        """Download the standard structure of the project."""
        request_url = self.project_url + '/structure'
        self._get_file(request_url, 'standard structure', output_file)

    def download_trajectory(self,
        output_file: 'File',
        frame_selection: Optional[str] = None,
        atom_selection: Optional[str] = None,
        format: Optional[str] = None
    ):
        """Download the main trajectory."""
        if [frame_selection, atom_selection, format] == [None, None, 'xtc']:
            # If we dont have a specific request, we can download the main trajectory
            # directly from the trajectory.xtc file so it is faster
            request_url = f'{self.project_url}/files/trajectory.xtc'
        else:
            # Set the base URL
            request_url = self.project_url + '/trajectory'
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
        print(f'Downloading main trajectory ({output_file.path})')
        # Create a temporal file to download the trajectory
        # Thus if the download is interrupted we will know the trajectory is incomplete
        incomplete_trajectory = output_file.get_prefixed_file(INCOMPLETE_PREFIX)
        # If we have a previous incomplete trajectory then remove it
        if incomplete_trajectory.exists: incomplete_trajectory.remove()
        try:
            response = workflow_urlopen(request_url, context=self.context)
            pbar = tqdm(unit='B', unit_scale=True, unit_divisor=1024,
                        miniters=1, desc=' Progress', leave=False)
            with open(incomplete_trajectory.path, 'wb') as file:
                while True:
                    chunk = response.read(CHUNK_SIZE)
                    if chunk: pbar.update(len(chunk))
                    else: break
                    file.write(chunk)
        except Exception as error:
            raise Exception(f'Something went wrong when downloading the main trajectory: {request_url} with error: {error}')
        # Once the trajectory is fully downloaded we change its filename
        incomplete_trajectory.rename_to(output_file)

    def get_inputs_file_source_filename(self) -> str:
        """Get the inputs file name with the source format."""
        return f'source_{self.database.alias}_{self.accession}_{DEFAULT_INPUTS_FILENAME}'

    def download_inputs_file(self, output_file: 'File'):
        """Download the inputs file."""
        request_url = self.project_url + '/inputs'
        # In case this is a json file we must specify the format in the query
        is_json = output_file.format == 'json'
        if is_json:
            request_url += '?format=json'
        self._get_file(request_url, 'inputs file', output_file)
        if is_json:
            file_content = load_json(output_file.path)
            save_json(file_content, output_file.path, indent=4)

    def download_analysis_data(self, analysis_type: str, output_file: 'File'):
        """Download analysis data."""
        request_url = f'{self.project_url}/analyses/{analysis_type}'
        self._get_file(request_url, f'{analysis_type} analysis data', output_file)
        # Format JSON if needed
        try:
            file_content = load_json(output_file.path)
            save_json(file_content, output_file.path, indent=4)
        except Exception:
            pass  # If not JSON, skip formatting


class Database:
    """Class to handle database operations."""
    def __init__(self, url: str, no_ssl_authentication: bool = False):
        """Initialize the database handler.

        Args:
            url: URL of the database (e.g. https://irb-dev.mddbr.eu/api/)
            no_ssl_authentication: If True, SSL certificates will not be authenticated.

        """
        if '://' not in url:
            raise InputError(f'Invalid database URL "{url}"')
        self.url = url
        # If the URL already includes /rest/... then clean this part away
        if '/rest' in self.url:
            self.url = self.url.split('/rest')[0] + '/'
        # Set an alias for this database
        self.alias = self.url.split('://')[1].split('.')[0].split('-')[0]
        # Set the context
        self.context = NO_SSL_CONTEXT if no_ssl_authentication else None

    def __str__(self) -> str:
        return f'< Database {self.url} >'

    def is_alive(self) -> bool:
        """Check if the database is alive.

        WARNING: Note that this function requires internet connection.
        WARNING: Do not run in by default.
        """
        try:
            response = workflow_urlopen(self.url, context=self.context)
            response.read(1)
            return True
        except urllib.error.HTTPError as error:
            # Server error
            if error.code == 503:
                warn('MDDB Service unavailable. Please try again later.')
                return False
            # Unknown HTTP error
            return False
        except urllib.error.URLError as error:
            # SSL error
            if 'SSL: CERTIFICATE_VERIFY_FAILED' in str(error):
                raise RemoteServiceError(f'Failed to verify SSL certificate from {self.url}\n' +
                    ' Use the "--ssleep" flag to avoid SSL authentication if you trust the source.')
            # Timeout error
            # The error variable as is do not work properly
            # Its 'errno' value is None (at least for tiemout errors)
            actual_error = error.args[0]
            if actual_error.errno == 110:
                warn('Timeout error when requesting MDposit. Is the node fallen?')
                return False
            # Unknown URL error
            return False
        except:
            # Unknown error
            return False

    def get_remote_project(self, accession: str) -> Remote:
        """Instantiate the remote project handler."""
        return Remote(self, accession, context=self.context)

    def get_reference_data(self, reference: str, id: str) -> Optional[dict]:
        """Check if the required sequence is already in the MDDB database."""
        # Make sure the database is alive (and thus the provided database URL is valid)
        # If not then we return None and allow the workflow to keep going
        # Probably if the reference can not be obtained the workflow will generate it again
        if not self.is_alive(): return None
        # Request the specific data
        request_url = f'{self.url}rest/v1/references/{reference}/{id}'
        try:
            with workflow_urlopen(request_url, context=self.context) as response:
                return json.loads(response.read().decode("utf-8", errors='ignore'))
        # Handle possible errors
        except urllib.error.HTTPError as error:
            # If the reference is not found in MDposit
            if error.code == 404: return None
            warn(f'Error when requesting MDposit: {request_url}')
            raise RuntimeError(f'Something went wrong with the MDposit request {request_url}')
        except Exception as error:
            raise RuntimeError(f'Something went wrong with the MDposit request {request_url} with error: {error}')
