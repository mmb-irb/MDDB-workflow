from os.path import exists, normpath
from shutil import which
from subprocess import run, Popen, PIPE
from select import select

from mddb_workflow.utils.auxiliar import load_yaml, InputError, RemoteServiceError
from mddb_workflow.utils.constants import LOADER_NODES_CONFIG_FILEPATH

# Se the nodeJS command
NODEJS_COMMAND = 'node'

# Tool to handle the the MDDB loader from the MDDB workflow
class Loader:
    
    def __init__(self, loader_directory : str, target_node : str):
        # Save input parameter(s)
        self.loader_directory = normpath(loader_directory)
        self.target_node = target_node
        # Make sure nodeJS is installed
        nodejs_executable = which(NODEJS_COMMAND)
        if not nodejs_executable:
            raise EnvironmentError('Trying to use the loader while Node JS is not installed.')
        # Check the loader nodes config file exists
        if not exists(LOADER_NODES_CONFIG_FILEPATH):
            raise InputError(f'Trying to use the loader while missing the loader nodes config file: {LOADER_NODES_CONFIG_FILEPATH}')
        # Make sure the loader directory looks correct by checking there is the index.js file
        self.loader_index = f'{self.loader_directory}/index.js'
        if not exists(self.loader_index):
            raise InputError(f'Loader directory {loader_directory} seems wrong. There is no index.js file in it.')
        # Load the config file
        nodes_config = load_yaml(LOADER_NODES_CONFIG_FILEPATH)
        # Get the config for the target node
        self.config = nodes_config.get(target_node, None)
        if self.config is None:
            raise InputError(f'There is no loader configuration for node "{target_node}". '
                f'Please either use one of the available nodes ({", ".join(nodes_config.keys())}) '
                f'or include your node in the loader config file {LOADER_NODES_CONFIG_FILEPATH}')
        # Set an internal variable to store the process holding the SSH connection
        self.ssh_process = None
        # Open the SSH connection
        if not self.open_ssh_connection():
            raise RemoteServiceError(f'Failed to open the SSH connection for node "{target_node}".\n'
                'Please make sure you have the SSH connection properly configured in your ~/.ssh/config')
        # Make sure the connection works
        if not self.has_database_access():
            raise RemoteServiceError(f'The loader failed to connect to node "{target_node}".\n'
                f'Please check the parameters are correct in the config file {LOADER_NODES_CONFIG_FILEPATH}\n'
                'Also make sure you have the SSH tunnel properly configured in your ~/.ssh/config\n'
                'There should be a line similar to "LocalForward <port> 127.0.0.1:27017"')

    # Make the SSH connection so we can access the node database
    # WARNING: Note that this function will only open a SSH connection
    # WARNING: the tunnel configuration must be in your ~/.ssh/config
    def open_ssh_connection(self) -> bool:
        # Get the hostname where we must aim at
        hostname = self.config.get('ssh_hostname', None)
        if hostname is None: raise InputError(f'Missing "ssh_hostname" for node {self.target_node} in {LOADER_NODES_CONFIG_FILEPATH}')
        # Run the SSH command and keep it alive
        print(f'Opening SSH connection to "{hostname}" for the loader...')
        # WARNING: stdin=PIPE prevents the terminal from becoming a zombi
        # The argument "-T" prevents the following error log:
        #   Pseudo-terminal will not be allocated because stdin is not a terminal.
        self.ssh_process = Popen([ "ssh", "-T", hostname ], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        # Wait for 30 seconds for the SSH response
        responses, _, _ = select([self.ssh_process.stdout, self.ssh_process.stderr], [], [], 30)
        if len(responses) == 0:
            print('  Exceeded timeout')
            return False
        # Get the first response we have
        first_response = responses[0]
        # If we have output then the SSH connection has worked
        if first_response is self.ssh_process.stdout:
            print('  The SSH connection is open')
            return True
        # If we have an error then show the problem
        if first_response is self.ssh_process.stderr:
            print('  Something went wrong')
            error_message = first_response.read1().decode()
            print(error_message)
            return False
        raise RuntimeError('Unexpected response in SSH connection')
            

    # Make sure the loader is connected and has access to the database
    def has_database_access(self) -> bool:
        print('Checking the loader has database access...')
        # Run the simplest and fastest loader command just to make sure it has database access
        check_process = run([NODEJS_COMMAND, self.loader_index, 'check'], stdout=PIPE, stderr=PIPE)
        check_output = check_process.stdout.decode().replace('\n','')
        if check_output:
            print(f'  Access confirmed -> The database has {check_output} projects')
            return True
        # If somthing goes wrong then report the problem
        print('  Something went wrong')
        check_errors = check_process.stderr.decode()
        print(check_errors)
        return False

    # Load a directory
    def load(self, target_directory: str, overwrite : bool = False) -> bool:
        print(f'Running loader to load directory {target_directory}')
        # We are not going to handle when the loader asks the user for input
        # Instead we will always use either the '-c' (conserve) or the '-o' (overwrite) arguments
        force_argument = '-o' if overwrite else '-c'
        # Run the load process
        load_process = run(
            [NODEJS_COMMAND, self.loader_index, 'load', target_directory, force_argument], 
            stderr=PIPE)
        # If there is any error then assume the load failed
        error_logs = load_process.stderr.decode()
        if error_logs:
            print('  Something went wrong')
            print(error_logs)
            return False
        return True
        