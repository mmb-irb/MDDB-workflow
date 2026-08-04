from os import environ
from os.path import exists, normpath
from shutil import which
from subprocess import run, Popen, PIPE
from select import select
from socket import create_connection
from threading import Event, Thread

from mddb_workflow.utils.auxiliar import load_yaml, InputError, RemoteServiceError, MISSING_VALUE
from mddb_workflow.utils.constants import LOADER_NODES_CONFIG_FILEPATH

# Set the time we wait for an SSH connection attempt before we give up
SSH_TIMEOUT = 30 # seconds

# Se the nodeJS command
NODEJS_COMMAND = 'node'

# The contents of the loader nodes config file must have the follwoing structure
# < node alias >:     < as in https://mdposit.mddbr.eu/api/rest/v1/nodes >
#   ssh_hostname:     < as in your own ~/.ssh/config file >
#   db_server:        < as in the loader's .env file of the corresponding node >
#   db_port:          < as in the loader's .env file of the corresponding node >
#   db_name:          < as in the loader's .env file of the corresponding node >
#   db_auth_user:     < as in the loader's .env file of the corresponding node >
#   db_auth_password: < as in the loader's .env file of the corresponding node >
#   db_auth_source:   < as in the loader's .env file of the corresponding node >

# Set which fields are actually loader environmental variables
LOADER_ENV_FIELDS = [
    'db_server',
    'db_port',
    'db_name',
    'db_auth_user',
    'db_auth_password',
    'db_authsource'
]

# Set all expected fields in the configuration file
CONFIG_ENV_FIELDS = [ 'ssh_hostname', *LOADER_ENV_FIELDS ]

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
        # Import all config variables and store them as instance variables as well
        for field in CONFIG_ENV_FIELDS:
            field_value = self.config.get(field, MISSING_VALUE)
            # If any of the expected fields is missing then complain
            # Note that a field may be 'None' and not missing
            if field_value is MISSING_VALUE:
                raise InputError(f'Missing "{field}" for node "{self.target_node}" in {LOADER_NODES_CONFIG_FILEPATH}')
            setattr(self, field, field_value)
        # Set an internal variable to store the process holding the SSH connection
        self.ssh_process = None
        # Open the SSH connection
        if not self.open_ssh_connection():
            raise RemoteServiceError(f'Failed to open the SSH connection for node "{target_node}".\n'
                'Please make sure you have the SSH connection properly configured in your ~/.ssh/config')
        # Make sure the connection works
        if not self.has_database_access():
            raise RemoteServiceError(f'The loader failed to access the database in node "{target_node}".\n'
                f'Please check the parameters are correct in the config file {LOADER_NODES_CONFIG_FILEPATH}\n'
                'Also make sure you have the SSH tunnel properly configured in your ~/.ssh/config\n'
                'There should be a line similar to "LocalForward <port> 127.0.0.1:27017"')

    # Make the SSH connection so we can access the node database
    # WARNING: Note that this function will only open a SSH connection
    # WARNING: the tunnel configuration must be in your ~/.ssh/config
    def open_ssh_connection(self) -> bool:
        # Run the SSH command and keep it alive
        print(f'Opening SSH connection to "{self.ssh_hostname}" for the loader...')
        # WARNING: stdin=PIPE prevents the terminal from becoming a zombi
        # The argument "-T" prevents the following error log:
        #   Pseudo-terminal will not be allocated because stdin is not a terminal.
        self.ssh_process = Popen([ "ssh", "-T", self.ssh_hostname ],
            stdin=PIPE, stdout=PIPE, stderr=PIPE)
        # Set a race between two independent processes to check if the SSH connection is open
        finished = Event()
        # WARNING: This has to be a list por the threads to edit it
        result = [None]
        # In one hand we check if the process returns output or error
        def check_ssh_response():
            # Watch for any writes in both stdout and stderr
            # The waiting stops as soon as any of the streams has content
            responses, _, _ = select([self.ssh_process.stdout, self.ssh_process.stderr], [], [], SSH_TIMEOUT)
            # If the other racing process finished already then stop here
            if finished.is_set(): return
            # If there are no responses then
            if len(responses) == 0: return
            # If the first response is output then we assume the connection is alive
            first_response = responses[0]
            if first_response is self.ssh_process.stdout:
                print('  The SSH connection is alive')
                result[0] = True
            # If the first response is error then we assume something went wrong
            elif first_response is self.ssh_process.stderr:
                error_message = first_response.read1().decode()
                print('  Something went wrong')
                print(error_message)
                result[0] = False
            # Set the race finished since this process got some response
            finished.set()
        # In the other hand we shcek if the SSH tunnel is open by trying to connect to the target port
        # Note that this is necessary for some hosts which do not return output on connect
        def check_port():
            # Do this while the other racing process did not finish yet
            while not finished.is_set():
                try:
                    with create_connection((self.db_server, self.db_port), timeout=1):
                        print('  The SSH tunnel is open')
                        result[0] = True
                        finished.set()
                        return
                except (ConnectionRefusedError, OSError):
                    finished.wait(0.5)
        # Run the two checks in paralel
        Thread(target=check_ssh_response, daemon=True).start()
        Thread(target=check_port, daemon=True).start()
        # Wait for any of them to finish
        finished.wait(timeout=SSH_TIMEOUT)
        # handle the results
        if result[0] is True: return True
        if result[0] is False: return False
        if result[0] is None:
            print('  Exceeded timeout')
            return False

    # Set the loader environment variables before calling the loader
    # Note that at this point we have already checked that the configuration fields exist
    def set_environment(self):
        for field in LOADER_ENV_FIELDS:
            caps_field = field.upper()
            field_value = getattr(self, field)
            if field_value is None: field_value = ""
            environ[caps_field] = field_value

    # Make sure the loader is connected and has access to the database
    def has_database_access(self) -> bool:
        print('Checking the loader has database access...')
        # Set the environmental variables before calling the loader
        self.set_environment()
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
        # Set the environmental variables before calling the loader
        self.set_environment()
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
        