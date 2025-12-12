from os import rename, mkdir
from os.path import exists
from shutil import rmtree
from inspect import getfullargspec
import time

from mddb_workflow.utils.arg_cksum import get_cksum_id
from mddb_workflow.utils.auxiliar import safe_getattr
from mddb_workflow.utils.constants import *
from mddb_workflow.utils.file import File
from mddb_workflow.utils.type_hints import *

# Set a special exception for missing task function arguments
# This is used for easy debug when a new functions is added wrongly
MISSING_ARGUMENT_EXCEPTION = Exception('Missing argument')

# Set a special exception for when a value is missing
MISSING_VALUE_EXCEPTION = Exception('Missing value')

# Name of the argument used by all functions to know where to write output
OUTPUT_DIRECTORY_ARG = 'output_directory'


class Task:
    """Descriptor class to handle a generic task.
    It implements lazy properties, caching and overwriting.

    Since its properties are static, results are stored in the parent object
    (MD/Project), or otherwise all MDs would share the same task values.
    """

    def __init__(self,
        flag: str,
        name: str,
        func: Callable,
        args: dict = {},
        output_filenames: dict[str | Callable[[], str]] = {},
        use_cache: bool = True,
        debug: bool = False,
    ):
        """Initialize the Task object.

        Args:
            flag (str):
                The task flag.
                This name is used by the include/exclude/overwrite arguments and to name the analysis output directory.

            name (str):
                The task user-friendly name is to be used in the logs.

            func (Callable):
                The task function.
                Function argument names must correspond with Project/MD property names.

            args (dict, optional):
                The task function additional arguments.
                Project/MD properties are automatically sent to the function as arguments.
                However some analyses have additional arguments (e.g. frames limit, cutoffs, etc.)

            output_filenames (dict, optional):
                The task output filenames.
                Paths will be set automatically relative to its project/MD.
                For those tasks which generate a directory with multiple outputs this is not necessary.
                This field may be passed in addition so output files are written outside the directory.
                This may come in handy when these output files are used later in this workflow.

            use_cache (bool, optional):
                Set if the returned output is to be cached.
                Note that argument values are always cached, this is not optional.

            debug (bool, optional):
                If the task is run in debug mode, producing more output logs. Defaults to False.

        """
        # Save input arguments
        self.flag = flag
        self.name = name
        self.func = func
        self.args = args
        self.output_filenames = output_filenames
        self.use_cache = use_cache
        self.debug = debug
        # Set the key used to store and retireve data in the parent and cache
        self.parent_output_key = f'_{self.flag}_task_output'
        self.parent_output_file_keys = {}
        for argument in self.output_filenames.keys():
            self.parent_output_file_keys[argument] = f'{self.flag}_task_{argument}'
        self.cache_arg_cksums_key = f'{self.flag}_task_arg_cksums'
        self.cache_output_key = f'{self.flag}_task_output'
        # Set keys for when a output-generator task has not output given the inputs
        self.cache_missing_output_dir_key = f'{self.flag}_missing_output_directory'
        self.cache_missing_output_file_keys = {}
        for argument in self.output_filenames.keys():
            self.cache_missing_output_file_keys[argument] = f'{self.flag}_missing_{argument}'
        # Para el arg_cksum
        self.__name__ = self.flag

    # Set internal functions to handle parent saved output
    # This output is not saved in the task itself, but in the parent, because the task is static
    def _get_parent_output(self, parent):
        return safe_getattr(parent, self.parent_output_key, MISSING_VALUE_EXCEPTION)
    def _set_parent_output(self, parent, new_output):
        return setattr(parent, self.parent_output_key, new_output)
    def _get_parent_output_file(self, parent, argument):
        parent_output_file_key = self.parent_output_file_keys[argument]
        return safe_getattr(parent, parent_output_file_key, MISSING_VALUE_EXCEPTION)
    def _set_parent_output_file(self, parent, argument, new_file: 'File'):
        parent_output_file_key = self.parent_output_file_keys[argument]
        return setattr(parent, parent_output_file_key, new_file)
    def get_output(self, parent):
        """Get the task output, running the task if necessary."""
        # If we already have a value stored from this run then return it
        output = self._get_parent_output(parent)
        if output != MISSING_VALUE_EXCEPTION: return output
        # Otherwise run the task and return the output
        return self(parent)
    output = property(get_output, None, None, "Task output (read only)")
    # Asking for the an output file implies running the Task, then returning the file
    # The argument must be specified, meaning the name of the output argument in the task function
    def get_output_file_getter(self, argument) -> Callable:
        def get_output_file(parent) -> 'File':
            # If we already have a file stored from this run then return it
            file = self._get_parent_output_file(parent, argument)
            if file != MISSING_VALUE_EXCEPTION: return file
            # Otherwise run the task and return the file
            self(parent)
            file = self._get_parent_output_file(parent, argument)
            if file != MISSING_VALUE_EXCEPTION: return file
            raise ValueError(f'Task {self.flag} has no output file after running')
        return get_output_file
    # output_file = property(get_output_file, None, None, "Task output file (read only)")

    def __repr__(self):
        """Represent the Task object with its flag."""
        return f'<Task ({self.flag})>'

    def __call__(self, parent: Union['Project', 'MD']):
        """Call the task, running it if necessary."""
        return self._run(parent)

    def _run (self, parent: Union['Project', 'MD']):
        # First of all check if this task has been already done in this very run
        # If so then return the stored vale
        output = self._get_parent_output(parent)
        if output != MISSING_VALUE_EXCEPTION: return output
        # Process the task function arguments
        processed_args = {}
        # Get the task function expected arguments
        specification = getfullargspec(self.func)
        expected_arguments = specification.args
        n_default_arguments = len(specification.defaults) if specification.defaults else 0
        # Find out which arguments are optional since they have default values
        default_arguments = set(expected_arguments[::-1][:n_default_arguments])
        # If output_filenames are among the expected arguments then set them here
        writes_output_files = len(self.output_filenames) > 0
        output_files = {}
        for argument, output_filename in self.output_filenames.items():
            # Make sure the declared output filename is expected
            if argument not in expected_arguments:
                raise RuntimeError(f'Unexpected argument "{argument}" in function "{self.func}"')
            # Set the actual output filename, which may be calculated through a function
            if type(output_filename) == str:
                output_filepath = parent.pathify(output_filename)
            elif callable(output_filename):
                # If it is a function then it must be fed with the parent
                # DANI: Solo hay un caso para esto: el topology_filepath del inpro
                # DANI: En este caso la topology ya viene con el path relativo a proyecto
                # DANI: No hay que hacerlo relativo a MD
                output_filepath = output_filename(parent)
            else: raise RuntimeError(f'Unexpected output filename type "{type(output_filename)}"')
            # Set the output file
            if type(output_filepath) == str:
                output_file = File(output_filepath)
            elif type(output_filepath) == Exception:
                output_file = output_filepath
            else: raise RuntimeError(f'Unexpected output filepath type "{type(output_filepath)}"')
            output_files[argument] = output_file
            self._set_parent_output_file(parent, argument, output_file)
            # Add it to the processed args
            processed_args[argument] = output_file
            # Remove the expected argument from the list
            expected_arguments.remove(argument)
        # If one of the expected arguments is the output_directory then set it here
        # We will set a new directory with the flag name of the task, in the correspoding path
        # Note that while the task is beeing done the output directory has a different name
        # Thus the directory is hidden and marked as incomplete
        # The final output directory is the one without the incomplete prefix
        writes_output_dir = OUTPUT_DIRECTORY_ARG in expected_arguments
        incomplete_output_directory = None
        final_output_directory = None
        if writes_output_dir:
            # Set the output directory path
            incomplete_output_directory = parent.pathify(INCOMPLETE_PREFIX + self.flag)
            final_output_directory = parent.pathify(self.flag)
            # Add it to the processed args
            processed_args[OUTPUT_DIRECTORY_ARG] = incomplete_output_directory
            # Remove the expected argument from the list
            expected_arguments.remove(OUTPUT_DIRECTORY_ARG)
        # Iterate the reamining expected arguments
        for arg in expected_arguments:
            # First find the argument among the parent properties
            arg_value = self._find_arg_value(arg, parent, default_arguments)
            if arg_value == MISSING_ARGUMENT_EXCEPTION: continue
            # Add the processed argument
            processed_args[arg] = arg_value
        # Check again if the task has output already
        # It may happend that some dependencies assign output on their own
        # e.g. charges, bonds
        # If so then return the stored vale
        output = self._get_parent_output(parent)
        if output != MISSING_VALUE_EXCEPTION: return output
        # Find if we have cached output
        if self.use_cache:
            output = parent.cache.retrieve(self.cache_output_key, MISSING_VALUE_EXCEPTION)
            self._set_parent_output(parent, output)
        # Check if this dependency is to be overwriten
        forced_overwrite = self.flag in parent.overwritables
        changed_inputs, had_cache, cache_cksums = self._get_changed_inputs(parent, processed_args)
        any_input_changed = len(changed_inputs) > 0
        # Update the cache inputs
        parent.cache.update(self.cache_arg_cksums_key, cache_cksums)
        # We must overwrite outputs either if inputs changed or if it was forced by the user
        must_overwrite = forced_overwrite or any_input_changed
        # Check if output already exists
        # If the final directory already exists then it means the task was started in a previous run
        existing_incomplete_output = writes_output_dir and exists(incomplete_output_directory)
        # If the final directory already exists then it means the task was done in a previous run
        existing_final_output = writes_output_dir and exists(final_output_directory)
        # Check if the output directory was missing after completing the task in the last run
        naturally_missing_output = writes_output_dir and \
            parent.cache.retrieve(self.cache_missing_output_dir_key, False)
        # If the output file already exists then it also means the task was done in a previous run
        any_existing_output_files = writes_output_files and \
            any( output_file == None or type(output_file) == Exception or output_file.exists
                 for output_file in output_files.values() )
        # Check if all output file exist or are naurally missing
        all_existing_or_naturally_missing_output_files = writes_output_files and \
            all( output_file == None or type(output_file) == Exception or output_file.exists
                 or parent.cache.retrieve(self.cache_missing_output_file_keys[argument], False)
                 for argument, output_file in output_files.items() )
        # If we already have a cached output result
        existing_output_data = output != MISSING_VALUE_EXCEPTION
        # If we must overwrite then purge previous outputs
        if must_overwrite:
            if existing_incomplete_output: rmtree(incomplete_output_directory)
            if existing_final_output: rmtree(final_output_directory)
            if any_existing_output_files:
                for output_file in output_files.values():
                    if output_file == None or type(output_file) == Exception: continue
                    if output_file.exists: output_file.remove()
            if existing_output_data: parent.cache.delete(self.cache_output_key)
        # If already existing output is not to be overwritten then check if it is already what we need
        else:
            # If output files/directories are expected then they must exist
            # If output data is expected then it must be cached
            satisfied_output = existing_output_data \
                and (not writes_output_dir or existing_final_output or naturally_missing_output) \
                and (not writes_output_files or all_existing_or_naturally_missing_output_files)
            if self.debug:
                print(f'existing output data: {existing_output_data}')
                print(f'existing output dir: {not writes_output_dir or existing_final_output or naturally_missing_output}')
                print(f'existing output files: {not writes_output_files or all_existing_or_naturally_missing_output_files}')
            # If we already have the expected output then we can skip the task at all
            if satisfied_output:
                print(f'{GREY_HEADER}-> Task {self.flag} ({self.name}) already completed{COLOR_END}')
                return output
        # If we are at this point then we are missing some output so we must proceed to run the task
        # Use the final output directory instead of the incomplete one if exists
        # Note that we must check if it exists again since it may have been deleted since the last check
        if writes_output_dir and exists(final_output_directory):
            processed_args[OUTPUT_DIRECTORY_ARG] = final_output_directory
        # Create the incomplete output directory, if necessary
        missing_incomplete_output = writes_output_dir \
            and not exists(incomplete_output_directory) \
            and not exists(final_output_directory)
        if missing_incomplete_output: mkdir(incomplete_output_directory)
        # Finally call the function
        print(f'{GREEN_HEADER}-> Running task {self.flag} ({self.name}){COLOR_END}')
        # If the task is to be run again because an inputs changed then let the user know
        if any_input_changed and had_cache and not forced_overwrite:
            changes = ''.join([ '\n   - ' + inp for inp in changed_inputs ])
            print(f'{GREEN_HEADER}   The task is run again since the following inputs changed:{changes}{COLOR_END}')
        # Save a few internal values the task although the task is static
        # We save it right before calling the function in case the function uses this task as input
        self.changed_inputs = changed_inputs
        self.cache_cksums = cache_cksums
        # Run the actual task
        start_time = time.time()
        output = self.func(**processed_args)
        end_time = time.time()
        print(f'   Task {self.flag} completed in {end_time - start_time:.2f} seconds{COLOR_END}')
        self._set_parent_output(parent, output)
        if self.use_cache:
            # Save output in cache unless it is marked to not save it
            # Note that all must be JSON serializable values
            parent.cache.update(self.cache_output_key, output)
        # Update the overwritables so this is not remade further in the same run
        parent.overwritables.discard(self.flag)
        # Check if, although it should have written output, it did not
        # Some tasks may no write output under certain conditions
        # e.g. protmap when there is no protein
        # In these cases the task must remember there is no output to not run again further
        if writes_output_dir and not exists(incomplete_output_directory):
            parent.cache.update(self.cache_missing_output_dir_key, True)
        for argument, output_file in output_files.items():
            if output_file == None or type(output_file) == Exception: continue
            if output_file.exists: continue
            cache_key = self.cache_missing_output_file_keys[argument]
            parent.cache.update(cache_key, True)
        # Change the incomplete directory name to its final name
        # We do not remove the directory if it is empty anymore
        # The empty directory stands as a proof that the task was run successfully
        # Thus its existance prevents the task to be run again further
        # Note that some tasks clean their own intermediate steps to save disk (e.g. inpro)
        if writes_output_dir and exists(incomplete_output_directory):
            rename(incomplete_output_directory, final_output_directory)
        return output

    def _find_arg_value(self, arg: str, parent: Union['Project', 'MD'], default_arguments: set):
        """Find argument values, thus running any dependency if necessary."""
        # Word 'task' is reserved for getting the task itself
        if arg == 'task': return self
        # Word 'self' is reserved for getting the caller Project/MD
        if arg == 'self': return parent
        # Check if the argument is an MD property
        arg_value = safe_getattr(parent, arg, MISSING_ARGUMENT_EXCEPTION)
        if arg_value != MISSING_ARGUMENT_EXCEPTION: return arg_value
        # If the parent is an MD then it may happen the property is from the Project
        # We can not use the 'isinstance' function here because we can not import the MD class
        if parent.__class__.__name__ == 'MD':
            arg_value = safe_getattr(parent.project, arg, MISSING_ARGUMENT_EXCEPTION)
            if arg_value != MISSING_ARGUMENT_EXCEPTION: return arg_value
        # If the property is missing then search among the additional arguments
        arg_value = self.args.get(arg, MISSING_ARGUMENT_EXCEPTION)
        if arg_value != MISSING_ARGUMENT_EXCEPTION: return arg_value
        # It may also happen that the argument has a default value
        # If this is the case then we can skip it
        if arg in default_arguments: return MISSING_ARGUMENT_EXCEPTION
        # NEVER FORGET: Function arguments must have the same name that the Project/MD property
        # If the argument is still missing then you programmed the function wrongly or...
        # You may have forgotten the additional argument in the task args
        raise RuntimeError(f'Function "{self.func.__name__}" from task "{self.flag}" expects argument "{arg}" but it is missing')

    def _get_changed_inputs(self,
        parent: Union['Project', 'MD'],
        processed_args: dict,
    ) -> tuple[list[str], bool, dict]:
        """Find out if input arguments changed regarding the last run.

        Get the list of inputs which have changed compared to a previous run.
        If had_cache is false then it means this is the first time the task is ever done.

        Warning:
            Always get changed inputs, since this function updates the cache.

        """
        # Get cache argument references
        cache_cksums = parent.cache.retrieve(self.cache_arg_cksums_key, MISSING_VALUE_EXCEPTION)
        had_cache = False if cache_cksums == MISSING_VALUE_EXCEPTION else True
        if cache_cksums == MISSING_VALUE_EXCEPTION: cache_cksums = {}
        # Check argument by argument
        # Keep a list with arguments which have changed
        unmatched_arguments = []
        for arg_name, arg_value in processed_args.items():
            # Skip the output directory argument and an output filename arguments
            # Changes in these arguments are not actual "input" changes
            if arg_name == OUTPUT_DIRECTORY_ARG: continue
            if arg_name in self.output_filenames: continue
            # Get the cksum from the new argument value
            new_cksum = get_cksum_id(arg_value)
            # Retrieve the cksum from the old argument value
            old_cksum = cache_cksums.get(arg_name, None)
            if self.debug: print(f'Task "{self.name}" -> argument "{arg_name}"\n' +
                f' new value: {arg_value}\n' +
                f' new value cksum: {new_cksum}\n' +
                f' old value cksum: {old_cksum}\n' +
                f' match: {new_cksum == old_cksum}')
            # Compare new and old cksums
            if new_cksum != old_cksum:
                # If we found a missmatch then add it to the list
                unmatched_arguments.append(arg_name)
                # Update the references
                cache_cksums[arg_name] = new_cksum
        return unmatched_arguments, had_cache, cache_cksums

    def prefill(self, parent: Union['Project', 'MD'], output, inputs):
        """Assign an output value to a task thus marking it as already run.
        This function is used in a very specific scenario:
        Tasks which generate no output files/directory inside the 'inpro' task.
        """
        # Save the task output in the parent instance
        self._set_parent_output(parent, output)
        # Update cache output unless it is marked to not save it
        if self.use_cache: parent.cache.update(self.cache_output_key, output)
        # Save input cksums to avoid repeating this task in future runs
        _, _, cache_cksums = self._get_changed_inputs(parent, inputs)
        parent.cache.update(self.cache_arg_cksums_key, cache_cksums)
