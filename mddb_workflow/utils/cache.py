from mddb_workflow.utils.auxiliar import load_json, save_json
from mddb_workflow.utils.arg_cksum import get_cksum_id
from mddb_workflow.utils.type_hints import *
from mddb_workflow.tools.get_inchi_keys import InChIKeyData

class Cache:
    """The cache is used to store data to be reused between different runs."""

    def __init__ (self, cache_file : 'File'):
        # Save the previous cache file
        self.file = cache_file
        # Set an auxiliar file which is used to not overwrite directly the previous cache
        self._aux_file = cache_file.get_prefixed_file('.aux')
        # Set a dict to store the actual cached data
        self.data = {}
        # Load data from the cache file, if it already exists
        self.load()
        # Save the entry for the first time
        self.save()

    def retrieve (self, key : str, fallback = None):
        """Read a value from the cache."""
        return self.data.get(key, fallback)

    def update (self, key : str, value):
        """Update the cache and save it to disk."""
        self.data[key] = value
        self.save()

    def delete (self, key : str):
        """Delete a value in the cache."""
        if key not in self.data: return
        del self.data[key]
        self.save()

    def reset (self):
        """Reset the cache. This is called when some input files are modified."""
        self.data = {}
        self.save()

    def load (self):
        """Load the cache to memory, as a dict."""
        # Load data from the cache file, if it already exists
        if not self.file.exists: return
        # Read the cache in disk
        self.data = load_json(self.file.path)
        if 'inchikeys_task_output' in self.data:
            self.data['inchikeys_task_output'] = InChIKeyData.load_cache(self.data['inchikeys_task_output'])


    def save (self):
        """Save the cache to disk, as a json file."""
        # Write entries to disk
        # Write to the auxiliar file, thus not overwritting the original cache yet
        save_json(self.data, self._aux_file.path, indent = 4)
        # If the new cache is successfully written then replace the old one
        self._aux_file.rename_to(self.file)

# Set a function to create a cached function
# The cached function will store input cksums and output in cache
# The cached function will return cached output when input cksums are identical
def get_cached_function (function : Callable, cache : Cache) -> Callable:
    cache_key = function.__name__
    def cached_function (*inputs):
        # Get all input cksums
        current_result = {}
        for i, input in enumerate(inputs, 1):
            input_key = f'input_{i}_cksum'
            input_cksum = get_cksum_id(input)
            current_result[input_key] = input_cksum
        # Find if there is a cached result with identical input keys
        cached_results = cache.retrieve(cache_key, [])
        for cached_result in cached_results:
            # Make sure every input cksum matches
            missmatch = False
            for input_key, new_input_cksum in current_result.items():
                cached_input_cksum = cached_result[input_key]
                if cached_input_cksum != new_input_cksum:
                    missmatch = True
                    break
            if missmatch: continue
            # If we make it this far with no missmatch then it means we have a previous result which matches
            # Just return the cached output
            return cached_result['output']
        # Otherwise we must run the function
        current_result['output'] = function(*inputs)
        # Save the result in the cache
        cached_results.append(current_result)
        cache.update(cache_key, cached_results)
        # And finally return the output
        return current_result['output']
    return cached_function