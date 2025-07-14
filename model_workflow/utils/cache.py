from model_workflow.utils.auxiliar import load_json, save_json
from model_workflow.utils.type_hints import *

# The cache is used to store data to be reused between different runs
class Cache:
    def __init__ (self, cache_file : 'File'):
        # Save the previous cache file
        self.file = cache_file
        # Set an auxiliar file which is used to not overwrite directly the previous cache
        self._aux_file = cache_file.get_prefixed_file('.aux')
        # Set a dict to store the actual cached data
        self.data = {}
        # Load data from the cache file, if it already exists
        if self.file.exists:
            # Read the cache in disk
            previous_data = load_json(self.file.path)
            # Inherit every field
            for field_name, field_value in previous_data.items():
                self.data[field_name] = field_value
        # Save the entry for the first time
        self.save()

    # Read a value from the cache
    def retrieve (self, key : str, fallback = None):
        return self.data.get(key, fallback)

    # Update the cache and save it to disk
    def update (self, key : str, value):
        self.data[key] = value
        self.save()

    # Delete a value in the cache
    def delete (self, key : str):
        if key not in self.data: return
        del self.data[key]
        self.save()

    # Reset the cache
    # This is called when some input files are modified
    def reset (self):
        self.data = {}
        self.save()

    # Save the cache to disk, as a json file
    def save (self):
        # Write entries to disk
        # Write to the auxiliar file, thus not overwritting the original cache yet
        save_json(self.data, self._aux_file.path, indent = 4)
        # If the new cache is successfully written then replace the old one
        self._aux_file.rename_to(self.file)