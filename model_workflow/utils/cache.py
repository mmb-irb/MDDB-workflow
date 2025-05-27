from model_workflow.utils.auxiliar import load_json, save_json
from model_workflow.utils.type_hints import *

# The cache is used to store data to be reused between different runs
class Cache:
    def __init__ (self, cache_file : 'File'):
        # Save the previous cache file
        self.file = cache_file
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

    # Reset the cache
    # This is called when some input files are modified
    def reset (self):
        self.data = {}
        self.save()

    # Save the cache to disk, as a json file
    def save (self):
        print(' Saving cache')
        # Write entries to disk
        save_json(self.data, self.file.path, indent = 4)