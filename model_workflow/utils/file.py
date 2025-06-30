from os import remove, symlink, rename, readlink
from os.path import exists, isabs, abspath, relpath, split, islink, normpath, getmtime, getsize
from time import strftime, gmtime
from shutil import copyfile
from typing import Optional

from model_workflow.utils.constants import EXTENSION_FORMATS, PYTRAJ_SUPPORTED_FORMATS, PYTRAJ_PARM_FORMAT
from model_workflow.utils.constants import DATE_STYLE, GLOBALS
from model_workflow.utils.auxiliar import InputError

LOCAL_PATH = '.'

# A file handler
# Absolute paths are used in runtime
# Relative paths are used to store paths
class File:
    def __init__ (self, relative_or_basolute_path : Optional[str]):
        # Declare all attributes as none by default
        self.absolute_path = self.relative_path = self.path = None
        self.basepath = self.filename = None
        self.extension = None
        self.extensionless_filename = None
        # If there is no path then leave everything as none
        if not relative_or_basolute_path:
            return
        # If input path is absolute
        if isabs(relative_or_basolute_path[0]):
            self.absolute_path = relative_or_basolute_path
            self.relative_path = relpath(self.absolute_path, LOCAL_PATH)
        # If it is relative
        else:
            self.relative_path = relative_or_basolute_path
            self.absolute_path = abspath(self.relative_path)
        # When simply a path is requested we return the relative path
        # Note that normalizing the path is essential to recognize same filepaths
        # Otherwise we could have './myfile' and 'myfile' considered as different filepaths
        self.path = normpath(self.relative_path)
        # Capture the filename and the basepath
        self.basepath, self.filename = split(self.path)
        # If the basepath is empty then it means the file is in the local directroy
        # WARNING: If the basepath is left empty an exists(basepath) would return false
        # WARNING: For this reason we must replace '' by '.'
        if not self.basepath:
            self.basepath = LOCAL_PATH
        # Set the file extension
        self.extension = self.filename.split('.')[-1]
        if self.extension == self.filename:
            self.extension = None
        # Set the extensionless filename
        self.extensionless_filename = self.filename
        if self.extension:
            extension_size = len(self.extension) + 1 # We include here the dot
            self.extensionless_filename = self.filename[-extension_size:]
        # Set internal values
        self._cksum = None

    # We must display the cksum here
    # Note that this is critical for the task args cksum when we handle lists of files
    # e.g. input_trajectory_files in process_input_files
    def __repr__ (self) -> str:
        if not self.filename:
            return '< No file >'
        return f'< File {self.cksum} >'

    def __str__ (self) -> str:
        return self.__repr__()

    def __hash__ (self) -> str:
        return hash(self.path) # Path is already normalized

    def __bool__ (self) -> bool:
        return bool(self.filename)

    def __eq__ (self, other : 'File') -> bool:
        if isinstance(other, self.__class__):
            return self.path == other.path # Paths are already normalized
        return False

    # Check if file exists
    def check_existence (self) -> bool:
        return exists(self.path)
    exists = property(check_existence, None, None, "Does the file exists? (read only)")

    # Get file format based on the extension
    # If the extension is not recognized then raise an error
    def get_format (self) -> Optional[str]:
        if not self.extension:
            return None
        extension_format = EXTENSION_FORMATS.get(self.extension, None)
        if not extension_format:
            raise InputError(f'Not recognized format extension "{self.extension}" from file "{self.filename}"')
        return extension_format
    format = property(get_format, None, None, "File standard format (read only)")

    # Get the file last modification time
    def get_mtime (self) -> str:
        raw_mtime = getmtime(self.path)
        return strftime(DATE_STYLE, gmtime(raw_mtime))
    mtime = property(get_mtime, None, None, "File last modification date (read only)")

    # Get the file size in bytes
    def get_size (self) -> str:
        return getsize(self.path)
    size = property(get_size, None, None, "File size in bytes (read only)")

    # Get a cksum code used to compare identical file content
    # DANI: This is provisional and it is not yet based in a cksum neither the file content
    def get_cksum (self) -> str:
        # If we already have an internal value then use it
        if self._cksum != None: return self._cksum
        # Calculate it otherwise
        if not self.exists: self._cksum = f'missing {self.path}'
        else: self._cksum = f'{self.path} -> {self.mtime} {self.size}'
        return self._cksum
    cksum = property(get_cksum, None, None, "Cksum code used to compare identical file content (read only)")

    # Set a couple of additional functions according to pytraj format requirements
    def is_pytraj_supported (self) -> bool:
        return self.format in PYTRAJ_SUPPORTED_FORMATS
    def get_pytraj_parm_format (self) -> Optional[str]:
        return PYTRAJ_PARM_FORMAT.get(self.format, None)

    # Remove the file
    def remove (self):
        remove(self.path)

    # Given a file who has non-standard extension of a supported format we set a symlink with the standard extension
    def get_standard_file (self) -> 'File':
        # If current file already has the extension then there is nothing to return
        if self.extension == self.format:
            return self
        # Set the filename with the standard extension and initiate the file
        standard_filename = self.extensionless_filename + '.' + self.format
        standard_file = File(standard_filename)
        # If standard file does not exist then set a symlink
        if not standard_file.exists:
            symlink(self.relative_path, standard_file.path)
        return standard_file

    # Get a prefixed file using this file name as the name base
    def get_prefixed_file (self, prefix : str) -> 'File':
        return File(f'{self.basepath}/{prefix}{self.filename}')
    
    # Get a file in the same path but whith a different name
    def get_neighbour_file (self, filename : str) -> 'File':
        return File(f'{self.basepath}/{filename}')

    # Get the symlink target of this file
    def get_symlink (self) -> Optional['File']:
        target_filepath = readlink(self.path)
        if not target_filepath:
            return None
        return File(self.basepath + '/' + target_filepath)

    # Set this file a symlink to another file
    def set_symlink_to (self, other_file : 'File'):
        # Check if symlinks are allowed
        no_symlinks = GLOBALS['no_symlinks']
        # If symlinks are now allowed then copy the file instead
        if no_symlinks:
            other_file.copy_to(self)
            return
        # Self file must not exist
        if self.exists:
            raise Exception('Cannot set a symlink from an already existing file: ' + str(self))
        # Note that symlink path must be relative to this file
        relative_path = relpath(other_file.path, self.basepath)
        # Set the symlink
        symlink(relative_path, self.path)

    # Check if a file is already a symlink
    def is_symlink (self) -> bool:
        return islink(self.path)

    # Copy a file to another
    def copy_to (self, other_file : 'File'):
        copyfile(self.path, other_file.path)

    # Rename a file to another
    def rename_to (self, other_file : 'File'):
        rename(self.path, other_file.path)