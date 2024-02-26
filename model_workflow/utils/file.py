from os import remove, symlink, rename
from os.path import exists, isabs, abspath, relpath, split
from typing import Optional

from model_workflow.utils.constants import EXTENSION_FORMATS, PYTRAJ_SUPPORTED_FORMATS, PYTRAJ_PARM_FORMAT
from model_workflow.utils.auxiliar import InputError

# A file handler
# Absolute paths are used in runtime
# Relative paths are used to store paths
class File:
    def __init__ (self, relative_or_basolute_path : Optional[str]):
        # Declare all attributes as none by default
        self.absolute_path = self.relative_path = self.path = None
        self.basepath = self.filename = None
        self.extension = self.format = None
        self.extensionless_filename = None
        # If there is no path then leave everythin as none
        if not relative_or_basolute_path:
            return
        # If input path is absolute
        if isabs(relative_or_basolute_path[0]):
            self.absolute_path = relative_or_basolute_path
            self.relative_path = relpath(self.absolute_path, '.')
        # If it is relative
        else:
            self.relative_path = relative_or_basolute_path
            self.absolute_path = abspath(self.relative_path)
        # When simply a path is requested we return the absolute path
        self.path = self.absolute_path
        # Capture the filename
        self.basepath, self.filename = split(self.path)
        # Set the file extension
        self.extension = self.filename.split('.')[-1]
        if self.extension == self.filename:
            self.extension = None
        # Set the extensionless filename
        self.extensionless_filename = self.filename
        if self.extension:
            extension_size = len(self.extension) + 1 # We include here the dot
            self.extensionless_filename = self.filename[-extension_size:]
        # Set the file format
        if self.extension:
            extension_format = EXTENSION_FORMATS.get(self.extension, None)
            if not extension_format:
                raise InputError('Not recognized format extension "' + self.extension + '" from file "' + self.filename + '"')
            self.format = extension_format
        # Set a couple of additional values according to pytraj format requirements
        self.is_pytraj_supported = self.format in PYTRAJ_SUPPORTED_FORMATS
        self.pytraj_parm_format = PYTRAJ_PARM_FORMAT.get(self.format, None)

    def __repr__ (self) -> str:
        if not self.filename:
            return '< No file >'
        return '< File ' + self.filename + ' >'

    def __str__ (self) -> str:
        return self.__repr__()

    def __bool__ (self) -> bool:
        return bool(self.filename)

    def __eq__ (self, other : 'File') -> bool:
        if isinstance(other, self.__class__):
            return self.path == other.path
        return False

    def check_existence (self) -> bool:
        return exists(self.absolute_path)
    exists = property(check_existence, None, None, "Does the file exists? (read only)")

    # Remove the file
    def remove (self):
        remove(self.absolute_path)

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
        return File(self.basepath + '/' + prefix + self.filename)

    # Set this file a symlink to another file
    def set_symlink_to (self, other_file : 'File'):
        # Self file must not exist
        if self.exists:
            raise Exception('Cannot set a symlink from an already existing file: ' + str(self))
        # Note that symlink path must be relative to this file
        relative_path = relpath(other_file.path, self.basepath)
        # Set the symlink
        symlink(relative_path, self.path)