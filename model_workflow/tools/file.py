from os.path import exists, isabs, abspath, relpath, split
from typing import Optional

from model_workflow.constants import EXTENSION_FORMATS, PYTRAJ_SUPPORTED_FORMATS, PYTRAJ_PARM_FORMAT

# A file handler
# Absolute paths are used in runtime
# Relative paths are used to store paths
class File:
    def __init__ (self, relative_or_basolute_path : Optional[str]):
        self.absolute_path = None
        self.relative_path = None
        self.filename = None
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
        self.filename = split(relative_or_basolute_path)[-1]
        # Set the file extension
        self.extension = self.filename.split('.')[-1]
        if self.extension == self.filename:
            self.extension = None
        # Set the file format
        self.format = None
        if self.extension:
            extension_format = EXTENSION_FORMATS.get(self.extension, None)
            if not extension_format:
                raise Exception('Not recognized format extension "' + self.extension + '" from file "' + self.filename + '"')
            self.format = extension_format
        # Set a couple of additional values according to pytraj format requirements
        self.is_pytraj_supported = self.format in PYTRAJ_SUPPORTED_FORMATS
        self.pytraj_parm_format = PYTRAJ_PARM_FORMAT.get(self.format, None)

    def __repr__ (self) -> str:
        return '< File ' + self.filename + ' >'

    def __str__ (self) -> str:
        return self.__repr__()

    def __bool__ (self) -> bool:
        return bool(self.filename)

    def check_existence (self) -> bool:
        return exists(self.absolute_path)
    exists = property(check_existence, None, None, "Does the file exists? (read only)")