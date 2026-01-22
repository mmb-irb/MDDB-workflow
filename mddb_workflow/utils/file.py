from os import remove, symlink, rename, readlink
from os.path import exists, isabs, abspath, relpath, split, islink, normpath, getmtime, getsize
from time import strftime, gmtime
from shutil import copyfile
from typing import Optional
from time import time

import xxhash

from mddb_workflow.utils.constants import EXTENSION_FORMATS, PYTRAJ_SUPPORTED_FORMATS, PYTRAJ_PARM_FORMAT
from mddb_workflow.utils.constants import DATE_STYLE, GLOBALS
from mddb_workflow.utils.auxiliar import InputError

LOCAL_PATH = '.'


class File:
    """File handler class.
    Absolute paths are used in runtime.
    Relative paths are used to store paths.
    """
    def __init__(self, relative_or_absolute_path: str):
        # If there is no path then complain
        if not relative_or_absolute_path:
            raise RuntimeError('Declared file with no path')
        # If input path is absolute
        if isabs(relative_or_absolute_path[0]):
            self.absolute_path = relative_or_absolute_path
            self.relative_path = relpath(self.absolute_path, LOCAL_PATH)
        # If it is relative
        else:
            self.relative_path = relative_or_absolute_path
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
        self.extensionless_filepath = self.path
        if self.extension:
            extension_size = len(self.extension) + 1  # We include here the dot
            self.extensionless_filename = self.filename[:-extension_size]
            self.extensionless_filepath = self.path[:-extension_size]
        # Set internal values
        self._cksum = None

    # We must display the cksum here
    # Note that this is critical for the task args cksum when we handle lists of files
    # e.g. input_trajectory_files in process_input_files
    def __repr__(self) -> str:
        return f'< File {self.path} >'

    def __str__(self) -> str:
        return self.__repr__()

    def __hash__(self) -> int:
        return hash(self.path)  # Path is already normalized

    def __bool__(self) -> bool:
        return bool(self.filename)

    def __eq__(self, other: 'File') -> bool:
        if isinstance(other, self.__class__):
            return self.path == other.path  # Paths are already normalized
        return False

    def check_existence(self) -> bool:
        """Check if file exists."""
        return exists(self.path)
    exists = property(check_existence, None, None, "Does the file exists? (read only)")

    def get_format(self) -> Optional[str]:
        """Get file format based on the extension.
        If the extension is not recognized then raise an error.
        """
        if not self.extension:
            return None
        extension_format = EXTENSION_FORMATS.get(self.extension, None)
        if not extension_format:
            raise InputError(f'Not recognized format extension "{self.extension}" from file "{self.filename}"')
        return extension_format
    format = property(get_format, None, None, "File standard format (read only)")

    def get_mtime(self) -> str:
        """Get the file last modification time."""
        raw_mtime = getmtime(self.path)
        return strftime(DATE_STYLE, gmtime(raw_mtime))
    mtime = property(get_mtime, None, None, "File last modification date (read only)")

    def get_size(self) -> str:
        """Get the file size in bytes."""
        return getsize(self.path)
    size = property(get_size, None, None, "File size in bytes (read only)")

    CKSUM_UNSAFE_SIZE_LIMIT = 1024 * 1024 * 100  # 100 MB
    def get_cksum(self, unsafe: bool = False, verbose: bool = False) -> str:
        """Get a cksum code used to compare identical file content.
        Use the unsafe argument to make it way faster for large files by reading them partially.
        """
        # If we already have a value then return it
        if self._cksum != None: return self._cksum
        # If the file does not exist then there is no cksum
        if not self.exists: return None
        # Set if the cksum will be unsafe
        # Note that files lighter than the size limit will always have safe cksums
        is_unsafe = unsafe and self.size > self.CKSUM_UNSAFE_SIZE_LIMIT
        # Calculate the xxhash of the whole file content
        # This should be the faster method available whcih still reads all content
        if verbose: start_time = time()
        hasher = xxhash.xxh64()
        with open(self.path, 'rb') as file:
            if is_unsafe: hasher.update(file.read(self.CKSUM_UNSAFE_SIZE_LIMIT))
            # DANI: This is not memory safe, a big file could consume all memory
            # DANI: We should iterate chunks but I was on a hurry
            else: hasher.update(file.read())
        final_xxhash = hasher.hexdigest()
        if verbose:
            end_time = time()
            total_time = end_time - start_time
            print(f'Got cksum for {self.path} ({self.size} Bytes) in {total_time:.2f} seconds -> {final_xxhash}')
        self._cksum = f'{self.size}-{"(UNSAFE)" if is_unsafe else ""}{final_xxhash}'
        return self._cksum

    # Set a couple of additional functions according to pytraj format requirements
    def is_pytraj_supported(self) -> bool:
        return self.format in PYTRAJ_SUPPORTED_FORMATS
    def get_pytraj_parm_format(self) -> Optional[str]:
        return PYTRAJ_PARM_FORMAT.get(self.format, None)

    def remove(self):
        """Remove the file."""
        remove(self.path)

    def get_standard_file(self) -> 'File':
        """Given a file who has non-standard extension of a supported format we set a symlink with the standard extension."""
        # If current file already has the extension then there is nothing to return
        if self.extension == self.format:
            return self
        return self.reformat(self.format)

    def reformat(self, new_extension: str) -> 'File':
        """Given a file and a new extension we set a symlink from a new file with that extension."""
        # Set the filename with the standard extension and initiate the file
        reformatted_filename = f'{self.extensionless_filepath}.{new_extension}'
        reformatted_file = File(reformatted_filename)
        # If standard file does not exist then set a symlink
        if not reformatted_file.exists:
            reformatted_file.set_symlink_to(self)
        return reformatted_file

    def get_prefixed_file(self, prefix: str) -> 'File':
        """Get a prefixed file using this file name as the name base."""
        return File(f'{self.basepath}/{prefix}{self.filename}')

    def get_neighbour_file(self, filename: str) -> 'File':
        """Get a file in the same path but with a different name."""
        return File(f'{self.basepath}/{filename}')

    def get_symlink(self) -> Optional['File']:
        """Get the symlink target of this file."""
        target_filepath = readlink(self.path)
        if not target_filepath:
            return None
        return File(self.basepath + '/' + target_filepath)

    def set_symlink_to(self, other_file: 'File', force: bool = False):
        """Set this file a symlink to another file.
        Use the "force" argument to delete an already existing file/symlink.
        """
        # Self file must not exist
        if self.exists or self.is_symlink():
            if force: self.remove()
            else: raise Exception(f'Cannot set a symlink from an already existing file or symlink: {self}')
        # Check if symlinks are allowed
        no_symlinks = GLOBALS['no_symlinks']
        # If symlinks are now allowed then copy the file instead
        if no_symlinks:
            other_file.copy_to(self)
            return
        # Note that symlink path must be relative to this file
        relative_path = relpath(other_file.path, self.basepath)
        # Set the symlink
        symlink(relative_path, self.path)

    def is_symlink(self) -> bool:
        """Check if a file is already a symlink."""
        return islink(self.path)

    def copy_to(self, other_file: 'File'):
        """Copy a file to another."""
        copyfile(self.path, other_file.path)

    def rename_to(self, other_file: 'File'):
        """Rename a file to another."""
        rename(self.path, other_file.path)
