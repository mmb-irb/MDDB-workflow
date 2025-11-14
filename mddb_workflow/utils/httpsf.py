#!/usr/bin/env python

import os
import sys
import requests
import stat      # For file mode constants
from fuse import FUSE, Operations  # FUSE bindings with fusepy


# Set the functions which are set by default when the file handler is passed to the fuse class
# These functions are better not overwritten in the file handler
# These functions were found experimentally but there may be more
default_fuse_functions = ['init', 'destroy', 'getattr', 'access']

# Set the functions which the fuse standards may expect to find in any compliant class and are not yet implemented
# These functions are set as whistleblowers
# These functions were found in the internet:
#     https://www.cs.hmc.edu/~geoff/classes/hmc.cs135.201109/homework/fuse/fuse_doc.html
#     https://www.cs.hmc.edu/~geoff/classes/hmc.cs135.201109/homework/fuse/fuse_doc.html#gotchas
#     https://libfuse.github.io/doxygen/structfuse__operations.html
not_implememnted_fuse_functions = ['fgetattr', 'readlink', 'opendir', 'readdir', 'mknod', 'mkdir', 'unlink',
    'rmdir', 'symlink', 'rename', 'link', 'chmod', 'chown', 'truncate', 'ftruncate', 'utimens', 'open', 'read',
    'write', 'statfs', 'release', 'releasedir', 'fsync', 'fsyncdir', 'flush', 'lock', 'bmap', 'setxattr',
    'getxattr', 'listxattr', 'removexattr', 'ioctl', 'poll', 'create', 'openat', 'write_buf', 'read_buf',
    'flock', 'fallocate', 'copy_file_range', 'lseek']

# FUSE implementation for a single file
class FileHandler(Operations):
    def __init__ (self, url : str):
        print("Initializing API Virtual File System...")
        self.url = url
        self._response = None
        # Set all missing fuse functions as whistleblowers
        # Most commands will return a clueless 'function not implemented' error when trying to access any of these functions
        # This prevents the error and intead tells the user which function is being called and not yet implemented
        # This double function is to avoid overwritting the same function every time
        def create_whistleblower (name : str):
            def whistleblower (*args):
                raise RuntimeError('FUSE function not implemented: ' + name)
            return whistleblower
        for function_name in not_implememnted_fuse_functions:
            if hasattr(self, function_name):
                continue
            setattr(self, function_name, create_whistleblower(function_name))

    def __call__ (self, function_name : str, *args, **kwargs):
        function = getattr(self, function_name)
        return function(*args)

    def _get_file_size(self):
        response = requests.head(self.url, allow_redirects=True)
        response.raise_for_status()  # Raise an exception for HTTP errors
        # Check if Content-Length header exists
        size = int(response.headers.get('Content-Length', 0))
        return size
        
    def _fetch_file_content(self):
        """Fetch content from remote API."""
        self._response = requests.get(self.url)
        self._response.raise_for_status()
    
    # Attribute keys
    # Got from https://github.com/skorokithakis/python-fuse-sample/blob/master/passthrough.py
    stat_keys = 'st_atime', 'st_ctime', 'st_mtime', 'st_gid', 'st_uid', 'st_mode', 'st_nlink', 'st_size'

    # This must be implemented, although there is a default function
    # FUSE relies in this function a lot and not implementing it properly leads to silent errors
    def getattr (self, path, fh=None):
        # Get the stat of current directory to get time and id values
        dot_st = os.lstat('.')
        # Example regular file:
        # (st_mode=33188, st_ino=13765084, st_dev=66306, st_nlink=1, st_uid=16556, st_gid=500, st_size=9410,
        # st_atime=1696945741, st_mtime=1696945600, st_ctime=1696945600)
        # Example fifo
        # (st_mode=4516, st_ino=13765935, st_dev=66306, st_nlink=1, st_uid=16556, st_gid=500, st_size=0,
        # st_atime=1697017261, st_mtime=1697017261, st_ctime=1697017261)
        st = dict((key, getattr(dot_st, key)) for key in self.stat_keys)
        # Set a size bigger than the expected response, no matter how much bigger
        # DANI: Intenté muchas formas de que la size fuese nula o infinita para que no haya limite, pero no pude
        # DANI: Prové a cambiar el st_mode para convertir el punto de montaje en ostras cosas, pero solo puede ser file
        #stat['st_size'] = sys.maxsize
        st['st_size'] = self._get_file_size()
        # Set a custom mode for FUSE to treat it as a file with infinite size
        # To see the actual meaningful value of st_mode you have to do 'oct(st_mode)'
        # e.g. regular file: 33188 -> 0o100644, fifo: 4516 -> 0o10644
        # Last three numbers stand for the permissions
        # The starting numbers stand for the file type
        # UNIX file type docs
        # https://man7.org/linux/man-pages/man7/inode.7.html (section 'The file type and mode')
        # Set the mode for a FIFO
        st['st_mode'] = 0o644 | stat.S_IFIFO
        return st

    # File methods
    # ============

    # Open the URL request
    def open (self, path : str, flags) -> int:
        # Lazy load content if not already fetched
        if self._response is None:
            self._fetch_file_content()
        # Return a file descriptor (dummy in this case)
        return 0
    
    # Read length is set automatically by FUSE
    # This length is calculated from the getattr 'st_size' value
    # However the length is not the exact value of st_size, but the first multiple of 4096 which exceeds its value
    # Also length has a maximum value: 131072, experimentally observed 
    # When the maximum value is not enought to cover the length the read function is called several times
    # Each time with different offset value
    def read (self, path : str, length : int, offset : int, fh : int) -> str:
        if self._response is None:
            self._fetch_file_content()
        return self._response.content[offset:offset + length]

    # Close the URL request
    def destroy (self, path : str):
        if self._response:
            self._response.close()

    # The following functions are not implemented
    # However they must be defined or FUSE complains when trying to read the file
    flush = None
    release = None
    lock = None
    # Note the FUSE will complain for several functions to be not implemented when trying to write in the file
    # User should never try to write this file so these cases are not handled
    
# Mount a fake file which returns an url content when read
def mount (url : str, mountpoint : str):
    # Create a fake file to be opened in case it does not exists yet
    handler = FileHandler(url)
    # Create the mount point in case it does not exist yet
    if not os.path.exists(mountpoint):
        open(mountpoint, 'w').close()
    # Mount the handler in the mount point
    FUSE(handler, mountpoint, nothreads=True, foreground=True)

if __name__ == '__main__':
    # Example usage
    url = 'https://cineca.mddbr.eu/api/rest/v1/projects/A0001/files/structure.pdb'
    mountpoint = '/home/rchaves/Downloads/structure.pdb'
    mount(url, mountpoint)