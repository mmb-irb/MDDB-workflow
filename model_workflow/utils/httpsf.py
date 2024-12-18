#!/usr/bin/env python

import os
import sys
import urllib.request

from fuse import FUSE 


# Set the functions which ar set by default when the file handler is passed to the fuse class
# These functions are better not overwritten in the file handler
# These functions were found experimentally but there may be more
default_fuse_functions = ['init', 'destroy', 'getattr', 'access']
# Set the functions which the fuse standards may expect to find in any compliant class and are not yet implemented
# These functions are set as whistleblowers
# These functions were found in the internet:
#     https://www.cs.hmc.edu/~geoff/classes/hmc.cs135.201109/homework/fuse/fuse_doc.html
#     https://www.cs.hmc.edu/~geoff/classes/hmc.cs135.201109/homework/fuse/fuse_doc.html#gotchas
#     https://libfuse.github.io/doxygen/structfuse__operations.html
not_implememnted_fuse_functions = ['fgetattr', 'readlink' 'opendir', 'readdir', 'mknod', 'mkdir', 'unlink',
    'rmdir', 'symlink', 'rename', 'link', 'chmod', 'chown', 'truncate', 'ftruncate', 'utimens', 'open', 'read',
    'write', 'statfs', 'release', 'releasedir', 'fsync', 'fsyncdir', 'flush', 'lock', 'bmap', 'setxattr',
    'getxattr', 'listxattr', 'removexattr', 'ioctl', 'poll', 'create', 'openat', 'write_buf', 'read_buf',
    'flock', 'fallocate', 'copy_file_range', 'lseek']

# FUSE implementation for a single file
class FileHandler ():
    def __init__ (self, url : str):
        self.url = url
        self._response = None
        # Set all missing fuse functions as whistleblowers
        # Most commands will return a clueless 'function not implemented' error when trying to access any of these functions
        # This prevents the error and intead tells the user which function is being called and not yet implemented
        # This double function is to avoid overwritting the same function every time
        def create_whistleblower (name : str):
            def whistleblower (*args):
                raise SystemExit('FUSE function not implemented: ' + name)
            return whistleblower
        for function_name in not_implememnted_fuse_functions:
            if hasattr(self, function_name):
                continue
            setattr(self, function_name, create_whistleblower(function_name))

    def __call__ (self, function_name : str, *args, **kwargs):
        function = getattr(self, function_name)
        return function(*args)

    # Attribute keys
    # Got from https://github.com/skorokithakis/python-fuse-sample/blob/master/passthrough.py
    stat_keys = 'st_atime', 'st_ctime', 'st_mtime', 'st_gid', 'st_uid', 'st_mode', 'st_nlink', 'st_size'

    # This must be implemented, although there is a default function
    # FUSE relies in this function a lot and not implementing it properly leads to silent errors
    def getattr (self, path, fh=None):
        # Get the stat of current directory to get time and id values
        st = os.lstat('.')
        # Example regular file:
        # (st_mode=33188, st_ino=13765084, st_dev=66306, st_nlink=1, st_uid=16556, st_gid=500, st_size=9410,
        # st_atime=1696945741, st_mtime=1696945600, st_ctime=1696945600)
        # Example fifo
        # (st_mode=4516, st_ino=13765935, st_dev=66306, st_nlink=1, st_uid=16556, st_gid=500, st_size=0,
        # st_atime=1697017261, st_mtime=1697017261, st_ctime=1697017261)
        stat = dict((key, getattr(st, key)) for key in self.stat_keys)
        # Set a size bigger than the expected response, no matter how much bigger
        # DANI: Intenté muchas formas de que la size fuese nula o infinita para que no haya limite, pero no pude
        # DANI: Prové a cambiar el st_mode para convertir el punto de montaje en ostras cosas, pero solo puede ser file
        stat['st_size'] = sys.maxsize
        # Set a custom mode for FUSE to treat it as a file with infinite size
        # To see the actual meaningful value of st_mode you have to do 'oct(st_mode)'
        # e.g. regular file: 33188 -> 0o100644, fifo: 4516 -> 0o10644
        # Last three numbers stand for the permissions
        # The starting numbers stand for the file type
        # UNIX file type docs
        # https://man7.org/linux/man-pages/man7/inode.7.html (section 'The file type and mode')
        # Set the mode in octal forma: regular file
        mode = '0100644'
        # Set the stat mode in decimal format, as it is expected
        stat['st_mode'] = int(mode, 8)
        return stat

    # File methods
    # ============

    # Open the URL request
    def open (self, path : str, flags) -> int:
        self._response = urllib.request.urlopen(self.url)
        return 0
    
    # Read length is set automatically by FUSE
    # This length is calculated from the getattr 'st_size' value
    # However the length is not the exact value of st_size, but the first multiple of 4096 which exceeds its value
    # Also length has a maximum value: 131072, experimentally observed 
    # When the maximum value is not enought to cover the length the read function is called several times
    # Each time with different offset value
    def read (self, path : str, length : int, offset : int, fh : int) -> str:
        return self._response.read(length)

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
    FUSE(handler, mountpoint, nothreads=False, foreground=False)