# This script is used to merge trajectory files and convert both topology and trajectory files formats
# Everything is carried by MDtraj

import os

from subprocess import run, PIPE, Popen

import mdtraj as mdt

# Multiple files may be selected with bash syntax (e.g. *.dcd)
# Tested supported input formats are .dcd
# Tested supported output formats are .xtc
def merge_dcd_files (
    input_filenames : list,
    output_filename : str
    ):

    # Run MDtraj
    logs = run([
        "mdconvert",
        "-o",
        output_filename,
        *input_filenames,
    ], stdout=PIPE).stdout.decode()

# Get the first frame from a dcd file
def get_dcd_first_frame (
    input_filename : str,
    output_filename : str):

    mdt.iterload(input_filename)