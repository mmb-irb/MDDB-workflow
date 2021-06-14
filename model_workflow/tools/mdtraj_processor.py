# This script is used to merge trajectory files and convert both topology and trajectory files formats
# Everything is carried by MDtraj

import os

from subprocess import run, PIPE, Popen

import mdtraj as mdt

# Multiple files may be selected with bash syntax (e.g. *.dcd)
# Tested supported input formats are .dcd
# Tested supported output formats are .xtc
def merge_and_convert_traj (
    input_filenames : list,
    output_filename : str
    ):

    print('Converting: ' + str(input_filenames) + ' -> ' + output_filename)

    # Run MDtraj
    run([
        "mdconvert",
        "-o",
        output_filename,
        *input_filenames,
    ])