# This script is used to get the last trajectory frame
import os
from subprocess import run, PIPE, Popen
from model_workflow.tools.get_frames_count import get_frames_count

# input_topology_filename - The name string of the input topology file (path)
# Tested supported formats are .pdb and .tpr
# input_trajectory_filename - The name string of the input trajectory file (path)
# Tested supported formats are .trr and .xtc
# first_frame_filename - The name string of the output first frame file (path)
# Tested supported format is .pdb
def get_last_frame (
    input_topology_filename : str,
    input_trajectory_filename : str,
    last_frame_filename : str
    ):

    # Stop here in case the file already exists
    if os.path.exists(last_frame_filename):
        return

    # Get the number of frames in trajectory
    frames_count = get_frames_count(input_topology_filename, input_trajectory_filename)

    # Run Gromacs
    p = Popen([
        "echo",
        "System",
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "trjconv",
        "-s",
        input_topology_filename,
        "-f",
        input_trajectory_filename,
        '-o',
        last_frame_filename,
        '-dump',
        str(frames_count - 1),
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE, stderr=PIPE).stdout.decode()
    p.stdout.close()