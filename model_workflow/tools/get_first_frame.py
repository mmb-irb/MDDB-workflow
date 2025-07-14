# This script is used to get the first trajectory frame
# This process is carried by Gromacs

from os.path import exists
from subprocess import run, PIPE, Popen

from model_workflow.utils.constants import GROMACS_EXECUTABLE
from model_workflow.utils.type_hints import *

# Get the trajectory first frame in PDB format
# DANI: No lo muevo a gmx spells porque allí ya hay un get_first_frame con otra finalidad
def get_first_frame (
    structure_file : 'File',
    trajectory_file : 'File',
    output_filepath : str
    ):

    # Run Gromacs
    p = Popen([
        "echo",
        "System",
    ], stdout=PIPE)
    process = run([
        GROMACS_EXECUTABLE,
        "trjconv",
        "-s",
        structure_file.path,
        "-f",
        trajectory_file.path,
        '-o',
        output_filepath,
        '-dump',
        '0',
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE, stderr=PIPE)
    logs = process.stdout.decode()
    p.stdout.close()

    # If output has not been generated then warn the user
    if not exists(output_filepath):
        print(logs)
        error_logs = process.stderr.decode()
        print(error_logs)
        raise SystemExit('Something went wrong with Gromacs')