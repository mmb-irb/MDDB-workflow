# Get frames count

from subprocess import run, PIPE
import re

# Get the trajectory frames number using a preinstalled bash command called 'xtc-length'
# This command prints trajectory metadata: frames count, time length and atoms count
# More info at https://github.com/jag1g13/xtc-length
# WARNING: Frames count is always correct but time length may be wrong so do not use it
def get_frames_count (input_trajectory_filename : str) -> int:
    
    # Run xtc-length
    logs = run([
        "xtc-length",
        input_trajectory_filename,
    ], stdout=PIPE).stdout.decode()

    # Get the frames number from the xtc-length logs
    match = re.search("Trajectory contains ([0-9.]+) frames", logs)
    if match:
        raw_frames = match.group(1)
        frames = int(raw_frames.replace('.',''))
        return frames

    # If nothing is found there must be an error
    # Print logs and stop here
    raise NameError('Error when counting frames')