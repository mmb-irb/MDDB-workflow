import glob
import os
from subprocess import run, PIPE

# Remove trash files
def remove_trash():

    # The subprocess.run command does not deal with the '*' bash syntax to grab multiple files
    # Instead, we have to do it this way
    gromacs_backups = glob.glob('#*')
    reduced_trajectories = glob.glob('f*.trajectory.xtc')

    # Remove gromacs backups
    if len(gromacs_backups) > 0:

        logs = run([
            "rm",
            *gromacs_backups,
        ], stdout=PIPE).stdout.decode()

    # Remove reduced trajectories
    if len(reduced_trajectories) > 0:

        logs = run([
            "rm",
            *reduced_trajectories,
        ], stdout=PIPE).stdout.decode()

    # Remove the 'restart' cmip file
    restart_filename = 'restart'
    if os.path.exists(restart_filename):

        logs = run([
            "rm",
            restart_filename,
        ], stdout=PIPE).stdout.decode()
