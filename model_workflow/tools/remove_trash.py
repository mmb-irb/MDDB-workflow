import glob
from subprocess import run, PIPE

# Remove gromacs back ups
def remove_trash():

    # The subprocess.run command does not deal with the '*' bash syntax to grab multiple files
    # Instead, we have to do it this way
    gromacs_backups = glob.glob('#*')
    reduced_trajectories = glob.glob('f*.trajectory.xtc')

    if len(gromacs_backups) > 0:

        logs = run([
            "rm",
            *gromacs_backups,
        ], stdout=PIPE).stdout.decode()

    if len(reduced_trajectories) > 0:

        logs = run([
            "rm",
            *reduced_trajectories,
        ], stdout=PIPE).stdout.decode()