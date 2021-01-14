import glob
from subprocess import run, PIPE

# Remove gromacs back ups
def remove_trash():

    # The subprocess.run command does not deal with the '*' bash syntax to grab multiple files
    # Instead, we have to do it this way
    trash_files = glob.glob('#*')

    if len(trash_files) > 0:

        logs = run([
            "rm",
            *trash_files,
        ], stdout=PIPE).stdout.decode()