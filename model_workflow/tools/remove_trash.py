from glob import glob
from os import remove
from os.path import exists

# Remove trash files
def remove_trash (directory : str):
    # Remove this time.txt file which is beeing created sometimes and I don't know where it comes from
    trash = [ directory + '/time.txt' ]
    # Find gromacs backups
    trash += glob(directory + '/#*')
    # Find reduced trajectories
    trash += glob(directory + '/f*.trajectory.xtc')
    # Remove each trash file
    for filepath in trash:
        if exists(filepath):
            remove(filepath)
