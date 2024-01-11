import glob
import os

# Remove trash files
def remove_trash (md_directory : str):
    trash = []
    # Find gromacs backups
    trash += glob.glob(md_directory + '/#*')
    # Find reduced trajectories
    trash += glob.glob(md_directory + '/f*.trajectory.xtc')
    # Remove each trash file
    for filepath in trash:
        os.remove(filepath)
