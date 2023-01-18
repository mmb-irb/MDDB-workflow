import glob
import os

# Remove trash files
def remove_trash():
    trash = []
    # Find gromacs backups
    trash += glob.glob('#*')
    # Find reduced trajectories
    trash += glob.glob('f*.trajectory.xtc')
    # Find the last frame file
    last_frame_filename = 'last_frame.pdb'
    if os.path.exists(last_frame_filename):
        trash.append(last_frame_filename)
    # Remove each trash file
    for filename in trash:
        os.remove(filename)
