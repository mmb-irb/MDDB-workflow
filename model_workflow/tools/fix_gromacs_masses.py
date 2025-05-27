from os import environ
from pathlib import Path

from model_workflow.utils.file import File

# Replace the gromacs masses file (atommass.dat) by our custom file
# This file includes all atoms also in all caps letters
# e.g. Zn and ZN
# This way we avoid gormacs not finding atoms because names are in caps

# LORE: This was done with a symlink before
# In the gromacs data directory from the conda enviornment we created a relative symlink to the workflow atommass.dat
# This worked fine until we had a very specific issue with a cluster (gpucluster at IRB Barcelona)
# In this cluster there was two different paths to reach the home directory: /home/username and /orozco/homes/username
# The python __file__ values were using the /orozco/homes/username root despite we were calling it from /home/username
# Thus relative paths were all passing through the root and the different number of jumps was breaking the paths
# I spent half a day trying to get __file__ values using the /home/username root and I did not succeed
# Then we decided to just write the masses file in the conda enviornment every time

# Set the path to our custom masses file
resources = str(Path(__file__).parent.parent / "resources")
custom_masses_file = File(resources + '/atommass.dat')

# Set the path to the gromacs masses file
gromacs_masses_directory = environ.get('CONDA_PREFIX','') + '/share/gromacs/top'
gromacs_masses_file = File(gromacs_masses_directory + '/atommass.dat')

# Set the backup filename
# We keep the original atom masses file just in case
backup_file = File(gromacs_masses_directory + '/backup_atommass.dat')

# Replace the original file by a symlink to our custom file if it is not done yet
def fix_gromacs_masses ():

    # Check if the backup file exists and, if not, then rename the current masses file as the backup
    if not backup_file.exists:
        gromacs_masses_file.rename_to(backup_file)

    # Always copy our custom masses file to the enviornment just in case it was modified since last time
    custom_masses_file.copy_to(gromacs_masses_file)