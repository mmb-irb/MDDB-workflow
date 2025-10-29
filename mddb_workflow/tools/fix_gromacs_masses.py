from os import environ
from pathlib import Path

from mddb_workflow.utils.file import File

# Replace the default gromacs masses file (atommass.dat) by our custom masses file
# This file includes relevant atoms also in all caps letters
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

# LORE: This was done by modifying the conda enviornment before
# However this lead to problems when containerizing the environment
# In one hand, there was no CONDA_PREFIX environmental variable
# But the real problem was that writting in the enviornment was not allowed
# Otherwise the container could not be shared between different user in a cluster

# Set the absolute path to the workflow resources directory
resources = str(Path(__file__).parent.parent / "resources")

# Replace the original file by a symlink to our custom file if it is not done yet
def fix_gromacs_masses ():

    # Set the source file
    # WARNING: Note that this must be done here, not outside the function
    # Otherwise a workflow called with a '-dir' parameter would have a wrong relative path
    source_custom_masses_file = File(resources + '/atommass.dat')

    # Set the path to a local copy of the workflow custom masses file
    # According to Justin Lemkul:
    # "All GROMACS programs will read relevant database files from the working directory 
    # before referencing them from $GMXLIB."
    local_custom_masses_file = File('atommass.dat')

    # Check if the backup file exists and, if not, then rename the current masses file as the backup
    if not local_custom_masses_file.exists:
        # If it does not exists then it means it is a symlink pointing to a not valid direction
        # This may happend when working in different machines
        # Simply remove the old symlink
        if local_custom_masses_file.is_symlink(): local_custom_masses_file.remove()
        local_custom_masses_file.set_symlink_to(source_custom_masses_file)