from os import environ
from pathlib import Path

from model_workflow.utils.file import File

# Replace the gromacs masses file (atommass.dat) by our custom file
# This file includes all atoms also in all caps letters
# e.g. Zn and ZN
# This way we avoid gormacs not finding atoms because names are in caps

# Set the path to our custom masses file
resources = str(Path(__file__).parent.parent / "resources")
custom_masses_file = File(resources + '/atommass.dat')

# Set the path to the gromacs masses file
gromacs_masses_directory = environ['CONDA_PREFIX'] + '/share/gromacs/top'
gromacs_masses_file = File(gromacs_masses_directory + '/atommass.dat')

# Set the backup filename
# We keep the original atom masses file just in case
backup_file = File(gromacs_masses_directory + '/backup_atommass.dat')

# Replace the original file by a symlink to our custom file if it is not done yet
def fix_gromacs_masses ():

    # Check if the gromacs masses file is already a symlink
    # If so, then it means we already did the fix so we can stop here
    if gromacs_masses_file.is_symlink():
        return

    # Copy the original file to make a backup
    gromacs_masses_file.rename_to(backup_file)

    # Now make a symlink to our custom masses file
    gromacs_masses_file.set_symlink_to(custom_masses_file)