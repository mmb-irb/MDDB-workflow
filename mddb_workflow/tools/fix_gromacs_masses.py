from mddb_workflow.utils.constants import GROMACS_CUSTOM_MASSES_FILEPATH
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

# Replace the original file by a symlink to our custom file if it is not done yet
def fix_gromacs_masses ():

    # Set the source file
    # WARNING: Note that this must be done here, not outside the function
    # Otherwise a workflow called with a '-dir' parameter would have a wrong relative path
    source_custom_masses_file = File(GROMACS_CUSTOM_MASSES_FILEPATH)

    # Set the path to a local copy of the workflow custom masses file
    # According to Justin Lemkul:
    # "All GROMACS programs will read relevant database files from the working directory
    # before referencing them from $GMXLIB."
    local_custom_masses_file = File('atommass.dat')

    # This was a symlink before
    # Make sure any reamining symlinks are removed
    if local_custom_masses_file.is_symlink(): local_custom_masses_file.remove()

    # Check if the backup file exists and, if not, then copy the reference
    if not local_custom_masses_file.exists:
        source_custom_masses_file.copy_to(local_custom_masses_file)

# Set the symbol in the atommass file representing 'any residue name'
ANY_RESIDUE_NAME = '???'
# Set a header for MWF extended masses
MWF_HEADER = '; MWF extension\n'

# Extend masses in the gromacs file
# New masses is a list of tuples with 3 values: residue name (optional), atom name, mass
def extend_gromacs_masses (new_masses : set[tuple]):

    # Get the local custom masses file
    local_custom_masses_file = File('atommass.dat')
    # If the file does not exist yet then create it
    if not local_custom_masses_file.exists:
        fix_gromacs_masses()
    
    # Read masses already listed in the file
    already_modified = False
    current_masses = {}
    with open(local_custom_masses_file.path, 'r') as file:
        for line in file.readlines():
            if line == MWF_HEADER: already_modified = True
            # Skip comment lines
            if line[0] == ';': continue
            # Skip empty lines
            if line == '\n': continue
            # Mine the mass
            residue_name, atom_name, mass = line.split()
            if residue_name == ANY_RESIDUE_NAME:
                residue_name = None
            current_masses[(residue_name, atom_name)] = mass
    
    # Add the new masses
    with open(local_custom_masses_file.path, 'a') as file:
        # Add a header for our own masses
        # DANI: Si alguien más modifica el atommass.dat esta sección se mezclará
        if not already_modified: file.write('; MWF extension\n')
        # Iterate new masses
        for new_mass in new_masses:
            # Get the new mass values
            residue_name, atom_name, mass = new_mass
            if not residue_name: residue_name = None
            # If we already have a mass value for this residue/atom combination then skip it
            atom_config = (residue_name, atom_name)
            if atom_config in current_masses: continue
            # Add the new mass
            if residue_name is None: residue_name = ANY_RESIDUE_NAME
            file.write(f'{residue_name} {atom_name} {mass}\n')