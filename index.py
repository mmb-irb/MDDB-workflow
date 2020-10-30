# This is the starter script

# Import python libraries
import os
from pathlib import Path
import urllib.request
import json

# Import tools
from vmd_processor import processor
from image_and_fit import image_and_fit
from topology_manager import TopologyReference
from topology_corrector import topology_corrector
from get_first_frame import get_first_frame
from get_summarized_trajectory import get_summarized_trajectory
from get_box_size import get_box_size
from get_atoms_count import get_atoms_count

# Import analyses
from rmsd_per_residue import rmsd_per_residue

# General inputs ----------------------------------------------------------------------------

# Main topology and trajectory output file names
# WARNING: Changing these names here means making changes in all the project
# i.e. update the api, client and already loaded projects in the database at least
topology_filename = "md.imaged.rot.dry.pdb"
trajectory_filename = "md.imaged.rot.xtc"

# Set the inputs filename
inputs_filename = "inputs.json"

# Set this workflow to skip steps where the ouput file already exist
# Change it to False if you want all steps to be done anyway
skip_repeats = True

# Set this workflow to display 3D molecular representations and analysis graphs
# Change it to False to skip renders making the workflow faster and memory cheaper
displayResults = True

# ------------------------------------------------------------------------------------------

# To analyze a local directory which is meant to include the topology and trajetory files
# 'directory' is the string path to the working directory. e.g. '/home/dbeltran/Desktop/my_directory/'
def analyze_directory (directory : str):
    # Move to the specified directory
    os.chdir(directory)
    # Run the analyses
    run_analyses()

# To download topology and trajectory files from an already uploaded project
# Also the directory where both the downloaded files and the produced files will be stored
# 'project' is a string code to the project to analyze. e.g. 'MCV1900002' or '5e95cc9466d570732286ef90'
# 'directory' is the string path to the working directory. e.g. '/home/dbeltran/Desktop/my_directory/'
def analyze_project (project : str, directory : str, url = 'https://bioexcel-cv19.bsc.es'):
    # Create the directory if it does not exists
    if os.path.exists(directory) == False :
        Path(directory).mkdir(parents=True, exist_ok=True)
    # Move to the specified directory
    os.chdir(directory)
    # Download the topology file if it does not exists
    if os.path.exists(topology_filename) == False :
        print('Downloading topology')
        topology_url = url + '/api/rest/current/projects/' + project + '/files/' + topology_filename
        urllib.request.urlretrieve (topology_url, topology_filename)
    # Download the trajectory file if it does not exists
    if os.path.exists(trajectory_filename) == False :
        print('Downloading trajectory')
        trajectory_url = url + '/api/rest/current/projects/' + project + '/files/' + trajectory_filename
        urllib.request.urlretrieve (trajectory_url, trajectory_filename)
    # Run the analyses
    run_analyses()

# Run all analyses with the provided topology and trajectory files
def run_analyses ():

    # Load the inputs file
    with open(inputs_filename, 'r') as file:
        inputs = json.load(file)

    # Get the input topology and trajectory filenames
    original_topology_filename = inputs['original_topology_filename']
    original_trajectory_filename = inputs['original_trajectory_filename']

    # Preprocessing ---------------------------------------------------------------------------------
    
    print('Preprocessing')

    # Process the topology and or trajectory files using VMD
    # Files are converted to supported formats and trajectory pieces are merged into a single file
    # In addition, some irregularities in the topology may be fixed by VMD
    # If the output topology and trajectory files already exists it is assumed they are already processed
    if required(topology_filename) or required(trajectory_filename):
        logs = processor(
            original_topology_filename,
            original_trajectory_filename,
            topology_filename,
            trajectory_filename,
        )

    # Image the trajectory if it is required
    # i.e. make the trajectory uniform avoiding atom jumps and making molecules to stay whole
    # Fit the trajectory by removing the translation and rotation if it is required
    preprocess_protocol = inputs['preprocess_protocol']
    if preprocess_protocol > 0:
        image_and_fit(topology_filename, trajectory_filename, trajectory_filename, preprocess_protocol)

    # Examine and correct the topology file using ProDy
    topology_corrector(topology_filename, topology_filename)
    
    # Create an object with the topology data in both ProDy and Pytraj formats
    # This object also include functions to convert residue numeration from one format to another
    topology_reference = TopologyReference(topology_filename)

    # Get the first trajectory frame
    first_frame = 'firstFrame.pdb'
    if required(first_frame):
        get_first_frame(topology_filename, topology_filename, first_frame)

    # Get the sumarized trajectory
    summarized_trajectory = 'md.imaged.rot.100.xtc'
    if required(summarized_trajectory):
        get_summarized_trajectory(topology_filename, trajectory_filename, summarized_trajectory)

    # Metadata mining --------------------------------------------------------------------------

    print('Mining metadata')

    # Find out the box size (x, y and z)
    (boxsizex, boxsizey, boxsizez) = get_box_size(topology_filename, trajectory_filename)

    # Count different type of atoms and residues
    (systats, protats, prot, dppc, sol, na, cl) = get_atoms_count(topology_filename)

    # Write the metadata file
    metadata = {
        'pdbId': inputs['pdbId'],
        'name': inputs['name'],
        'unit': inputs['unit'],
        'description': inputs['description'],
        'authors': inputs['authors'],
        'program': inputs['program'],
        'version': inputs['version'],
        'license': inputs['license'],
        'linkcense': inputs['linkcense'],
        'citation': inputs['citation'],
        'length': inputs['length'],
        'timestep': inputs['timestep'],
        'snapshots': inputs['snapshots'],
        'frequency': inputs['frequency'],
        'ff': inputs['ff'],
        'temp': inputs['temp'],
        'wat': inputs['wat'],
        'boxType': inputs['boxtype'],
        'boxSizeX': boxsizex,
        'boxSizeY': boxsizey,
        'boxSizeZ': boxsizez,
        'ensemble': inputs['ensemble'],
        'pcoupling': inputs['pcoupling'],
        'membrane': inputs['membrane'],
        'systats': systats,
        'protats': protats,
        'prot': prot,
        'dppc': dppc,
        'sol': sol,
        'na': na,
        'cl': cl,
        'ligands': inputs['ligands'],
        'domains': inputs['domains'],
        'interfaces': inputs['interfaces'],
        'chainnames': inputs['chainnames'],
    }
    metadata_filename = 'metadata.json'
    with open(metadata_filename, 'w') as file:
        json.dump(metadata, file)

    # Analyses ---------------------------------------------------------------------------------

    print('Running analyses')

    # Set the RMSd per resiude analysis file name
    rmsd_perres_analysis = 'md.rmsd.perres.xvg'
    if required(rmsd_perres_analysis):
        rmsd_per_residue(topology_filename, topology_filename, rmsd_perres_analysis, topology_reference)

    print('Done!')

# Set a function to check if a process must be run (True) or skipped (False)
# i.e. check if the output file already exists and reapeated analyses must be skipped
def required (analysis_filename : str):
    if os.path.exists(analysis_filename) and skip_repeats:
        return False
    return True