# This is the starter script

# Import python libraries
import os
from pathlib import Path
import urllib.request
import json

# Import external analysis tools
import pytraj as pt

# Import local tools
from vmd_processor import processor
from image_and_fit import image_and_fit
from topology_manager import TopologyReference
from topology_corrector import topology_corrector
from get_first_frame import get_first_frame
from get_backbone import get_backbone
from get_summarized_trajectory import get_summarized_trajectory
from get_frames_count import get_frames_count
from get_box_size import get_box_size
from get_atoms_count import get_atoms_count

# Import local analyses
from generic_analyses import rmsd, rmsf, rgyr
from pca import pca
from rmsd_per_residue import rmsd_per_residue
from rmsd_pairwise import rmsd_pairwise
from hydrogen_bonds import hydrogen_bonds
from energies import energies
from pockets import pockets

# General inputs ----------------------------------------------------------------------------

# Main topology and trajectory output file names
# WARNING: Changing these names here means making changes in all the project
# i.e. update the api, client and already loaded projects in the database at least
topology_filename = "md.imaged.rot.dry.pdb"
trajectory_filename = "md.imaged.rot.xtc"

# Set the inputs filename
# This file may be easily generated using the 'input_setter' notebook in the 'dev' directory
inputs_filename = "inputs.json"

# Set this workflow to skip steps where the ouput file already exist
# Change it to False if you want all steps to be done anyway
skip_repeats = True

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
    # It is used further in the RMSd analysis
    first_frame_filename = 'firstFrame.pdb'
    if required(first_frame_filename):
        get_first_frame(topology_filename, topology_filename, first_frame_filename)

    # Get the backbone structure
    # It is further loaded to database and used to represent PCA projections
    backbone_filename = 'backbone.pdb'
    if required(backbone_filename):
        get_backbone(topology_reference, backbone_filename)

    # Get the sumarized trajectory
    # It is used further in some trajectory analyses
    summarized_trajectory = 'md.imaged.rot.100.xtc'
    if required(summarized_trajectory):
        get_summarized_trajectory(topology_filename, trajectory_filename, summarized_trajectory)

    # Interfaces definition --------------------------------------------------------------------

    print('Processing interfaces')

    # Read the defined interfaces from the inputs file
    interfaces = inputs['interfaces']

    # Find out all residues and interface residues for each interface 'agent'
    # Interface residues are defined as by a cuttoff distance in Angstroms
    cutoff_distance = 5
    for interface in interfaces:
        # residues_1 is the list of all residues in agent_1
        interface['residues_1'] = topology_reference.topology_selection(
            interface['agent_1']
        )
        # residues_2 is the list of all residues in agent_2
        interface['residues_2'] = topology_reference.topology_selection(
            interface['agent_2']
        )
        # interface_1 is the list of residues from agent_1 which are close to the agent_2
        interface['interface_1'] = topology_reference.topology_selection(
            interface['agent_1'] +
            ' and same residue as exwithin ' +
            str(cutoff_distance) + 
            ' of ' + 
            interface['agent_2'])
        # interface_2 is the list of residues from agent_2 which are close to the agent_1
        interface['interface_2'] = topology_reference.topology_selection(
            interface['agent_2'] +
            ' and same residue as exwithin ' +
            str(cutoff_distance) + 
            ' of ' + 
            interface['agent_1'])

        # Set a paralel pytraj interfaces dict
        # These values are used along the workflow but not added to metadata
        # Translate all selections to pytraj residue notation
        interface.update(
            {
                'pt_residues_1': list(map(topology_reference.source2pytraj, interface['residues_1'])),
                'pt_residues_2': list(map(topology_reference.source2pytraj, interface['residues_2'])),
                'pt_interface_1': list(map(topology_reference.source2pytraj, interface['interface_1'])),
                'pt_interface_2': list(map(topology_reference.source2pytraj, interface['interface_2'])),
            }
        )

        print('1 -> ' + str(interface['interface_1'] + interface['interface_2']))

    # Metadata mining --------------------------------------------------------------------------

    print('Mining metadata')

    # Count the number of snapshots
    snapshots = get_frames_count(trajectory_filename)

    # Find out the box size (x, y and z)
    (boxsizex, boxsizey, boxsizez) = get_box_size(topology_filename, trajectory_filename)

    # Count different type of atoms and residues
    (systats, protats, prot, dppc, sol, na, cl) = get_atoms_count(topology_filename)

    # Extract some additional metadata from the inputs file which is required further
    ligands = inputs['ligands']

    # Set the metadata interfaces
    metadata_interfaces = [{
        'name': interface['name'],
        'interface': str(interface['interface_1'] + interface['interface_2']),
    } for interface in interfaces]

    # Write the metadata file
    # Metadata keys must be in caps, as they are in the client
    metadata = {
        'PDBID': inputs['pdbId'],
        'NAME': inputs['name'],
        'UNIT': inputs['unit'],
        'DESCRIPTION': inputs['description'],
        'AUTHORS': inputs['authors'],
        'PROGRAM': inputs['program'],
        'VERSION': inputs['version'],
        'LICENSE': inputs['license'],
        'LINKCENSE': inputs['linkcense'],
        'CITATION': inputs['citation'],
        'LENGTH': inputs['length'],
        'TIMESTEP': inputs['timestep'],
        'SNAPSHOTS': snapshots,
        'FREQUENCY': inputs['frequency'],
        'FF': inputs['ff'],
        'TEMP': inputs['temp'],
        'WAT': inputs['wat'],
        'BOXTYPE': inputs['boxtype'],
        'BOXSIZEX': boxsizex,
        'BOXSIZEY': boxsizey,
        'BOXSIZEZ': boxsizez,
        'ENSEMBLE': inputs['ensemble'],
        'PCOUPLING': inputs['pcoupling'],
        'MEMBRANE': inputs['membrane'],
        'SYSTATS': systats,
        'PROTATS': protats,
        'PROT': prot,
        'DPPC': dppc,
        'SOL': sol,
        'NA': na,
        'CL': cl,
        'LIGANDS': inputs['ligands'],
        'DOMAINS': inputs['domains'],
        'INTERFACES': metadata_interfaces,
        'CHAINNAMES': inputs['chainnames'],
    }
    metadata_filename = 'metadata.json'
    with open(metadata_filename, 'w') as file:
        json.dump(metadata, file)

    # Analyses ---------------------------------------------------------------------------------

    print('Running analyses')

    # Set the RMSd analysis file name and run the analysis
    rmsd_analysis = 'md.rmsd.xvg'
    if required(rmsd_analysis):
        rmsd(first_frame_filename, trajectory_filename, rmsd_analysis)

    # Set the fluctuation analysis file name and run the analysis
    rmsf_analysis = 'md.rmsf.xvg'
    if required(rmsf_analysis):
        rmsf(topology_filename, trajectory_filename, rmsf_analysis)

    # Set the RMSd per resiude analysis file name and run the analysis
    rgyr_analysis = 'md.rgyr.xvg'
    if required(rgyr_analysis):
        rgyr(topology_filename, trajectory_filename, rgyr_analysis)

    # Set the pca output filenames and run the analysis
    # WARNING: This analysis will generate several output files
    # Files 'average.pdb' and 'covar.log' are generated by the PCA but never used
    eigenvalues_filename = 'pca.eigenval.xvg'
    eigenvectors_filename = 'eigenvec.trr'
    if required(eigenvalues_filename) or required(eigenvectors_filename):
        pca(topology_filename, trajectory_filename, eigenvalues_filename, eigenvectors_filename, snapshots)

    # Set the pytraj trayectory, which is further used in all pytraj analyses
    pt_trajectory = pt.iterload(trajectory_filename, topology_filename)
    reduced_pt_trajectory = pt_trajectory[0:2000:10]

    # Set the RMSd per resiude analysis file name and run the analysis
    rmsd_perres_analysis = 'md.rmsd.perres.xvg'
    if required(rmsd_perres_analysis):
        rmsd_per_residue(reduced_pt_trajectory, rmsd_perres_analysis, topology_reference)

    # Set the RMSd pairwise analysis file name and run the analysis
    rmsd_pairwise_analysis = 'md.rmsd.pairwise.json'
    if required(rmsd_pairwise_analysis):
        rmsd_pairwise(reduced_pt_trajectory, rmsd_pairwise_analysis, interfaces)

    # Set the hydrogen bonds analysis file name and run the analysis
    hbonds_analysis = 'md.hbonds.json'
    if required(hbonds_analysis) and len(interfaces) > 0:
        hydrogen_bonds(reduced_pt_trajectory, hbonds_analysis, topology_reference, interfaces)

    # Set the energies analysis filename and run the analysis
    energies_analysis = 'md.energies.json'
    if required(energies_analysis) and len(ligands) > 0:
        energies(topology_filename, trajectory_filename, energies_analysis, topology_reference, snapshots, ligands)

    # Set the pockets analysis filename and run the analysis
    pockets_analysis = 'md.pockets.json'
    if required(pockets_analysis):
        pockets(topology_filename, trajectory_filename, pockets_analysis, topology_reference, snapshots)

    print('Done!')

# Set a function to check if a process must be run (True) or skipped (False)
# i.e. check if the output file already exists and reapeated analyses must be skipped
def required (analysis_filename : str):
    if os.path.exists(analysis_filename) and skip_repeats:
        return False
    return True