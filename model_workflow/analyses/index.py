# This is the starter script
# Import python libraries
import os
import math
from pathlib import Path
import urllib.request
import json

# Import external analysis tools
import pytraj as pt

# Import local tools
from model_workflow.tools.vmd_processor import vmd_processor
from model_workflow.tools.image_and_fit import image_and_fit
from model_workflow.tools.topology_manager import TopologyReference
from model_workflow.tools.topology_corrector import topology_corrector
from model_workflow.tools.get_first_frame import get_first_frame
from model_workflow.tools.get_backbone import get_backbone
from model_workflow.tools.get_average import get_average
from model_workflow.tools.get_summarized_trajectory import get_summarized_trajectory
from model_workflow.tools.get_frames_count import get_frames_count
from model_workflow.tools.get_box_size import get_box_size
from model_workflow.tools.get_atoms_count import get_atoms_count
from model_workflow.tools.remove_trash import remove_trash

# Import local analyses
from model_workflow.analyses.rmsds import rmsds
from model_workflow.analyses.generic_analyses import rmsf, rgyr
from model_workflow.analyses.pca import pca
from model_workflow.analyses.pca_contacts import pca_contacts
from model_workflow.analyses.rmsd_per_residue import rmsd_per_residue
from model_workflow.analyses.rmsd_pairwise import rmsd_pairwise
from model_workflow.analyses.distance_per_residue import distance_per_residue
from model_workflow.analyses.hydrogen_bonds import hydrogen_bonds
from model_workflow.analyses.energies import energies
from model_workflow.analyses.pockets import pockets

# Set this workflow to skip steps where the ouput file already exist
# Change it to False if you want all steps to be done anyway
skip_repeats = True

# ------------------------------------------------------------------------------------------

# Set a function to check if a process must be run (True) or skipped (False)
# i.e. check if the output file already exists and reapeated analyses must be skipped


def required(analysis_filename: str, skip_repeats=True):
    if os.path.exists(analysis_filename) and skip_repeats:
        return False
    return True

# ------------------------------------------------------------------------------------------


# The OUTPUT filenames. They are not customizable since the loader uses these names
# Thus, these are also the MoDEL names of the files in the database
topology_filename = "md.imaged.rot.dry.pdb"
trajectory_filename = "md.imaged.rot.xtc"
# trajectory_filename = "md.trr"

current_directory = os.getcwd()


# The inputs filename is input and output at the same time
# Input filename if it is already in the directory
# Output filename if you download it from MoDEL
def start(
        directory: str = current_directory,
        project: str = None,
        url: str = 'https://bioexcel-cv19-dev.bsc.es',
        inputs_filename: str = 'inputs.json'):
    # Create the directory if it does not exists
    if not os.path.exists(directory):
        Path(directory).mkdir(parents=True, exist_ok=True)
    # Move to the specified directory
    os.chdir(directory)
    if project:
        # Download the topology file if it does not exists
        if not os.path.exists(topology_filename):
            print('Downloading topology')
            topology_url = url + '/api/rest/current/projects/' + \
                project + '/files/' + topology_filename
            urllib.request.urlretrieve(topology_url, topology_filename)
        # Download the trajectory file if it does not exists
        if not os.path.exists(trajectory_filename):
            print('Downloading trajectory')
            trajectory_url = url + '/api/rest/current/projects/' + \
                project + '/files/' + trajectory_filename
            urllib.request.urlretrieve(trajectory_url, trajectory_filename)
        # Download the inputs json file if it does not exists
        if not os.path.exists(inputs_filename):
            print('Downloading inputs')
            inputs_url = url + '/api/rest/current/projects/' + \
                project + '/inputs/'
            urllib.request.urlretrieve(inputs_url, inputs_filename)

# Run all analyses with the provided topology and trajectory files


def analysis_prep(
        topology_filename="md.imaged.rot.dry.pdb",
        trajectory_filename="md.imaged.rot.xtc",
        inputs_filename="inputs.json",
        interface_cutoff_distance=5):
    # Load the inputs file
    with open(inputs_filename, 'r') as file:
        inputs = json.load(file)

    # Set a function to retrieve 'inputs' values and handle missing keys
    def getInput(input: str):
        return inputs.get(input, None)

    # Get the input topology and trajectory filenames
    original_topology_filename = getInput('original_topology_filename')
    original_trajectory_filename = getInput('original_trajectory_filename')

    # Preprocessing ---------------------------------------------------------------------------------

    print('Preprocessing')

    # Process the topology and or trajectory files using VMD
    # Files are converted to supported formats and trajectory pieces are merged into a single file
    # In addition, some irregularities in the topology may be fixed by VMD
    # If the output topology and trajectory files already exists it is assumed they are already processed
    if required(topology_filename) or required(trajectory_filename):
        logs = vmd_processor(
            original_topology_filename,
            original_trajectory_filename,
            topology_filename,
            trajectory_filename,
        )

    # Image the trajectory if it is required
    # i.e. make the trajectory uniform avoiding atom jumps and making molecules to stay whole
    # Fit the trajectory by removing the translation and rotation if it is required
    preprocess_protocol = getInput('preprocess_protocol')
    if preprocess_protocol > 0:
        image_and_fit(topology_filename, trajectory_filename,
                      trajectory_filename, preprocess_protocol)

    # Examine and correct the topology file using ProDy
    topology_corrector(topology_filename, topology_filename)

    # Create an object with the topology data in both ProDy and Pytraj formats
    # This object also include functions to convert residue numeration from one format to another
    topology_reference = TopologyReference(topology_filename)

    # Pytraj setup ---------------------------------------------------------------------------------

    # Set the pytraj trayectory, which is further used in all pytraj analyses
    pt_trajectory = pt.iterload(trajectory_filename, topology_filename)
    trajectory_frames = pt_trajectory.n_frames

    # Set a reduced trajectory used for heavy analyses
    reduced_pt_trajectory = None
    # First, set the maximum number of frames for the reduced trajectory
    # WARNING: The final number of frames in some analyses may be +1
    reduced_trajectory_frames = 200
    # If the current trajectory has already less frames than the maximum then use it also as reduced
    if trajectory_frames < reduced_trajectory_frames:
        reduced_pt_trajectory = pt_trajectory
        # Add a step value which will be required later
        reduced_pt_trajectory.step = 1
    # Otherwise, create a reduced trajectory with as much frames as specified above
    # These frames are picked along the trajectory
    else:
        # Calculate how many frames we must jump between each reduced frame to never exceed the limit
        # The '- 1' is because the first frame is 0 (you have to do the math to understand)
        step = math.floor(trajectory_frames / (reduced_trajectory_frames - 1))
        reduced_pt_trajectory = pt_trajectory[0:trajectory_frames:step]
        # DANI: hay que chequear, porque sis siempre son 201 frames el -1 de arriba no tiene sentido
        #print(reduced_pt_trajectory.n_frames)
        # Add the step value to the reduced trajectory explicitly. It will be required later
        reduced_pt_trajectory.step = step

    # Additional refinement --------------------------------------------------------------------

    # Get the first trajectory frame
    # It is used further in the RMSd analysis
    first_frame_filename = 'firstFrame.pdb'
    if required(first_frame_filename):
        get_first_frame(topology_filename, topology_filename,
                        first_frame_filename)

    # Get the backbone structure
    # It is further loaded to database and used to represent PCA projections
    backbone_filename = 'backbone.pdb'
    if required(backbone_filename):
        get_backbone(topology_reference, backbone_filename)

    # Get the average structure in frame format
    # It is further loaded to database and used to represent pockets
    average_structure_filename = 'average.pdb'
    if required(average_structure_filename):
        get_average(pt_trajectory, average_structure_filename)

    # Get the average structure in frame format
    # It is further loaded to database and used to represent pockets
    average_frame_filename = 'average.xtc'
    if required(average_frame_filename):
        get_average(pt_trajectory, average_frame_filename)

    # Interactions setup --------------------------------------------------------------------

    print('Processing interactions')

    # Read the defined interactions from the inputs file
    interactions = getInput('interactions')

    if interactions:
    # Find out all residues and interface residues for each interaction 'agent'
    # Interface residues are defined as by a cuttoff distance in Angstroms
        cutoff_distance = 5
        for interaction in interactions:
            # residues_1 is the list of all residues in the first agent
            interaction['residues_1'] = topology_reference.topology_selection(
                interaction['selection_1']
            )
            # residues_2 is the list of all residues in the second agent
            interaction['residues_2'] = topology_reference.topology_selection(
                interaction['selection_2']
            )
            # interface_1 is the list of residues from the agent 1 which are close to the agent 2
            interaction['interface_1'] = topology_reference.topology_selection(
                '(' + interaction['selection_1'] +
                ') and same residue as exwithin ' +
                str(cutoff_distance) +
                ' of (' +
                interaction['selection_2'] + ')')
            # interface_2 is the list of residues from agent 2 which are close to the agent 1
            interaction['interface_2'] = topology_reference.topology_selection(
                '(' + interaction['selection_2'] +
                ') and same residue as exwithin ' +
                str(cutoff_distance) +
                ' of (' +
                interaction['selection_1'] + ')')

            # Translate all residues selections to pytraj notation
            # These values are used along the workflow but not added to metadata
            interaction.update(
                {
                    'pt_residues_1': list(map(topology_reference.source2pytraj, interaction['residues_1'])),
                    'pt_residues_2': list(map(topology_reference.source2pytraj, interaction['residues_2'])),
                    'pt_interface_1': list(map(topology_reference.source2pytraj, interaction['interface_1'])),
                    'pt_interface_2': list(map(topology_reference.source2pytraj, interaction['interface_2'])),
                }
            )

            print(
                '1 -> ' + str(interaction['interface_1'] + interaction['interface_2']))

    else:
        interactions = []

    # Metadata mining --------------------------------------------------------------------------

    print('Mining metadata')

    # Count the number of snapshots
    snapshots = get_frames_count(topology_filename, trajectory_filename)

    # Find out the box size (x, y and z)
    (boxsizex, boxsizey, boxsizez) = get_box_size(
        topology_filename, trajectory_filename)

    # Count different type of atoms and residues
    (systats, protats, prot, dppc, sol, na,
     cl) = get_atoms_count(topology_filename)

    # Extract some additional metadata from the inputs file which is required further
    ligands = getInput('ligands')

    if not ligands:
        ligands = []

    # Set the metadata interactions
    metadata_interactions = [{
        'name': interaction['name'],
        'agent_1': interaction['agent_1'],
        'agent_2': interaction['agent_2'],
        'selection_1': interaction['selection_1'],
        'selection_2': interaction['selection_2'],
        'residues_1': [ str(residue) for residue in interaction['residues_1'] ],
        'residues_2': [ str(residue) for residue in interaction['residues_2'] ],
        'interface_1': [ str(residue) for residue in interaction['interface_1'] ],
        'interface_2': [ str(residue) for residue in interaction['interface_2'] ],
    } for interaction in interactions]

    # Write the metadata file
    # Metadata keys must be in CAPS, as they are in the client
    metadata = {
        'PDBID': getInput('pdbId'),
        'NAME': getInput('name'),
        'UNIT': getInput('unit'),
        'DESCRIPTION': getInput('description'),
        'AUTHORS': getInput('authors'),
        'GROUPS': getInput('groups'),
        'PROGRAM': getInput('program'),
        'VERSION': getInput('version'),
        'LICENSE': getInput('license'),
        'LINKCENSE': getInput('linkcense'),
        'CITATION': getInput('citation'),
        'LENGTH': getInput('length'),
        'TIMESTEP': getInput('timestep'),
        'SNAPSHOTS': snapshots,
        'FREQUENCY': getInput('frequency'),
        'FF': getInput('ff'),
        'TEMP': getInput('temp'),
        'WAT': getInput('wat'),
        'BOXTYPE': getInput('boxtype'),
        'BOXSIZEX': boxsizex,
        'BOXSIZEY': boxsizey,
        'BOXSIZEZ': boxsizez,
        'ENSEMBLE': getInput('ensemble'),
        'PCOUPLING': getInput('pcoupling'),
        'MEMBRANE': getInput('membrane'),
        'SYSTATS': systats,
        'PROTATS': protats,
        'PROT': prot,
        'DPPC': dppc,
        'SOL': sol,
        'NA': na,
        'CL': cl,
        'LIGANDS': ligands,
        'DOMAINS': getInput('domains'),
        'INTERACTIONS': metadata_interactions,
        'CHAINNAMES': getInput('chainnames'),
    }
    metadata_filename = 'metadata.json'
    with open(metadata_filename, 'w') as file:
        json.dump(metadata, file)

    return topology_reference, pt_trajectory, reduced_pt_trajectory, interactions, ligands, snapshots


# All analyses ---------------------------------------------------------------------------------
def run_analyses(
        topology_filename="md.imaged.rot.dry.pdb",
        trajectory_filename="md.imaged.rot.xtc"):

    topology_reference, pt_trajectory, reduced_pt_trajectory, interactions, ligands, snapshots = analysis_prep(
        topology_filename,
        trajectory_filename)

    print('Running analyses')

    # Run the RMSD analyses
    rmsd_references = [first_frame_filename, average_structure_filename]
    rmsds_analysis = 'md.rmsds.json'
    if required(rmsds_analysis):
        print('- RMSds analysis')
        rmsds(trajectory_filename, rmsds_analysis, rmsd_references)

    # Set the fluctuation analysis file name and run the analysis
    rmsf_analysis = 'md.rmsf.xvg'
    if required(rmsf_analysis):
        print('- Fluctuation')
        rmsf(topology_filename, trajectory_filename, rmsf_analysis)

    # Set the RMSd per resiude analysis file name and run the analysis
    rgyr_analysis = 'md.rgyr.xvg'
    if required(rgyr_analysis):
        print('- Radius of gyration')
        rgyr(topology_filename, trajectory_filename, rgyr_analysis)

    # Set the pca output filenames and run the analysis
    # WARNING: This analysis will generate several output files
    # Files 'average.pdb' and 'covar.log' are generated by the PCA but never used
    eigenvalues_filename = 'pca.eigenval.xvg'
    eigenvectors_filename = 'eigenvec.trr'
    if required(eigenvalues_filename) or required(eigenvectors_filename):
        print('- PCA')
        pca(topology_filename, trajectory_filename,
            eigenvalues_filename, eigenvectors_filename, snapshots)

    contacts_pca_filename = 'contacts_PCA.json'
    # DANI: Intenta usar mucha memoria, hay que revisar
    # DANI: De momento me lo salto
    if required(contacts_pca_filename) and False:
        print('- PCA on residue contacts')
        pca_contacts(
            trajectory_filename,
            topology_filename,
            interactions,
            contacts_pca_filename)

    # Set the RMSd per resiude analysis file name and run the analysis
    # WARNING: This analysis is fast enought to use the full trajectory instead of the reduced one
    # WARNING: However, the output file size depends on the trajectory size. It may be pretty big
    rmsd_perres_analysis = 'md.rmsd.perres.json'
    if required(rmsd_perres_analysis):
        print('- RMSd per residue')
        rmsd_per_residue(reduced_pt_trajectory,
                         rmsd_perres_analysis, topology_reference)

    # Set the RMSd pairwise analysis file name and run the analysis
    # WARNING: This analysis is fast enought to use the full trajectory instead of the reduced one
    # WARNING: However, the output file size depends on the trajectory size exponentially. It may be huge
    rmsd_pairwise_analysis = 'md.rmsd.pairwise.json'
    if required(rmsd_pairwise_analysis):
        print('- RMSd pairwise')
        rmsd_pairwise(reduced_pt_trajectory,
                      rmsd_pairwise_analysis, interactions)

    # Set the distance per residue analysis file name and run the analysis
    # WARNING: This analysis is not fast enought to use the full trajectory. It would take a while
    distance_perres_analysis = 'md.dist.perres.json'
    if required(distance_perres_analysis) and len(interactions) > 0:
        print('- Distance per residue')
        distance_per_residue(reduced_pt_trajectory,
                             distance_perres_analysis, interactions)

    # Set the hydrogen bonds analysis file name and run the analysis
    # WARNING: This analysis is fast enought to use the full trajectory instead of the reduced one
    # WARNING: However, the output file size depends on the trajectory
    # WARNING: Files have no limit, but analyses must be no heavier than 16Mb in BSON format
    # WARNING: In case of large surface interaction the output analysis may be larger than the limit
    # DANI: Esto no puede quedar así
    # DANI: Me sabe muy mal perder resolución con este análisis, porque en cáculo es muy rápido
    # DANI: Hay que crear un sistema de carga en mongo alternativo para análisis pesados
    hbonds_analysis = 'md.hbonds.json'
    if required(hbonds_analysis) and len(interactions) > 0:
        print('- Hydrogen bonds')
        hydrogen_bonds(reduced_pt_trajectory, hbonds_analysis,
                       topology_reference, interactions)

    # Set the energies analysis filename and run the analysis
    energies_analysis = 'md.energies.json'
    if required(energies_analysis) and len(ligands) > 0:
        print('- Energies')
        energies(topology_filename, trajectory_filename,
                 energies_analysis, topology_reference, snapshots, ligands)

    # Set the pockets analysis filename and run the analysis
    pockets_analysis = 'md.pockets.json'
    if required(pockets_analysis):
        print('- Pockets')
        pockets(topology_filename, trajectory_filename,
                pockets_analysis, topology_reference, snapshots)

    # Remove gromacs backups
    remove_trash()
    print('Done!')
