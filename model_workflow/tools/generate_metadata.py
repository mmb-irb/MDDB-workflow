from model_workflow.tools.get_box_size import get_box_size
from model_workflow.tools.get_atoms_count import get_atoms_count

import json
from pathlib import Path

# Generate a JSON file with all the metadata
def generate_metadata (
    input_topology_filename : str,
    input_trajectory_filename : str,
    inputs_filename : str,
    snapshots : int,
    residues_map : dict,
    output_metadata_filename : str):

    # Set a function to retrieve 'inputs' values and handle missing keys
    inputs = None
    with open(inputs_filename, 'r') as file:
        inputs = json.load(file)
    def getInput(input: str):
        return inputs.get(input, None)

    # Calculate the frequency
    # Divide the simulation time, which is in nanoseconds (ns) by the number of frames
    # Multiply it by 1000 since we want the frequency in picoseconds (ps)
    length = getInput('length')
    frequency = (length / snapshots) * 1000

    # Find out the box size (x, y and z)
    (boxsizex, boxsizey, boxsizez) = get_box_size(
        input_topology_filename, input_trajectory_filename)

    # Count different type of atoms and residues
    (systats, protats, prot, dppc, sol, na,
     cl) = get_atoms_count(input_topology_filename)

    # Extract some additional metadata from the inputs file which is required further
    ligands = getInput('ligands')
    if not ligands:
        ligands = []

    # Get the references from the residues map
    references = residues_map['references'] if residues_map else []

    # Write the metadata file
    # Metadata keys must be in CAPS, as they are in the client
    metadata = {
        'PDBIDS': getInput('pdbIds'),
        'NAME': getInput('name'),
        'UNIT': getInput('unit'),
        'DESCRIPTION': getInput('description'),
        'AUTHORS': getInput('authors'),
        'GROUPS': getInput('groups'),
        'CONTACT': getInput('contact'),
        'PROGRAM': getInput('program'),
        'VERSION': getInput('version'),
        'METHOD': getInput('method'),
        'LICENSE': getInput('license'),
        'LINKCENSE': getInput('linkcense'),
        'CITATION': getInput('citation'),
        'THANKS': getInput('thanks'),
        'LENGTH': length,
        'TIMESTEP': getInput('timestep'),
        'SNAPSHOTS': snapshots,
        'FREQUENCY': frequency,
        'FF': getInput('ff'),
        'TEMP': getInput('temp'),
        'WAT': getInput('wat'),
        'BOXTYPE': getInput('boxtype'),
        'BOXSIZEX': boxsizex,
        'BOXSIZEY': boxsizey,
        'BOXSIZEZ': boxsizez,
        'ENSEMBLE': getInput('ensemble'),
        'SYSTATS': systats,
        'PROTATS': protats,
        'PROT': prot,
        'DPPC': dppc,
        'SOL': sol,
        'NA': na,
        'CL': cl,
        'LIGANDS': ligands,
        'CUSTOMS': getInput('customs'),
        'INTERACTIONS': getInput('interactions'),
        'REFERENCES': references,
        'CHAINNAMES': getInput('chainnames'),
        'MEMBRANES': getInput('membranes'),
        'LINKS': getInput('links'),
    }
    metadata_filename = 'metadata.json'
    with open(metadata_filename, 'w') as file:
        json.dump(metadata, file)
