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
    output_metadata_filename : str,
    register : dict
    ):

    # Set a function to retrieve 'inputs' values and handle missing keys
    inputs = None
    with open(inputs_filename, 'r') as file:
        inputs = json.load(file)
    def get_input(input: str):
        return inputs.get(input, None)

    # Calculate the frequency
    # Divide the simulation time, which is in nanoseconds (ns) by the number of frames
    # Multiply it by 1000 since we want the frequency in picoseconds (ps)
    length = get_input('length')
    frequency = (length / snapshots) * 1000 if length != None else None

    # Find out the box size (x, y and z)
    (boxsizex, boxsizey, boxsizez) = get_box_size(
        input_topology_filename, input_trajectory_filename)

    # Count different type of atoms and residues
    (systats, protats, prot, dppc, sol, na,
     cl) = get_atoms_count(input_topology_filename)

    # Extract some additional metadata from the inputs file which is required further
    ligands = get_input('ligands')
    if not ligands:
        ligands = []

    # Get the references from the residues map
    references = residues_map['references'] if residues_map else []

    # Make the forcefields a list in case it is a single string
    forcefields = get_input('ff')
    if type(forcefields) == str:
        forcefields = [forcefields]

    # Collections must be null in case there are not collections
    collections = get_input('collections')
    if not collections:
        collections = None

    # Write the metadata file
    # Metadata keys must be in CAPS, as they are in the client
    metadata = {
        'PDBIDS': get_input('pdbIds'),
        'NAME': get_input('name'),
        'COLLECTIONS': collections,
        'DESCRIPTION': get_input('description'),
        'AUTHORS': get_input('authors'),
        'GROUPS': get_input('groups'),
        'CONTACT': get_input('contact'),
        'PROGRAM': get_input('program'),
        'VERSION': get_input('version'),
        'TYPE': get_input('type'),
        'METHOD': get_input('method'),
        'LICENSE': get_input('license'),
        'LINKCENSE': get_input('linkcense'),
        'CITATION': get_input('citation'),
        'THANKS': get_input('thanks'),
        'LENGTH': length,
        'TIMESTEP': get_input('timestep'),
        'SNAPSHOTS': snapshots,
        'FREQUENCY': frequency,
        'FF': forcefields,
        'TEMP': get_input('temp'),
        'WAT': get_input('wat'),
        'BOXTYPE': get_input('boxtype'),
        'BOXSIZEX': boxsizex,
        'BOXSIZEY': boxsizey,
        'BOXSIZEZ': boxsizez,
        'ENSEMBLE': get_input('ensemble'),
        'SYSTATS': systats,
        'PROTATS': protats,
        'PROT': prot,
        'DPPC': dppc,
        'SOL': sol,
        'NA': na,
        'CL': cl,
        'LIGANDS': ligands,
        'CUSTOMS': get_input('customs'),
        'INTERACTIONS': get_input('interactions'),
        'FORCED_REFERENCES': get_input('forced_references'),
        'REFERENCES': references,
        'CHAINNAMES': get_input('chainnames'),
        'MEMBRANES': get_input('membranes'),
        'LINKS': get_input('links'),
        'ORIENTATION': get_input('orientation'),
        'WARNINGS': register['warnings'],
        # Collection specifics
        'CV19_UNIT': get_input('cv19_unit')
    }
    metadata_filename = 'metadata.json'
    with open(metadata_filename, 'w') as file:
        json.dump(metadata, file)
