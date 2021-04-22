from model_workflow.tools.get_box_size import get_box_size
from model_workflow.tools.get_atoms_count import get_atoms_count
from model_workflow.tools.generate_map import generate_map

import json
from pathlib import Path

# Set the path to stored topology references which are required for this analysis
resources = str(Path(__file__).parent.parent / "utils" / "resources")
toporefs_source = resources + '/toporefs.json'

# Generate a JSON file with all the metadata
def generate_metadata (
    input_topology_filename : str,
    input_trajectory_filename : str,
    inputs_filename : str,
    processed_interactions : list,
    snapshots : int,
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
    } for interaction in processed_interactions]

    # Get the input topology reference
    toporefs = getInput('toporefs')

    # Load topology references stored in this workflow
    stored_toporefs = None
    with open(toporefs_source, 'r') as file:
        stored_toporefs = json.load(file)

    # For each input topology reference, align the reference sequence with the topology sequence
    # Each metadata 'toporef' must include the input toporef name and the align map
    metadata_toporefs = []
    for toporef in toporefs:
        # Find the corresponding stored topology reference and get its residues sequence
        name = toporef['name']
        reference = next(st for st in stored_toporefs if st['name'] == name)
        reference_sequence = reference['sequence']
        reference_map = generate_map(reference_sequence, input_topology_filename)
        metadata_toporefs.append({
            'name': name,
            'map': reference_map,
        })

    # Write the metadata file
    # Metadata keys must be in CAPS, as they are in the client
    metadata = {
        'PDBID': getInput('pdbId'),
        'NAME': getInput('name'),
        'UNIT': getInput('unit'),
        'DESCRIPTION': getInput('description'),
        'AUTHORS': getInput('authors'),
        'GROUPS': getInput('groups'),
        'CONTACT': getInput('contact'),
        'PROGRAM': getInput('program'),
        'VERSION': getInput('version'),
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
        'TOPOREFS': metadata_toporefs,
        'CUSTOMS': getInput('customs'),
        'INTERACTIONS': metadata_interactions,
        'CHAINNAMES': getInput('chainnames'),
    }
    metadata_filename = 'metadata.json'
    with open(metadata_filename, 'w') as file:
        json.dump(metadata, file)
