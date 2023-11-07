from model_workflow.tools.get_box_size import get_box_size
from model_workflow.tools.get_atoms_count import get_atoms_count
from model_workflow.tools.generate_map import get_sequence_metadata

import json
from pathlib import Path

# Generate a JSON file with all project metadata
def generate_project_metadata (
    input_structure_filename : str,
    input_trajectory_filename : str,
    inputs_filename : str,
    structure : 'Structure',
    residues_map : dict,
    interactions : list,
    register : dict,
    output_metadata_filename : str
    ):

    # Set a function to retrieve 'inputs' values and handle missing keys
    inputs = None
    with open(inputs_filename, 'r') as file:
        inputs = json.load(file)
    def get_input(input: str):
        return inputs.get(input, None)

    # Find out the box size (x, y and z)
    (boxsizex, boxsizey, boxsizez) = get_box_size(
        input_structure_filename, input_trajectory_filename)

    # Count different type of atoms and residues
    (systats, protats, prot, dppc, sol, na,
     cl) = get_atoms_count(input_structure_filename)

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

    # Metadata interactions are simply the input interactions
    # Final interactions are used only to check which interactions were discarded
    # Thus failed interactions are removed from metadata
    metadata_interactions = get_input('interactions')
    final_interaction_names = [ interaction['name'] for interaction in interactions ]
    metadata_interactions = [ interaction for interaction in metadata_interactions if interaction['name'] in final_interaction_names ]

    # Get additional metadata related to the aminoacids sequence
    sequence_metadata = get_sequence_metadata(structure, residues_map)

    # Find the PTMs
    # Save only their names for now
    # DANI: Esto es temporal y de momento solo busca ser un parámetro de facil query
    # DANI: Cuando esté más maduro también almacenaremos residuo afectado, como mínimo
    ptms = structure.find_ptms()
    ptm_names = list(set([ ptm['name'] for ptm in ptms ]))

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
        'FRAMESTEP': get_input('framestep'),
        'TIMESTEP': get_input('timestep'),
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
        'INTERACTIONS': metadata_interactions,
        'PBC_SELECTION': get_input('pbc_selection'),
        'FORCED_REFERENCES': get_input('forced_references'),
        'REFERENCES': references,
        'SEQUENCES': sequence_metadata['sequences'],
        'DOMAINS': sequence_metadata['domains'],
        'PTM': ptm_names,
        'MULTIMERIC' : get_input('multimeric'),
        'CHAINNAMES': get_input('chainnames'),
        'MEMBRANES': get_input('membranes'),
        'LINKS': get_input('links'),
        'ORIENTATION': get_input('orientation'),
        'WARNINGS': register.warnings,
        # Collection specifics
        'CV19_UNIT': get_input('cv19_unit'),
        'CV19_STARTCONF': get_input('cv19_startconf'),
        'CV19_ABS': get_input('cv19_abs'),
        'CV19_NANOBS': get_input('cv19_nanobs'),
        'CV19_VARIANT': sequence_metadata['cv19_variant']
    }
    
    # Write metadata to a file
    with open(output_metadata_filename, 'w') as file:
        json.dump(metadata, file)

# Generate a JSON file with MD metadata
def generate_md_metadata (
    md_inputs : dict,
    structure : 'Structure',
    snapshots : int,
    register : dict,
    output_metadata_filename : str
    ):

    # Write the metadata file
    metadata = {
        'name': md_inputs['name'],
        'frames': snapshots,
        'atoms': len(structure.atoms), # Should be always the same but we better have explicit confirmation
        'warnings': register.warnings,
    }

    # Write metadata to a file
    with open(output_metadata_filename, 'w') as file:
        json.dump(metadata, file)