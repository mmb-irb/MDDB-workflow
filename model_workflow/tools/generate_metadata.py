from model_workflow.tools.get_box_size import get_box_size
from model_workflow.tools.get_atoms_count import get_atoms_count
from model_workflow.tools.generate_map import get_sequence_metadata
from model_workflow.utils.auxiliar import InputError, save_json

from pathlib import Path
from typing import Callable

# Generate a JSON file with all project metadata
def generate_project_metadata (
    input_structure_filename : str,
    input_trajectory_filename : str,
    get_input : Callable,
    structure : 'Structure',
    residue_map : dict,
    interactions : list,
    register : dict,
    output_metadata_filename : str
    ):

    # Find out the box size (x, y and z)
    (boxsizex, boxsizey, boxsizez) = get_box_size(
        input_structure_filename, input_trajectory_filename)

    # Count different type of atoms and residues
    (systats, protats, prot, dppc, sol, na,
     cl) = get_atoms_count(input_structure_filename)

    # Get protein references from the residues map
    # Get ligand references from the residues map
    protein_references = []
    ligand_references = []
    references = residue_map['references']
    if references and len(references) > 0:
        for ref, ref_type in zip(references, residue_map['reference_types']):
            if ref_type == 'protein':
                protein_references.append(ref)
            elif ref_type == 'ligand':
                ligand_references.append(ref)

    # Make the forcefields a list in case it is a single string
    forcefields = get_input('ff')
    if type(forcefields) == str:
        forcefields = [forcefields]

    # Collections must be null in case there are not collections
    collections = get_input('collections')
    if not collections:
        collections = []

    # Metadata interactions are simply the input interactions
    # Final interactions are used only to check which interactions were discarded
    # Thus failed interactions are removed from metadata
    final_interaction_names = [ interaction['name'] for interaction in interactions ]
    final_metadata_interactions = []
    input_interactions = get_input('interactions')
    if input_interactions and len(input_interactions) > 0:
        for input_interaction in input_interactions:
            if input_interaction['name'] in final_interaction_names:
                final_metadata_interactions.append(input_interaction)

    # Get additional metadata related to the aminoacids sequence
    sequence_metadata = get_sequence_metadata(structure, residue_map)

    # Find the PTMs
    # Save only their names for now
    # DANI: Esto es temporal y de momento solo busca ser un parámetro de facil query
    # DANI: Cuando esté más maduro también almacenaremos residuo afectado, como mínimo
    ptms = structure.find_ptms()
    ptm_names = list(set([ ptm['name'] for ptm in ptms ]))

    # Check chainnames to actually exist in the structure
    structure_chains = set([ chain.name for chain in structure.chains ])
    chainnames = get_input('chainnames')
    if chainnames:
        for chain in chainnames.keys():
            if chain not in structure_chains:
                raise InputError(f'Chain {chain} from chainnames does not exist in the structure')

    # Write the metadata file
    # Metadata keys must be in CAPS, as they are in the client
    metadata = {
        'NAME': get_input('name'),
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
        'LINKS': get_input('links'),
        'PDBIDS': get_input('pdbIds'),
        'FORCED_REFERENCES': get_input('forced_references'),
        'REFERENCES': protein_references,
        'INPUT_LIGANDS': get_input('ligands'),
        'LIGANDS': ligand_references,
        'SEQUENCES': sequence_metadata['sequences'],
        'DOMAINS': sequence_metadata['domains'],
        'FRAMESTEP': get_input('framestep'),
        'TIMESTEP': get_input('timestep'),
        'TEMP': get_input('temp'),
        'ENSEMBLE': get_input('ensemble'),
        'FF': forcefields,
        'WAT': get_input('wat'),
        'BOXTYPE': get_input('boxtype'),
        'BOXSIZEX': boxsizex,
        'BOXSIZEY': boxsizey,
        'BOXSIZEZ': boxsizez,
        'SYSTATS': systats,
        'PROTATS': protats,
        'PROT': prot,
        'DPPC': dppc,
        'SOL': sol,
        'NA': na,
        'CL': cl,
        'INTERACTIONS': final_metadata_interactions,
        'PBC_SELECTION': get_input('pbc_selection'),
        'CHAINNAMES': chainnames,
        'MEMBRANES': get_input('membranes'),
        'CUSTOMS': get_input('customs'),
        'ORIENTATION': get_input('orientation'),
        'PTM': ptm_names,
        'MULTIMERIC' : get_input('multimeric'),
        'COLLECTIONS': collections,
        'WARNINGS': register.warnings,
    }

    # Add collection specific fields
    if 'cv19' in collections:
        metadata = { **metadata,
            'CV19_UNIT': get_input('cv19_unit'),
            'CV19_STARTCONF': get_input('cv19_startconf'),
            'CV19_ABS': get_input('cv19_abs'),
            'CV19_NANOBS': get_input('cv19_nanobs'),
            'CV19_VARIANT': sequence_metadata['cv19_variant']
        }
    
    # Write metadata to a file
    save_json(metadata, output_metadata_filename)

# Generate a JSON file with MD metadata
def generate_md_metadata (
    md_inputs : dict,
    structure : 'Structure',
    snapshots : int,
    reference_frame : int,
    register : dict,
    output_metadata_filename : str
    ):

    # Write the metadata file
    metadata = {
        'name': md_inputs['name'],
        'frames': snapshots,
        'atoms': len(structure.atoms), # Should be always the same but we better have explicit confirmation
        'refframe': reference_frame,
        'warnings': register.warnings,
    }

    # Write metadata to a file
    save_json(metadata, output_metadata_filename)
