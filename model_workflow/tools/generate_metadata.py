from model_workflow.tools.get_box_size import get_box_size
from model_workflow.tools.get_atoms_count import get_atoms_count
from model_workflow.tools.generate_map import get_sequence_metadata
from model_workflow.utils.auxiliar import InputError, save_json
from model_workflow.utils.constants import MD_DIRECTORY

from pathlib import Path
from typing import Callable

# Generate a JSON file with all project metadata
def generate_project_metadata (
    input_structure_filename : str,
    input_trajectory_filename : str,
    get_input : Callable,
    structure : 'Structure',
    residue_map : dict,
    protein_references_file : 'File',
    register : dict,
    output_metadata_filename : str,
    ligand_customized_names : str
    ):

    print('-> Generating project metadata')

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

    # Get ligand names if any
    input_ligands = get_input('ligands')
    if len(ligand_customized_names) == 0:
        ligand_customized_names = None

    # Make the forcefields a list in case it is a single string
    forcefields = get_input('ff')
    if type(forcefields) == str:
        forcefields = [forcefields]

    # Collections must be null in case there are not collections
    collections = get_input('collections')
    if not collections:
        collections = []

    # Get additional metadata related to the aminoacids sequence
    sequence_metadata = get_sequence_metadata(structure, protein_references_file, residue_map)

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

    # Get the MD type
    md_type = get_input('type')
    # In case this is an ensemble and not a time related trajectory and not an ensemble, the framestep may be missing
    framestep = None if md_type == 'ensemble' else get_input('framestep')

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
        'TYPE': md_type,
        'METHOD': get_input('method'),
        'LICENSE': get_input('license'),
        'LINKCENSE': get_input('linkcense'),
        'CITATION': get_input('citation'),
        'THANKS': get_input('thanks'),
        'LINKS': get_input('links'),
        'PDBIDS': get_input('pdbIds'),
        'FORCED_REFERENCES': get_input('forced_references'),
        'REFERENCES': protein_references,
        'INPUT_LIGANDS': input_ligands,
        'LIGANDS': ligand_references,
        'LIGANDNAMES': ligand_customized_names,
        'SEQUENCES': sequence_metadata['sequences'],
        'DOMAINS': sequence_metadata['domains'],
        'FRAMESTEP': framestep,
        'TIMESTEP': get_input('timestep'),
        'TEMP': get_input('temp'),
        'ENSEMBLE': get_input('ensemble'),
        'FF': forcefields,
        'WAT': get_input('wat'),
        'BOXTYPE': get_input('boxtype'),
        'SYSTATS': systats,
        'PROTATS': protats,
        'PROT': prot,
        'DPPC': dppc,
        'SOL': sol,
        'NA': na,
        'CL': cl,
        'INTERACTIONS': get_input('interactions'),
        'PBC_SELECTION': get_input('pbc_selection'),
        'CHAINNAMES': chainnames,
        'MEMBRANES': get_input('membranes'),
        'CUSTOMS': get_input('customs'),
        'ORIENTATION': get_input('orientation', optional=True),
        'PTM': ptm_names,
        'MULTIMERIC' : get_input('multimeric'),
        'COLLECTIONS': collections,
        'WARNINGS': register.warnings,
    }
    # Add boxsizes only if any of them is 0
    if boxsizex > 0 and boxsizey > 0 and boxsizez > 0:
        metadata['BOXSIZEX'] = boxsizex
        metadata['BOXSIZEY'] = boxsizey
        metadata['BOXSIZEZ'] = boxsizez
    # Add collection specific fields
    if 'cv19' in collections:
        cv19_unit = get_input('cv19_unit', optional=True)
        cv19_startconf = get_input('cv19_startconf', optional=True)
        cv19_abs = get_input('cv19_abs', optional=True)
        cv19_nanobs = get_input('cv19_nanobs', optional=True)
        cv19_variant = sequence_metadata['cv19_variant']

        if cv19_unit is not None:
            metadata['CV19_UNIT'] = cv19_unit

        if cv19_startconf is not None:
            metadata['CV19_STARTCONF'] = cv19_startconf

        if cv19_abs is not None:
            metadata['CV19_ABS'] = cv19_abs

        if cv19_nanobs is not None:
            metadata['CV19_NANOBS'] = cv19_nanobs

        if cv19_variant is not None:
            metadata['CV19_VARIANT'] = cv19_variant
    
    # Write metadata to a file
    save_json(metadata, output_metadata_filename)

metadata_fields = set([ 'NAME', 'DESCRIPTION', 'AUTHORS', 'GROUPS', 'CONTACT', 'PROGRAM', 'VERSION',
    'TYPE', 'METHOD', 'LICENSE', 'LINKCENSE', 'CITATION', 'THANKS', 'LINKS', 'PDBIDS', 'FORCED_REFERENCES', 
    'REFERENCES', 'INPUT_LIGANDS', 'LIGANDS', 'LIGANDNAMES', 'SEQUENCES', 'DOMAINS', 'FRAMESTEP', 'TIMESTEP',
    'TEMP', 'ENSEMBLE', 'FF', 'WAT', 'BOXTYPE', 'SYSTATS', 'PROTATS', 'PROT', 'DPPC', 'SOL', 'NA', 'CL',
    'INTERACTIONS', 'PBC_SELECTION', 'CHAINNAMES', 'MEMBRANES', 'CUSTOMS', 'ORIENTATION', 'PTM', 
    'MULTIMERIC', 'COLLECTIONS', 'WARNINGS', 'BOXSIZEX', 'BOXSIZEY', 'BOXSIZEZ', 'CV19_UNIT', 'CV19_STARTCONF',
    'CV19_ABS', 'CV19_NANOBS', 'CV19_VARIANT'
])

# Generate a JSON file with MD metadata
def generate_md_metadata (
    md_inputs : dict,
    structure : 'Structure',
    snapshots : int,
    reference_frame : int,
    register : dict,
    output_metadata_filename : str
    ):

    print('-> Generating MD metadata')

    # Mine name and directory from MD inputs
    name = md_inputs.get('name', None)
    directory = md_inputs.get(MD_DIRECTORY, None)
    
    # Write the metadata file
    md_metadata = {
        'name': name,
        'frames': snapshots,
        'atoms': len(structure.atoms), # Should be always the same but we better have explicit confirmation
        'refframe': reference_frame,
        'warnings': register.warnings,
    }

    # Get other MD inputs than the name and the directory
    other_md_inputs = { k: v for k, v in md_inputs.items() }
    # Remove name from MD inputs to not further overwrite project metadata
    if name:
        del other_md_inputs['name']
    # Remove the directory name form MD inputs since it is not to be uploaded to the database
    if directory:
        del other_md_inputs[MD_DIRECTORY]

    # Inherit all metadata fields
    metadata = {}
    for field in metadata_fields:
        input_field = field.lower()
        field_value = other_md_inputs.get(input_field, None)
        if field_value:
            metadata[field] = field_value

    # Add the matadata field only if there is at least one value
    if len(metadata) > 0:
        md_metadata['metadata'] = metadata

    # Write metadata to a file
    save_json(md_metadata, output_metadata_filename)
