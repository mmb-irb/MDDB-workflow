from mddb_workflow.tools.get_box_size import get_box_size
from mddb_workflow.tools.get_atoms_count import get_atoms_count
from mddb_workflow.tools.generate_map import get_sequence_metadata
from mddb_workflow.utils.auxiliar import InputError, save_json
from mddb_workflow.utils.constants import MD_DIRECTORY
from mddb_workflow.utils.type_hints import *

# Input fields + interaction type
METADATA_INTERACTION_FIELDS = { "name", "agent_1", "agent_2", "selection_1", "selection_2", "type" }


def prepare_project_metadata (
    structure_file : 'File',
    trajectory_file : 'File',
    output_file : 'File',
    structure : 'Structure',
    residue_map : dict,
    protein_references_file : 'File',
    pdb_ids : list[str],
    ligand_references : dict,
    input_protein_references : list[str] | dict,
    input_ligands : list[dict],
    interactions : list[dict],
    warnings : dict,
    # Set all inputs to be loaded as they are
    input_force_fields : list[str],
    input_collections : list[str],
    input_chain_names : list[str],
    input_type : str,
    input_framestep : float,
    input_name : str,
    input_description : str,
    input_authors : list[str],
    input_groups : list[str],
    input_contact : str,
    input_program : str,
    input_version : str,
    input_method : str,
    input_license : str,
    input_linkcense : str,
    input_citation : str,
    input_thanks : str,
    input_links : list[dict],
    input_timestep : float,
    input_temperature : float,
    input_ensemble : str,
    input_water : str,
    input_boxtype : str,
    input_pbc_selection : str,
    input_cg_selection : str,
    input_customs : list[dict],
    input_orientation : list[float],
    input_multimeric : list[str],
    # Additional topic-specific inputs
    input_cv19_unit : str,
    input_cv19_startconf : str,
    input_cv19_abs : bool,
    input_cv19_nanobs : bool,
    ):
    """ Prepare a JSON file with all project metadata. """

    # Find out the box size (x, y and z)
    (boxsizex, boxsizey, boxsizez) = get_box_size(
        structure_file.path, trajectory_file.path)

    # Count different types of atoms and residues
    (system_atoms, system_residues, protein_atoms, protein_residues,
    dna_atoms, dna_residues, rna_atoms, rna_residues, lipid_atoms, lipid_residues,
    carbohydrates_atoms, carbohydrates_residues, solvent_atoms, solvent_residues,
    counter_cations, counter_anions, counter_ions, non_counter_ions, other_atoms) = get_atoms_count(structure)

    # Get protein references from the residues map
    # Get ligand references from the residues map
    protein_references = []
    ligand_references = []
    inchikey_references = []
    references = residue_map['references']
    if references and len(references) > 0:
        for ref, ref_type in zip(references, residue_map['reference_types']):
            if ref_type == 'protein':
                protein_references.append(ref)
            elif ref_type == 'ligand':
                ligand_references.append(ref)
            elif ref_type == 'inchikey':
                inchikey_references.append(ref)

    # Get ligand names if any
    forced_ligand_names = {
        lig['name']: lig['forced_name'] for lig in ligand_references if lig.get('forced_name', False) }
    if len(forced_ligand_names) == 0:
        forced_ligand_names = None

    # Make the forcefields a list in case it is a single string
    forcefields = input_force_fields
    if type(forcefields) == str:
        forcefields = [forcefields]

    # Collections must be null in case there are not collections
    collections = input_collections
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
    chainnames = input_chain_names
    if chainnames:
        for chain in chainnames.keys():
            if chain not in structure_chains:
                raise InputError(f'Chain {chain} from chainnames does not exist in the structure')

    # Get the MD type
    md_type = input_type
    # In case this is an ensemble and not a time related trajectory and not an ensemble, the framestep may be missing
    framestep = None if md_type == 'ensemble' else input_framestep

    # Metadata interactions are input interactions and the interaction types combined
    # Thus we take the processed interactions and remove the field we are not interested in
    metadata_interactions = []
    if interactions is not None:
        for interaction in interactions:
            metadata_interaction = { k: v for k, v in interaction.items() if k in METADATA_INTERACTION_FIELDS }
            metadata_interactions.append(metadata_interaction)

    # Make sure links are correct
    links = input_links
    if links != None:
        if type(links) != list: links = [ links ]
        for link in input_links:
            if type(link) != dict: raise InputError('Links must be a list of objects')
            if link.get('name', None) == None: raise InputError('Links must have a name')
            if link.get('url', None) == None: raise InputError('Links must have a URL')

    # Write the metadata file
    # Metadata keys must be in CAPS, as they are in the client
    metadata = {
        'NAME': input_name,
        'DESCRIPTION': input_description,
        'AUTHORS': input_authors,
        'GROUPS': input_groups,
        'CONTACT': input_contact,
        'PROGRAM': input_program,
        'VERSION': input_version,
        'TYPE': md_type,
        'METHOD': input_method,
        'LICENSE': input_license,
        'LINKCENSE': input_linkcense,
        'CITATION': input_citation,
        'THANKS': input_thanks,
        'LINKS': input_links,
        'PDBIDS': pdb_ids,
        'FORCED_REFERENCES': input_protein_references,
        'REFERENCES': protein_references,
        'INPUT_LIGANDS': input_ligands,
        # TODO: Ligands are now inchikeys only, remove after checking it does not break the client removing this
        'LIGANDS': [],
        'LIGANDNAMES': forced_ligand_names,
        'INCHIKEYS': inchikey_references,
        'PROTSEQ': sequence_metadata['protein_sequences'],
        'NUCLSEQ': sequence_metadata['nucleic_sequences'],
        'DOMAINS': sequence_metadata['domains'],
        'FRAMESTEP': framestep,
        'TIMESTEP': input_timestep,
        'TEMP': input_temperature,
        'ENSEMBLE': input_ensemble,
        'FF': forcefields,
        'WAT': input_water,
        'BOXTYPE': input_boxtype,
        'SYSTATS': system_atoms,
        'SYSTRES': system_residues,
        'PROTATS': protein_atoms,
        'PROTRES': protein_residues,
        'DNAATS': dna_atoms,
        'DNARES': dna_residues,
        'RNAATS': rna_atoms,
        'RNARES': rna_residues,
        'LIPIATS': lipid_atoms,
        'LIPIRES': lipid_residues,
        'CARBATS': carbohydrates_atoms,
        'CARBRES': carbohydrates_residues,
        'SOLVATS': solvent_atoms,
        'SOLVRES': solvent_residues,
        'COUNCAT': counter_cations,
        'COUNANI': counter_anions,
        'COUNION': counter_ions,
        'NOCNION': non_counter_ions,
        'OTHRATS': other_atoms,
        'INTERACTIONS': metadata_interactions,
        'PBC_SELECTION': input_pbc_selection,
        'CG_SELECTION': input_cg_selection,
        'CHAINNAMES': chainnames,
        'CUSTOMS': input_customs,
        'ORIENTATION': input_orientation,
        'PTM': ptm_names,
        'MULTIMERIC': input_multimeric,
        'COLLECTIONS': collections,
        'WARNINGS': warnings,
    }
    # Add boxsizes only if any of them is 0
    if boxsizex > 0 and boxsizey > 0 and boxsizez > 0:
        metadata['BOXSIZEX'] = boxsizex
        metadata['BOXSIZEY'] = boxsizey
        metadata['BOXSIZEZ'] = boxsizez
    # Add collection specific fields
    if 'cv19' in collections:
        cv19_unit = input_cv19_unit
        cv19_startconf = input_cv19_startconf
        cv19_abs = input_cv19_abs
        cv19_nanobs = input_cv19_nanobs
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
    save_json(metadata, output_file.path)

metadata_fields = set([ 'NAME', 'DESCRIPTION', 'AUTHORS', 'GROUPS', 'CONTACT', 'PROGRAM', 'VERSION',
    'TYPE', 'METHOD', 'LICENSE', 'LINKCENSE', 'CITATION', 'THANKS', 'LINKS', 'PDBIDS', 'FORCED_REFERENCES',
    'REFERENCES', 'INPUT_LIGANDS', 'LIGANDS', 'LIGANDNAMES', 'PROTSEQ', 'NUCLSEQ', 'DOMAINS', 'FRAMESTEP', 'TIMESTEP',
    'TEMP', 'ENSEMBLE', 'FF', 'WAT', 'BOXTYPE', 'SYSTATS', 'PROTATS', 'PROT', 'DPPC', 'SOL', 'NA', 'CL',
    'INTERACTIONS', 'PBC_SELECTION', 'CHAINNAMES', 'MEMBRANES', 'CUSTOMS', 'ORIENTATION', 'PTM',
    'MULTIMERIC', 'COLLECTIONS', 'WARNINGS', 'BOXSIZEX', 'BOXSIZEY', 'BOXSIZEZ', 'CV19_UNIT', 'CV19_STARTCONF',
    'CV19_ABS', 'CV19_NANOBS', 'CV19_VARIANT'
])

def generate_md_metadata (
    md_inputs : dict,
    structure : 'Structure',
    snapshots : int,
    reference_frame : int,
    warnings : dict,
    output_file : 'File'
    ):
    """Produce the MD metadata file to be uploaded to the database."""

    # Mine name and directory from MD inputs
    name = md_inputs.get('name', None)
    directory = md_inputs.get(MD_DIRECTORY, None)

    # Write the metadata file
    md_metadata = {
        'name': name,
        'frames': snapshots,
        'atoms': len(structure.atoms), # Should be always the same but we better have explicit confirmation
        'refframe': reference_frame,
        'warnings': warnings,
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
    save_json(md_metadata, output_file.path)
