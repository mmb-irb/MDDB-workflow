from packaging.version import Version
from os import mkdir, replace
from os.path import exists
from re import sub
from time import time

from mddb_workflow.utils.auxiliar import load_json, save_json, warn, InputError
from mddb_workflow.utils.constants import PROTEIN_REFERENCE_VERSION, PROTEIN_REFERENCES_FILENAME
from mddb_workflow.utils.constants import INCHI_REFERENCE_VERSION, INCHIKEY_REFERENCES_FILENAME
from mddb_workflow.utils.constants import PDB_REFERENCE_VERSION, PDB_REFERENCES_FILENAME
from mddb_workflow.utils.database import Database, get_available_nodes
from mddb_workflow.utils.loader import Loader

from mddb_workflow.tools.generate_map import get_uniprot_reference
from mddb_workflow.tools.generate_pdb_references import get_pdb_reference


# Set the available reference types and their configurations
REFERENCE_TYPE_CONFIGURATIONS = {
    'uniprot': {
        'version': PROTEIN_REFERENCE_VERSION,
        'endpoint': 'proteins',
        'id_key': 'uniprot',
        'maker': get_uniprot_reference,
        'output': PROTEIN_REFERENCES_FILENAME,
    },
    # 'inchi': {
    #     'version': INCHI_REFERENCE_VERSION,
    #     'endpoint': 'inchikeys',
    #     'id_key': 'inchi',
    #     'maker': None, # DANI: Falta una función que solo requiera del inchi
    #     'output': INCHIKEY_REFERENCES_FILENAME,
    # },
    'pdb': {
        'version': PDB_REFERENCE_VERSION,
        'endpoint': 'pdbs',
        'id_key': 'id',
        'maker': get_pdb_reference,
        'output': PDB_REFERENCES_FILENAME,
    },
}
# Set a flag to update them all
ALL_FLAG = 'all'
# Collect all requestable available reference types
AVAILABLE_REFERENCE_TYPES = list(REFERENCE_TYPE_CONFIGURATIONS.keys())
# Set every how many updated references we upload to the database
LOADER_BATCH = 100



def update_references (
    node_url_or_alias : str,
    reference_type : str,
    loader_directory : str | None = None):

    # If all nodes are requested then simply call this same function with every node
    if node_url_or_alias == ALL_FLAG:
        # Get the available nodes
        available_nodes = get_available_nodes()
        for node_alias in available_nodes.keys():
            update_references(node_alias, reference_type)
        return

    # If all reference types are requested then simply call this same function with every type
    if reference_type == ALL_FLAG:
        for current_reference_type in REFERENCE_TYPE_CONFIGURATIONS.keys():
            update_references(node_url_or_alias, current_reference_type)
        return

    # Make sure the type is known
    if reference_type not in REFERENCE_TYPE_CONFIGURATIONS:
        raise InputError(f'Unknown reference type "{reference_type}". '
            f'Please select any of the following: {", ".join(AVAILABLE_REFERENCE_TYPES)}')

    print(f'Updating references of type "{reference_type}" from "{node_url_or_alias}"')
    reference_config = REFERENCE_TYPE_CONFIGURATIONS[reference_type]

    # Instantiate the database handler
    database = Database(node_url_or_alias)

    # If we are to upload the new references to the corresponding database then prepare the loader
    loader = Loader(loader_directory, node_url_or_alias) if loader_directory else None

    # Paginate to get all available references and their versions
    references = database.get_all_refences_data(reference_config['endpoint'])

    # Get the references which are not updated according to this workflow latest versions
    updated_version = Version(reference_config['version'])
    id_key = reference_config['id_key']
    outdated_reference_ids = [
        ref[id_key] for ref in references
        if Version(ref.get('version', '0.0.0')) < updated_version
    ]
    remote_outdated_count = len(outdated_reference_ids)
    print(f'  Found {remote_outdated_count} remote outdated references out of {len(references)}')

    # Make a directory for the current url or alias
    directory = sub('https?://', '', node_url_or_alias.replace('/api/', ''))
    if not exists(directory): mkdir(directory)

    # Set the filepaths where output is to be written
    output_references_filename = reference_config['output']
    output_references_filepath = f"{directory}/{output_references_filename}"
    provisional_output_references_filepath = f"{directory}/provisional_{output_references_filename}"

    # Load already updated references in this directory, if any
    # Thus there is no need to update them again
    print(f'  Results will be saved to {output_references_filepath}')
    updated_references = []
    already_updated_references_count = 0
    if exists(output_references_filepath):
        updated_references = load_json(output_references_filepath)
        already_updated_references_count = len(updated_references)
        print(f'  There are {already_updated_references_count} already updated references in {output_references_filepath}')

    # Set a function to upload the already updated references
    def load_batch():
        nonlocal updated_references
        print(f'Loading batch of {len(updated_references)} updated references')
        # Upload the already updated references
        if loader.load(directory):
            # After the load we make a cleanup both in memory and disk
            # 1) Rename the updated references file to save apart
            # Thus new updated references will be saved in a new fresh JSON
            # Otherwise the rewrite operations become expensive
            # Use timestamps for the batched files to avoid data loss
            timestamp = round(time())
            batch_backup = f'{directory}/uploaded_time{timestamp}_{output_references_filename}'
            replace(output_references_filepath, batch_backup)
            # 2) Cleanup the already uploaded references from the updated references list
            updated_references = []
        # If something went wrong with the upload then stop here
        else: raise RuntimeError('Something went wrong when uploading updated references')

    # Substract their ids from the outdated ids list
    locally_updated_references = set([ ref[id_key] for ref in updated_references ])
    outdated_reference_ids = [
        ref_id for ref_id in outdated_reference_ids
        if ref_id not in locally_updated_references
    ]
    actual_outdated_count = len(outdated_reference_ids)
    print(f'  There are {actual_outdated_count} outdated references left after using the locally updated references')

    # Keep track of the every reference we fail to update
    # Note that the PDB may return a random 500 (internal server error) sometimes
    # If we fail to update some references then just trying again may solve the issues
    failed_references_count = 0

    # Remake every outdated reference
    reference_maker = reference_config['maker']
    for o, outdated_reference_id in enumerate(outdated_reference_ids, already_updated_references_count + 1):
        # If the amount of updated references is enough for the loader batch then upload them
        if loader and len(updated_references) >= LOADER_BATCH:
            load_batch()
        # Update the next reference
        print(f'  Remaking {outdated_reference_id} ({o}/{remote_outdated_count})')
        # Wrap the updater in a try/except so we are resilient when having response issues
        try:
            new_reference = reference_maker(outdated_reference_id)
        except:
            warn(f'Something went wrong while updating {outdated_reference_id} -> Skipped')
            failed_references_count += 1
            continue
        updated_references.append(new_reference)
        # Write the updated references to disk after every update
        # Note that there may be a lot of references so the full update may take much time
        # We save the progress constantly so we don't have to start from scratch if something goes wrong
        # Note that we do not write directly to the output file in case the user interrupts the process
        # Instead, we first write to a provisional file and then replace the original file with it
        # Otherwise we could loose all the already updated references if we are not lucky
        save_json(updated_references, provisional_output_references_filepath)
        replace(provisional_output_references_filepath, output_references_filepath)

    # Load the last remaining batch of simulations
    if loader and len(updated_references) > 0:
        load_batch()

    # Warn the user about failed updates
    if failed_references_count > 0:
        warn(f' Failed to update {failed_references_count} out of {remote_outdated_count} references. '
            'This may be caused by network problems so please run the updater again.')