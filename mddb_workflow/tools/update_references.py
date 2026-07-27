from packaging.version import Version
from os.path import exists

from mddb_workflow.utils.auxiliar import load_json, save_json
from mddb_workflow.utils.constants import PROTEIN_REFERENCE_VERSION, INCHI_REFERENCE_VERSION, PDB_REFERENCE_VERSION
from mddb_workflow.utils.constants import PROTEIN_REFERENCES_FILENAME, INCHIKEY_REFERENCES_FILENAME, PDB_REFERENCES_FILENAME
from mddb_workflow.utils.database import Database

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
ALL_TYPES = 'all'
# Collect all requestable available reference types
AVAILABLE_REFERENCE_TYPES = [ *list(REFERENCE_TYPE_CONFIGURATIONS.keys()), ALL_TYPES ]

def update_references (database_url_or_alias : str, reference_type : str, ssleep : bool = False):

    # If all reference types are requested then simply call this same function with every type
    if reference_type == ALL_TYPES:
        for current_reference_type in REFERENCE_TYPE_CONFIGURATIONS.keys():
            update_references(database_url_or_alias, current_reference_type, ssleep)
        return

    # Make sure the type is known
    if reference_type not in REFERENCE_TYPE_CONFIGURATIONS:
        print(f'Unknown reference type "{reference_type}". Please select one of the following: {", ".join(AVAILABLE_REFERENCE_TYPES)}')
        return

    print(f'Updating references of type "{reference_type}"')
    reference_config = REFERENCE_TYPE_CONFIGURATIONS[reference_type]

    # Instantiate the database handler
    database = Database(database_url_or_alias, ssleep)

    # Paginate to get all available references and their versions
    references = database.get_all_refences_data(reference_config['endpoint'])

    # Get the references which are not updated according to this workflow latest versions
    updated_version = Version(reference_config['version'])
    id_key = reference_config['id_key']
    outdated_reference_ids = [
        ref[id_key] for ref in references
        if Version(ref.get('version', '0.0.0')) < updated_version
    ]
    outdated_count = len(outdated_reference_ids)
    print(f'  Found {outdated_count} remote outdated references out of {len(references)}')

    # Load already updated references in this directory, if any
    # Thus there is no need to update them again
    output_references_filepath = reference_config['output']
    updated_references = []
    if exists(output_references_filepath):
        load_json(output_references_filepath)
        print(f'  There are {len(updated_references)} locally already updated references')

    # Substract their ids from the outdated ids list
    locally_updated_references = set([ ref[id_key] for ref in updated_references ])
    outdated_reference_ids = [
        ref_id for ref_id in outdated_reference_ids
        if ref_id not in locally_updated_references
    ]
    outdated_count = len(outdated_reference_ids)
    print(f'  There are {outdated_count} outdated references left after using the locally updated references')

    # Remake every outdated reference
    reference_maker = reference_config['maker']
    for o, outdated_reference_id in enumerate(outdated_reference_ids, 1):
        print(f'  Remaking {outdated_reference_id} ({o}/{outdated_count})')
        new_reference = reference_maker(outdated_reference_id)
        updated_references.append(new_reference)
        # Write the updated references to disk after every update
        # Note that there may be a lot of references so the ful update may take much time
        # We save the progress constantly so we don't have to start from scratch if somethig goes wrong
        save_json(updated_references, output_references_filepath)