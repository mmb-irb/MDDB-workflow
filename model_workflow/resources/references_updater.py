import urllib.request
import json
from os.path import exists

from model_workflow.tools.generate_map import get_uniprot_reference
from model_workflow.utils.auxiliar import save_json

# Set the database API URL
database_api_url = 'https://mdposit-dev.mddbr.eu/api/rest/v1'

# IMPORTANT: Custom fields from old references are also rescued
custom_fields = ['entropies', 'epitopes', 'custom_domains']

# 1. Get all uniprot ids from project options
#    Get also all references in the database just to count the number of orphan references
# 2. Remake each reference object from uniprot fresh data
# 3. Create a huge references.json with the new references
#    This references.json is to be used by the loader after manually deleting all current references
#    It may be appended to any project, since the project will be not modified
def updater ():
    # In case there is a references.json in the current directory already abort
    references_filename = 'references.json'
    if exists(references_filename):
        raise SystemExit('File ' + references_filename + ' already exists')
    # Request the options
    options_url = database_api_url + '/projects/options?projection=metadata.REFERENCES'
    with urllib.request.urlopen(options_url) as response:
        options = json.loads(response.read().decode("utf-8"))
    project_ids = set(options['metadata.REFERENCES'].keys())
    # Remove null uniprot values
    project_ids -= set(['null', 'noref'])
    project_ids_count = len(project_ids)
    print('There are ' + str(project_ids_count) + ' different uniprot ids among projects')
    # Request all references
    references_url = database_api_url + '/references'
    with urllib.request.urlopen(references_url) as response:
        reference_ids = set(json.loads(response.read().decode("utf-8")))
    print('There are ' + str(len(reference_ids)) + ' references in the database')
    # Show the orphan uniprot ids
    orphan_uniprots = project_ids.difference(reference_ids)
    orphan_uniprots_count = len(orphan_uniprots)
    if orphan_uniprots_count > 0:
        print('We have ' + str(orphan_uniprots_count) + ' orphan project uniprot ids: ' + ', '.join(orphan_uniprots))
        print(' WARNING: This may be a problem')
    # Show the orphan references
    orphan_references = reference_ids.difference(project_ids)
    orphan_references_count = len(orphan_references)
    if orphan_references_count > 0:
        print('We have ' + str(orphan_references_count) + ' orphan references: ' + ', '.join(orphan_references))
        print(' -> They will be excluded from the new generated references.json')
    # Now build references for each uniprot from scratch
    new_references = []
    final_project = str(project_ids_count)
    for n, uniprot_id in enumerate(project_ids, 1):
        print("Building reference " + str(n) + '/' + final_project, end='\r')
        new_reference = get_uniprot_reference(uniprot_id)
        # Download the old reference and rescue custom fields from it
        old_reference_url = database_api_url + '/references/' + uniprot_id
        with urllib.request.urlopen(old_reference_url) as response:
            old_reference = json.loads(response.read().decode("utf-8"))
        for custom_field in custom_fields:
            old_value = old_reference.get(custom_field, None)
            if old_value:
                new_reference[custom_field] = old_value
        new_references.append(new_reference)
    # Write references to a new file
    save_json(new_references, references_filename)

updater()