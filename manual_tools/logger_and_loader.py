import urllib.request
import json
from os.path import exists
from subprocess import run, PIPE

# Set the file where project accessions are listed
source_file = 'pending.txt'
# Set database URL (or API URL)
database_url = 'https://mmb-dev.mddbr.eu/api'

# Mine accessions in the source file
projects = []
with open(source_file, 'r') as file:
    projects = file.readlines()
    projects = [ project.replace('\n', '') for project in projects ]

# Iterate projects
for accession in projects:
    project_url = f'{database_url}/rest/current/projects/{accession}'
    parsed_response = None
    try:
        response = urllib.request.urlopen(project_url)
        parsed_response = json.loads(response.read())
    except:
        raise Exception('EFE')

    # Log the project title
    title = parsed_response["metadata"]["NAME"]
    print(f'{accession} -> {title}')

    # Log the directory and if its exists or not
    directory_name = title.replace('EGFR - ','').replace(' ', '_')
    print(f'  {directory_name} -> {exists(directory_name)}')

    # Upload the directory
    process = run([
        "node",
        "/home/dbeltran/libraries/loader/index.js",
        "load",
        directory_name,
        "-a",
        accession,
    ], stdout=PIPE)
    logs = process.stdout.decode()
