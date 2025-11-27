import os
import json
import urllib
from math import ceil
from mddb_workflow.console import DEFAULT_API_URL

# The URL to the API from our node
default_url = DEFAULT_API_URL + 'rest/current'


def query_api(url: str, verbose=False) -> dict:
    """Parse the URL in case it contains any HTTP control characters
    Replace white spaces by the corresponding percent notation character.
    """
    parsed_url = url.replace(" ", "%20")
    if verbose:
        print(f'Querying API URL: {parsed_url}')
    with urllib.request.urlopen(parsed_url) as response:
        return json.loads(response.read().decode("utf-8"))


def count_projects(query_url):
    """Count the number of projects matching the query URL."""
    # Ask which projects have a tpr file
    response = query_api(query_url)
    n_projects = response['filteredCount']
    print(f'We found {n_projects} projects')
    return n_projects


def paginate_projects(selection_query, limit=100, api_url=default_url):
    """Paginate through all projects matching the selection query."""
    # Ask which projects have a tpr file
    query_url = api_url + f'/projects?query={{{selection_query}}}'
    n_projects = count_projects(query_url)
    # Calculate the expected number of pages
    pages = ceil(n_projects / limit)
    # Set a list to store all the mined accession values
    accessions = []
    projects_list = []
    # Iterate over pages
    # Note that pages are 1-base numerated, and NOT 0-based
    for page in range(1, pages + 1):
        print(f'Requesting page {page} / {pages}', end='\r')
        # Set the URL for the projects endpoint
        # This time include both limit and page parameters
        paginated_url = f'{query_url}&limit={limit}&page={page}'
        # Query the API
        response = query_api(paginated_url)
        # Mine target data
        projects = response['projects']
        projects_list.extend(projects)
        project_accessions = [project['accession'] for project in projects]
        accessions += project_accessions

    print(f'We have mined {len(accessions)} project accessions')
    return accessions, projects_list


def download_files(accessions: list, file: str, output_folder: str, api_url=default_url):
    """Download specific files from a list of project accessions."""
    topologies = {}
    os.makedirs(output_folder, exist_ok=True)
    for a, accession in enumerate(accessions, 1):
        file_path = f'{output_folder}/{accession}_{file}'
        if not os.path.exists(file_path):
            print(f'Parsing {file} from accession {accession} ({a}/{len(accessions)})          ', end='\r')
            if file == 'topology.json':
                file_url = f'{api_url}/projects/{accession}/topology'
            else:
                file_url = f'{api_url}/projects/{accession}/files/{file}'
            try:
                urllib.request.urlretrieve(file_url, file_path)
            except Exception:
                print(f'error: Failed to download {file_url}')
                continue
        with open(file_path) as f:
            if file == 'topology.json':
                topology = json.load(f)
                topologies[accession] = topology
            else:
                topologies[accession] = f.read()
    return topologies

# Example of usage
# Find all projects with a topology.top file
# accessions = paginate_projects('"files.name":"topology.top"')
# Download a specific file from all those projects
# download_files(accessions, 'topology.top', './topologies')
