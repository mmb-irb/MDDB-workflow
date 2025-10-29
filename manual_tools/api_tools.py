import json, urllib
from math import ceil
from mddb_workflow.console import DEFAULT_API_URL

# The URL to the API from our node
base_url = DEFAULT_API_URL+'rest/current'

def query_api (url : str, verbose=False) -> dict:
    """ Parse the URL in case it contains any HTTP control characters
    Replace white spaces by the corresponding percent notation character. """
    parsed_url = url.replace(" ", "%20")
    if verbose:
        print(f'Querying API URL: {parsed_url}')
    with urllib.request.urlopen(parsed_url) as response:
        return json.loads(response.read().decode("utf-8"))
    
def count_projects(query_url):
    # Ask which projects have a tpr file
    response = query_api(query_url)
    n_projects = response['filteredCount']
    print(f'We found {n_projects} projects')
    return n_projects

def paginate_projects(selection_query, limit = 100):
    # Ask which projects have a tpr file
    query_url = base_url + f'/projects?query={{{selection_query}}}'
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
        project_accessions = [ project['accession'] for project in projects ]
        accessions += project_accessions
        
    print(f'We have mined {len(accessions)} project accessions')
    return accessions, projects_list

def download_files(accessions: list, file: str, output_folder: str):
    # Store all JSON topologies in a list
    topologies = {}
    # Iterate the retrieved projects
    for a, accession in enumerate(accessions, 1):
        # Skip already parsed topologies
        if accession in topologies: continue
        print(f'Parsing topology from accession {accession} ({a}/{len(accessions)})          ', end='\r')
        # Download the TPR file
        topology_file_url = f'{base_url}/projects/{accession}/files/{file}'
        try:
            urllib.request.urlretrieve(topology_file_url, f'{output_folder}/{accession}_{file}')
        except:
            topologies[accession] = { 'error': f'Failed to download {topology_file_url}' }
            continue
    return topologies   

# Example of usage
# Find all projects with a topology.top file
# accessions = paginate_projects('"files.name":"topology.top"')
# Download a specific file from all those projects
# download_files(accessions, 'topology.top', './topologies')
