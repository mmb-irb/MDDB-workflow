from packaging.version import Version
from os import mkdir, listdir
from os.path import exists, isfile
from re import sub, compile, IGNORECASE
from shutil import rmtree

from mddb_workflow.utils.auxiliar import load_json, save_json, load_yaml, InputError
from mddb_workflow.utils.constants import PROJECT_METADATA_VERSION, TOPOLOGY_VERSION
from mddb_workflow.utils.constants import OUTPUT_METADATA_FILENAME, STANDARD_TOPOLOGY_FILENAME
from mddb_workflow.utils.constants import RED_HEADER, GREEN_HEADER, COLOR_END
from mddb_workflow.utils.constants import CHANNELS_ANALYSIS_VERSION, CLUSTERS_ANALYSIS_VERSION
from mddb_workflow.utils.constants import ENERGIES_ANALYSIS_VERSION, HBONDS_ANALYSIS_VERSION
from mddb_workflow.utils.constants import INTERACTIONS_ANALYSIS_VERSION, LIPID_INTERACTIONS_ANALYSIS_VERSION, LIPID_ORDER_ANALYSIS_VERSION
from mddb_workflow.utils.constants import RGYR_ANALYSIS_VERSION, RMSDS_ANALYSIS_VERSION, RMSF_ANALYSIS_VERSION
from mddb_workflow.utils.constants import SASA_ANALYSIS_VERSION, THICKNESS_ANALYSIS_VERSION, TMSCORES_ANALYSIS_VERSION
from mddb_workflow.utils.database import Database, get_available_nodes
from mddb_workflow.utils.loader import Loader
from mddb_workflow.utils.tasks import Task
from mddb_workflow.mwf import workflow, md_requestables

# Set the current version of every versioned analysis, by its workflow task flag
# Analyses with no version here are considered to be version 0.0.0
ANALYSIS_VERSIONS = {
    'channels': CHANNELS_ANALYSIS_VERSION,
    'clusters': CLUSTERS_ANALYSIS_VERSION,
    'energies': ENERGIES_ANALYSIS_VERSION,
    'hbonds': HBONDS_ANALYSIS_VERSION,
    'inter': INTERACTIONS_ANALYSIS_VERSION,
    'linter': LIPID_INTERACTIONS_ANALYSIS_VERSION,
    'lorder': LIPID_ORDER_ANALYSIS_VERSION,
    'rgyr': RGYR_ANALYSIS_VERSION,
    'rmsds': RMSDS_ANALYSIS_VERSION,
    'rmsf': RMSF_ANALYSIS_VERSION,
    'sas': SASA_ANALYSIS_VERSION,
    'thickness': THICKNESS_ANALYSIS_VERSION,
    'tmscore': TMSCORES_ANALYSIS_VERSION,
}
# Analysis output files are named after the task flag (e.g. mda.rgyr.json) and the loader names the analysis after the file
# Set here the legacy analyses whose output filename (mda.xxxx.json) does not match the task flag
ANALYSIS_FILENAME_EXCEPTIONS = {
    'dihedrals': 'dihenergies',
    'dist': 'dist_perres',
    'inter': 'interactions',
    'linter': 'lipid_inter',
    'lorder': 'lipid_order',
    'pairwise': 'rmsd_pairwise',
    'perres': 'rmsd_perres',
    'rmsf': 'fluctuation',
    'sas': 'sasa',
    'tmscore': 'tmscores',
}
# Set the flags of all MD analysis tasks in the workflow
# The interactions processing is not in the analyses module but it is versioned and uploaded as an analysis
ANALYSIS_TASKS = [ 'inter' ] + [ flag for flag, task in md_requestables.items()
    if isinstance(task, Task) and task.func.__module__.startswith('mddb_workflow.analyses.') ]
# Set the files to be loaded for every project task
PROJECT_TASK_FILES = {
    'pmeta': [ OUTPUT_METADATA_FILENAME ],
    'stopology': [ f'stopology/{STANDARD_TOPOLOGY_FILENAME}' ],
}
# Set the flags of all tasks whose version is tracked, so they may be updated
UPDATABLE_TASKS = list(PROJECT_TASK_FILES) + ANALYSIS_TASKS
# Set the label of missing analyses in the outdated versions
MISSING_VERSION = 'missing'
# Set a flag to update them all
ALL_FLAG = 'all'
# Set the name of the file where we keep track of the already updated projects
UPDATED_PROJECTS_FILENAME = 'updated_projects.json'
# Set the pattern used by the loader to find the inputs file in the project directory
LOADER_INPUTS_FILE_PATTERN = compile(r'inputs.(yaml|yml|json)$', IGNORECASE)


def get_analysis_database_names (flag : str) -> list[str]:
    """Get the possible names of an analysis in the database, given its task flag.
    If several names are returned then the first name available in the database is to be used.
    Note that analyses split by interaction (or cluster) have their version in the numbered analyses (e.g. energies-00)
    The overview analysis (e.g. energies) has no version, but it is the only one available in some legacy projects
    """
    filename_name = ANALYSIS_FILENAME_EXCEPTIONS.get(flag, flag)
    # The loader replaces every '_' by '-' in the analysis name
    name = filename_name.replace('_', '-')
    return [ f'{name}-00', name ]


def find_outdated_tasks (versions : dict, include : list[str] | None = None) -> dict[str, tuple[list[str], str]]:
    """Given the remote versions of a project, get the flags of the workflow tasks to be run again.
    Every flag comes with the outdated versions found in the database and the current version.
    Note that a task is outdated if it is outdated in at least one MD.
    Analyses which are missing in the database are not considered, since they may not apply to the project.
    However, analyses which are explicitly included are considered outdated when missing.
    """
    outdated_tasks = {}
    # Check project metadata and topology
    raw_metadata_version = versions.get('metadata', None)
    metadata_version = Version(raw_metadata_version) if raw_metadata_version else Version('0.0.0')
    if metadata_version < Version(PROJECT_METADATA_VERSION):
        outdated_tasks['pmeta'] = ([ str(metadata_version) ], PROJECT_METADATA_VERSION)
    raw_topology_version = versions.get('topology', None)
    topology_version = Version(raw_topology_version) if raw_topology_version else Version('0.0.0')
    if topology_version < Version(TOPOLOGY_VERSION):
        outdated_tasks['stopology'] = ([ str(topology_version) ], TOPOLOGY_VERSION)
    # Check the analyses of every MD
    for flag in ANALYSIS_TASKS:
        names = get_analysis_database_names(flag)
        updated_version = ANALYSIS_VERSIONS.get(flag, '0.0.0')
        outdated_versions = []
        for md in versions.get('mds', []):
            md_analyses = md.get('analyses', {})
            name = next((name for name in names if name in md_analyses), None)
            if name is None:
                if include and flag in include and MISSING_VERSION not in outdated_versions:
                    outdated_versions.append(MISSING_VERSION)
                continue
            raw_analysis_version = md_analyses[name]
            analysis_version = Version(raw_analysis_version) if raw_analysis_version else Version('0.0.0')
            if analysis_version < Version(updated_version) and str(analysis_version) not in outdated_versions:
                outdated_versions.append(str(analysis_version))
        if len(outdated_versions) > 0:
            outdated_tasks[flag] = (outdated_versions, updated_version)
    return outdated_tasks


def filter_tasks (
    outdated_tasks : dict[str, tuple[list[str], str]],
    include : list[str] | None = None,
    exclude : list[str] | None = None) -> dict[str, tuple[list[str], str]]:
    """Keep only the outdated tasks which are included (if any include is passed) and not excluded."""
    return { task: versions for task, versions in outdated_tasks.items()
        if (not include or task in include) and (not exclude or task not in exclude) }


def get_task_files (task : str) -> list[str]:
    """Get the files generated by a task, as paths or wildcards relative to the project directory.
    These are passed to the loader so only these files are uploaded.
    """
    project_files = PROJECT_TASK_FILES.get(task, None)
    if project_files: return project_files
    # Analyses write their outputs in a directory named after the task flag in every MD directory
    return [ f'*/{task}/*' ]


def check_inputs_accession (project_directory : str, accession : str):
    """Make sure the inputs file to be read by the loader has the expected accession.
    Otherwise the loader would create a new project instead of updating the current one.
    """
    inputs_filenames = [ filename for filename in listdir(project_directory)
        if isfile(f'{project_directory}/{filename}') and LOADER_INPUTS_FILE_PATTERN.search(filename) ]
    if len(inputs_filenames) == 0:
        raise RuntimeError(f'No inputs file was found in {project_directory}')
    for inputs_filename in inputs_filenames:
        inputs_accession = load_yaml(f'{project_directory}/{inputs_filename}').get('accession', None)
        if inputs_accession != accession:
            raise RuntimeError(f'Inputs file {inputs_filename} in {project_directory} has accession'
                f' "{inputs_accession}" while "{accession}" was expected')


def update_projects (
    node_url_or_alias : str,
    loader_directory : str | None = None,
    accessions : list[str] | None = None,
    query : str | None = None,
    trust : list[str] | bool = [],
    mercy : list[str] | bool = [],
    faith : bool = False,
    include : list[str] | None = None,
    exclude : list[str] | None = None):

    # If all nodes are requested then simply call this same function with every node
    if node_url_or_alias == ALL_FLAG:
        # Get the available nodes
        available_nodes = get_available_nodes()
        for node_alias in available_nodes.keys():
            update_projects(node_alias, loader_directory, accessions, query, trust, mercy, faith, include, exclude)
        return

    print(f'----- Updating projects from "{node_url_or_alias}" -----')

    # Instantiate the database handler
    database = Database(node_url_or_alias)

    # If we are to upload the updated projects to the corresponding database then prepare the loader
    loader = Loader(loader_directory, node_url_or_alias) if loader_directory else None

    # Get all available projects unless specific projects were requested
    # If a search query was passed then get only the matching projects
    if accessions and query:
        raise InputError('Specific projects and a search query can not be requested at the same time')
    if accessions:
        accessions_count = len(accessions)
    else:
        accessions_count, accessions = database.get_all_project_accessions(search=query)
    print(f'  Checking {accessions_count} projects' + (f' matching the query "{query}"' if query else ''))
    if include: print(f'  Only these tasks are checked: {", ".join(include)}')
    if exclude: print(f'  These tasks are not checked: {", ".join(exclude)}')

    # Make a directory for the current url or alias
    directory = sub('https?://', '', node_url_or_alias.replace('/api/', '')).strip('/')
    if not exists(directory): mkdir(directory)

    # Load the already updated projects in this directory, if any
    # For every accession we keep track of which task versions were updated and which of them were loaded
    # Thus there is no need to update them again and, if they were not loaded yet, we can load them in a later run
    # Note that the database is always checked, so tasks updated or loaded to an older version will be updated again
    updated_projects_filepath = f'{directory}/{UPDATED_PROJECTS_FILENAME}'
    updated_projects = load_json(updated_projects_filepath) if exists(updated_projects_filepath) else {}
    if type(updated_projects) is not dict or any(type(record['updated']) is not dict for record in updated_projects.values()):
        raise ValueError(f'{updated_projects_filepath} has a legacy format.'
            ' Now it must be a dict with the updated and loaded task versions of every accession. Please remove or convert it.')
    if len(updated_projects) > 0:
        print(f'  There are {len(updated_projects)} already updated projects in {updated_projects_filepath}')

    # Iterate projects
    for a, accession in enumerate(accessions, 1):
        project_directory = f'{directory}/{accession}'
        # Find which tasks are outdated in the database
        versions = database.get_project_versions(accession)
        # Tasks out of the include/exclude selection are ignored, even if they are outdated
        outdated_tasks = filter_tasks(find_outdated_tasks(versions, include), include, exclude)
        if len(outdated_tasks) == 0:
            print(f'--- {accession} ({a}/{accessions_count}) is up to date in the database ---')
            continue
        project_record = updated_projects.get(accession, { 'updated': {}, 'loaded': {} })
        # If a task was already loaded with the current version then the database should not be outdated
        inconsistent_tasks = [ task for task, (_, new_version) in outdated_tasks.items()
            if project_record['loaded'].get(task, None) == new_version ]
        if len(inconsistent_tasks) > 0:
            raise RuntimeError(f'Project {accession} is outdated in the database but it was already loaded'
                f' with the current version according to {updated_projects_filepath} -> {", ".join(inconsistent_tasks)}')
        # Update the outdated tasks unless they were already updated to the current version in a previous run
        tasks_to_update = [ task for task, (_, new_version) in outdated_tasks.items()
            if project_record['updated'].get(task, None) != new_version ]
        if len(tasks_to_update) > 0:
            print(f'--- Updating {accession} ({a}/{accessions_count}) ---')
            for task in tasks_to_update:
                old_versions, new_version = outdated_tasks[task]
                print(f'  {task}: {RED_HEADER}{", ".join(old_versions)}{COLOR_END} -> {GREEN_HEADER}{new_version}{COLOR_END}')
            # Tasks updated to an older version in a previous run have stale outputs which must be overwritten
            # Otherwise the workflow would find them and skip the task
            stale_tasks = [ task for task in tasks_to_update if task in project_record['updated'] ]
            # Run the workflow only for the outdated tasks
            if not exists(project_directory): mkdir(project_directory)
            workflow(
                project_parameters={
                    'directory': project_directory,
                    'accession': accession,
                    'database_url': database.url,
                    'trust': trust,
                    'mercy': mercy,
                    'faith': faith,
                },
                include=tasks_to_update,
                overwrite=stale_tasks,
            )
            # Keep track of the updated task versions
            for task in tasks_to_update:
                project_record['updated'][task] = outdated_tasks[task][1]
            updated_projects[accession] = project_record
            save_json(updated_projects, updated_projects_filepath)
        else:
            print(f'--- {accession} ({a}/{accessions_count}) is up to date locally ({", ".join(outdated_tasks)}) but not in the database ---')
        # Upload only the outdated tasks outputs, even if other dependencies were generated as well
        if not loader: continue
        tasks_to_load = list(outdated_tasks)
        print(f'--- Loading {accession} ({a}/{accessions_count}) -> {", ".join(tasks_to_load)} ---')
        if not exists(project_directory):
            raise RuntimeError(f'Project {accession} was updated but not loaded and its directory {project_directory} is missing')
        files_to_load = [ file for task in tasks_to_load for file in get_task_files(task) ]
        check_inputs_accession(project_directory, accession)
        if not loader.load(project_directory, overwrite=True, include=files_to_load):
            raise RuntimeError(f'Something went wrong when uploading {accession}')
        # Make sure the remote project is now up to date
        still_outdated_tasks = filter_tasks(find_outdated_tasks(database.get_project_versions(accession), include), include, exclude)
        if len(still_outdated_tasks) > 0:
            raise RuntimeError(f'Project {accession} is still outdated after uploading -> {", ".join(still_outdated_tasks)}')
        # Keep track of the loaded task versions
        for task in tasks_to_load:
            project_record['loaded'][task] = outdated_tasks[task][1]
        save_json(updated_projects, updated_projects_filepath)
        # Once uploaded, remove the project directory since it may be heavy (e.g. trajectories)
        rmtree(project_directory)
