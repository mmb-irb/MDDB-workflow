"""Read the provenance of the simulations from the RO-Crates written by the BioBB.

The provenance may come in several pieces: one crate per workflow and one crate per step.
Pieces are joined through the files: a file written by a step and read by another step has the same cksum.
Then we walk back from the files of every MD to the steps which produced them.
Steps needed by all MDs go to the project metadata and the rest go to the metadata of each MD.
"""
from glob import glob
from os.path import exists, relpath

from mddb_workflow.utils.auxiliar import InputError, load_json
from mddb_workflow.utils.file import File
from mddb_workflow.utils.quickhash import QUICKHASH_SAMPLE_SIZE
from mddb_workflow.utils.type_hints import *

# Name of every RO-Crate file
CRATE_FILENAME = 'ro-crate-metadata.json'

# Step parameters about the machine where the step was run, which tell nothing about the simulation
EXECUTION_PARAMETERS = {
    'binary_path', 'sandbox_path', 'remove_tmp', 'restart',
    'container_image', 'container_path', 'container_shell_path', 'container_volume_path',
    'num_threads', 'num_threads_mpi', 'num_threads_omp', 'num_threads_omp_pme',
    'use_gpu', 'gpu_id', 'gpu_tasks',
}

# Steps skipped by the BioBB restart property are recorded afterwards from their outputs
# They have no real times and their description starts like this
RECORDED_STEP_DESCRIPTION = 'Skipped by the restart property'

# Fields of the machine where a step ran which we do not keep
IGNORED_ENVIRONMENT_FIELDS = {'slurm_job_id'}


def find_crate_files(provenance_directory: str) -> list[File]:
    """Find all RO-Crate files in the provenance directory."""
    if not exists(provenance_directory): return []
    filepaths = glob(f'{provenance_directory}/**/{CRATE_FILENAME}', recursive=True)
    return [File(filepath) for filepath in sorted(filepaths)]


def get_provenance(
    provenance_crate_files: list['File'],
    provenance_md_files: dict[str, dict],
) -> Optional[dict]:
    """Read the provenance RO-Crates and split the steps into the common ones and the ones of every MD."""
    if len(provenance_crate_files) == 0:
        print(' There is no provenance')
        return None
    steps, files = read_crates(provenance_crate_files)
    md_step_ids = find_md_steps(steps, provenance_md_files)
    # Write the steps of every part along with the files they use
    # Files of the common steps are written only in the project part
    project_step_ids = [ step_id for step_id in steps
        if all(step_id in step_ids for step_ids in md_step_ids.values()) ]
    project_part = get_part(project_step_ids, steps, files)
    project_files = set([ file['id'] for file in project_part['files'] ])
    mds_part = {}
    for md_key, md_files in provenance_md_files.items():
        step_ids = [ step_id for step_id in steps
            if step_id in md_step_ids[md_key] and step_id not in project_step_ids ]
        md_part = get_part(step_ids, steps, files, skipped_files=project_files)
        # Say which files are the ones used by the MD in the workflow
        md_part['md_files'] = {
            'topology': md_files['topology'].get_cksum(unsafe=True),
            'structure': md_files['structure'].get_cksum(unsafe=True),
            'trajectories': [ file.get_cksum(unsafe=True) for file in md_files['trajectories'] ],
        }
        mds_part[md_key] = md_part
    return { 'project': project_part, 'mds': mds_part }


def get_part(step_ids: list[str], steps: dict, files: dict, skipped_files: set = set()) -> dict:
    """Get the steps and the files they use, except for the skipped files."""
    file_ids = []
    for step_id in step_ids:
        step = steps[step_id]
        for file_id in step['inputs'] + step['outputs']:
            if file_id in skipped_files or file_id in file_ids: continue
            file_ids.append(file_id)
    return {
        'steps': [ steps[step_id] for step_id in step_ids ],
        'files': [ files[file_id] for file_id in file_ids ],
    }


def read_crates(crate_files: list['File']) -> tuple[dict, dict]:
    """Read the steps and files in all RO-Crates.
    Files are identified by their cksum so the steps of different crates are linked.
    The same step may be in several crates (the one of the step and the one of its workflow) so it is read once.
    """
    steps = {}
    files = {}
    # Read the workflow crates first, since they also say the workflow every step belongs to
    # Crates of every step are in a deeper directory than the crate of their workflow
    sorted_crate_files = sorted(crate_files, key=lambda file: (file.path.count('/'), file.path))
    for crate_file in sorted_crate_files:
        crate = load_json(crate_file.path)
        graph = { entity['@id']: entity for entity in crate['@graph'] }
        # Get the workflow name, if this is the crate of a workflow
        workflow_run = graph.get('#workflow-run', None)
        workflow = workflow_run['name'].removeprefix('Run of ') if workflow_run else None
        for entity in graph.values():
            if entity.get('@type') != 'CreateAction': continue
            if entity['@id'] == '#workflow-run': continue
            step_id = entity['@id'].removeprefix('#')
            if step_id in steps: continue
            steps[step_id] = read_step(entity, graph, workflow, files, crate_file)
    return steps, files


def read_step(action: dict, graph: dict, workflow: Optional[str], files: dict, crate_file: 'File') -> dict:
    """Read a step and add its files to the files dict."""
    instrument = graph[action['instrument']['@id']]
    # Steps skipped by the restart property are recorded from their outputs, so their times are not real
    is_recorded = action.get('description', '').startswith(RECORDED_STEP_DESCRIPTION)
    status = 'recorded' if is_recorded else get_status(action['actionStatus'])
    # Read the inputs: files and parameters
    inputs = []
    parameters = {}
    for reference in as_list(action.get('object', None)):
        entity = graph[reference['@id']]
        entity_type = entity['@type']
        if entity_type == 'File':
            inputs.append(read_file(entity, graph, files, crate_file))
        elif entity_type == 'PropertyValue':
            if entity['name'] in EXECUTION_PARAMETERS: continue
            parameters[entity['name']] = read_parameter(entity, graph)
        else: raise ValueError(f'Unexpected input type "{entity_type}" in step {action["@id"]} from {crate_file.path}')
    # Read the outputs, which are always files
    outputs = []
    for reference in as_list(action.get('result', None)):
        entity = graph[reference['@id']]
        if entity['@type'] != 'File':
            raise ValueError(f'Unexpected output type "{entity["@type"]}" in step {action["@id"]} from {crate_file.path}')
        outputs.append(read_file(entity, graph, files, crate_file))
    step = {
        'id': action['@id'].removeprefix('#'),
        'tool': instrument['name'],
        'status': status,
        'inputs': inputs,
        'outputs': outputs,
        'parameters': parameters,
    }
    # Add the optional fields only when we have them
    if workflow: step['workflow'] = workflow
    version = instrument.get('softwareVersion', None)
    if version: step['version'] = version
    url = instrument['@id']
    if url.startswith('http'): step['url'] = url
    if not is_recorded:
        if 'startTime' in action: step['start'] = action['startTime']
        if 'endTime' in action: step['end'] = action['endTime']
    # Read who ran the step and where
    environment = read_environment(action, graph)
    username = environment.pop('username', None)
    user = read_agent(action, graph) or username
    if user: step['user'] = user
    if environment: step['environment'] = environment
    return step


def read_agent(action: dict, graph: dict) -> Optional[str]:
    """Get the name of the person who ran the step, if it is declared."""
    reference = action.get('agent', None)
    if reference is None: return None
    return graph[reference['@id']]['name']


def read_environment(action: dict, graph: dict) -> dict:
    """Get the machine where the step ran (hostname, cpus, cluster, etc.).
    Old BioBB versions only write here a link to their documentation, which tells nothing.
    """
    reference = action.get('executionEnvironment', None)
    if not isinstance(reference, dict): return {}
    environment = read_parameter(graph[reference['@id']], graph)
    return { key: value for key, value in environment.items() if key not in IGNORED_ENVIRONMENT_FIELDS }


def get_status(action_status: str) -> str:
    """Get a short status from the schema.org action status.
    e.g. http://schema.org/CompletedActionStatus -> completed
    """
    return action_status.split('/')[-1].removesuffix('ActionStatus').lower()


def read_parameter(entity: dict, graph: dict):
    """Read the value of a parameter.
    Parameters with several values (e.g. a parsed mdp file) refer to one parameter per value.
    """
    references = entity.get('valueReference', None)
    if references is None: return entity.get('value', None)
    return { graph[reference['@id']]['name']: graph[reference['@id']].get('value', None)
        for reference in as_list(references) }


def read_file(entity: dict, graph: dict, files: dict, crate_file: 'File') -> str:
    """Read a file, add it to the files dict and return its id: our cksum."""
    quickhash_reference = entity.get('quickhash', None)
    if quickhash_reference is None:
        raise InputError(f'File "{entity["name"]}" has no quickhash in {crate_file.path}\n'
            ' The provenance must come from a BioBB version which writes the quickhash of every file.')
    quickhash = graph[quickhash_reference['@id']]
    if quickhash['algo'] != 'xxh64-sample' or quickhash['version'] != 1 or quickhash['sample_size'] != QUICKHASH_SAMPLE_SIZE:
        raise InputError(f'Not supported quickhash of file "{entity["name"]}" in {crate_file.path}: {quickhash}')
    # Build the same cksum we get with File.get_cksum(unsafe=True)
    size = int(entity['contentSize'])
    is_unsafe = size >= 3 * QUICKHASH_SAMPLE_SIZE
    cksum = f'{size}-{"(UNSAFE)" if is_unsafe else ""}{quickhash["value"]}'
    # The same content may have several names (e.g. a renamed file)
    file = files.get(cksum, None)
    if file is None:
        file = { 'id': cksum, 'names': [], 'size': size }
        file_format = entity.get('encodingFormat', None)
        if file_format: file['format'] = file_format['@id'] if type(file_format) == dict else file_format
        files[cksum] = file
    if entity['name'] not in file['names']: file['names'].append(entity['name'])
    return cksum


def find_md_steps(steps: dict, provenance_md_files: dict[str, dict]) -> dict[str, set]:
    """Find the steps of every MD and mark the ones which are hidden by default:
    - Redundant steps: all their outputs are identical to some input, so they add nothing.
    - Branch steps: they do not lead to any MD file.
    """
    # Find redundant steps
    redundant_step_ids = set()
    for step_id, step in steps.items():
        if len(step['outputs']) == 0: continue
        if all(output in step['inputs'] for output in step['outputs']):
            redundant_step_ids.add(step_id)
            step['redundant'] = True
    # Find which steps produce every file, apart from the redundant ones
    # Otherwise a redundant step would be the producer of its own input
    producers = {}
    for step_id, step in steps.items():
        if step_id in redundant_step_ids: continue
        for output in step['outputs']:
            producers.setdefault(output, []).append(step_id)
    # Walk back from the files of every MD to the steps which produced them
    md_step_ids = {}
    for md_key, md_files in provenance_md_files.items():
        md_file_list = [ md_files['topology'], md_files['structure'], *md_files['trajectories'] ]
        for md_file in md_file_list:
            if md_file.get_cksum(unsafe=True) in producers: continue
            raise InputError(f'File {md_file.path} of MD {md_key} was not written by any step in the provenance.\n'
                ' Use the very same files written by the provenance steps.')
        start = [ md_file.get_cksum(unsafe=True) for md_file in md_file_list ]
        md_step_ids[md_key] = walk_back(start, steps, producers)
    # Add the steps which do not lead to any MD file to the MDs whose files they read
    # Mark them as branches, unless they are already marked as redundant
    md_keys = set(md_step_ids)
    reached_step_ids = set().union(*md_step_ids.values())
    other_step_ids = [ step_id for step_id in steps if step_id not in reached_step_ids ]
    for step_id in other_step_ids:
        if step_id not in redundant_step_ids: steps[step_id]['branch'] = True
    for step_id, step_md_keys in assign_other_steps(other_step_ids, steps, md_step_ids, md_keys).items():
        for md_key in step_md_keys: md_step_ids[md_key].add(step_id)
    return md_step_ids


def walk_back(start: list[str], steps: dict, producers: dict) -> set:
    """Find the steps which lead to the start files.
    When a file was written by several steps one is enough: the one already in the graph, if any.
    Otherwise, the one with less new inputs, so we do not drag the steps of other MDs.
    """
    step_ids = set()
    seen = set(start)
    pending = list(start)
    while len(pending) > 0:
        # Files with a single producer go first, so we know which steps are in the graph before choosing
        pending.sort(key=lambda file_id: len(producers.get(file_id, [])))
        file_id = pending.pop(0)
        candidates = producers.get(file_id, [])
        if len(candidates) == 0: continue
        if any(candidate in step_ids for candidate in candidates): continue
        step_id = min(candidates, key=lambda candidate: sum(i not in seen for i in steps[candidate]['inputs']))
        step_ids.add(step_id)
        for input_id in steps[step_id]['inputs']:
            if input_id in seen: continue
            seen.add(input_id)
            pending.append(input_id)
    return step_ids


def assign_other_steps(step_ids: list[str], steps: dict, md_step_ids: dict[str, set], md_keys: set) -> dict[str, set]:
    """Find the MDs of the steps which do not lead to any MD file.
    A step belongs to the MDs which have all its inputs.
    Files read by no MD step (e.g. the first inputs) belong to all MDs.
    """
    # Find the MDs of every file used by the MD steps
    file_md_keys = {}
    for md_key, md_steps in md_step_ids.items():
        for step_id in md_steps:
            step = steps[step_id]
            for file_id in step['inputs'] + step['outputs']:
                file_md_keys.setdefault(file_id, set()).add(md_key)
    # Some of these steps may read the outputs of others, so wait for those to be assigned first
    def waits_for_others(step_id: str, pending: list[str]) -> bool:
        inputs = steps[step_id]['inputs']
        return any(input_id in steps[other_id]['outputs']
            for other_id in pending if other_id != step_id for input_id in inputs)
    assigned = {}
    pending = list(step_ids)
    while len(pending) > 0:
        ready = [ step_id for step_id in pending if not waits_for_others(step_id, pending) ]
        # Steps reading each other outputs in a loop belong to all MDs
        if len(ready) == 0: ready = pending
        for step_id in ready:
            input_md_keys = [ file_md_keys.get(input_id, md_keys) for input_id in steps[step_id]['inputs'] ]
            step_md_keys = set.intersection(md_keys, *input_md_keys)
            # A step reading the files of several MDs but none of all of them belongs to all those MDs
            if len(step_md_keys) == 0: step_md_keys = set.union(*input_md_keys)
            assigned[step_id] = step_md_keys
            for output_id in steps[step_id]['outputs']:
                file_md_keys.setdefault(output_id, set()).update(step_md_keys)
        pending = [ step_id for step_id in pending if step_id not in assigned ]
    return assigned


def as_list(value) -> list:
    """Make a list of a JSON-LD value, which may be missing, a single item or a list."""
    if value is None: return []
    if type(value) == list: return value
    return [value]
