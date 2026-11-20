
import os
import json
import time
import signal
import hashlib
import multiprocessing
from pathlib import Path
from shutil import rmtree

from urllib.request import urlopen
from urllib.parse import urlencode
from urllib.error import HTTPError

from mddb_workflow.utils.auxiliar import load_json, save_json, retry_request, warn, get_available_cpus
from mddb_workflow.utils.auxiliar import RemoteServiceError, DatabaseUnavailableError
from mddb_workflow.utils.interproscan import run_local_interproscan
from mddb_workflow.utils.interproscan import take_directory_lock, stop_on_termination, raise_on_termination
from mddb_workflow.utils.interproscan import MAX_ATTEMPTS, RUN_TIMEOUT
from mddb_workflow.utils.type_hints import *

# Set the URL of the EBI InterProScan 6 service
INTERPROSCAN_URL = 'https://www.ebi.ac.uk/Tools/services/rest/iprscan6'

# Set the time to wait for the InterProScan jobs (seconds)
# AGUS: a veces ha llegado a tardar ~6 minutos que es excesivo, creo que  minutos es suficiente tiempo de espera
INTERPROSCAN_TIMEOUT = 300

# Set analysis version
# 0.2: InterProScan 6 results (library labels and some match fields changed from InterProScan 5)
CHAINS_VERSION = '0.2'


@retry_request
def request_interproscan(sequence: str) -> str:
    """Get the sequence and name of the chain in the structure and request the InterProScan."""
    # Set the request URL
    request_url = f'{INTERPROSCAN_URL}/run'
    # Set the POST data
    data = urlencode({
        'email': 'daniel.beltran@irbbarcelona.org',
        'title': 'Chain X',
        'sequence': f'>chain X\n{sequence}'
    }).encode()
    parsed_response = None
    try:
        with urlopen(request_url, data=data) as response:
            parsed_response = response.read().decode("utf-8")
    except HTTPError as error:
        print('Error:', error.read().decode())
        if error.code == 404:
            print(' Not found')
            return None
        else:
            raise error
    return parsed_response


def check_interproscan_overload(error: HTTPError, message: str):
    """Stop with a clear message when InterProScan is overloaded.
    InterProScan answers 400 when its own job system gets too many requests from all users.
    """
    if error.code == 400 and 'rate_limit_exceeded' in message:
        raise RemoteServiceError('InterProScan is overloaded right now. Please try again in a while.')


@retry_request
def check_interproscan_status(jobid: str) -> str:
    """Check the status of the InterProScan job."""
    # Set the request URL
    request_url = f'{INTERPROSCAN_URL}/status/{jobid}'
    try:
        with urlopen(request_url) as response:
            return response.read().decode("utf-8")
    # Show what InterProScan says and let the error go up so server errors are retried
    except HTTPError as error:
        message = error.read().decode()
        print(message)
        check_interproscan_overload(error, message)
        raise


@retry_request
def check_interproscan_result(jobid: str) -> dict:
    """Obtain the result of the InterProScan job in json format."""
    # Set the request URL
    request_url = f'{INTERPROSCAN_URL}/result/{jobid}/json'
    try:
        with urlopen(request_url) as response:
            return json.loads(response.read().decode("utf-8"))
    # Show what InterProScan says and let the error go up so server errors are retried
    except HTTPError as error:
        message = error.read().decode()
        print(message)
        check_interproscan_overload(error, message)
        raise


# Set the expected ligand data fields
CHAIN_DATA_FIELDS = set(['sequence', 'interproscan'])


def import_chains(chains_references_file: 'File') -> dict:
    """Import the chains data from a file if exists."""
    # Read the file
    imported_chains = load_json(chains_references_file.path)
    # Format data as the process expects to find it
    for imported_chain in imported_chains:
        for expected_field in CHAIN_DATA_FIELDS:
            if expected_field not in imported_chain:
                imported_chain[expected_field] = None
    return imported_chains


def clean_interproscan_result(interproscan_result: dict) -> dict:
    """Remove versions and pathways so Mongo don't get confused when they change."""
    interproscan_result.pop('interproscan-version', None)
    interproscan_result.pop('interpro-version', None)
    # RUBEN: creo que results siempre tiene un solo elemento, pero por si acaso iteramos
    for result in interproscan_result['results']:
        for match in result['matches']:
            if match['signature'].get('entry') is not None:
                match['signature']['entry'].pop('pathwayXRefs', None)
    return interproscan_result


def save_interproscan_result(sequence: str, interproscan_result: dict, chains_data: list, chains_references_file: 'File'):
    """Add the InterProScan result of a sequence to its chain data and save all chains data.
    Every result is saved as soon as it is done, so it is kept if the run is interrupted.
    """
    chain_data = next(data for data in chains_data if data['sequence'] == sequence)
    chain_data['version'] = CHAINS_VERSION
    chain_data['interproscan'] = clean_interproscan_result(interproscan_result)
    save_json(chains_data, chains_references_file.path)


def prepare_chain_references(
    structure: 'Structure',
    chains_references_file: 'File',
    database: 'Database',
    local_interproscan: bool = False,
):
    """Define the main function that will be called from the main script.
    This function will get the parsed chains from the structure and request the InterProScan service
    to obtain the data for each chain.
    If 'local_interproscan' is True, InterProScan runs locally instead of using the remote EBI service.
    """
    # Obtain the protein parsed chains from the structure
    protein_parsed_chains = structure.get_parsed_chains(only_protein=True)

    # Get unique sequences
    protein_sequences = set([chain['sequence'] for chain in protein_parsed_chains])

    print(f' Found {len(protein_parsed_chains)} protein chains with {len(protein_sequences)} unique sequences')

    # Save data from all chains to be saved in a file
    chains_data = []
    # Load the chains file if exists already
    if chains_references_file.exists:
        chains_data += import_chains(chains_references_file)

    # Save the jobids of every call to InterProScan
    interproscan_jobids = {}
    # Save the sequences to be run in the local InterProScan
    local_sequences = []

    # Iterate protein sequences
    for sequence in protein_sequences:
        # Check if the chain data already exists in the chains file
        chain_data = next((data for data in chains_data if data['sequence'] == sequence), None)
        # If we have no previous chain data then check if the sequence is already in the MDDB database
        if chain_data is None:
            try:
                chain_data = database.get_reference_data('chains', sequence)
            # If the database is not available then the chain data will be generated again
            except DatabaseUnavailableError:
                chain_data = None
            if chain_data is not None:
                chains_data.append(chain_data)
                # Save the chains data at this point
                # This may seem redundant since data will be not loaded further in the database
                # However, saving the chain in the local backup file is useufl to further run without internet connection
                save_json(chains_data, chains_references_file.path)
        # If we still have no chain data then create a new chain data dict
        # Set an object with the results of every call to InterProScan
        if chain_data is None:
            chain_data = {
                'sequence': sequence,
                'interproscan': None
            }
            chains_data.append(chain_data)
        # If chain data is missing any analysis then send a job
        # Request the InterProScan service (sequence length must be at least 3)
        # Keep the returned job ids to check the status and get the results later
        # Old references from the database (e.g. version 0.1 of short peptides) may have no 'interproscan' field at all
        if chain_data.get('interproscan') is None and len(chain_data['sequence']) > 2:
            if local_interproscan:
                local_sequences.append(sequence)
                continue
            interproscan_jobid = request_interproscan(sequence)
            interproscan_jobids[sequence] = interproscan_jobid

    # Run the local InterProScan for all missing sequences
    if len(local_sequences) > 0:
        print(f' Running InterProScan locally for {len(local_sequences)} sequences. It takes 20 seconds to 2 minutes per sequence...')
        for sequence in local_sequences:
            save_interproscan_result(sequence, run_local_interproscan(sequence), chains_data, chains_references_file)
        print(' Protein chains data obtained')
        return

    # Get the pending interpsocan jobids
    pending_jobids = list(interproscan_jobids.values())

    # If we already have the results of all the chains then we can skip the next steps
    if len(pending_jobids) == 0:
        print(' All reference chains are already in the backup file')
        # DANI: No es necesario devolver los datos, no se usa en el workflow
        # return chains_data
    # RUBEN: Separated functions so can be used in the references updater
    get_interproscan_results(pending_jobids, interproscan_jobids, chains_data, chains_references_file)


def get_interproscan_results(
    pending_jobids: list,
    interproscan_jobids: dict,
    chains_data: list,
    chains_references_file: 'File',
) -> None:
    start_time = time.time()
    # Jobs which fail in InterProScan are sent again, since failures are often temporary problems on their side
    MAX_RESUBMISSIONS = 2
    resubmissions = {}
    # Iterate over the jobids to check the status and get the results
    # If the status is 'FINISHED' then we can get the results and eliminate the jobid from the list
    # until there are no more jobids in either list
    while len(pending_jobids) >= 1:
        # Stop with an error so the task is not taken as done and the missing chains are requested in the next run
        # Note that chains already obtained are saved, so they are not requested again
        if time.time() - start_time > INTERPROSCAN_TIMEOUT:
            raise RemoteServiceError(f'{len(pending_jobids)} InterProScan jobs did not finish in '
                f'{INTERPROSCAN_TIMEOUT // 60} minutes. Run the workflow again later, or run InterProScan locally (-lip).')
        time.sleep(3)  # Wait for 3 seconds
        print(f' We are still waiting for {len(pending_jobids)} jobs to finish', end='\r')
        for sequence, interproscan_jobid in interproscan_jobids.items():
            # If the jobid is already processed then skip it
            if interproscan_jobid not in pending_jobids:
                continue
            status = check_interproscan_status(interproscan_jobid)
            # We only know four possible status for interproscan jobs, but it could be more
            if status == 'RUNNING' or status == 'PENDING' or status == 'QUEUED':
                continue
            # If the job failed then send it again
            # The wait time starts again so the new job has the full time to finish
            if status in ('ERROR', 'FAILURE') and resubmissions.get(sequence, 0) < MAX_RESUBMISSIONS:
                resubmissions[sequence] = resubmissions.get(sequence, 0) + 1
                new_jobid = request_interproscan(sequence)
                print(f' InterProScan job {interproscan_jobid} ended with status {status}.' +
                    f' Sending it again as {new_jobid} ({resubmissions[sequence]}/{MAX_RESUBMISSIONS})')
                interproscan_jobids[sequence] = new_jobid
                pending_jobids[pending_jobids.index(interproscan_jobid)] = new_jobid
                start_time = time.time()
                continue
            if status != 'FINISHED':
                raise RemoteServiceError(f'InterProScan job {interproscan_jobid} ended with status {status}.' +
                    f' Check the details in {INTERPROSCAN_URL}/result/{interproscan_jobid}/error')
            # Retrive the results from InterProScan and add them to the corresponding chain data
            interproscan_result = check_interproscan_result(interproscan_jobid)
            save_interproscan_result(sequence, interproscan_result, chains_data, chains_references_file)
            # Remove the jobid from the queue list
            pending_jobids.remove(interproscan_jobid)

    print(' Protein chains data obtained              ')


# Chain references are remade with the local InterProScan, which takes from 20 seconds to 2 minutes per sequence
# Every new reference is saved in a directory shared by all nodes and processes, in a file per sequence
# Thus every sequence is run once, although it is in several nodes or several processes share the work
CHAIN_REFERENCES_DIRECTORY = 'chain_references'
# A process claims a sequence before running it, so no other process runs it at the same time
# A claim of a process in another host is taken as abandoned when it is older than the longest possible run
CHAIN_CLAIM_EXPIRATION = MAX_ATTEMPTS * RUN_TIMEOUT + 600
# Set the result of a sequence which is being run by another process
CLAIMED = 'claimed'


def get_chain_reference_filepath(sequence: str) -> Path:
    """Get the file where the new reference of a sequence is saved."""
    return Path(CHAIN_REFERENCES_DIRECTORY) / f'{hashlib.md5(sequence.encode()).hexdigest()}.json'


def make_chain_reference(sequence: str) -> str | None:
    """Run the local InterProScan for a sequence in a worker process and save its new reference in the shared directory.
    Return the error, if any, or CLAIMED if another process is running it.
    """
    os.makedirs(CHAIN_REFERENCES_DIRECTORY, exist_ok=True)
    filepath = get_chain_reference_filepath(sequence)
    claim_directory = filepath.with_suffix('.claim')
    if take_directory_lock(claim_directory, CHAIN_CLAIM_EXPIRATION) is not None: return CLAIMED
    try:
        # Another process may have finished it while we were claiming it
        if filepath.exists(): return None
        # Sequences shorter than 3 residues are not run, as in the chains task
        interproscan = clean_interproscan_result(run_local_interproscan(sequence)) if len(sequence) > 2 else None
        reference = { 'sequence': sequence, 'version': CHAINS_VERSION, 'interproscan': interproscan }
        # Write it with a temporary name so an interrupted write is never taken as a reference
        incomplete_filepath = filepath.with_suffix(f'.{os.getpid()}.part')
        save_json(reference, str(incomplete_filepath))
        os.replace(incomplete_filepath, filepath)
    except Exception as error:
        return f'{type(error).__name__}: {error}'
    finally:
        rmtree(claim_directory, ignore_errors=True)
    return None


def assign_cpus(cpu_sets: 'multiprocessing.Queue'):
    """Restrict a worker process to its own CPUs, so the InterProScan runs at once do not compete for them.
    Workers are terminated when the pool is stopped, so they stop their InterProScan run and release its claim.
    """
    cpus = cpu_sets.get()
    # With only a number of CPUs, the run is told to use this many, as if it was a SLURM job of its own
    if isinstance(cpus, int): os.environ['SLURM_CPUS_PER_TASK'] = str(cpus)
    else: os.sched_setaffinity(0, cpus)
    signal.signal(signal.SIGTERM, stop_on_termination)


def make_chain_references(sequences: list[str], workers: int = 1):
    """Make the new references which are missing in the shared directory, with several InterProScan runs at once.
    Runs are done in worker processes, which are forked so they keep the decision about the Matches API made by this process.
    The CPUs of this process (or of its SLURM job) are split among the runs.
    """
    missing_sequences = [sequence for sequence in sequences if not get_chain_reference_filepath(sequence).exists()]
    print(f'  Running InterProScan for {len(missing_sequences)} sequences with no new reference yet, {workers} at once')
    if len(missing_sequences) == 0: return
    cpus = sorted(os.sched_getaffinity(0))
    available_cpus = get_available_cpus()
    cpus_per_worker = max(1, available_cpus // workers)
    context = multiprocessing.get_context('fork')
    cpu_sets = context.Queue()
    for worker in range(workers):
        # SLURM does not always pin the job to its CPUs, and then the CPUs it may run on are shared with other jobs
        # Then workers are not pinned to any of them, but every run uses only its share of the CPUs of the job
        if len(cpus) > available_cpus: cpu_sets.put(cpus_per_worker)
        else: cpu_sets.put(cpus[worker * cpus_per_worker:(worker + 1) * cpus_per_worker] or cpus)
    # When terminated (e.g. by scancel), stop the InterProScan runs and release their claims
    # Note that the pool stops its workers when leaving, also if this process is interrupted
    failed_count = claimed_count = 0
    with raise_on_termination(), context.Pool(workers, initializer=assign_cpus, initargs=(cpu_sets,)) as pool:
        errors = pool.imap(make_chain_reference, missing_sequences)
        # Show the progress of the runs and warn about the failed ones as they end
        for s, (sequence, error) in enumerate(zip(missing_sequences, errors), 1):
            if error == CLAIMED:
                claimed_count += 1
            elif error:
                failed_count += 1
                warn(f'InterProScan failed for sequence {sequence[:20]}...: {error}')
            print(f'  InterProScan done for {s}/{len(missing_sequences)} sequences '
                f'({failed_count} failed, {claimed_count} run by another process)')
