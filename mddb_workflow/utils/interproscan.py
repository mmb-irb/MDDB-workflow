"""Run InterProScan locally, instead of the remote EBI service.

InterProScan 6 is a Nextflow pipeline whose analyses run in containers (Apptainer, Singularity or Docker).
It is the same engine as the remote EBI service (iprscan6), so local results have the same JSON format.
Licensed analyses (Phobius, SignalP) and the machine learning model InterPro-N, which EBI runs, are not run locally.
Known sequences take their precalculated matches from the InterPro Matches API, which InterProScan queries,
but only when the API serves the same InterPro release as the local data, so all results come from one release.
Otherwise (e.g. after a new InterPro release, or with no internet) all matches are calculated locally.

The pipeline, the member databases data (about 35 GB) and the container images are set up once with:
    python -m mddb_workflow.utils.interproscan --setup
"""
import os
import sys
import json
import time
import signal
import socket
from pathlib import Path
from functools import cache
from contextlib import contextmanager
from shutil import which, rmtree, copyfileobj
from subprocess import run, Popen, PIPE, TimeoutExpired, CompletedProcess
from tempfile import TemporaryDirectory, mkdtemp
from urllib.request import urlopen

from mddb_workflow.utils.auxiliar import ToolError, warn, get_available_cpus
from mddb_workflow.utils.constants import RESOURCES_DIRECTORY_PATH


INTERPROSCAN_VERSION = '6.0.3.0'
# Pin the InterPro data release so local results do not change with every new release
INTERPRO_VERSION = '111.0'
INTERPROSCAN_DIRECTORY = Path(RESOURCES_DIRECTORY_PATH) / 'interproscan'
PIPELINE_DIRECTORY = INTERPROSCAN_DIRECTORY / f'interproscan6-{INTERPROSCAN_VERSION}'
PIPELINE_URL = f'https://github.com/ebi-pf-team/interproscan6/archive/refs/tags/{INTERPROSCAN_VERSION}.tar.gz'
DATA_DIRECTORY = INTERPROSCAN_DIRECTORY / f'data-{INTERPRO_VERSION}'
MATCHES_API_URL = 'https://www.ebi.ac.uk/interpro/matches/api'
IMAGES_DIRECTORY = INTERPROSCAN_DIRECTORY / 'images'
# Nextflow and Apptainer keep their own files in the home directory by default
# Keep them in the installation instead, so runs find them with no internet and home quotas are not filled
NEXTFLOW_HOME = INTERPROSCAN_DIRECTORY / 'nextflow'
APPTAINER_CACHE = INTERPROSCAN_DIRECTORY / '.apptainer'
# This file is written once the setup is complete
SETUP_MARKER = INTERPROSCAN_DIRECTORY / f'.setup-{INTERPROSCAN_VERSION}-{INTERPRO_VERSION}'
# The setup lock is a directory, since creating a directory is atomic also in network filesystems (e.g. NFS)
SETUP_LOCK = INTERPROSCAN_DIRECTORY / '.setup.lock'
# Set the command to set up the local installation, to be shown in errors
SETUP_COMMAND = 'python -m mddb_workflow.utils.interproscan --setup'
# Container runtimes supported by InterProScan, by order of preference (same names as their Nextflow profiles)
CONTAINER_RUNTIMES = ['apptainer', 'singularity', 'docker']
# Number of times a local run is tried for every sequence
# About 1 run in 20 fails in its first seconds by a race inside Nextflow, so trying again is cheap
MAX_ATTEMPTS = 5
# Time limit of a local run (seconds), so a run which hangs is stopped and tried again
RUN_TIMEOUT = 3600
# Ubiquitin, used to run every analysis once during the setup
SETUP_SEQUENCE = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'


def find_container_runtime() -> str | None:
    """Find the first available container runtime."""
    return next((runtime for runtime in CONTAINER_RUNTIMES if which(runtime)), None)


def check_requirements():
    """Make sure Nextflow and a container runtime are available, or explain how to get them."""
    if not which('nextflow'):
        raise ToolError('Cannot find nextflow, which is required by InterProScan 6 (Nextflow 25.10.4 or newer). '
            'Install it with conda install -c bioconda nextflow')
    if not find_container_runtime():
        raise ToolError('Cannot find a container runtime, which InterProScan 6 needs to run its analyses. '
            f'Make sure one of {", ".join(CONTAINER_RUNTIMES)} is in the PATH')


def local_interproscan_exists() -> bool:
    """Check if the local InterProScan is set up, so runs need no download.
    Check it before launching many workflow runs at once (e.g. a SLURM array) with '-lip'.
    """
    return SETUP_MARKER.exists()


def get_matches_api_release() -> str | None:
    """Get the InterPro release served by the Matches API, or None if it is not reachable."""
    try:
        with urlopen(f'{MATCHES_API_URL}/info', timeout=10) as response:
            return json.loads(response.read().decode()).get('release')
    except Exception:
        return None


@cache
def prepare_local_interproscan() -> bool:
    """Make sure the local InterProScan is set up and can run, or explain how to fix it.
    Then check if the precalculated matches of the Matches API may be used, and return it.
    They are used only when the API serves the same InterPro release as the local data,
    since InterProScan does not check it, and the results of a project would otherwise mix two releases.
    It is done once per process, so all chains of a project or a references update are run the same way.
    """
    if not local_interproscan_exists():
        raise ToolError(f'The local InterProScan is not set up. Set it up once with "{SETUP_COMMAND}" '
            '(about 37 GB of disk)')
    check_requirements()
    release = get_matches_api_release()
    if release == INTERPRO_VERSION: return True
    if release is None:
        warn('The InterPro Matches API is not reachable, so all InterProScan matches are calculated locally')
    else:
        warn(f'The InterPro Matches API serves InterPro {release} while the local InterProScan data are InterPro '
            f'{INTERPRO_VERSION}, so all InterProScan matches are calculated locally. To use precalculated matches '
            'again, set INTERPRO_VERSION to the new release in mddb_workflow/utils/interproscan.py and run the setup again')
    return False


def get_last_lines(text: str, count: int) -> str:
    """Get the last lines of a text (e.g. logs), to be shown in errors."""
    return '\n'.join(text.strip().split('\n')[-count:])


def run_process(command: list[str], cwd: Path, env: dict, timeout: int | None = None) -> CompletedProcess:
    """Run a command in its own process group, so its whole process tree (e.g. containers) is stopped
    when it takes longer than the timeout or when this process is interrupted.
    """
    def stop_process_group(signal_number: int):
        try:
            os.killpg(process.pid, signal_number)
        except ProcessLookupError:
            pass
    process = Popen(command, cwd=cwd, env=env, stdout=PIPE, stderr=PIPE, start_new_session=True)
    try:
        stdout, stderr = process.communicate(timeout=timeout)
    except BaseException as error:
        # Give the processes some time to end on their own before killing them
        stop_process_group(signal.SIGTERM)
        try:
            stdout, stderr = process.communicate(timeout=60)
        except TimeoutExpired:
            stop_process_group(signal.SIGKILL)
            stdout, stderr = process.communicate()
        if not isinstance(error, TimeoutExpired): raise
        raise ToolError(f'Local InterProScan did not finish in {timeout // 60} minutes:\n'
            + get_last_lines(stdout.decode() + stderr.decode(), 30))
    return CompletedProcess(command, process.returncode, stdout, stderr)


def run_pipeline(fasta_filepath: Path, run_directory: Path, data_directory: Path,
    use_matches_api: bool, check_interpro_version: bool = False, timeout: int | None = RUN_TIMEOUT) -> dict:
    """Run the InterProScan pipeline in a directory where Nextflow writes its logs and intermediate files.
    Analyses are run like the remote service does by default, with GO terms and pathways.
    Return the parsed JSON output.
    """
    container_runtime = find_container_runtime()
    command = ['nextflow', 'run', str(PIPELINE_DIRECTORY), '-profile', container_runtime,
        '-work-dir', str(run_directory / 'work'),
        '--input', str(fasta_filepath), '--datadir', str(data_directory), '--interpro', INTERPRO_VERSION,
        '--formats', 'json', '--outdir', str(run_directory), '--outprefix', 'chain',
        '--goterms', '--pathways', '--max-workers', str(get_available_cpus())]
    # The version check needs internet, so it is done only in the setup
    if not check_interpro_version: command.append('--skip-interpro-version-check')
    # Use the same Matches API whose release was checked
    if use_matches_api: command += ['--matches-api-url', MATCHES_API_URL]
    else: command.append('--no-matches-api')
    # Docker containers run as root by default, which leaves files the run directory can not be removed with
    if container_runtime == 'docker':
        docker_config_filepath = run_directory / 'docker.config'
        docker_config_filepath.write_text(f"docker.runOptions = '-u {os.getuid()}:{os.getgid()}'\n")
        command += ['-c', str(docker_config_filepath)]
    # Keep the files of Nextflow and the container images in the installation
    # Note that Apptainer and Singularity read the variables with their own names
    environment = { **os.environ, 'NXF_HOME': str(NEXTFLOW_HOME), 'NXF_ANSI_LOG': 'false',
        'NXF_APPTAINER_CACHEDIR': str(IMAGES_DIRECTORY), 'NXF_SINGULARITY_CACHEDIR': str(IMAGES_DIRECTORY),
        'APPTAINER_CACHEDIR': str(APPTAINER_CACHE), 'SINGULARITY_CACHEDIR': str(APPTAINER_CACHE) }
    # Keep Nextflow from looking for updates online
    if not check_interpro_version: environment['NXF_OFFLINE'] = 'true'
    process = run_process(command, cwd=run_directory, env=environment, timeout=timeout)
    output_filepath = run_directory / 'chain.json'
    if process.returncode != 0 or not output_filepath.exists():
        # The Nextflow log has the details of the error (e.g. the stack trace), but it is removed with the run directory
        nextflow_log_filepath = run_directory / '.nextflow.log'
        nextflow_log = ('\nEnd of the Nextflow log:\n' + get_last_lines(nextflow_log_filepath.read_text(), 40)
            if nextflow_log_filepath.exists() else '')
        raise ToolError('Local InterProScan failed:\n'
            + get_last_lines(process.stdout.decode() + process.stderr.decode(), 30) + nextflow_log)
    with open(output_filepath) as file:
        return json.load(file)


def take_directory_lock(lock_directory: Path, max_age: float | None = None) -> list[str] | None:
    """Take a lock which is a directory, since creating a directory is atomic also in network filesystems (e.g. NFS).
    The host, process and time of the owner are written inside, so a lock left by a process which was killed is found:
    in this host because its process is not running, and in any host because it is older than 'max_age' (seconds).
    Return None if the lock was taken, or the owner of the lock otherwise (empty if unknown).
    """
    for _ in range(2):
        try:
            os.mkdir(lock_directory)
        except FileExistsError:
            # Note that the owner may be missing if the lock is being taken right now
            owner_filepath = lock_directory / 'owner'
            owner = owner_filepath.read_text().split() if owner_filepath.exists() else []
            # Signal 0 only checks if the process is running (a process of another user can not be signaled, but it runs)
            killed = False
            if len(owner) == 3 and owner[0] == socket.gethostname():
                try:
                    os.kill(int(owner[1]), 0)
                except ProcessLookupError:
                    killed = True
                except PermissionError:
                    pass
            expired = len(owner) == 3 and max_age is not None and time.time() - float(owner[2]) > max_age
            if not killed and not expired: return owner
            warn(f'Removing a lock of a process which is not running anymore: {lock_directory}')
            rmtree(lock_directory, ignore_errors=True)
            continue
        (lock_directory / 'owner').write_text(f'{socket.gethostname()} {os.getpid()} {time.time()}\n')
        return None
    return []


def stop_on_termination(signal_number: int, _frame):
    """Stop with an exception when terminated (e.g. by scancel), so locks are released and InterProScan is stopped."""
    raise SystemExit(f'Terminated by signal {signal_number}')


@contextmanager
def raise_on_termination():
    """Stop with an exception when terminated while in this context (signal handlers can be set in the main thread only)."""
    try:
        previous_handler = signal.signal(signal.SIGTERM, stop_on_termination)
    except ValueError:
        previous_handler = None
    try:
        yield
    finally:
        if previous_handler is not None: signal.signal(signal.SIGTERM, previous_handler)


def setup_local_interproscan():
    """Download the InterProScan pipeline, the member databases data and the container images, if not done yet.
    Data (about 35 GB) and images (about 1.5 GB) are downloaded by InterProScan itself, running every analysis once.
    It takes about 40 minutes, mostly downloading from EBI.
    Data are downloaded in a temporary directory which is renamed at the end,
    so an interrupted setup never leaves incomplete data behind.
    """
    if local_interproscan_exists():
        print(f'InterProScan {INTERPROSCAN_VERSION} with InterPro {INTERPRO_VERSION} is already set up at {INTERPROSCAN_DIRECTORY}')
        return
    check_requirements()
    os.makedirs(INTERPROSCAN_DIRECTORY, exist_ok=True)
    # Make sure no other setup is running at the same time, since they would download the same files
    # A lock left by a setup which was killed in this host (e.g. by SLURM) is removed
    owner = take_directory_lock(SETUP_LOCK)
    if owner is not None:
        owner_label = f' in host {owner[0]} (process {owner[1]})' if len(owner) == 3 else ''
        raise ToolError(f'Another InterProScan setup is running{owner_label}. '
            f'If it is not, remove {SETUP_LOCK} and try again.')
    # Release the lock also when the setup is terminated
    with raise_on_termination():
        try:
            # Another setup may have finished while we were taking the lock
            if local_interproscan_exists():
                print(f'InterProScan {INTERPROSCAN_VERSION} with InterPro {INTERPRO_VERSION} is already set up at {INTERPROSCAN_DIRECTORY}')
                return
            if not (PIPELINE_DIRECTORY / 'main.nf').exists():
                print(f'Downloading InterProScan {INTERPROSCAN_VERSION}...')
                tarball = INTERPROSCAN_DIRECTORY / f'interproscan6-{INTERPROSCAN_VERSION}.tar.gz'
                incomplete_directory = Path(mkdtemp(prefix='.incomplete_', dir=INTERPROSCAN_DIRECTORY))
                # The tarball is removed at the end, also when the download is interrupted
                try:
                    with urlopen(PIPELINE_URL) as response, open(tarball, 'wb') as file:
                        copyfileobj(response, file, 2**24)
                    process = run(['tar', '-xzf', str(tarball), '-C', str(incomplete_directory)], stdout=PIPE, stderr=PIPE)
                    if process.returncode != 0:
                        raise ToolError(f'Failed to extract {tarball}:\n' + process.stderr.decode())
                    os.rename(incomplete_directory / PIPELINE_DIRECTORY.name, PIPELINE_DIRECTORY)
                finally:
                    rmtree(incomplete_directory, ignore_errors=True)
                    tarball.unlink(missing_ok=True)
            # Download data with a temporary name
            data_directory = DATA_DIRECTORY
            if not DATA_DIRECTORY.exists():
                data_directory = INTERPROSCAN_DIRECTORY / f'.incomplete_{DATA_DIRECTORY.name}'
                # If a previous setup was interrupted then its data are downloaded again from the start
                # InterProScan takes any database directory as complete, so a partly extracted one would be used
                if data_directory.exists():
                    print('Removing the data of an interrupted setup, which are downloaded again from the start...')
                    rmtree(data_directory)
            print(f'Downloading the InterPro {INTERPRO_VERSION} data and the container images, '
                'and running every analysis once. This takes about 40 minutes...')
            with TemporaryDirectory(prefix='mwf_interproscan_setup_') as run_directory:
                fasta_filepath = Path(run_directory) / 'setup.fasta'
                fasta_filepath.write_text(f'>setup\n{SETUP_SEQUENCE}\n')
                # Do not use the Matches API so every analysis is run and every container image is pulled
                # There is no time limit since downloads may be slow
                run_pipeline(fasta_filepath, Path(run_directory), data_directory,
                    use_matches_api=False, check_interpro_version=True, timeout=None)
            if data_directory != DATA_DIRECTORY: os.rename(data_directory, DATA_DIRECTORY)
            # The Apptainer cache is used only to pull images, which are already in the images directory
            rmtree(APPTAINER_CACHE, ignore_errors=True)
            SETUP_MARKER.touch()
        finally:
            rmtree(SETUP_LOCK, ignore_errors=True)
    print(f'InterProScan {INTERPROSCAN_VERSION} with InterPro {INTERPRO_VERSION} is set up at {INTERPROSCAN_DIRECTORY}')


def run_local_interproscan(sequence: str) -> dict:
    """Run InterProScan locally for a protein sequence.
    Return the result in the JSON format returned by the remote EBI service.
    Sequences are run one by one, as the remote service does, since some e-values depend on the number
    of sequences in the run (e.g. PIRSR, searched with hmmsearch).
    Thus the result of a sequence does not depend on which other chains were run with it.
    A sequence takes about 20 seconds when its matches are precalculated, and about 1.5 minutes otherwise, whatever its length.
    Precalculated matches of known sequences are taken from the InterPro Matches API
    only when it serves the same InterPro release as the local data.
    A failed run is tried again, as jobs which fail in the remote service are sent again,
    since in our tests 5 runs in about 110 failed in their first seconds with a Nextflow 'Buffer underflow' error
    (several InterProScan subworkflows split the same FASTA at once and Nextflow shares the split between them),
    and the run tried again worked.
    """
    use_matches_api = prepare_local_interproscan()
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            with TemporaryDirectory(prefix='mwf_interproscan_') as run_directory:
                # Use the same header the remote service gets, since InterProScan reports it as cross-reference
                fasta_filepath = Path(run_directory) / 'chain.fasta'
                fasta_filepath.write_text(f'>chain X\n{sequence}\n')
                output = run_pipeline(fasta_filepath, Path(run_directory), DATA_DIRECTORY, use_matches_api)
            break
        except ToolError as error:
            if attempt == MAX_ATTEMPTS: raise
            warn(f'{error}\nRunning InterProScan again ({attempt + 1}/{MAX_ATTEMPTS})')
    if len(output['results']) != 1:
        raise ToolError(f'Local InterProScan returned {len(output["results"])} results for sequence {sequence}')
    return output


if __name__ == '__main__':
    if '--setup' in sys.argv[1:]:
        setup_local_interproscan()
    else:
        print(f'Local InterProScan {INTERPROSCAN_VERSION} with InterPro {INTERPRO_VERSION}: '
            + ('set up' if local_interproscan_exists() else 'not set up'))
        print(f'InterPro Matches API: InterPro {get_matches_api_release() or "not reachable"}')
        print(f'Set it up in {INTERPROSCAN_DIRECTORY} with: {SETUP_COMMAND}')
