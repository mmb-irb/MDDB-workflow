"""Tests for InterProScan in the chains task: the remote EBI service and the local InterProScan 6 (-lip).

InterProScan itself is not run here: Nextflow is faked with outputs in the InterProScan 6 JSON format.
"""
import os
import io
import json
import time
import socket
import hashlib
import tarfile
import subprocess
from pathlib import Path
from subprocess import CompletedProcess
from urllib.error import URLError

import pytest

import mddb_workflow.utils.interproscan as interproscan
import mddb_workflow.tools.chains as chains
from mddb_workflow.utils.file import File
from mddb_workflow.utils.auxiliar import ToolError, DatabaseUnavailableError, load_json


def build_match(accession: str, entry_accession: str | None) -> dict:
    """Build a match in the InterProScan 6 JSON format, as both the remote service and the local run return it."""
    entry = None
    if entry_accession:
        entry = {'accession': entry_accession, 'name': 'Name', 'description': 'Description', 'type': 'Domain',
            'goXRefs': [], 'pathwayXRefs': [{'name': 'Pathway', 'databaseName': 'Reactome', 'id': 'R-HSA-1'}]}
    return {
        'signature': {'accession': accession, 'name': accession, 'description': None, 'type': 'Family',
            'signatureLibraryRelease': {'library': 'Pfam', 'version': '39.0'}, 'entry': entry},
        'locations': [{'start': 1, 'end': 3, 'representative': True,
            'location-fragments': [{'start': 1, 'end': 3, 'dc-status': 'CONTINUOUS'}]}],
        'model-ac': accession,
        'source': 'Pfam',
    }


def option(command: list[str], name: str) -> str:
    return command[command.index(name) + 1]


def fake_nextflow(command: list[str], **kwargs) -> CompletedProcess:
    """Write what InterProScan 6 writes: the data it downloads and the JSON output for the input FASTA."""
    assert command[0] == 'nextflow'
    # Data are downloaded only when missing
    databases = Path(option(command, '--datadir')) / 'interpro' / interproscan.INTERPRO_VERSION / 'databases.json'
    if not databases.exists():
        databases.parent.mkdir(parents=True)
        databases.write_text('{}')
    lines = Path(option(command, '--input')).read_text().split('\n')
    results = [{
        'sequence': sequence,
        'md5': hashlib.md5(sequence.encode()).hexdigest(),
        'matches': [build_match('PF00001', 'IPR000001'), build_match('PF00002', None)],
        'xref': [{'name': header[1:], 'id': header[1:].split()[0]}],
    } for header, sequence in zip(lines[0::2], lines[1::2])]
    output_filepath = Path(option(command, '--outdir')) / f'{option(command, "--outprefix")}.json'
    output_filepath.write_text(json.dumps({'interproscan-version': interproscan.INTERPROSCAN_VERSION,
        'interpro-version': interproscan.INTERPRO_VERSION, 'results': results}))
    return CompletedProcess(command, 0, b'', b'')


def fake_which(available: list[str]):
    return lambda executable: f'/usr/bin/{executable}' if executable in available else None


def record_runs(monkeypatch, run=fake_nextflow) -> list[tuple[list[str], dict, str]]:
    """Fake Nextflow with 'run' and record the command, the arguments and the input FASTA of every run."""
    runs = []
    def recording_run(command, **kwargs):
        runs.append((command, kwargs, Path(option(command, '--input')).read_text()))
        return run(command, **kwargs)
    monkeypatch.setattr(interproscan, 'run_process', recording_run)
    return runs


@pytest.fixture(autouse=True)
def forget_matches_api_decision():
    """The local InterProScan is checked once per process, so forget it between tests."""
    interproscan.prepare_local_interproscan.cache_clear()


@pytest.fixture
def resources(monkeypatch, tmp_path) -> Path:
    """Use a temporary resources directory, with Nextflow and Apptainer available and Nextflow faked."""
    directory = tmp_path / 'interproscan'
    monkeypatch.setattr(interproscan, 'INTERPROSCAN_DIRECTORY', directory)
    monkeypatch.setattr(interproscan, 'PIPELINE_DIRECTORY', directory / 'interproscan6')
    monkeypatch.setattr(interproscan, 'DATA_DIRECTORY', directory / 'data')
    monkeypatch.setattr(interproscan, 'IMAGES_DIRECTORY', directory / 'images')
    monkeypatch.setattr(interproscan, 'SETUP_MARKER', directory / '.setup')
    monkeypatch.setattr(interproscan, 'SETUP_LOCK', directory / '.setup.lock')
    monkeypatch.setattr(interproscan, 'NEXTFLOW_HOME', directory / 'nextflow')
    monkeypatch.setattr(interproscan, 'APPTAINER_CACHE', directory / '.apptainer')
    monkeypatch.setattr(interproscan, 'which', fake_which(['nextflow', 'apptainer']))
    monkeypatch.setattr(interproscan, 'run_process', fake_nextflow)
    # The Matches API serves the same InterPro release as the local data
    monkeypatch.setattr(interproscan, 'get_matches_api_release', lambda: interproscan.INTERPRO_VERSION)
    return directory


@pytest.fixture
def installed(resources) -> Path:
    """Pretend the local InterProScan is set up."""
    (resources / 'data').mkdir(parents=True)
    interproscan.SETUP_MARKER.touch()
    return resources


@pytest.mark.unit_int
def test_missing_setup(resources, monkeypatch):
    """A clear error tells how to set up InterProScan when it is not, without waiting for the Matches API."""
    def no_api_check():
        raise AssertionError('The Matches API must not be checked')
    monkeypatch.setattr(interproscan, 'get_matches_api_release', no_api_check)
    with pytest.raises(ToolError, match='--setup'):
        interproscan.run_local_interproscan('MKV')


@pytest.mark.unit_int
@pytest.mark.parametrize('available, error', [(['apptainer'], 'nextflow'), (['nextflow'], 'container runtime')])
def test_missing_requirements(installed, monkeypatch, available, error):
    """A clear error tells Nextflow or the container runtime is missing, instead of an InterProScan failure."""
    monkeypatch.setattr(interproscan, 'which', fake_which(available))
    with pytest.raises(ToolError, match=error):
        interproscan.run_local_interproscan('MKV')


@pytest.mark.unit_int
def test_results_have_the_remote_format(installed, monkeypatch):
    """Every sequence is run alone and gets its own result, as the remote service runs and returns it."""
    runs = record_runs(monkeypatch)
    results = {sequence: interproscan.run_local_interproscan(sequence) for sequence in ['MKVLAAGIV', 'GSHMSELLK']}
    # Some e-values depend on the number of sequences in the run, so every sequence is run alone
    # Sequences get the same header as in the remote service
    assert [fasta for _, _, fasta in runs] == ['>chain X\nMKVLAAGIV\n', '>chain X\nGSHMSELLK\n']
    command, kwargs, _ = runs[0]
    # Analyses are run with GO terms and pathways like the remote service, with the pinned InterPro data
    assert {'--goterms', '--pathways', '--skip-interpro-version-check'} <= set(command)
    assert option(command, '-profile') == 'apptainer'
    assert option(command, '--formats') == 'json'
    assert option(command, '--interpro') == interproscan.INTERPRO_VERSION
    assert option(command, '--datadir') == str(interproscan.DATA_DIRECTORY)
    # Precalculated matches are taken from the Matches API, since it serves the same InterPro release
    assert '--no-matches-api' not in command
    assert option(command, '--matches-api-url') == interproscan.MATCHES_API_URL
    # Nextflow works offline, with its files and the images in the workflow resources, and writes its files in the run directory
    assert kwargs['env']['NXF_OFFLINE'] == 'true'
    assert kwargs['env']['NXF_HOME'] == str(interproscan.NEXTFLOW_HOME)
    assert kwargs['env']['NXF_APPTAINER_CACHEDIR'] == str(interproscan.IMAGES_DIRECTORY)
    assert kwargs['env']['APPTAINER_CACHEDIR'] == str(interproscan.APPTAINER_CACHE)
    assert kwargs['timeout'] == interproscan.RUN_TIMEOUT
    # Containers are not run with Docker, so there is no Docker configuration
    assert '-c' not in command
    assert option(command, '--outdir') == str(kwargs['cwd'])
    for sequence, result in results.items():
        assert result['interpro-version'] == interproscan.INTERPRO_VERSION
        assert len(result['results']) == 1
        assert result['results'][0]['sequence'] == sequence
        assert result['results'][0]['xref'] == [{'name': 'chain X', 'id': 'chain'}]


@pytest.mark.unit_int
@pytest.mark.parametrize('release, uses_api, warning', [
    (interproscan.INTERPRO_VERSION, True, None),
    ('999.0', False, f'serves InterPro 999.0 while the local InterProScan data are InterPro {interproscan.INTERPRO_VERSION}'),
    (None, False, 'not reachable'),
])
def test_matches_api_release_check(installed, monkeypatch, capsys, release, uses_api, warning):
    """Precalculated matches are used only when they come from the same InterPro release as the local data."""
    checks = []
    def recording_release():
        checks.append(True)
        return release
    monkeypatch.setattr(interproscan, 'get_matches_api_release', recording_release)
    runs = record_runs(monkeypatch)
    for sequence in ['MKVLAAGIV', 'GSHMSELLK']:
        interproscan.run_local_interproscan(sequence)
    # The release is checked once, so all sequences are run the same way
    assert len(checks) == 1
    assert [('--no-matches-api' not in command) for command, _, _ in runs] == [uses_api, uses_api]
    output = capsys.readouterr().out
    if warning: assert warning in output
    else: assert 'Matches API' not in output


@pytest.mark.unit_int
def test_matches_api_release(monkeypatch):
    """The release is read from the API information, and an unreachable API has no release."""
    urls = []
    def fake_urlopen(url, timeout=None):
        urls.append(url)
        return FakeResponse(json.dumps({'api': '0.6.0', 'release': '111.0', 'release_date': '2026-11-05'}))
    monkeypatch.setattr(interproscan, 'urlopen', fake_urlopen)
    assert interproscan.get_matches_api_release() == '111.0'
    assert urls == [f'{interproscan.MATCHES_API_URL}/info']
    def unreachable(url, timeout=None):
        raise URLError('Name or service not known')
    monkeypatch.setattr(interproscan, 'urlopen', unreachable)
    assert interproscan.get_matches_api_release() is None


@pytest.mark.unit_int
def test_container_runtime_preference(installed, monkeypatch):
    """Apptainer is used first, then Singularity, then Docker."""
    monkeypatch.setattr(interproscan, 'which', fake_which(['nextflow', 'docker', 'singularity']))
    assert interproscan.find_container_runtime() == 'singularity'
    monkeypatch.setattr(interproscan, 'which', fake_which(['nextflow', 'docker']))
    assert interproscan.find_container_runtime() == 'docker'


def failing_nextflow(command: list[str], **kwargs) -> CompletedProcess:
    """Fail as InterProScan did once in our tests, writing the error in the Nextflow log."""
    (Path(kwargs['cwd']) / '.nextflow.log').write_text('INFO Launching\nERROR ~ Buffer underflow.\n')
    return CompletedProcess(command, 1, b'Launching...\n', b'Process PANTHER failed\n')


@pytest.mark.unit_int
def test_failure_shows_the_logs(installed, monkeypatch):
    """A failed run is tried again a few times, and then it is an error which shows the end of the Nextflow logs."""
    runs = record_runs(monkeypatch, failing_nextflow)
    with pytest.raises(ToolError, match='(?s)Process PANTHER failed.*Buffer underflow'):
        interproscan.run_local_interproscan('MKVLAAGIV')
    assert len(runs) == interproscan.MAX_ATTEMPTS == 5


@pytest.mark.unit_int
def test_failed_run_is_tried_again(installed, monkeypatch):
    """A run which fails now and then (a race inside Nextflow) works the second time."""
    def failing_once(command, **kwargs):
        if len(runs) == 1: return failing_nextflow(command, **kwargs)
        return fake_nextflow(command, **kwargs)
    runs = record_runs(monkeypatch, failing_once)
    result = interproscan.run_local_interproscan('MKVLAAGIV')
    assert len(runs) == 2
    assert result['results'][0]['sequence'] == 'MKVLAAGIV'


class FakeStructure:
    def __init__(self, sequences: list[str]):
        self.sequences = sequences

    def get_parsed_chains(self, only_protein: bool = False) -> list[dict]:
        return [{'sequence': sequence} for sequence in self.sequences]


class FakeDatabase:
    def get_reference_data(self, reference: str, sequence: str):
        raise DatabaseUnavailableError('No database')


@pytest.mark.unit_int
@pytest.mark.parametrize('api_release_matches', [True, False])
def test_chains_task_runs_local_interproscan(installed, monkeypatch, tmp_path, api_release_matches):
    """Missing chains are run locally and saved as remote results, and the remote service is not used."""
    def no_remote_request(sequence):
        raise AssertionError('The remote InterProScan must not be requested')
    monkeypatch.setattr(chains, 'request_interproscan', no_remote_request)
    runs = record_runs(monkeypatch)
    # The Matches API release is checked once per task, so all chains are done the same way
    checks = []
    def recording_release():
        checks.append(True)
        return interproscan.INTERPRO_VERSION if api_release_matches else '999.0'
    monkeypatch.setattr(interproscan, 'get_matches_api_release', recording_release)
    # One chain is already in the references file from a previous run
    previous_chain = {'sequence': 'GSHMSELLK', 'interproscan': {'results': []}, 'version': '0.1'}
    (tmp_path / 'chains.json').write_text(json.dumps([previous_chain]))
    chains_file = File(str(tmp_path / 'chains.json'))
    structure = FakeStructure(['MKVLAAGIV', 'GSHMSELLK', 'WWWWWWW', 'MK'])
    chains.prepare_chain_references(structure, chains_file, FakeDatabase(), local_interproscan=True)
    # Sequences shorter than 3 residues are never run, as with the remote service
    assert sorted(fasta for _, _, fasta in runs) == ['>chain X\nMKVLAAGIV\n', '>chain X\nWWWWWWW\n']
    assert [('--no-matches-api' not in command) for command, _, _ in runs] == [api_release_matches] * 2
    assert len(checks) == 1
    chains_data = {chain['sequence']: chain for chain in load_json(chains_file.path)}
    assert chains_data['GSHMSELLK'] == previous_chain
    assert chains_data['MK']['interproscan'] is None
    for sequence in ['MKVLAAGIV', 'WWWWWWW']:
        chain_data = chains_data[sequence]
        assert chain_data['version'] == chains.CHAINS_VERSION == '0.2'
        # Versions and pathways are removed, as for remote results
        assert set(chain_data['interproscan']) == {'results'}
        matches = chain_data['interproscan']['results'][0]['matches']
        assert [match['signature']['accession'] for match in matches] == ['PF00001', 'PF00002']
        assert 'pathwayXRefs' not in matches[0]['signature']['entry']
        assert 'goXRefs' in matches[0]['signature']['entry']


@pytest.mark.unit_int
def test_chains_task_keeps_finished_chains(installed, monkeypatch, tmp_path):
    """Every chain is saved as soon as it is run, so an interrupted task does not run it again."""
    local_runs = []
    def failing_second_run(sequence):
        local_runs.append(sequence)
        if len(local_runs) == 2: raise ToolError('Killed')
        return interproscan.run_local_interproscan(sequence)
    monkeypatch.setattr(chains, 'run_local_interproscan', failing_second_run)
    chains_file = File(str(tmp_path / 'chains.json'))
    with pytest.raises(ToolError):
        chains.prepare_chain_references(FakeStructure(['MKVLAAGIV', 'WWWWWWW']), chains_file, FakeDatabase(),
            local_interproscan=True)
    saved = {chain['sequence']: chain['interproscan'] for chain in load_json(chains_file.path)}
    assert saved[local_runs[0]] is not None
    assert [sequence for sequence, result in saved.items() if result is None] == [local_runs[1]]


class FakeResponse:
    def __init__(self, content: str):
        self.content = content

    def read(self) -> bytes:
        return self.content.encode()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


@pytest.mark.unit_int
def test_remote_requests_use_interproscan_6(monkeypatch):
    """The remote service is the EBI InterProScan 6 one."""
    urls = []
    def fake_urlopen(url, data=None):
        urls.append(url)
        return FakeResponse('iprscan6-R00000000-000000-0000-00000000-p1m' if url.endswith('/run') else 'FINISHED')
    monkeypatch.setattr(chains, 'urlopen', fake_urlopen)
    jobid = chains.request_interproscan('MKVLAAGIV')
    chains.check_interproscan_status(jobid)
    assert urls == ['https://www.ebi.ac.uk/Tools/services/rest/iprscan6/run',
        f'https://www.ebi.ac.uk/Tools/services/rest/iprscan6/status/{jobid}']


def no_download(url: str):
    raise AssertionError(f'Unexpected download of {url}')


def download_fake_pipeline(url: str) -> io.BytesIO:
    """Serve the tarball GitHub serves for an InterProScan release."""
    assert url == interproscan.PIPELINE_URL
    script = b'workflow {}\n'
    script_info = tarfile.TarInfo(f'{interproscan.PIPELINE_DIRECTORY.name}/main.nf')
    script_info.size = len(script)
    tarball = io.BytesIO()
    with tarfile.open(fileobj=tarball, mode='w:gz') as tar:
        tar.addfile(script_info, io.BytesIO(script))
    tarball.seek(0)
    return tarball


@pytest.mark.unit_int
def test_setup(resources, monkeypatch):
    """The setup downloads the pipeline and runs it once, with every analysis, to download data and images."""
    def setup_run(command, **kwargs):
        # The lock is a directory with its owner inside, while the setup runs
        assert (interproscan.SETUP_LOCK / 'owner').read_text().split()[:2] == [socket.gethostname(), str(os.getpid())]
        # Apptainer uses its cache to pull images
        Path(kwargs['env']['APPTAINER_CACHEDIR']).mkdir(parents=True)
        return fake_nextflow(command, **kwargs)
    runs = record_runs(monkeypatch, setup_run)
    monkeypatch.setattr(interproscan, 'urlopen', download_fake_pipeline)
    assert not interproscan.local_interproscan_exists()
    interproscan.setup_local_interproscan()
    assert interproscan.local_interproscan_exists()
    assert (interproscan.PIPELINE_DIRECTORY / 'main.nf').exists()
    assert (interproscan.DATA_DIRECTORY / 'interpro' / interproscan.INTERPRO_VERSION / 'databases.json').exists()
    # Every analysis is run, the pinned InterPro version is checked, and Nextflow may go online
    command, kwargs, _ = runs[0]
    assert '--no-matches-api' in command and '--skip-interpro-version-check' not in command
    assert 'NXF_OFFLINE' not in kwargs['env']
    # Downloads may be slow, so the setup has no time limit
    assert kwargs['timeout'] is None
    # Only the installation is left: no tarball, no lock, no Apptainer cache and no incomplete directories
    assert sorted(path.name for path in resources.iterdir()) == sorted(['.setup', 'interproscan6', 'data'])
    # A second setup does nothing
    monkeypatch.setattr(interproscan, 'urlopen', no_download)
    interproscan.setup_local_interproscan()
    assert len(runs) == 1


@pytest.mark.unit_int
def test_setup_failure_leaves_no_data(resources, monkeypatch):
    """When the setup fails, data stay out of place and a new setup starts the download again."""
    monkeypatch.setattr(interproscan, 'urlopen', download_fake_pipeline)
    def failing_run(command, **kwargs):
        fake_nextflow(command, **kwargs)
        return CompletedProcess(command, 1, b'', b'Download failed\n')
    monkeypatch.setattr(interproscan, 'run_process', failing_run)
    with pytest.raises(ToolError, match='Download failed'):
        interproscan.setup_local_interproscan()
    assert not interproscan.local_interproscan_exists()
    assert not interproscan.DATA_DIRECTORY.exists()
    assert not (resources / '.setup.lock').exists()
    # A new setup removes the incomplete data and works
    monkeypatch.setattr(interproscan, 'run_process', fake_nextflow)
    interproscan.setup_local_interproscan()
    assert interproscan.local_interproscan_exists()
    assert not any(path.name.startswith('.incomplete_') for path in resources.iterdir())


def take_lock(lock_directory: Path, host: str, pid: int, age: float = 0):
    """Take a lock as another process would have taken it 'age' seconds ago."""
    lock_directory.mkdir(parents=True)
    (lock_directory / 'owner').write_text(f'{host} {pid} {time.time() - age}\n')


def finished_process_id() -> int:
    """Get the id of a process which is not running anymore."""
    process = subprocess.Popen(['true'])
    process.wait()
    return process.pid


@pytest.mark.unit_int
@pytest.mark.parametrize('host, running', [('this', True), ('another-host', True), ('this', False)])
def test_setup_lock(resources, monkeypatch, host, running):
    """Two setups at the same time would download the same files, so the second one stops.
    A lock left by a setup which was killed in this host is removed, but a lock from another host can not be checked.
    """
    monkeypatch.setattr(interproscan, 'urlopen', download_fake_pipeline)
    take_lock(interproscan.SETUP_LOCK, socket.gethostname() if host == 'this' else host,
        os.getpid() if running else finished_process_id())
    if running:
        with pytest.raises(ToolError, match='Another InterProScan setup is running in host'):
            interproscan.setup_local_interproscan()
        assert interproscan.SETUP_LOCK.is_dir()
    else:
        interproscan.setup_local_interproscan()
        assert interproscan.local_interproscan_exists()
        assert not interproscan.SETUP_LOCK.exists()


@pytest.mark.unit_int
def test_setup_finished_while_waiting(resources, monkeypatch):
    """A setup which finished while this one was taking the lock is not done again."""
    take_directory_lock = interproscan.take_directory_lock
    def lock_after_another_setup(lock_directory, max_age=None):
        interproscan.SETUP_MARKER.touch()
        return take_directory_lock(lock_directory, max_age)
    monkeypatch.setattr(interproscan, 'take_directory_lock', lock_after_another_setup)
    monkeypatch.setattr(interproscan, 'urlopen', no_download)
    interproscan.setup_local_interproscan()
    assert not interproscan.SETUP_LOCK.exists()


@pytest.mark.unit_int
def test_run_process_timeout(tmp_path):
    """A run which takes too long is stopped with its whole process tree, and it is an error."""
    command = ['bash', '-c', 'sleep 60 & echo $! > child.pid; echo started; wait']
    with pytest.raises(ToolError, match='did not finish in 0 minutes'):
        interproscan.run_process(command, cwd=tmp_path, env=dict(os.environ), timeout=3)
    child = int((tmp_path / 'child.pid').read_text())
    # The child may take a moment to end after the signal
    # Note that an ended process may still be listed for a moment, as zombie (Z) or while it is removed (X)
    # Reading the status of a process which is being removed fails with ProcessLookupError instead
    def child_ended() -> bool:
        try:
            return open(f'/proc/{child}/stat').read().split()[2] in ('Z', 'X')
        except (FileNotFoundError, ProcessLookupError):
            return True
    for _ in range(50):
        if child_ended(): break
        time.sleep(0.1)
    else:
        pytest.fail(f'Process {child} is still running')


@pytest.mark.unit_int
@pytest.mark.parametrize('slurm_cpus', [None, '1'])
def test_runs_use_the_cpus_of_the_job(installed, monkeypatch, slurm_cpus):
    """Runs use the CPUs of the process, or only those requested by its SLURM job, which may not be pinned to them."""
    monkeypatch.delenv('SLURM_CPUS_ON_NODE', raising=False)
    if slurm_cpus: monkeypatch.setenv('SLURM_CPUS_PER_TASK', slurm_cpus)
    else: monkeypatch.delenv('SLURM_CPUS_PER_TASK', raising=False)
    runs = record_runs(monkeypatch)
    interproscan.run_local_interproscan('MKVLAAGIV')
    command, _, _ = runs[0]
    assert option(command, '--max-workers') == (slurm_cpus or str(len(os.sched_getaffinity(0))))


@pytest.mark.unit_int
def test_docker_runs_as_the_user(installed, monkeypatch):
    """Docker containers run as the user, so the run directory can be removed."""
    monkeypatch.setattr(interproscan, 'which', fake_which(['nextflow', 'docker']))
    configs = []
    def recording_run(command, **kwargs):
        configs.append(Path(option(command, '-c')).read_text())
        return fake_nextflow(command, **kwargs)
    monkeypatch.setattr(interproscan, 'run_process', recording_run)
    interproscan.run_local_interproscan('MKVLAAGIV')
    assert configs == [f"docker.runOptions = '-u {os.getuid()}:{os.getgid()}'\n"]


@pytest.mark.unit_int
def test_remote_timeout_stops_the_task(monkeypatch, tmp_path):
    """Remote jobs which do not finish in time stop the task, so it is not taken as done with missing chains."""
    monkeypatch.setattr(chains, 'INTERPROSCAN_TIMEOUT', -1)
    monkeypatch.setattr(chains, 'check_interproscan_status', lambda jobid: 'QUEUED')
    chains_data = [{'sequence': 'MKVLAAGIV', 'interproscan': None}]
    with pytest.raises(chains.RemoteServiceError, match='did not finish'):
        chains.get_interproscan_results(['job1'], {'MKVLAAGIV': 'job1'}, chains_data, File(str(tmp_path / 'chains.json')))


@pytest.mark.unit_int
def test_cleanup_with_no_entry_field():
    """Matches with no InterPro entry field are kept as they are."""
    match = build_match('PF00001', None)
    del match['signature']['entry']
    result = chains.clean_interproscan_result({'interproscan-version': '6.0.3.0', 'interpro-version': '111.0',
        'results': [{'matches': [match, build_match('PF00002', 'IPR000002')]}]})
    assert set(result) == {'results'}
    assert 'pathwayXRefs' not in result['results'][0]['matches'][1]['signature']['entry']
