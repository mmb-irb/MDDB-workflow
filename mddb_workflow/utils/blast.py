import os
from pathlib import Path
from subprocess import run, PIPE
from shutil import which
import xmltodict
from mddb_workflow.utils.auxiliar import ToolError, retry_request
from mddb_workflow.utils.constants import RESOURCES_DIRECTORY_PATH
from Bio.Blast import NCBIWWW


BLASTP_EXECUTABLE = which('blastp')
UPDATE_BLASTDB_EXECUTABLE = which('update_blastdb.pl')
BLASTDB_DIRECTORY = Path(RESOURCES_DIRECTORY_PATH) / 'blastdb'
BLASTDB_NAME = 'swissprot'
BLASTDB_DIR = str(Path(BLASTDB_DIRECTORY) / BLASTDB_NAME)


def local_blastdb_exists() -> bool:
    """Check if the local Swiss-Prot BLAST database has already been downloaded."""
    # Single-volume databases are named '<name>.pin' while multi-volume ones are named '<name>.00.pin', etc.
    if os.path.exists(f'{BLASTDB_DIR}.pin'):
        return True
    return os.path.exists(f'{BLASTDB_DIR}.00.pin')


def update_local_blastdb():
    """Download (or update, if already downloaded) the pre-formatted Swiss-Prot BLAST database from NCBI.
    This is the same database used by the remote NCBIWWW.qblast counterpart, so results should match
    as long as this local copy is kept reasonably up to date.
    Relies on the 'update_blastdb.pl' script shipped with the BLAST+ suite.
    """
    if not UPDATE_BLASTDB_EXECUTABLE:
        raise ToolError('Cannot find the BLAST+ update_blastdb.pl script. Is BLAST+ installed? '
                        'Add it to the PATH or install it with conda install -c bioconda blast')
    os.makedirs(BLASTDB_DIRECTORY, exist_ok=True)
    print(f'Downloading/updating the local {BLASTDB_NAME} BLAST database at {BLASTDB_DIRECTORY}...')
    process = run([UPDATE_BLASTDB_EXECUTABLE, '--decompress', '--quiet', BLASTDB_NAME],
        cwd=BLASTDB_DIRECTORY, stdout=PIPE, stderr=PIPE)
    if process.returncode != 0:
        raise ToolError(f'Failed to download/update the local {BLASTDB_NAME} BLAST database:\n'
            + process.stderr.decode())


def local_blastp(sequence: str) -> str:
    """Given an amino acids sequence, run a local blastp against the Swiss-Prot database.
    Returns the raw BLAST XML output (outfmt 5), exactly as the remote NCBIWWW.qblast counterpart does,
    so callers may keep using the very same parsing logic regardless of which one was used.
    Downloads the local database the first time it is needed, if it is missing.
    """
    if not BLASTP_EXECUTABLE:
        raise ToolError('Cannot find blastp. Is BLAST+ installed? '
                        'Add it to the PATH or install it with conda install -c bioconda blast')
    if not local_blastdb_exists():
        update_local_blastdb()
    process = run([BLASTP_EXECUTABLE, '-db', BLASTDB_DIR, '-outfmt', '5'],
        input=sequence.encode(), stdout=PIPE, stderr=PIPE)
    if process.returncode != 0:
        raise ToolError('blastp failed:\n' + process.stderr.decode())
    return process.stdout.decode()


@retry_request
def remote_blastp(sequence: str) -> str:
    """Given an amino acids sequence, run a remote blastp against the Swiss-Prot database."""
    result = NCBIWWW.qblast(
        program="blastp",
        database="swissprot",  # UniProtKB / Swiss-Prot
        sequence=sequence,
    )
    return result.read()


def run_blastp(sequence: str, local: bool = False) -> str | None:
    """Given an aminoacids sequence, return a list of uniprot ids.
    Note that we are blasting against UniProtKB / Swiss-Prot so results will always be valid UniProt accessions.
    WARNING: This always means results will correspond to curated entries only.
    If your sequence is from an exotic organism the result may be not from it but from other more studied organism.
    Since this function may take some time we always cache the result.
    """
    print(f'Throwing blast for sequence {sequence}. This may take some time...')
    xml_result = local_blastp(sequence) if local else remote_blastp(sequence)
    parsed_result = xmltodict.parse(xml_result)
    hits = parsed_result['BlastOutput']['BlastOutput_iterations']['Iteration']['Iteration_hits']
    # When there is no result return None
    # Note that this is possible although hardly unprobable
    if not hits:
        return None
    # Get the first result only
    # Note that when there is only one result the Hit isnot an list, but the hit itself
    results = hits['Hit']
    if type(results) is list: first_result = results[0]
    elif type(results) is dict: first_result = results
    else: raise RuntimeError('Invalid hit format')
    # Return the accession
    # DANI: Si algun día tienes problemas porque te falta el '.1' al final del accession puedes sacarlo de Hit_id
    accession = first_result['Hit_accession']
    print('Result: ' + accession)
    return accession
