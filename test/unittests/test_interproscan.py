"""Tests for the EBI InterProScan service in the chains task."""
import pytest

import mddb_workflow.tools.chains as chains
from mddb_workflow.utils.file import File


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
