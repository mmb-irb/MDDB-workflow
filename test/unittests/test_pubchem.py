"""Tests for the PubChem requests."""

import socket
import pytest
from unittest.mock import MagicMock
from rdkit import Chem
import mddb_workflow.utils.pubchem as pubchem
import mddb_workflow.tools.get_ligands as get_ligands
from mddb_workflow.core.dataclass import InChIKeyData
from mddb_workflow.utils.auxiliar import RequestError
from mddb_workflow.utils.cache import Cache
from mddb_workflow.utils.file import File

PHENOL_INCHI = 'InChI=1S/C6H6O/c7-6-4-2-1-3-5-6/h1-5,7H'
PHENOL_CID = '996'


def fake_response(status_code: int, body: dict) -> MagicMock:
    """Fake a PubChem response with the given status and JSON body."""
    response = MagicMock()
    response.status_code = status_code
    response.headers = {}
    response.json.return_value = body
    response.text = str(body)
    return response


TIMEOUT_RESPONSE = fake_response(504, {'Fault': {'Code': 'PUGREST.Timeout', 'Message': 'Request timed out'}})


@pytest.fixture
def fake_post(monkeypatch):
    """Replace the PubChem POST requests and count them."""
    pubchem.pubchem_standardization.cache_clear()
    post = MagicMock()
    monkeypatch.setattr(pubchem, '_make_post_request', post)
    yield post
    pubchem.pubchem_standardization.cache_clear()


@pytest.mark.unit_int
def test_standardization_timeout_gives_no_result(fake_post):
    """Test that a PubChem standardization timeout is not retried and gives no result.

    Some molecules always time out in PubChem after ~28 s (e.g. the ligand in PDB 1NZ7).
    The timeout raised an error which made the whole ligands task fail.
    """
    fake_post.return_value = TIMEOUT_RESPONSE
    assert pubchem.pubchem_standardization(PHENOL_INCHI) is None
    assert fake_post.call_count == 1


@pytest.mark.unit_int
def test_other_gateway_timeouts_still_fail(fake_post):
    """Test that a 504 which is not the PubChem timeout fault still raises an error."""
    fake_post.return_value = fake_response(504, {})
    fake_post.return_value.json.side_effect = ValueError('Not JSON')
    with pytest.raises(RequestError):
        pubchem.pubchem_standardization(PHENOL_INCHI)


@pytest.mark.unit_int
def test_ligand_matching_continues_after_standardization_timeout(fake_post, monkeypatch, tmp_path):
    """Test that the ligands are still matched to PDB-derived ligands when the standardization times out."""
    fake_post.return_value = TIMEOUT_RESPONSE
    # Nothing is found in the database or by InChIKey, so the standardization is reached
    monkeypatch.setattr(get_ligands, 'get_database_inchikey_reference', lambda *args: None)
    monkeypatch.setattr(get_ligands, 'inchikey_2_pubchem', lambda *args: None)
    monkeypatch.setattr(get_ligands, 'pdbs_2_pubchems',
        lambda *args: [{'pubchem': PHENOL_CID, 'inchi': PHENOL_INCHI, 'pdb': True}])
    monkeypatch.setattr(get_ligands, 'obtain_pubchem_data_from_input', lambda *args: {})
    inchikey = Chem.MolToInchiKey(Chem.MolFromInchi(PHENOL_INCHI))
    inchikeys = {inchikey: InChIKeyData(inchi=PHENOL_INCHI, resindices=[0], fragments=[], resnames={'IPH'},
        molname='IPH', moltype='residue', classification=set(), frag_len=1, references={})}
    cache = Cache(File(str(tmp_path / '.mwf_cache.json')))
    old_socket_timeout = socket.getdefaulttimeout()
    try:
        references = get_ligands.generate_ligand_references(None, cache, [], ['1ABC'], inchikeys)
    finally:
        socket.setdefaulttimeout(old_socket_timeout)
    assert fake_post.call_count == 1
    assert references[inchikey]['pubchem'] == PHENOL_CID
