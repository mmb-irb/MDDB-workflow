"""Tests for the protein mapping: the check of forced references and the UniProt requests."""

import io
import json
import urllib.request
from urllib.error import HTTPError
import pytest
from mddb_workflow.tools.generate_map import ProteinMapper, get_uniprot_reference, read_uniprot_response
from mddb_workflow.utils.auxiliar import InputError, RequestError


def make_mapper(forced_references: list | dict, matched_uniprot_ids: list[str]) -> ProteinMapper:
    """Build a mapper with the given forced references and one chain matched to each uniprot id.
    Skip the constructor, which needs a structure and connects to the internet.
    """
    mapper = ProteinMapper.__new__(ProteinMapper)
    mapper.input_protein_references = forced_references
    mapper.strict_references = type(forced_references) is dict
    mapper.protein_parsed_chains = [{'match': {'ref': {'uniprot': uniprot_id}}} for uniprot_id in matched_uniprot_ids]
    return mapper


def test_forced_reference_matched():
    make_mapper(['P12345'], ['P12345'])._check_forced_references_matched()


def test_forced_reference_matched_by_isoform():
    """A chain may align better to an isoform of the forced reference (e.g. IDRome A6/NG/B0/215_276)."""
    make_mapper(['A6NGB0'], ['A6NGB0-1'])._check_forced_references_matched()
    make_mapper({'A': 'A6NGB0'}, ['A6NGB0-1'])._check_forced_references_matched()


def test_forced_reference_not_matched():
    """A different protein, e.g. a blast hit of a paralog, is still not a match of the forced reference."""
    with pytest.raises(InputError):
        make_mapper(['A0A0B4J2F2'], ['P57059'])._check_forced_references_matched()


def fake_uniprot(monkeypatch, responses: list) -> list:
    """Make UniProt answer the given HTTP error codes or JSON data, in order, and skip the waits.
    Return the list of waits, which is filled as the request is tried again.
    """
    waits = []
    monkeypatch.setattr('mddb_workflow.utils.auxiliar.time.sleep', waits.append)
    def urlopen(request):
        response = responses.pop(0)
        if type(response) is int:
            raise HTTPError(request.full_url, response, 'Error', {}, None)
        return io.BytesIO(json.dumps(response).encode())
    monkeypatch.setattr(urllib.request, 'urlopen', urlopen)
    return waits


def test_uniprot_request_retried_on_server_error(monkeypatch):
    """EBI answers 503 for a few minutes now and then (e.g. 6 minutes in IDRome batch 02)."""
    waits = fake_uniprot(monkeypatch, [503, 502, {'accession': 'P12345'}])
    request = urllib.request.Request('https://www.ebi.ac.uk/proteins/api/proteins/P12345')
    assert read_uniprot_response(request) == {'accession': 'P12345'}
    assert waits == [15, 30]


def test_uniprot_request_gives_up(monkeypatch):
    """After about 8 minutes of server errors the mapping fails as before."""
    responses = [503] * 6
    waits = fake_uniprot(monkeypatch, responses)
    with pytest.raises(RequestError):
        get_uniprot_reference('P12345')
    assert waits == [15, 30, 60, 120, 240]
    assert responses == []


def test_uniprot_entry_not_found(monkeypatch):
    """A missing entry is not requested again."""
    waits = fake_uniprot(monkeypatch, [404, {'accession': 'P12345'}])
    assert get_uniprot_reference('P12345') is None
    assert waits == []
