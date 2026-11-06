"""Tests for the check of forced references in the protein mapping."""

import pytest
from mddb_workflow.tools.generate_map import ProteinMapper
from mddb_workflow.utils.auxiliar import InputError


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
