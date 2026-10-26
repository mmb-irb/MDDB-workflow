"""Outdated tasks of projects, found from their versions in the database."""
import pytest
from mddb_workflow.utils.constants import PROJECT_METADATA_VERSION, TOPOLOGY_VERSION, THICKNESS_ANALYSIS_VERSION
from mddb_workflow.tools.membrane_mapping import MEMBRANE_MAPPING_VERSION
from mddb_workflow.tools.update_projects import find_outdated_tasks, get_task_files, MISSING_VERSION

# Versions of an up to date project with a membrane mapping and a thickness analysis, as returned by the database API
UP_TO_DATE = {
    'metadata': PROJECT_METADATA_VERSION,
    'topology': TOPOLOGY_VERSION,
    'analyses': {'mem-map': MEMBRANE_MAPPING_VERSION},
    'mds': [{'analyses': {'thickness': THICKNESS_ANALYSIS_VERSION}}],
}


def test_up_to_date():
    assert find_outdated_tasks(UP_TO_DATE) == {}


@pytest.mark.parametrize('stored, outdated', [('0.2.0', '0.2.0'), (None, '0.0.0')])
def test_outdated_membrane_mapping(stored, outdated):
    """The membrane mapping is a project analysis (with no MD), uploaded from the project directory."""
    versions = {**UP_TO_DATE, 'analyses': {'mem-map': stored}}
    assert find_outdated_tasks(versions) == {'memmap': ([outdated], MEMBRANE_MAPPING_VERSION)}
    assert get_task_files('memmap') == ['mda.mem_map.json']


def test_missing_membrane_mapping():
    """Projects with no membrane have no membrane mapping, so it is only updated when included."""
    versions = {**UP_TO_DATE, 'analyses': {}}
    assert find_outdated_tasks(versions) == {}
    assert find_outdated_tasks(versions, include=['memmap']) == {'memmap': ([MISSING_VERSION], MEMBRANE_MAPPING_VERSION)}


def test_outdated_md_analysis():
    versions = {**UP_TO_DATE, 'mds': [{'analyses': {'thickness': '0.1.0'}}, {'analyses': {'thickness': None}}]}
    assert find_outdated_tasks(versions) == {'thickness': (['0.1.0', '0.0.0'], THICKNESS_ANALYSIS_VERSION)}
