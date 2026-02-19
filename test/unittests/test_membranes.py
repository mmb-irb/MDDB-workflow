from mddb_workflow.mwf import Project, MD, input_files
from mddb_workflow.tools.membrane_mapping import generate_membrane_mapping
from mddb_workflow.utils.auxiliar import save_json
from mddb_workflow.utils.mda_spells import get_mda_universe
from mddb_workflow.utils.cache import Cache
from mddb_workflow.utils.file import File
import MDAnalysis as mda
import pytest


def generate_test_caches(cache_str: str, test_file: str, keep_refbonds: list = False):
    """Generate test cases by filtering the cache data."""
    cache = Cache(File(cache_str))
    keep = ['inchikeys_task_output', 'lipmap_task_output']
    if keep_refbonds:
        keep.append('refbonds_task_output')
    to_delete = [k for k in cache.data.keys() if k not in keep]
    for k in to_delete:
        cache.delete(k)
    save_json(cache.data, test_file, indent=None, separators=(',', ': '))


def download_project(project: 'Project'):
    """Download the project files."""
    project.overwritables = {'itopology'}
    input_files['itopology'](project)
    md: MD = project.mds[0]
    md.overwritables = {'istructure', 'itrajectory'}
    input_files['istructure'](md)
    input_files['itrajectory'](md)


class TestMembraneMapping:
    """Unit tests for the membrane mapping tool.

    - A01IP: base case
    - A01J5: glucolipids
    - A02F9: only lipids
    - OTRMG: no topology/charges
    - cg_test: coarse-grained system
    - cg_test_04: coarse-grained system with changing membranes

    """
    @pytest.fixture(scope='class', params=['A01IP', 'A01J5', 'A02F9', 'OTRMG', 'cg_test'])
    def test_accession(self, request):
        """Fixture to provide different test accessions."""
        return request.param

    def test_generate_membrane_mapping(self, project: 'Project', test_accession, test_proj_dir, test_data_dir):
        """Test the generate_membrane_mapping function."""
        download_project(project)
        cache = Cache(File(f'{test_data_dir}/input/caches/{test_accession}_memmap.json'))
        inchikeys = cache.data['inchikeys_task_output']
        lipid_references = cache.data['lipmap_task_output']
        top = project.input_topology_file.absolute_path
        trajectory_file = project.mds[0].input_trajectory_files[0]
        structure_file = project.mds[0].input_structure_file
        if top.endswith('.top'):
            universe = mda.Universe(top, trajectory_file.absolute_path, topology_format='ITP')
        elif top.endswith('.json'):
            universe = get_mda_universe(structure_file, trajectory_file, cache.data['refbonds_task_output'], charges=None)
        else:
            universe = mda.Universe(all_coordinates=structure_file.absolute_path, topology=top)
        if test_accession == 'cg_test':
            universe.atoms.types = ['Cg'] * len(universe.atoms)
        output_file = File(f'{test_proj_dir}_membrane_mapping_output.json')
        membrane_map = generate_membrane_mapping(inchikeys, lipid_references, structure_file, universe, output_file)
        assert membrane_map['n_mems'] == 1
        top_size = len(membrane_map['mems']['0']['leaflets']['top'])
        bot_size = len(membrane_map['mems']['0']['leaflets']['bot'])
        # Check the leaflets have similar size with 5% tolerance
        assert abs(top_size - bot_size) / max(top_size, bot_size) < 0.20
