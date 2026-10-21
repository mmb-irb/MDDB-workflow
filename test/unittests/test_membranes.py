from mddb_workflow.mwf import Project, MD, input_files
from mddb_workflow.tools.membrane_mapping import generate_membrane_mapping, find_leaflets, BILAYER, VESICLE, MICELLE
from mddb_workflow.utils.auxiliar import save_json
from mddb_workflow.utils.mda_spells import get_mda_universe
from mddb_workflow.utils.cache import Cache
from mddb_workflow.utils.file import File
import MDAnalysis as mda
import numpy as np
import pathlib
import pytest
import shutil

# Set up paths
data_dir = pathlib.Path(__file__).parent.parent / 'data'
test_dir = data_dir / 'output/test_membranes'


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
    input_files['itopology'](project)
    md: MD = project.mds[0]
    input_files['istructure'](md)
    input_files['itrajectory'](md)


class TestMembraneMapping:
    """Unit tests for the membrane mapping tool.

    - A01IP: base case
    - A01J5: glucolipids
    - A02F9: only lipids
    - OTRMG: no topology/charges, this should fail with the current implementation, leaflets are not similar size when they should
    - cg_test: coarse-grained system
    - cg_test_04: coarse-grained system with changing membranes
    - A01J7: united atoms
    - A05SA: LNP-endosome fusion stalk vesicle, leaflets merge at the default LeafletFinder cutoff in the first frame

    """
    @pytest.fixture(scope='class', params=['A01IP', 'A01J5', 'A02F9', 'OTRMG', 'cg_test', 'A01J7', 'A05SA'])
    def test_accession(self, request):
        """Fixture to provide different test accessions."""
        return request.param

    def test_generate_membrane_mapping(self, project: 'Project', test_accession, test_data_dir):
        """Test the generate_membrane_mapping function."""
        test_dir.mkdir(parents=True, exist_ok=True)
        download_project(project)
        cache_src = f'{test_data_dir}/input/caches/{test_accession}_memmap.json'
        cache_cpy = f'{test_dir}/{test_accession}_cache.json'
        cache_file = File(cache_cpy)
        if cache_file.exists: cache_file.remove()
        shutil.copy(cache_src, cache_cpy)
        cache = Cache(File(cache_cpy))
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
            universe = mda.Universe(top, trajectory_file.absolute_path)
        cg_selection = universe.select_atoms('resname DPPC').residues.resindices if test_accession in ['cg_test', 'A01J7'] else []
        output_file = File((test_dir / f'{test_accession}_membrane_mapping_output.json').as_posix())
        if output_file.exists: output_file.remove()
        membrane_map = generate_membrane_mapping(
            inchikeys,
            lipid_references,
            universe,
            cg_selection,
            output_file
        )
        assert membrane_map['n_mems'] == 1
        membrane = membrane_map['mems']['0']
        assert membrane['type'] == (VESICLE if test_accession == 'A05SA' else BILAYER)
        # Compare number of lipids (one polar atom per lipid) instead of number of atoms
        sizes = [len(polar_atoms) for polar_atoms in membrane['polar_atoms'].values()]
        assert len(sizes) == 2
        # Check the leaflets have similar size, the fusion stalk in A05SA makes its leaflets asymmetric
        tolerance = 0.30 if test_accession == 'A05SA' else 0.20
        assert abs(sizes[0] - sizes[1]) / max(sizes) < tolerance


def fibonacci_sphere(n_points: int, radius: float) -> np.ndarray:
    """Get evenly distributed points on a sphere."""
    indices = np.arange(n_points) + 0.5
    phi = np.arccos(1 - 2 * indices / n_points)
    theta = np.pi * (1 + 5 ** 0.5) * indices
    return radius * np.column_stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)])


def sphere_points(radius: float, spacing: float = 8) -> np.ndarray:
    """Get points on a sphere separated approximately by the given spacing."""
    return fibonacci_sphere(int(4 * np.pi * radius ** 2 / spacing ** 2), radius)


def build_lipids(heads: np.ndarray, tail_directions: np.ndarray, box: list[float]) -> tuple['mda.Universe', 'mda.AtomGroup']:
    """Build a universe of 3 bead lipids (head and two tail beads) and return it with its head atoms."""
    n_lipids = len(heads)
    tail_directions = tail_directions / np.linalg.norm(tail_directions, axis=1)[:, None]
    universe = mda.Universe.empty(n_lipids * 3, n_residues=n_lipids, atom_resindex=np.repeat(np.arange(n_lipids), 3), trajectory=True)
    universe.add_TopologyAttr('names', ['P', 'C1', 'C2'] * n_lipids)
    universe.add_TopologyAttr('resnames', ['LIP'] * n_lipids)
    positions = np.stack([heads, heads + 7 * tail_directions, heads + 14 * tail_directions], axis=1).reshape(-1, 3)
    universe.atoms.positions = positions + np.array(box) / 2
    universe.dimensions = list(box) + [90, 90, 90]
    return universe, universe.select_atoms('name P')


class TestAggregateClassification:
    """Classification of amphiphilic aggregates on synthetic systems."""

    def test_micelle(self):
        heads = sphere_points(20)
        universe, headgroups = build_lipids(heads, -heads, [150, 150, 150])
        membrane_map = find_leaflets(universe, headgroups)
        assert membrane_map['n_mems'] == 1
        assert membrane_map['mems']['0']['type'] == MICELLE
        assert len(membrane_map['mems']['0']['polar_atoms']['outer']) == len(heads)

    def test_vesicle(self):
        outer, inner = sphere_points(75), sphere_points(35)
        heads = np.concatenate([outer, inner])
        universe, headgroups = build_lipids(heads, np.concatenate([-outer, inner]), [250, 250, 250])
        membrane_map = find_leaflets(universe, headgroups)
        assert membrane_map['n_mems'] == 1
        membrane = membrane_map['mems']['0']
        assert membrane['type'] == VESICLE
        assert len(membrane['polar_atoms']['outer']) == len(outer)
        assert len(membrane['polar_atoms']['inner']) == len(inner)

    def test_bilayer(self):
        grid = np.array([[x, y] for x in np.arange(-56, 64, 8) for y in np.arange(-56, 64, 8)], dtype=float)
        bot = np.column_stack([grid, np.full(len(grid), -20)])
        top = np.column_stack([grid, np.full(len(grid), 20)])
        directions = np.concatenate([np.tile([0, 0, 1], (len(grid), 1)), np.tile([0, 0, -1], (len(grid), 1))])
        universe, headgroups = build_lipids(np.concatenate([bot, top]), directions, [120, 120, 100])
        membrane_map = find_leaflets(universe, headgroups)
        assert membrane_map['n_mems'] == 1
        membrane = membrane_map['mems']['0']
        assert membrane['type'] == BILAYER
        # Top leaflet is the one further along the normal
        top_z = universe.atoms[membrane['polar_atoms']['top']].positions[:, 2].mean()
        bot_z = universe.atoms[membrane['polar_atoms']['bot']].positions[:, 2].mean()
        assert top_z > bot_z

    def test_merged_leaflets_retry(self):
        """Leaflets 13 Å apart are merged at the default 15 Å cutoff and must be split with a smaller one."""
        outer, inner = sphere_points(53), sphere_points(40)
        heads = np.concatenate([outer, inner])
        universe, headgroups = build_lipids(heads, np.concatenate([-outer, inner]), [200, 200, 200])
        membrane_map = find_leaflets(universe, headgroups)
        assert membrane_map['n_mems'] == 1
        assert membrane_map['mems']['0']['type'] == VESICLE
