from mddb_workflow.mwf import Project, MD, input_files
from mddb_workflow.tools.membrane_mapping import generate_membrane_mapping, find_leaflets, BILAYER, VESICLE, MICELLE, MIDPLANE
from mddb_workflow.utils.auxiliar import save_json
from mddb_workflow.utils.mda_spells import get_mda_universe, get_mda_universe_from_stopology
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


# Expected aggregate type for the test accessions which are not flat bilayers
EXPECTED_TYPES = {'A05SA': VESICLE, 'A05RV': MICELLE, 'A05Q3': VESICLE}
# Expected number of unassigned lipids for the test accessions where it is known
EXPECTED_UNASSIGNED = {'A01JD': 1}


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
    - A05SA: LNP-endosome fusion stalk vesicle, leaflets merge at the default LeafletFinder cutoff in the first frame, MC3 and MC3H forced as lipids
    - A05RV: empty lipid nanoparticle micelle, no topology/charges, MC3H and DPG2 forced as lipids
    - A05Q3: LNP-endosome fusion stalk vesicle, no topology/charges, MC3T forced as lipid
    - A01JD: united atoms channel with a lipid inside the pore which must remain unassigned

    """
    @pytest.fixture(scope='class', params=['A01IP', 'A01J5', 'A02F9', 'OTRMG', 'cg_test', 'A05SA', 'A05RV', 'A05Q3', 'A01JD'])
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
        elif top.endswith('.json') and 'refbonds_task_output' in cache.data:
            universe = get_mda_universe(structure_file, trajectory_file, cache.data['refbonds_task_output'], charges=None)
        elif top.endswith('.json'):
            # Bonds are read from the standard topology, too big to be cached for large systems
            universe = get_mda_universe_from_stopology(top, trajectory_file.absolute_path)
        else:
            universe = mda.Universe(top, trajectory_file.absolute_path)
        cg_selection = universe.select_atoms('resname DPPC').residues.resindices if test_accession in ['cg_test', 'A01JD'] else []
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
        assert membrane['type'] == EXPECTED_TYPES.get(test_accession, BILAYER)
        if test_accession in EXPECTED_UNASSIGNED:
            unassigned_lipids = len(set(universe.atoms[membrane_map['no_mem_lipid']].resindices))
            assert unassigned_lipids == EXPECTED_UNASSIGNED[test_accession]
        # Compare number of lipids (one polar atom per lipid) instead of number of atoms
        if membrane['type'] == MICELLE:
            assert list(membrane['polar_atoms']) == ['outer']
            return
        leaflet_names = ['outer', 'inner'] if membrane['type'] == VESICLE else ['top', 'bot']
        sizes = [len(membrane['polar_atoms'][name]) for name in leaflet_names]
        # Check the leaflets have similar size, the fusion stalk vesicles have asymmetric leaflets
        tolerance = 0.30 if membrane['type'] == VESICLE else 0.20
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

    def test_core_lipids_join_micelle(self):
        """Lipids in the core of a micelle, too far from the polar atoms of the shell, join it by contact."""
        shell = sphere_points(20, spacing=6)
        core = np.array([[0, 0, 0], [0, 0, 4], [0, 4, 0]], dtype=float)
        heads = np.concatenate([shell, core])
        directions = np.concatenate([-shell, np.tile([1, 0, 0], (len(core), 1))])
        universe, headgroups = build_lipids(heads, directions, [150, 150, 150])
        membrane_map = find_leaflets(universe, headgroups)
        assert membrane_map['n_mems'] == 1
        assert len(membrane_map['mems']['0']['polar_atoms']['outer']) == len(heads)
        assert membrane_map['no_mem_lipid'] == []

    def test_midplane_lipid_stays_unassigned(self):
        """A lipid in the bilayer midplane touching both leaflets alike is not assigned to any of them."""
        grid = np.array([[x, y] for x in np.arange(-56, 64, 8) for y in np.arange(-56, 64, 8)], dtype=float)
        bot = np.column_stack([grid, np.full(len(grid), -19)])
        top = np.column_stack([grid, np.full(len(grid), 19)])
        midplane = np.array([[4, 4, 0]], dtype=float)
        directions = np.concatenate([np.tile([0, 0, 1], (len(grid), 1)), np.tile([0, 0, -1], (len(grid), 1)), [[1, 0, 0]]])
        universe, headgroups = build_lipids(np.concatenate([bot, top, midplane]), directions, [120, 120, 100])
        membrane_map = find_leaflets(universe, headgroups)
        assert membrane_map['n_mems'] == 1
        assert membrane_map['mems']['0']['type'] == BILAYER
        assert len(membrane_map['no_mem_lipid']) == 3  # The 3 beads of the midplane lipid

    def test_vesicle_midplane(self):
        """A lipid between both leaflets of a vesicle touching them alike is added to its midplane leaflet."""
        # Tail ends are at 61 Å (outer) and 49 Å (inner) from the center, set the same tail end spacing in both
        outer, inner = sphere_points(75), sphere_points(35, spacing=8 * 61 / 75 * 35 / 49)
        midplane = np.array([[54.5, 0, 0]], dtype=float)
        heads = np.concatenate([outer, inner, midplane])
        directions = np.concatenate([-outer, inner, [[0, 1, 0]]])
        universe, headgroups = build_lipids(heads, directions, [250, 250, 250])
        membrane_map = find_leaflets(universe, headgroups)
        membrane = membrane_map['mems']['0']
        assert membrane['type'] == VESICLE
        assert len(membrane['polar_atoms'][MIDPLANE]) == 1
        assert membrane_map['no_mem_lipid'] == []

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


class TestLipidReferences:
    """Lipid references of residues classified as fatty which are not in lipid databases."""

    @pytest.fixture(autouse=True)
    def no_databases(self, monkeypatch):
        """Make every residue missing in the MDDB, SwissLipids and LIPID MAPS databases."""
        import mddb_workflow.tools.get_lipids as get_lipids
        monkeypatch.setattr(get_lipids, 'get_database_inchikey_reference', lambda *args: None)
        monkeypatch.setattr(get_lipids, 'is_in_swisslipids', lambda *args, **kwargs: None)
        monkeypatch.setattr(get_lipids, 'is_in_LIPID_MAPS', lambda *args, **kwargs: None)

    @staticmethod
    def fatty_inchikeys():
        from mddb_workflow.core.dataclass import InChIKeyData
        return {'NRLNQCOGCKAESA-KWXKLSQISA-O': InChIKeyData(
            inchi='InChI=1S/test', resindices=[0, 1], resnames={'MC3H'}, molname='MC3H', classification={'fatty'})}

    def test_fatty_not_lipid_fails(self):
        from mddb_workflow.tools.get_lipids import generate_lipid_references
        from mddb_workflow.utils.auxiliar import InputError
        with pytest.raises(InputError, match='forced_lipids'):
            generate_lipid_references(self.fatty_inchikeys(), cache=None)

    def test_fatty_not_lipid_mercy(self):
        from mddb_workflow.tools.get_lipids import generate_lipid_references
        assert generate_lipid_references(self.fatty_inchikeys(), cache=None, mercy=['lipids']) == {}

    def test_forced_lipids(self):
        from mddb_workflow.tools.get_lipids import generate_lipid_references
        references = generate_lipid_references(self.fatty_inchikeys(), cache=None, input_forced_lipids=['MC3H'])
        assert references['NRLNQCOGCKAESA-KWXKLSQISA-O']['forced'] is True
