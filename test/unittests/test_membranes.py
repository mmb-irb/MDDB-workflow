from mddb_workflow.mwf import Project, MD, input_files
from mddb_workflow.core.dataclass import InChIKeyData
from mddb_workflow.tools.membrane_mapping import generate_membrane_mapping, find_leaflets, BILAYER, VESICLE, MICELLE, MIDPLANE
from mddb_workflow.utils.auxiliar import save_json, load_json
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


def shift_and_wrap(universe: 'mda.Universe', shift: list[float]):
    """Translate all atoms and put the residues back in the box, so lipids crossing the box boundary are whole."""
    universe.atoms.translate(shift)
    universe.atoms.wrap(compound='residues', center='cog')


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

    @pytest.mark.parametrize('shift', [0, 30])
    def test_bilayer_in_short_box(self, shift):
        """Leaflets 40 Å apart in a 70 Å high box are 30 Å apart back to back through the box boundary.
        The leaflets facing each other must be paired, also when the bilayer crosses the box boundary.
        """
        grid = np.array([[x, y] for x in np.arange(-56, 64, 8) for y in np.arange(-56, 64, 8)], dtype=float)
        bot = np.column_stack([grid, np.full(len(grid), -20)])
        top = np.column_stack([grid, np.full(len(grid), 20)])
        directions = np.concatenate([np.tile([0, 0, 1], (len(grid), 1)), np.tile([0, 0, -1], (len(grid), 1))])
        universe, headgroups = build_lipids(np.concatenate([bot, top]), directions, [120, 120, 70])
        shift_and_wrap(universe, [0, 0, shift])
        membrane_map = find_leaflets(universe, headgroups)
        assert membrane_map['n_mems'] == 1
        membrane = membrane_map['mems']['0']
        assert membrane['type'] == BILAYER
        assert set(membrane['polar_atoms']['top']) == set(headgroups[len(grid):].indices)

    def test_vesicle_across_box_boundary(self):
        """A vesicle split by the box boundary is a single vesicle."""
        outer, inner = sphere_points(75), sphere_points(35)
        universe, headgroups = build_lipids(np.concatenate([outer, inner]), np.concatenate([-outer, inner]), [250, 250, 250])
        shift_and_wrap(universe, [125, 0, 0])
        membrane_map = find_leaflets(universe, headgroups)
        assert membrane_map['n_mems'] == 1
        membrane = membrane_map['mems']['0']
        assert membrane['type'] == VESICLE
        assert len(membrane['polar_atoms']['outer']) == len(outer)
        assert len(membrane['polar_atoms']['inner']) == len(inner)

    def test_vesicle_with_lumen_micelle(self):
        """The inner leaflet of a vesicle merges at 15 Å with a micelle in its lumen and they are split at 12 Å.
        The outer leaflet must wait for the inner leaflet instead of becoming a micelle.
        """
        outer, inner, lumen = sphere_points(75), sphere_points(35), sphere_points(22)
        heads = np.concatenate([outer, inner, lumen])
        universe, headgroups = build_lipids(heads, np.concatenate([-outer, inner, -lumen]), [250, 250, 250])
        membrane_map = find_leaflets(universe, headgroups)
        assert [mem['type'] for mem in membrane_map['mems'].values()] == [VESICLE, MICELLE]
        vesicle, micelle = membrane_map['mems']['0'], membrane_map['mems']['1']
        assert len(vesicle['polar_atoms']['outer']) == len(outer)
        assert len(vesicle['polar_atoms']['inner']) == len(inner)
        assert len(micelle['polar_atoms']['outer']) == len(lumen)
        assert membrane_map['no_mem_lipid'] == []

    def test_rough_bilayer_patch_is_not_micelle(self):
        """Leaflets of a small and rough bilayer patch are not flat, but their lipids are not oriented to the center either.
        They must not be micelles.
        """
        grid = np.array([[x, y] for x in np.arange(-20, 28, 8) for y in np.arange(-20, 28, 8)], dtype=float)
        roughness = np.where((grid[:, 0] + grid[:, 1]) % 16 == 0, 4.5, -4.5)
        bot = np.column_stack([grid, roughness - 20])
        top = np.column_stack([grid, roughness + 20])
        directions = np.concatenate([np.tile([0, 0, 1], (len(grid), 1)), np.tile([0, 0, -1], (len(grid), 1))])
        universe, headgroups = build_lipids(np.concatenate([bot, top]), directions, [150, 150, 150])
        membrane_map = find_leaflets(universe, headgroups)
        assert all(mem['type'] != MICELLE for mem in membrane_map['mems'].values())

# Molecules to test the lipid classification, with their expected LIPID MAPS category
LIPIDS = {
    'POPC': (r'CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCC/C=C\CCCCCCCC', 'GP'),
    'DSPC': ('CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCCCCCC', 'GP'),
    'cholesterol': ('CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C', 'ST'),
    'palmitoyl sphingomyelin': ('CCCCCCCCCCCCC/C=C/[C@@H](O)[C@H](COP(=O)([O-])OCC[N+](C)(C)C)NC(=O)CCCCCCCCCCCCCCC', 'SP'),
    'triglyceride': (r'CCCCCCCCCCCCCCCC(=O)OCC(OC(=O)CCCCCCC/C=C\CCCCCCCC)COC(=O)CCCCCCCCCCCCCCC', 'GL'),
    'MC3': (r'CCCCC/C=C\C/C=C\CCCCCCCCC(CCCCCCCC/C=C\C/C=C\CCCCC)OC(=O)CCCN(C)C', 'FA'),
    'MC3H': (r'CCCCC/C=C\C/C=C\CCCCCCCCC(CCCCCCCC/C=C\C/C=C\CCCCC)OC(=O)CCC[NH+](C)C', 'FA'),
}
NON_LIPIDS = {
    'ibuprofen': 'CC(C)Cc1ccc(C(C)C(=O)O)cc1',
    'alanine': 'C[C@H](N)C(=O)O',
    'ATP': 'Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O',
    'glucose': 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O',
    # Steroid hormones and drugs have no sterol side chain and they do not form membranes
    'testosterone': 'C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@@]43C)[C@@H]1CC[C@@H]2O',
    # Too short chain for a fatty acyl
    'octanoic acid': 'CCCCCCCC(=O)O',
}


def smiles_to_inchi(smiles: str) -> str:
    """Get the InChI of a SMILES."""
    from rdkit import Chem
    return Chem.MolToInchi(Chem.MolFromSmiles(smiles))


class TestLipidClassification:
    """Lipid recognition from the structure with SMARTS of the LIPID MAPS categories."""

    @pytest.mark.parametrize('name', LIPIDS)
    def test_lipids(self, name):
        from mddb_workflow.tools.lipid_classification import classify_lipid
        smiles, category = LIPIDS[name]
        assert classify_lipid(smiles_to_inchi(smiles)) == category

    @pytest.mark.parametrize('name', NON_LIPIDS)
    def test_non_lipids(self, name):
        from mddb_workflow.tools.lipid_classification import classify_lipid
        assert classify_lipid(smiles_to_inchi(NON_LIPIDS[name])) is None

    def test_invalid_inchi(self):
        from mddb_workflow.tools.lipid_classification import classify_lipid
        assert classify_lipid('InChI=1S/test') is None

    @pytest.mark.parametrize('test_accession', ['A01IP', 'A01J5', 'A02F9', 'OTRMG', 'A05SA', 'A05RV', 'A05Q3'])
    def test_cached_lipid_references(self, test_accession, test_data_dir):
        """The lipids of the membrane mapping test caches, found in lipid databases or forced, are recognised."""
        from mddb_workflow.tools.lipid_classification import classify_lipid
        # Read the cache file directly since the Cache rewrites it
        cache_data = load_json(f'{test_data_dir}/input/caches/{test_accession}_memmap.json')
        lipid_references = cache_data['lipmap_task_output']
        for inchikey, data in cache_data['inchikeys_task_output'].items():
            is_lipid = classify_lipid(data['inchi']) is not None
            assert is_lipid == (inchikey in lipid_references), f'{inchikey} {data["resnames"]}'


def inchikey_data(smiles: str, resname: str, classification: str, n_molecules: int = 1) -> tuple[str, InChIKeyData]:
    """Get the InChI key and InChIKeyData of a molecule repeated in consecutive residues."""
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToInchiKey(mol), InChIKeyData(inchi=Chem.MolToInchi(mol), resindices=list(range(n_molecules)),
        resnames={resname}, molname=resname, classification={classification})


class TestLipidReferences:
    """Lipid references of residues recognised as lipids, with the lipid databases monkeypatched."""

    @pytest.fixture(autouse=True)
    def databases(self, monkeypatch):
        """Make every residue missing in the MDDB database and record the SwissLipids and LIPID MAPS requests.
        Only POPC is found in SwissLipids. SwissLipids fails while swisslipids_down is set.
        """
        import mddb_workflow.tools.get_lipids as get_lipids
        from mddb_workflow.utils.auxiliar import RequestError
        self.popc_key, _ = inchikey_data(LIPIDS['POPC'][0], 'POPC', 'fatty')
        self.requests = set()
        self.swisslipids_down = False

        def is_in_swisslipids(inchikey, only_first_layer=False):
            if self.swisslipids_down: raise RequestError('SwissLipids is down')
            self.requests.add(inchikey)
            return {'entity_id': 'SLM:000000001'} if inchikey == self.popc_key else None

        def is_in_LIPID_MAPS(inchikey, only_first_layer=False):
            self.requests.add(inchikey)
            return False

        monkeypatch.setattr(get_lipids, 'get_database_inchikey_reference', lambda *args: None)
        monkeypatch.setattr(get_lipids, 'is_in_swisslipids', is_in_swisslipids)
        monkeypatch.setattr(get_lipids, 'is_in_LIPID_MAPS', is_in_LIPID_MAPS)

    def test_lipid_references(self):
        """Lipids are recognised from the structure and the databases are only requested for them."""
        from mddb_workflow.tools.get_lipids import generate_lipid_references
        inchikeys = dict([
            inchikey_data(LIPIDS['POPC'][0], 'POPC', 'fatty', n_molecules=100),
            inchikey_data(LIPIDS['MC3H'][0], 'MC3H', 'fatty', n_molecules=100),
            inchikey_data(NON_LIPIDS['ibuprofen'], 'IBU', 'other'),
        ])
        popc_key, mc3h_key, ibuprofen_key = inchikeys
        references = generate_lipid_references(inchikeys, cache=None)
        assert set(references) == {popc_key, mc3h_key}
        assert references[popc_key]['category'] == 'GP'
        assert references[popc_key]['swisslipids'] == {'entity_id': 'SLM:000000001'}
        # Synthetic lipids are lipids without database references
        assert references[mc3h_key]['category'] == 'FA'
        assert references[mc3h_key]['swisslipids'] is None and references[mc3h_key]['lipidmaps'] is None
        # The connection test is not a lipid request
        assert self.requests - {'test'} == {popc_key, mc3h_key}

    def test_fatty_not_lipid_warns(self):
        """A few fatty residues which are not lipids (e.g. ligands with alkyl chains) are not lipids."""
        from mddb_workflow.tools.get_lipids import generate_lipid_references
        inchikeys = dict([inchikey_data(NON_LIPIDS['octanoic acid'], 'OCA', 'fatty', n_molecules=2)])
        assert generate_lipid_references(inchikeys, cache=None) == {}
        assert self.requests == set()

    def test_abundant_fatty_not_lipid_fails(self):
        """Fatty residues which are not lipids, but enough to be part of a membrane, make the workflow fail."""
        from mddb_workflow.tools.get_lipids import generate_lipid_references
        from mddb_workflow.utils.auxiliar import InputError
        inchikeys = dict([inchikey_data(NON_LIPIDS['octanoic acid'], 'OCA', 'fatty', n_molecules=30)])
        with pytest.raises(InputError, match='forced_lipids'):
            generate_lipid_references(inchikeys, cache=None)

    def test_abundant_fatty_not_lipid_mercy(self):
        from mddb_workflow.tools.get_lipids import generate_lipid_references
        inchikeys = dict([inchikey_data(NON_LIPIDS['octanoic acid'], 'OCA', 'fatty', n_molecules=30)])
        assert generate_lipid_references(inchikeys, cache=None, mercy=['lipids']) == {}

    def test_forced_lipids(self):
        from mddb_workflow.tools.get_lipids import generate_lipid_references
        inchikeys = dict([inchikey_data(NON_LIPIDS['octanoic acid'], 'OCA', 'fatty', n_molecules=30)])
        oca_key = next(iter(inchikeys))
        references = generate_lipid_references(inchikeys, cache=None, input_forced_lipids=['OCA'])
        assert references[oca_key]['forced'] is True
        assert references[oca_key]['category'] is None
        assert self.requests == set()

    def test_inchikey_references_lipid_data(self, tmp_path):
        """InChI key references tell lipids, with their category, from other molecules, also without database references."""
        from mddb_workflow.tools.get_lipids import generate_lipid_references
        from mddb_workflow.tools.get_inchi_keys import generate_inchi_references
        inchikeys = dict([
            inchikey_data(LIPIDS['POPC'][0], 'POPC', 'fatty', n_molecules=100),
            inchikey_data(LIPIDS['MC3H'][0], 'MC3H', 'fatty', n_molecules=100),
            inchikey_data(NON_LIPIDS['octanoic acid'], 'OCA', 'fatty', n_molecules=30),
            inchikey_data(NON_LIPIDS['ibuprofen'], 'IBU', 'other'),
        ])
        popc_key, mc3h_key, oca_key, ibuprofen_key = inchikeys
        lipid_references = generate_lipid_references(inchikeys, cache=None, input_forced_lipids=['OCA'])
        output_file = File(str(tmp_path / 'inchikey_references.json'))
        generate_inchi_references(inchikeys, lipid_references, {}, {'no_mem_lipid': []}, None, output_file)
        references = {reference['inchikey']: reference for reference in load_json(output_file.path)}
        assert references[popc_key]['lipid'] == {'category': 'GP', 'forced': False}
        assert references[mc3h_key]['lipid'] == {'category': 'FA', 'forced': False}
        assert references[mc3h_key]['swisslipids'] is None and references[mc3h_key]['lipidmaps'] is None
        assert references[oca_key]['lipid'] == {'category': None, 'forced': True}
        assert references[ibuprofen_key]['lipid'] is None

    def test_remade_inchi_reference_lipid_data(self, monkeypatch):
        """InChI references remade by mwf upref also have the lipid data."""
        import mddb_workflow.tools.update_references as update_references
        monkeypatch.setattr(update_references, 'generate_ligand_references', lambda *args, **kwargs: {})
        reference = update_references.get_inchi_reference(smiles_to_inchi(LIPIDS['MC3H'][0]), cache=None)
        assert reference['lipid'] == {'category': 'FA', 'forced': False}
        reference = update_references.get_inchi_reference(smiles_to_inchi(NON_LIPIDS['ibuprofen']), cache=None)
        assert reference['lipid'] is None

    def test_references_without_lipid_databases(self):
        """Lipid references made while SwissLipids is down are outdated, so they are not reused from the database.
        SwissLipids is checked again for the next lipids, since it may be a transient failure.
        """
        from mddb_workflow.tools.get_lipids import generate_lipid_references
        from mddb_workflow.tools.get_inchi_keys import make_inchikey_reference
        from mddb_workflow.utils.constants import LIPID_REFERENCE_VERSION, INCHI_REFERENCE_VERSION
        from packaging.version import Version
        inchikeys = dict([inchikey_data(LIPIDS['POPC'][0], 'POPC', 'fatty', n_molecules=100)])
        self.swisslipids_down = True
        lipid_reference = generate_lipid_references(inchikeys, cache=None)[self.popc_key]
        assert lipid_reference['swisslipids'] is None
        assert Version(lipid_reference['version']) < Version(LIPID_REFERENCE_VERSION)
        reference = make_inchikey_reference(self.popc_key, inchikeys[self.popc_key].inchi, lipid_reference, None)
        assert Version(reference['version']) < Version(INCHI_REFERENCE_VERSION)
        # SwissLipids answers again
        self.swisslipids_down = False
        lipid_reference = generate_lipid_references(inchikeys, cache=None)[self.popc_key]
        assert lipid_reference['swisslipids'] == {'entity_id': 'SLM:000000001'}
        assert lipid_reference['version'] == LIPID_REFERENCE_VERSION
        reference = make_inchikey_reference(self.popc_key, inchikeys[self.popc_key].inchi, lipid_reference, None)
        assert reference['version'] == INCHI_REFERENCE_VERSION

    def test_remade_inchi_reference_without_lipid_databases(self, monkeypatch):
        """mwf upref fails to remake lipid references while SwissLipids is down and it succeeds once it answers again."""
        import mddb_workflow.tools.update_references as update_references
        from mddb_workflow.utils.auxiliar import RemoteServiceError
        monkeypatch.setattr(update_references, 'generate_ligand_references', lambda *args, **kwargs: {})
        popc_inchi = smiles_to_inchi(LIPIDS['POPC'][0])
        self.swisslipids_down = True
        with pytest.raises(RemoteServiceError):
            update_references.get_inchi_reference(popc_inchi, cache=None)
        self.swisslipids_down = False
        reference = update_references.get_inchi_reference(popc_inchi, cache=None)
        assert reference['swisslipids'] == {'entity_id': 'SLM:000000001'}

    def test_membrane_lipids_are_not_ligands(self, tmp_path):
        """Lipids with ligand references are ligands only when most of them are out of membranes."""
        from mddb_workflow.tools.get_inchi_keys import generate_inchi_references
        from mddb_workflow.utils.selections import Selection

        class Structure:
            """Structure with 3 atoms per residue."""
            def select_residue_indices(self, residue_indices):
                return Selection([3 * index + atom for index in residue_indices for atom in range(3)])

        # A few MC3 lipids in a membrane and many POPC lipids out of membranes
        mc3_key, mc3_data = inchikey_data(LIPIDS['MC3'][0], 'MC3', 'fatty', n_molecules=2)
        popc_key, popc_data = inchikey_data(LIPIDS['POPC'][0], 'POPC', 'fatty', n_molecules=10)
        popc_data.resindices = list(range(2, 12))
        inchikeys = {mc3_key: mc3_data, popc_key: popc_data}
        lipid_references = {key: {'category': 'GP', 'version': '0.1.0'} for key in inchikeys}
        ligand_references = {key: {'formula': 'C42H82NO8P'} for key in inchikeys}
        membrane_map = {'no_mem_lipid': Structure().select_residue_indices(popc_data.resindices).atom_indices}
        inchikey_map = generate_inchi_references(inchikeys, lipid_references, ligand_references, membrane_map,
            Structure(), File(str(tmp_path / 'inchikey_references.json')))
        is_ligand = {data['inchikey']: data['is_ligand'] for data in inchikey_map}
        assert is_ligand == {mc3_key: False, popc_key: True}


class TestUpdateReferences:
    """Remaking database references with mwf upref, with the database and the reference maker monkeypatched."""

    @pytest.fixture(autouse=True)
    def database(self, monkeypatch, tmp_path):
        """Set a database with two InChI references at the latest version and record the remade references."""
        import mddb_workflow.tools.update_references as update_references
        from mddb_workflow.utils.constants import INCHI_REFERENCE_VERSION
        self.update_references = update_references
        self.remade = []

        class FakeDatabase:
            def __init__(self, url_or_alias): pass
            def get_all_refences_data(self, reference_type):
                return [{'inchi': inchi, 'version': INCHI_REFERENCE_VERSION} for inchi in ['A', 'B']]

        def get_inchi_reference(inchi, cache):
            self.remade.append(inchi)
            return {'inchi': inchi, 'remade': True}

        monkeypatch.setattr(update_references, 'Database', FakeDatabase)
        monkeypatch.setattr(update_references, 'get_inchi_reference', get_inchi_reference)
        monkeypatch.chdir(tmp_path)
        self.directory = tmp_path / 'mddb.test'
        self.directory.mkdir()
        # A reference already remade in a previous run
        self.output = self.directory / 'inchikey_references.json'
        save_json([{'inchi': 'A', 'remade': False}], str(self.output))

    def test_update_skips_latest_and_locally_remade(self):
        self.update_references.update_references('https://mddb.test/api/', 'inchi')
        assert self.remade == []
        assert load_json(str(self.output)) == [{'inchi': 'A', 'remade': False}]

    def test_forced_update_remakes_all(self):
        """Forced updates remake all references, also those already remade in the local output file."""
        self.update_references.update_references('https://mddb.test/api/', 'inchi', force=True)
        assert self.remade == ['A', 'B']
        assert load_json(str(self.output)) == [{'inchi': 'A', 'remade': True}, {'inchi': 'B', 'remade': True}]
        # The previous output file is kept apart
        previous = [path for path in self.directory.iterdir() if path.name.endswith('_inchikey_references.json')]
        assert len(previous) == 1 and load_json(str(previous[0])) == [{'inchi': 'A', 'remade': False}]
