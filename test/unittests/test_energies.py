"""Energies analysis (per-atom CMIP interaction energies) on a tiny generated system.

The system has 37 atoms and two interactions, with no downloads:
- salt bridge: methylammonium (agent 1) and acetate (agent 2). They are at contact in frames 1 and 3, and 40 Å apart in
  frame 2, which is beyond FAR_AGENTS_CUTOFF, so CMIP is skipped and the energies of that frame are 0.
- peptide: Ace-Ala (agent 1) and Nme (agent 2), covalently bonded. The two atoms of the bond between agents (a strong
  bond) get no energies (None).

Refactors must reproduce the golden output (atol 1e-6, one OpenMP thread). Regenerate it only in the commit of an
intended numerical change: `python test/unittests/test_energies.py`.
The tests marked xfail(strict=True) are known bugs. They pass, and so fail here, once the bug is fixed: remove the mark then.
"""
import os
import re
import sys
import json
import contextlib
import io
from pathlib import Path
import numpy as np
import pytest
import mdtraj

# Run as a script, import the workflow of this repository, not another installed copy
if __name__ == '__main__':
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from mddb_workflow.analyses.energies import energies, mine_cmip_output, compute_new_grid, format_data
from mddb_workflow.analyses.energies import set_cmip_elements, name_terminal_residues, OUTPUT_MINING_EXCEPTION
from mddb_workflow.analyses.energies import cap_grid_points, mine_cmip_atom_energies, GRID_POINTS_LIMIT
from mddb_workflow.analyses.energies import get_cmip_energies, write_cmip_vdw_table, zero_clashes
from mddb_workflow.analyses import energies as energies_module
from mddb_workflow.utils.constants import CMIP_VDW_SOURCE
from mddb_workflow.utils.constants import GLOBALS, ENERGIES_ANALYSIS_VERSION
from mddb_workflow.utils.structures import Structure, Atom, Residue, Chain
from mddb_workflow.utils.selections import Selection
from mddb_workflow.utils.file import File

pytestmark = pytest.mark.CI

GOLDEN_FILEPATH = Path(__file__).parent / 'test_energies_golden.json'
GRIDFREE_GOLDEN_FILEPATH = Path(__file__).parent / 'test_energies_gridfree_golden.json'

# Chain, residue name, residue number, atom name, element and charge of every atom
ATOMS = [
    ('A', 'MAM', 1, 'C1', 'C', 0.01), ('A', 'MAM', 1, 'N1', 'N', -0.30),
    ('A', 'MAM', 1, 'H1', 'H', 0.10), ('A', 'MAM', 1, 'H2', 'H', 0.10), ('A', 'MAM', 1, 'H3', 'H', 0.10),
    ('A', 'MAM', 1, 'HN1', 'H', 0.33), ('A', 'MAM', 1, 'HN2', 'H', 0.33), ('A', 'MAM', 1, 'HN3', 'H', 0.33),
    ('B', 'ACT', 1, 'C1', 'C', -0.30), ('B', 'ACT', 1, 'C2', 'C', 0.75), ('B', 'ACT', 1, 'O1', 'O', -0.80),
    ('B', 'ACT', 1, 'O2', 'O', -0.80), ('B', 'ACT', 1, 'H1', 'H', 0.05), ('B', 'ACT', 1, 'H2', 'H', 0.05),
    ('B', 'ACT', 1, 'H3', 'H', 0.05),
    ('C', 'ACE', 1, 'CH3', 'C', -0.3662), ('C', 'ACE', 1, 'C', 'C', 0.5972), ('C', 'ACE', 1, 'O', 'O', -0.5679),
    ('C', 'ACE', 1, 'H1', 'H', 0.1123), ('C', 'ACE', 1, 'H2', 'H', 0.1123), ('C', 'ACE', 1, 'H3', 'H', 0.1123),
    ('C', 'ALA', 2, 'N', 'N', -0.4157), ('C', 'ALA', 2, 'H', 'H', 0.2719), ('C', 'ALA', 2, 'CA', 'C', 0.0337),
    ('C', 'ALA', 2, 'HA', 'H', 0.0823), ('C', 'ALA', 2, 'CB', 'C', -0.1825), ('C', 'ALA', 2, 'HB1', 'H', 0.0603),
    ('C', 'ALA', 2, 'HB2', 'H', 0.0603), ('C', 'ALA', 2, 'HB3', 'H', 0.0603), ('C', 'ALA', 2, 'C', 'C', 0.5973),
    ('C', 'ALA', 2, 'O', 'O', -0.5679),
    ('C', 'NME', 3, 'N', 'N', -0.4157), ('C', 'NME', 3, 'H', 'H', 0.2719), ('C', 'NME', 3, 'C', 'C', -0.1490),
    ('C', 'NME', 3, 'H1', 'H', 0.0976), ('C', 'NME', 3, 'H2', 'H', 0.0976), ('C', 'NME', 3, 'H3', 'H', 0.0976),
]
BONDS = [
    (0, 1), (0, 2), (0, 3), (0, 4), (1, 5), (1, 6), (1, 7),
    (8, 9), (8, 12), (8, 13), (8, 14), (9, 10), (9, 11),
    (15, 16), (15, 18), (15, 19), (15, 20), (16, 17), (16, 21),
    (21, 22), (21, 23), (23, 24), (23, 25), (25, 26), (25, 27), (25, 28), (23, 29), (29, 30), (29, 31),
    (31, 32), (31, 33), (33, 34), (33, 35), (33, 36),
]
MAM, ACT, ACE_ALA, NME = list(range(0, 8)), list(range(8, 15)), list(range(15, 31)), list(range(31, 37))
# The peptide bond between Ala (agent 1) and Nme (agent 2)
STRONG_BOND = [29, 31]
# Salt bridge: the acetate O1 is 1.75 Å from the methylammonium HN1, along the N1-HN1 bond (RDKit + MMFF geometries)
FRAME = np.array([
    (-0.749, 0.080, -0.013), (0.718, -0.077, 0.012), (-1.147, -0.619, -0.751), (-0.970, 1.112, -0.291),
    (-1.126, -0.147, 0.986), (0.988, -1.035, 0.270), (1.152, 0.560, 0.694), (1.133, 0.125, -0.907),
    (1.558, -4.890, -0.053), (1.779, -3.842, 1.024), (1.448, -2.666, 0.708), (2.268, -4.285, 2.101),
    (2.511, -5.362, -0.312), (1.134, -4.447, -0.959), (0.869, -5.657, 0.313),
    (2.603, 21.539, 0.683), (1.733, 20.542, -0.031), (1.720, 20.433, -1.254), (3.604, 21.522, 0.242),
    (2.694, 21.300, 1.746), (2.168, 22.536, 0.576), (0.921, 19.803, 0.795), (0.738, 20.163, 1.723),
    (-0.031, 18.823, 0.265), (0.406, 18.354, -0.625), (-0.297, 17.745, 1.308), (-0.748, 18.168, 2.213),
    (0.635, 17.250, 1.602), (-0.980, 16.982, 0.920), (-1.345, 19.510, -0.180), (-2.455, 19.155, 0.210),
    (-1.166, 20.512, -1.114), (-0.245, 20.582, -1.547), (-2.282, 21.151, -1.758), (-2.662, 20.488, -2.540),
    (-3.076, 21.351, -1.033), (-1.937, 22.088, -2.202),
])


def shift_acetate(shift: tuple) -> np.ndarray:
    """Get the frame with the acetate moved by the given vector (Å)."""
    frame = FRAME.copy()
    frame[ACT] += shift
    return frame


# Contact, far apart (beyond FAR_AGENTS_CUTOFF), and contact again a bit further
FRAMES = np.array([FRAME, shift_acetate((40.0, 0.0, 0.0)), shift_acetate((0.3, -0.2, 0.4))])
# The acetate O1 at 0.53 Å from the methylammonium HN1: a clash
CLASH_FRAMES = np.array([shift_acetate((-0.322, 1.142, -0.307))])
INTERACTIONS = [
    dict(name='salt bridge', agent_1='methylammonium', agent_2='acetate', atom_indices_1=MAM, atom_indices_2=ACT,
         strong_bonds=[]),
    dict(name='peptide', agent_1='Ace-Ala', agent_2='Nme', atom_indices_1=ACE_ALA, atom_indices_2=NME,
         strong_bonds=[STRONG_BOND]),
]
AGENT_KEYS = ['vdw', 'es', 'both', 'ivdw', 'ies', 'iboth', 'fvdw', 'fes', 'fboth']


def build_structure(elements: dict = {}) -> Structure:
    """Build the structure of the tiny system, with its bonds, and optionally other elements for some atoms."""
    atoms, residues, chains = [], [], []
    for index, (chain_name, residue_name, residue_number, atom_name, element, _) in enumerate(ATOMS):
        element = elements.get(index, element)
        atoms.append(Atom(name=atom_name, element=element, coords=tuple(FRAME[index])))
        if not residues or (chains[-1].name, residues[-1].number) != (chain_name, residue_number):
            if not chains or chains[-1].name != chain_name:
                chains.append(Chain(name=chain_name))
                chains[-1].residue_indices = []
            residues.append(Residue(name=residue_name, number=residue_number, icode=''))
            residues[-1].atom_indices = []
            chains[-1].residue_indices = chains[-1].residue_indices + [len(residues) - 1]
        residues[-1].atom_indices = residues[-1].atom_indices + [index]
    structure = Structure(atoms=atoms, residues=residues, chains=chains)
    bonds = [[] for _ in atoms]
    for a, b in BONDS:
        bonds[a].append(b)
        bonds[b].append(a)
    structure.bonds = bonds
    return structure


def run_energies(work_directory: Path, frames: np.ndarray, interactions: list[dict], verbose: bool = False,
                 elements: dict = {}, engine: str = 'cmip') -> tuple[Path, str]:
    """Run the energies analysis with one OpenMP thread and return its output directory and its logs.

    CMIP leaves fort.* files in the current directory, so the analysis runs from the work directory.
    """
    trajectory_filepath = str(work_directory / 'trajectory.dcd')
    with mdtraj.formats.DCDTrajectoryFile(trajectory_filepath, 'w') as writer:
        writer.write(frames.astype(np.float32))
    output_directory = work_directory / 'energies'
    previous_directory, previous_threads = os.getcwd(), os.environ.get('OMP_NUM_THREADS')
    previous_working_directory = GLOBALS['working_directory']
    logs = io.StringIO()
    try:
        os.chdir(work_directory)
        os.environ['OMP_NUM_THREADS'] = '1'
        GLOBALS['working_directory'] = str(work_directory)
        with contextlib.redirect_stdout(logs):
            energies(
                trajectory_file=File(trajectory_filepath),
                output_directory=str(output_directory),
                structure=build_structure(elements),
                interactions=[dict(interaction) for interaction in interactions],
                charges=[atom[5] for atom in ATOMS],
                cg_selection=Selection([]),
                dummy_selection=Selection([]),
                snapshots=len(frames),
                verbose=verbose,
                energies_engine=engine,
            )
    finally:
        os.chdir(previous_directory)
        GLOBALS['working_directory'] = previous_working_directory
        if previous_threads is None:
            os.environ.pop('OMP_NUM_THREADS', None)
        else:
            os.environ['OMP_NUM_THREADS'] = previous_threads
    return output_directory, logs.getvalue()


def load_outputs(output_directory: Path) -> dict:
    """Load the output of every interaction, by analysis name."""
    summary = json.load(open(output_directory / 'mda.energies.json'))
    return {entry['analysis']: json.load(open(output_directory / f'mda.energies_{entry["analysis"][-2:]}.json'))
            for entry in summary}


@pytest.fixture(scope='module')
def run(tmp_path_factory) -> tuple[Path, dict]:
    """Run the analysis once on the three frames."""
    output_directory, _ = run_energies(tmp_path_factory.mktemp('energies'), FRAMES, INTERACTIONS)
    return output_directory, load_outputs(output_directory)


@pytest.fixture(scope='module')
def clash_run(tmp_path_factory) -> tuple[dict, str]:
    """Run the analysis once on the salt bridge with a clash."""
    output_directory, logs = run_energies(tmp_path_factory.mktemp('clash'), CLASH_FRAMES, INTERACTIONS[:1])
    return load_outputs(output_directory), logs


def assert_same_values(values: list, expected: list, label: str):
    """Assert two lists of energies are equal (atol 1e-6), with None in the same places."""
    assert len(values) == len(expected), label
    assert [v is None for v in values] == [e is None for e in expected], label
    numbers = np.array([v for v in values if v is not None], dtype=float)
    expected_numbers = np.array([e for e in expected if e is not None], dtype=float)
    assert np.allclose(numbers, expected_numbers, rtol=0, atol=1e-6), label


# ----- Whole analysis -----

def test_golden(run):
    """The per-atom output is the stored golden output."""
    _, outputs = run
    golden = json.load(open(GOLDEN_FILEPATH))
    assert sorted(outputs) == sorted(golden)
    for analysis, output in outputs.items():
        assert output['name'] == golden[analysis]['name']
        for agent in ('agent1', 'agent2'):
            assert output[agent]['step'] == golden[analysis][agent]['step']
            for key in AGENT_KEYS:
                assert_same_values(output[agent][key], golden[analysis][agent][key], f'{analysis} {agent} {key}')


def test_schema(run):
    """The summary lists one output per interaction, and each output has one value per agent atom."""
    output_directory, outputs = run
    summary = json.load(open(output_directory / 'mda.energies.json'))
    assert summary == [{'name': 'salt bridge', 'analysis': 'energies-00'}, {'name': 'peptide', 'analysis': 'energies-01'}]
    for interaction, analysis in zip(INTERACTIONS, ('energies-00', 'energies-01')):
        output = outputs[analysis]
        assert sorted(output) == ['agent1', 'agent2', 'name', 'version']
        assert output['name'] == interaction['name']
        assert output['version'] == ENERGIES_ANALYSIS_VERSION
        for agent, atom_indices in (('agent1', interaction['atom_indices_1']), ('agent2', interaction['atom_indices_2'])):
            assert sorted(output[agent]) == sorted(AGENT_KEYS + ['step', 'zeroed'])
            assert output[agent]['step'] == 1
            assert output[agent]['zeroed'] == [0, 0, 0]
            for key in AGENT_KEYS:
                assert len(output[agent][key]) == len(atom_indices)


def test_salt_bridge_attracts(run):
    """The electrostatic energy of the salt bridge is negative for both agents, in both frames at contact."""
    salt_bridge = run[1]['energies-00']
    for agent in ('agent1', 'agent2'):
        assert sum(salt_bridge[agent]['ies']) < -3
        assert sum(salt_bridge[agent]['fes']) < -3


def test_far_frame_is_zero(run):
    """The frame with the agents far apart adds 0: the average of three frames is the sum of the first and last over 3.

    With 3 frames, the first and last 20% are one frame each.
    """
    salt_bridge = run[1]['energies-00']
    for agent in ('agent1', 'agent2'):
        for energy in ('vdw', 'es', 'both'):
            first, last = np.array(salt_bridge[agent][f'i{energy}']), np.array(salt_bridge[agent][f'f{energy}'])
            assert np.allclose(salt_bridge[agent][energy], (first + last) / 3, rtol=0, atol=1e-12)
            assert not np.allclose(first, last)


def test_strong_bond_atoms_are_none(run):
    """The atoms of a covalent bond between agents get no energies, and every other atom does."""
    peptide = run[1]['energies-01']
    for agent, atom_indices in (('agent1', ACE_ALA), ('agent2', NME)):
        for key in AGENT_KEYS:
            nones = [atom_indices[i] for i, value in enumerate(peptide[agent][key]) if value is None]
            assert nones == [index for index in atom_indices if index in STRONG_BOND]


def test_clash_atoms_are_zeroed(clash_run):
    """Atoms whose energy is above 100 kcal/mol are set to 0, and counted in the output."""
    outputs, _ = clash_run
    salt_bridge = outputs['energies-00']
    hydrogen = MAM.index(5)
    oxygen = ACT.index(10)
    assert salt_bridge['agent1']['both'][hydrogen] == 0
    assert salt_bridge['agent2']['both'][oxygen] == 0
    # The output has the number of atoms set to 0 in each frame
    assert salt_bridge['agent1']['zeroed'][0] >= 1 and salt_bridge['agent2']['zeroed'][0] >= 1


def test_cmip_runs_in_the_output_directory(tmp_path):
    """CMIP runs in the output directory: a fort.* file of the caller survives, and CMIP leaves no fort.* files."""
    (tmp_path / 'fort.99').write_text('not a CMIP file')
    output_directory, _ = run_energies(tmp_path, CLASH_FRAMES, INTERACTIONS[:1])
    assert (tmp_path / 'fort.99').read_text() == 'not a CMIP file'
    assert not list(output_directory.glob('fort.*'))
    assert (output_directory / 'vdwprm').exists()


def test_interactions_over_the_atom_limit(tmp_path, monkeypatch):
    """Interactions over the CMIP atom limit are skipped, even all of them, and the input interactions do not change."""
    monkeypatch.setattr(energies_module, 'CMIP_ATOM_LIMIT', 7)
    interactions = [dict(interaction) for interaction in INTERACTIONS]
    output_directory, logs = run_energies(tmp_path, FRAMES[:1], interactions)
    assert 'atom count limit' in logs
    assert not (output_directory / 'mda.energies.json').exists()
    assert interactions == INTERACTIONS


def test_missing_strong_bonds(tmp_path):
    """Interactions with no strong bonds key are calculated as having no strong bonds."""
    interaction = {key: value for key, value in INTERACTIONS[0].items() if key != 'strong_bonds'}
    output_directory, _ = run_energies(tmp_path, FRAMES, [interaction])
    golden = json.load(open(GOLDEN_FILEPATH))['energies-00']
    output = load_outputs(output_directory)['energies-00']
    for agent in ('agent1', 'agent2'):
        for key in AGENT_KEYS:
            assert_same_values(output[agent][key], golden[agent][key], f'{agent} {key}')


def test_verbose_with_strong_bonds(tmp_path):
    """The verbose totals skip the atoms of strong bonds, which have no energies."""
    _, logs = run_energies(tmp_path, FRAMES[:1], INTERACTIONS[1:], verbose=True)
    assert 'Total energies for Ace-Ala' in logs


def test_failed_cmip_does_not_reuse_previous_output(tmp_path):
    """If CMIP fails to write its output, the output of a previous run is not read as if it was new."""
    for filename in ('vdwprm', 'cmip.in', 'guest.pdb', 'host.pdb'):
        (tmp_path / filename).write_text('not a CMIP input\n')
    (tmp_path / 'cmip_output.pdb').write_text(
        'ATOM      1  H   NME C   3      0.2719    -0.6282E-01     -3.398         -3.461    \n')
    with pytest.raises(RuntimeError):
        get_cmip_energies(File(str(tmp_path / 'cmip.in')), File(str(tmp_path / 'guest.pdb')),
                          File(str(tmp_path / 'host.pdb')), str(tmp_path), File(str(tmp_path / 'restart')))


def test_elements_missing_from_the_cmip_table(tmp_path):
    """Elements missing from the CMIP table get their AmberTools parameters, placed before the one-letter types.

    Without them CMIP would silently use the carbon parameters for calcium.
    """
    energies_by_element = {}
    for element in ('Ca', 'C'):
        (tmp_path / element).mkdir()
        output_directory, logs = run_energies(tmp_path / element, FRAMES[:1], INTERACTIONS[:1], elements={10: element})
        energies_by_element[element] = load_outputs(output_directory)['energies-00']['agent2']['vdw'][ACT.index(10)]
        table = (output_directory / 'vdwprm').read_text().splitlines()
        if element == 'Ca':
            assert table[1] == 'Ca 0.46 3.05 40.08  20  0.000'
            assert re.search(r'WARNING.*Elements Ca are not in the CMIP van der Waals table', logs)
        else:
            assert table[1].startswith('Cl ')
    assert energies_by_element['Ca'] != energies_by_element['C']


def test_elements_with_no_parameters(tmp_path):
    """Interactions with elements which have no van der Waals parameters for CMIP are skipped with a warning."""
    output_directory, logs = run_energies(tmp_path, FRAMES[:1], INTERACTIONS[:1], elements={10: 'Si'})
    assert re.search(r'WARNING: .*salt bridge has atoms of elements with no van der Waals parameters for CMIP \(Si\)', logs)
    assert not (output_directory / 'mda.energies.json').exists()


def test_selenium_takes_the_sulfur_type(tmp_path):
    """Selenium takes the CMIP type of sulfur, with a warning, as CMIP gave it by prefix up to version 1.1.0."""
    outputs = {}
    for element in ('Se', 'S'):
        (tmp_path / element).mkdir()
        output_directory, logs = run_energies(tmp_path / element, FRAMES[:1], INTERACTIONS[:1], elements={10: element})
        outputs[element] = load_outputs(output_directory)
        if element == 'Se':
            assert re.search(r'WARNING: .*1 Se atoms take the CMIP type of S', logs)
    assert outputs['Se'] == outputs['S']


def test_write_cmip_vdw_table(tmp_path):
    """Without extra elements the table is the CMIP one; extra elements go first, with sigma = 2 R*/2 / 2^(1/6)."""
    write_cmip_vdw_table(str(tmp_path / 'vdwprm'), set())
    assert (tmp_path / 'vdwprm').read_text() == open(CMIP_VDW_SOURCE).read()
    write_cmip_vdw_table(str(tmp_path / 'vdwprm'), {'Ca', 'Fe'})
    lines = (tmp_path / 'vdwprm').read_text().splitlines()
    assert lines[1:3] == ['Ca 0.46 3.05 40.08  20  0.000', 'Fe 0.01 2.41 55.84  26  0.000']
    assert lines[3:] == open(CMIP_VDW_SOURCE).read().splitlines()[1:]


# ----- Grid-free engine -----

@pytest.fixture(scope='module')
def gridfree_run(tmp_path_factory) -> tuple[Path, dict]:
    """Run the analysis once on the three frames with the grid-free engine."""
    output_directory, _ = run_energies(tmp_path_factory.mktemp('gridfree'), FRAMES, INTERACTIONS, engine='gridfree')
    return output_directory, load_outputs(output_directory)


def test_gridfree_golden(gridfree_run):
    """The grid-free per-atom output is the stored golden output."""
    _, outputs = gridfree_run
    golden = json.load(open(GRIDFREE_GOLDEN_FILEPATH))
    for analysis, output in outputs.items():
        for agent in ('agent1', 'agent2'):
            for key in AGENT_KEYS:
                assert_same_values(output[agent][key], golden[analysis][agent][key], f'{analysis} {agent} {key}')


def direct_pair_energies(frame: np.ndarray, guest: list[int], host: list[int], excluded: set = set()) -> list:
    """Compute the grid-free energies of the guest atoms with a plain double loop, as an independent check."""
    structure = build_structure()
    set_cmip_elements(structure)
    table = {line[:3].strip(): (float(line[3:7]), float(line[8:12])) for line in open(CMIP_VDW_SOURCE).read().splitlines()[1:]}
    energies = []
    for i in guest:
        if i in excluded:
            energies.append(None)
            continue
        vdw = es = 0.0
        for j in host:
            if j in excluded:
                continue
            r = float(np.linalg.norm(frame[i] - frame[j]))
            b = 78.4 + 8.5525
            dielectric = -8.5525 + b / (1 + 7.7839 * np.exp(-0.003627 * b * r))
            es += 332.0637 * ATOMS[i][5] * ATOMS[j][5] / (dielectric * r)
            (epsilon_i, sigma_i), (epsilon_j, sigma_j) = table[structure.atoms[i].element], table[structure.atoms[j].element]
            sigma, epsilon = 0.8 * (sigma_i + sigma_j) / 2, (epsilon_i * epsilon_j) ** 0.5
            if r <= 8:
                vdw += 4 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)
        energies.append((vdw, es))
    return energies


@pytest.mark.parametrize('analysis, guest, host, excluded', [
    ('energies-00', MAM, ACT, set()),
    ('energies-01', NME, ACE_ALA, set(STRONG_BOND)),
], ids=['salt bridge', 'peptide'])
def test_gridfree_matches_direct_sum(gridfree_run, analysis, guest, host, excluded):
    """The grid-free energies of the first frame are the plain sums over the other agent's atoms."""
    _, outputs = gridfree_run
    agent = 'agent1' if guest[0] < host[0] else 'agent2'
    for value, vdw, es in zip(direct_pair_energies(FRAMES[0], guest, host, excluded),
                              outputs[analysis][agent]['ivdw'], outputs[analysis][agent]['ies']):
        if value is None:
            assert vdw is None and es is None
        else:
            assert np.allclose(value, (vdw, es), rtol=0, atol=1e-9)


def test_gridfree_directions_agree(gridfree_run):
    """Both agents get the full energy of each pair, so both directions give the same totals."""
    _, outputs = gridfree_run
    for analysis in outputs:
        for key in ('ivdw', 'ies', 'fvdw', 'fes'):
            totals = [sum(v for v in outputs[analysis][agent][key] if v is not None) for agent in ('agent1', 'agent2')]
            assert abs(totals[0] - totals[1]) < 1e-9


def test_gridfree_far_frame_and_salt_bridge(gridfree_run):
    """The far frame adds 0 and the salt bridge attracts, as with CMIP."""
    _, outputs = gridfree_run
    salt_bridge = outputs['energies-00']
    for agent in ('agent1', 'agent2'):
        first, last = np.array(salt_bridge[agent]['ies']), np.array(salt_bridge[agent]['fes'])
        assert np.allclose(salt_bridge[agent]['es'], (first + last) / 3, rtol=0, atol=1e-12)
        assert sum(salt_bridge[agent]['ies']) < -3


def test_unknown_engine(tmp_path):
    """An unknown engine stops the analysis."""
    with pytest.raises(ValueError):
        run_energies(tmp_path, FRAMES[:1], INTERACTIONS[:1], engine='unknown')


# ----- Backup -----

class Interrupted(Exception):
    """Stops the analysis in the middle of a run, as a kill would."""


def watch_frame_energies(monkeypatch, stop_after: int | None = None) -> list:
    """Count the interaction frames the analysis calculates, and optionally stop it after this number of them."""
    calls = []
    calculate = energies_module.get_interaction_frame_energies

    def watched(*args, **kwargs):
        if stop_after is not None and len(calls) == stop_after:
            raise Interrupted()
        calls.append(1)
        return calculate(*args, **kwargs)
    monkeypatch.setattr(energies_module, 'get_interaction_frame_energies', watched)
    return calls


def interrupt(tmp_path: Path, monkeypatch, frames_done: int) -> Path:
    """Run the grid-free analysis and stop it after some frames; return its backup."""
    with monkeypatch.context() as patch:
        watch_frame_energies(patch, stop_after=frames_done * len(INTERACTIONS))
        with pytest.raises(Interrupted):
            run_energies(tmp_path, FRAMES, INTERACTIONS, engine='gridfree')
    return tmp_path / 'energies' / energies_module.ENERGIES_BACKUP_FILENAME


@pytest.mark.parametrize('frames_done', [1, 2])
@pytest.mark.parametrize('cut_line', [False, True], ids=['complete lines', 'cut last line'])
def test_backup_resumes(tmp_path, monkeypatch, gridfree_run, frames_done, cut_line):
    """An interrupted analysis resumes from its backup and gives the output of an uninterrupted run."""
    backup = interrupt(tmp_path, monkeypatch, frames_done)
    lines = backup.read_bytes().splitlines(keepends=True)
    # A header line and one line per frame done
    assert len(lines) == 1 + frames_done
    if cut_line:
        # The next frame was being written when the analysis was killed
        with open(backup, 'ab') as file:
            file.write(lines[-1][:len(lines[-1]) // 2])
    calls = watch_frame_energies(monkeypatch)
    output_directory, logs = run_energies(tmp_path, FRAMES, INTERACTIONS, engine='gridfree')
    assert f'Recovering {frames_done} frames' in logs
    assert len(calls) == (len(FRAMES) - frames_done) * len(INTERACTIONS)
    assert load_outputs(output_directory) == gridfree_run[1]
    assert not backup.exists()


def test_backup_of_another_setup_is_not_reused(tmp_path, monkeypatch, gridfree_run):
    """A backup written with another engine is discarded and every frame is calculated again."""
    backup = interrupt(tmp_path, monkeypatch, 2)
    lines = backup.read_bytes().splitlines(keepends=True)
    header = json.loads(lines[0])
    header['engine'] = 'cmip'
    backup.write_bytes(json.dumps(header).encode() + b'\n' + b''.join(lines[1:]))
    calls = watch_frame_energies(monkeypatch)
    output_directory, logs = run_energies(tmp_path, FRAMES, INTERACTIONS, engine='gridfree')
    assert 'another setup' in logs
    assert len(calls) == len(FRAMES) * len(INTERACTIONS)
    assert load_outputs(output_directory) == gridfree_run[1]


def test_coordinate_frames_match_pdb_frames(tmp_path, monkeypatch):
    """The grid-free frames, read from the trajectory, have exactly the coordinates of the PDB frames CMIP uses."""
    trajectory_filepath = str(tmp_path / 'trajectory.dcd')
    # Coordinates with more than 3 decimals, so the PDB rounds them
    frames = FRAMES + np.array([0.00041, -0.00062, 0.00051])
    with mdtraj.formats.DCDTrajectoryFile(trajectory_filepath, 'w') as writer:
        writer.write(frames.astype(np.float32))
    monkeypatch.chdir(tmp_path)
    monkeypatch.setitem(GLOBALS, 'working_directory', str(tmp_path))
    structure = build_structure()
    pdb_frames, *pdb_counts = energies_module.get_pdb_frame_structures(structure, trajectory_filepath, len(frames), 100)
    pdb_coordinates = [[atom.coords for atom in frame.atoms] for frame in pdb_frames]
    frame_structures, *counts = energies_module.get_coordinate_frame_structures(
        structure, trajectory_filepath, len(frames), 100)
    coordinates = [[atom.coords for atom in frame.atoms] for frame in frame_structures]
    assert counts == pdb_counts
    assert len(coordinates) == len(frames)
    assert coordinates == pdb_coordinates
    assert not np.allclose(np.array(coordinates), frames, rtol=0, atol=1e-5)


def test_old_backup_is_removed(tmp_path, gridfree_run):
    """A backup in the format of version 1.1.0 or older is removed, not reused."""
    old_backup = tmp_path / 'energies' / energies_module.OLD_ENERGIES_BACKUP_FILENAME
    old_backup.parent.mkdir()
    old_backup.write_text(json.dumps([[], []]))
    output_directory, logs = run_energies(tmp_path, FRAMES, INTERACTIONS, engine='gridfree')
    assert 'older version' in logs
    assert not old_backup.exists()
    assert load_outputs(output_directory) == gridfree_run[1]


# ----- Helpers -----

CMIP_LOGS = """
 AUTOMATIC GRID
 --------------
 Fill:    0.800
Automatic Grid
 Grid Size:      10.000  10.000  10.000  90.000  90.000  90.000
 Grid density:    20   20   20
 Grid units:      0.500   0.500   0.500
 Grid origin:    -4.997  -4.962  -4.960
 Grid center:     0.003  -10.038   0.040
       9261 Grid points
"""


def test_mine_cmip_output():
    """The automatic grid center, density and units are read from the CMIP logs."""
    assert mine_cmip_output(CMIP_LOGS) == ((0.003, -10.038, 0.04), (20.0, 20.0, 20.0), (0.5, 0.5, 0.5))


def test_mine_cmip_output_fails():
    """Logs without the grid return the mining exception."""
    logs = CMIP_LOGS.replace('Grid units', 'Grid unknown')
    assert mine_cmip_output(logs) is OUTPUT_MINING_EXCEPTION


@pytest.mark.parametrize('line, expected', [
    ('ATOM      1  H   NME C   3      0.2719    -0.6282E-01     -3.398         -3.461    \n', (-0.06282, -3.398, -3.461)),
    ('ATOM      2  C   NME C   3     -0.1490    -0.2134         0.8978         0.6844    \n', (-0.2134, 0.8978, 0.6844)),
    ('ATOM      6  HN1 MAM A   1      0.3300     0.4512E+05     -12.34          0.4511E+05\n', (45120, -12.34, 45110)),
], ids=['E notation', 'decimal', 'above 100 kcal/mol'])
def test_mine_cmip_atom_energies(line, expected):
    """The vdw, es and both energies are the last three numbers of a CMIP output line."""
    assert mine_cmip_atom_energies(line) == expected


def test_zero_clashes():
    """Atoms above 100 kcal/mol are set to 0 and counted; large negative energies and None are kept."""
    energies = [(1.0, 2.0, 3.0), None, (200.0, -10.0, 190.0), (-10.0, -100.0, -110.0)]
    assert zero_clashes(energies) == ([(1.0, 2.0, 3.0), None, (0, 0, 0), (-10.0, -100.0, -110.0)], 1)


def test_cap_grid_points_below_limit():
    """Grids below the limit of points are not changed."""
    assert cap_grid_points([100, 200, 300], 0.5) == ([100, 200, 300], 0.5)


def test_set_cmip_elements():
    """Hydrogens are typed by their bonded heavy atom: H (carbon), HO (oxygen) and HN (nitrogen or sulfur)."""
    structure = build_structure()
    set_cmip_elements(structure)
    elements = [atom.element for atom in structure.atoms]
    assert elements[:8] == ['C', 'N', 'H', 'H', 'H', 'HN', 'HN', 'HN']
    assert elements[21:23] == ['N', 'HN']


def test_set_cmip_elements_oxygen_and_sulfur():
    """Hydrogens on oxygen are HO and hydrogens on sulfur are HN."""
    structure = Structure(atoms=[Atom('O', 'O'), Atom('HO', 'H'), Atom('S', 'S'), Atom('HS', 'H')])
    structure.bonds = [[1], [0], [3], [2]]
    set_cmip_elements(structure)
    assert [atom.element for atom in structure.atoms] == ['O', 'HO', 'S', 'HN']


def test_set_cmip_elements_two_bonds():
    """A hydrogen must have exactly one bond."""
    structure = Structure(atoms=[Atom('P', 'P'), Atom('H1', 'H'), Atom('C', 'C')])
    structure.bonds = [[1], [0, 2], [1]]
    with pytest.raises(ValueError):
        set_cmip_elements(structure)


def test_set_cmip_elements_other_elements(capsys):
    """Hydrogens bonded to atoms other than C, N, O or S are typed as hydrogens bonded to carbon, with a warning."""
    structure = Structure(atoms=[Atom('P', 'P'), Atom('H1', 'H'), Atom('B', 'B'), Atom('H2', 'H')])
    structure.bonds = [[1], [0], [3], [2]]
    set_cmip_elements(structure)
    assert [atom.element for atom in structure.atoms] == ['P', 'H', 'B', 'H']
    assert re.search(r'WARNING.*1 bonded to B, 1 bonded to P', capsys.readouterr().out)


def make_chain_structure(residue_names: list[str]) -> Structure:
    """Build a structure with one chain of one-atom residues."""
    atoms = [Atom('X', 'C') for _ in residue_names]
    residues = [Residue(name=name, number=n + 1, icode='') for n, name in enumerate(residue_names)]
    for n, residue in enumerate(residues):
        residue.atom_indices = [n]
    chain = Chain(name='A')
    chain.residue_indices = list(range(len(residues)))
    return Structure(atoms=atoms, residues=residues, chains=[chain])


@pytest.mark.parametrize('names, expected', [
    (['ALA', 'GLY', 'SER'], ['ALAN', 'GLY', 'SERC']),
    (['ALAN', 'GLY', 'SERC'], ['ALAN', 'GLY', 'SERC']),
    (['ACE', 'ALA', 'NME'], ['ACE', 'ALA', 'NME']),
    (['DA', 'DC', 'DG'], ['DA5', 'DC', 'DG3']),
    (['LIG'], ['LIG']),
])
def test_name_terminal_residues(names, expected):
    """The first and last residues of protein and nucleic chains get terminal names, unless they already have them."""
    structure = make_chain_structure(names)
    name_terminal_residues(structure)
    assert [residue.name for residue in structure.residues] == expected


def test_format_data():
    """Atom averages over all frames, the first 20% and the last 20%, skipping None."""
    frames = [
        [(1.0, 2.0, 3.0), None, (0.0, 0.0, 0.0)],
        [(3.0, 4.0, 7.0), None, (1.0, 1.0, 2.0)],
        [(5.0, 6.0, 11.0), None, None],
        [(7.0, 8.0, 15.0), None, (3.0, 3.0, 6.0)],
        [(9.0, 10.0, 19.0), None, (4.0, 4.0, 8.0)],
    ]
    output = format_data(frames, 7)
    assert output == {
        'vdw': [5.0, None, 2.0], 'es': [6.0, None, 2.0], 'both': [11.0, None, 4.0],
        'ivdw': [1.0, None, 0.0], 'ies': [2.0, None, 0.0], 'iboth': [3.0, None, 0.0],
        'fvdw': [9.0, None, 4.0], 'fes': [10.0, None, 4.0], 'fboth': [19.0, None, 8.0],
        'step': 7,
    }


@pytest.mark.parametrize('frame_count, edge_frames', [(1, 1), (2, 1), (3, 1), (7, 1), (8, 2), (12, 2), (13, 3)])
def test_format_data_edges(frame_count, edge_frames):
    """The first and last 20% of frames, rounded and at least one frame."""
    frames = [[(float(f), 0.0, 0.0)] for f in range(frame_count)]
    output = format_data(frames, 1)
    assert output['ivdw'] == [np.mean(range(edge_frames))]
    assert output['fvdw'] == [np.mean(range(frame_count - edge_frames, frame_count))]


def test_format_data_backup_lists():
    """Frames restored from the JSON backup have lists instead of tuples."""
    output = format_data([[[1.0, 2.0, 3.0], None]], 1)
    assert output['both'] == [3.0, None]


# ----- Bugs B1 to B8 of energies version 1.1.0 -----

def read_cmip_pdb(filepath: Path) -> list[tuple[tuple, str]]:
    """Read the atoms of a CMIP input PDB: (residue name, residue number, atom name) and the element tail."""
    return [((line[17:20], int(line[22:26]), line[12:16].strip()), line[64:].strip()) for line in open(filepath)]


def test_dummies_keep_their_element(run):
    """The partner's atoms are dummies in the host file, but they keep their element so the cavity keeps its shape."""
    output_directory, _ = run
    elements = {(residue[:3], number, name): element for _, residue, number, name, element, _ in ATOMS}
    for agent in ('methylammonium', 'acetate', 'Ace-Ala', 'Nme'):
        for key, tail in read_cmip_pdb(output_directory / f'{agent}_host.cmip.pdb'):
            if tail.startswith('X') and elements[key] != 'H':
                assert tail == 'X' + elements[key], (agent, key, tail)


@pytest.mark.parametrize('analysis, key', [
    ('energies-00', 'ies'),
    ('energies-00', 'fes'),
    ('energies-01', 'es'),
], ids=['salt bridge frame 1', 'salt bridge frame 3', 'peptide'])
def test_both_directions_agree(run, analysis, key):
    """The electrostatic interaction energy is the same with either agent as the guest, to 1%.

    It needs the same low-dielectric cavity in both directions: all atoms of both agents, with their volume.
    As of ENERGIES_ANALYSIS_VERSION 1.1.0: 11.7% and 22.9% for the salt bridge, 21.6% for the peptide (B1).
    With dummies keeping their element: 0.02% and 0.10%, and 6.3% for the peptide, whose host files lacked the other
    agent's strong bond atom (B8).
    """
    _, outputs = run
    agent1 = sum(value for value in outputs[analysis]['agent1'][key] if value is not None)
    agent2 = sum(value for value in outputs[analysis]['agent2'][key] if value is not None)
    assert abs(agent1 - agent2) <= 0.01 * abs(agent1 + agent2) / 2, (agent1, agent2)


def test_new_grid_contains_both_grids():
    """The grid for both agents contains the automatic grid of each agent (center ± density × units / 2).

    The analysis writes the new density with 0.5 Å between grid points.
    """
    grid_1 = ((0.0, 0.0, 0.0), (80.0, 80.0, 80.0), (0.5, 0.5, 0.5))
    grid_2 = ((25.0, 0.0, 0.0), (20.0, 20.0, 20.0), (0.5, 0.5, 0.5))
    center, density = compute_new_grid(*grid_1, *grid_2)
    low = np.array(center) - np.array(density) * 0.5 / 2
    high = np.array(center) + np.array(density) * 0.5 / 2
    for grid_center, grid_density, grid_units in (grid_1, grid_2):
        half_size = np.array(grid_density) * np.array(grid_units) / 2
        assert np.all(low <= np.array(grid_center) - half_size)
        assert np.all(high >= np.array(grid_center) + half_size)


def test_cap_grid_points_keeps_resolution():
    """A grid above the limit of points keeps its box and is only coarsened as much as needed.

    The capped grid has at most GRID_POINTS_LIMIT points and at least half of them.
    As of ENERGIES_ANALYSIS_VERSION 1.1.0, 501³ points (126 M) were cut to 260³ (18 M), a 0.968 Å grid. Now 402³ (65.0 M),
    a 0.624 Å grid.
    """
    density, unit_size = cap_grid_points([500, 500, 500], 0.5)
    points = np.prod(np.array(density) + 1)
    assert np.all(np.array(density) * unit_size >= 500 * 0.5)
    assert GRID_POINTS_LIMIT / 2 <= points <= GRID_POINTS_LIMIT


def test_clash_warns(clash_run):
    """Zeroing the energies of clashing atoms is reported with a warning naming the interaction."""
    _, logs = clash_run
    assert re.search(r'WARNING.*salt bridge', logs)


def write_golden():
    """Run the analysis on the tiny system with each engine and store its output as the golden output."""
    import tempfile
    for engine, golden_filepath in (('cmip', GOLDEN_FILEPATH), ('gridfree', GRIDFREE_GOLDEN_FILEPATH)):
        with tempfile.TemporaryDirectory() as work_directory:
            output_directory, _ = run_energies(Path(work_directory), FRAMES, INTERACTIONS, engine=engine)
            outputs = load_outputs(output_directory)
        with open(golden_filepath, 'w') as file:
            json.dump(outputs, file, indent=1)
        print(f'Golden output written to {golden_filepath}')


if __name__ == '__main__':
    sys.exit(write_golden())
