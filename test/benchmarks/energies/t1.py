"""T1 regression: the workflow's energies analysis on one frame of each T1 system, compared with an earlier run.

Each system is one frame run through `mddb_workflow.analyses.energies.energies` as the workflow runs it, so refactors
and fixes are measured on the code users run. The frames are fetched and built the first time they are needed
(t1_inputs.py, downloads/t1/<system>.npz). A run's per-atom output (vdw, es and both for each agent; NaN where the
analysis gives None) is stored under a label in downloads/t1/results/<label>/. Git ignores downloads/.

Run in mwf_env, from the repository root:
    python -m test.benchmarks.energies t1 --label <label>       # before a change: store the results under a label
    python -m test.benchmarks.energies t1 --baseline <label>    # after it: compare with the label (atol 1e-6)
    python -m test.benchmarks.energies t1 --engine gridfree ... # the same with the grid-free engine (CMIP by default)
The run uses one OpenMP thread (--threads) so CMIP results are reproducible.
"""
import os
import io
import re
import sys
import json
import time
import shutil
import argparse
import tempfile
import contextlib
from pathlib import Path
import numpy as np
from mddb_workflow.analyses.energies import DEFAULT_ENERGIES_ENGINE, ENERGIES_ENGINES
from test.benchmarks.energies.t1_inputs import get_input_filepath

LABELS_DIRECTORY = Path(__file__).parent / 'downloads' / 't1' / 'results'
SYSTEMS = ['A01TW', 'A0221', '1VFB']
ENERGY_KEYS = ['vdw', 'es', 'both']
ATOL = 1e-6


def build_structure(data):
    """Build a workflow structure, with bonds, from the arrays of a T1 input file."""
    from mddb_workflow.utils.structures import Structure, Atom, Residue, Chain
    atoms, residues, chains = [], [], {}
    residue_by_index = {}
    for index in range(len(data['names'])):
        atoms.append(Atom(name=str(data['names'][index]), element=str(data['elements'][index]),
                          coords=tuple(float(c) for c in data['coords'][index])))
        residue_index = int(data['resindex'][index])
        if residue_index not in residue_by_index:
            residue = Residue(name=str(data['resnames'][index]), number=int(data['resnums'][index]),
                              icode=str(data['icodes'][index]))
            residue.atom_indices = []
            residue_by_index[residue_index] = len(residues)
            residues.append(residue)
            chain_name = str(data['chains'][index])
            if chain_name not in chains:
                chains[chain_name] = Chain(name=chain_name)
                chains[chain_name].residue_indices = []
            chains[chain_name].residue_indices = chains[chain_name].residue_indices + [len(residues) - 1]
        residue = residues[residue_by_index[residue_index]]
        residue.atom_indices = residue.atom_indices + [index]
    structure = Structure(atoms=atoms, residues=residues, chains=list(chains.values()))
    bonds = [[] for _ in atoms]
    for a, b in data['bonds']:
        bonds[int(a)].append(int(b))
        bonds[int(b)].append(int(a))
    structure.bonds = bonds
    return structure


def read_grid_box(cmip_inputs_filepath: str) -> tuple[np.ndarray, np.ndarray]:
    """Get the lower and upper corners (Å) of the grid written in a CMIP inputs file."""
    values = dict(re.findall(r'(\w+)=([-\d.]+)', open(cmip_inputs_filepath).read()))
    center = np.array([float(values[f'cen{axis}']) for axis in 'xyz'])
    size = np.array([float(values[f'dim{axis}']) * float(values[f'int{axis}']) for axis in 'xyz'])
    return center - size / 2, center + size / 2


def count_atoms_outside(cmip_pdb_filepath: str, low: np.ndarray, high: np.ndarray) -> int:
    """Count the atoms of a CMIP input PDB which are outside a grid box."""
    coords = np.array([[float(line[30:38]), float(line[38:46]), float(line[46:54])] for line in open(cmip_pdb_filepath)])
    return int((~np.all((coords >= low) & (coords <= high), axis=1)).sum())


def run_system(system: str, threads: int, keep: bool = False, engine: str = None) -> dict:
    """Run the energies analysis on one T1 system and return its per-atom energies and measures."""
    # Read every array once: an NpzFile decompresses the whole array again at each access
    with np.load(get_input_filepath(system)) as source:
        data = {key: source[key] for key in source.files}
    interaction = json.loads(str(data['interaction']))
    return run_workflow_energies(data, interaction, threads, keep, system, engine)


def run_workflow_energies(data, interaction: dict, threads: int, keep: bool = False, label: str = '',
                          engine: str = None) -> dict:
    """Run the workflow's energies analysis on one frame and return its per-atom energies and measures.

    Data has the arrays of a T1 input file (coords, atoms, residues, chains, charges, bonds, agent1, agent2).
    """
    import mdtraj
    from mddb_workflow.analyses.energies import energies
    from mddb_workflow.utils.constants import GLOBALS
    from mddb_workflow.utils.selections import Selection
    from mddb_workflow.utils.file import File
    interaction = dict(interaction)
    interaction['atom_indices_1'] = [int(i) for i in data['agent1']]
    interaction['atom_indices_2'] = [int(i) for i in data['agent2']]
    # CMIP truncates file paths longer than 150 characters, so the work directory must have a short path
    work_directory = tempfile.mkdtemp(prefix='mwf_t1_')
    trajectory_filepath = f'{work_directory}/trajectory.dcd'
    with mdtraj.formats.DCDTrajectoryFile(trajectory_filepath, 'w') as writer:
        writer.write(data['coords'][np.newaxis].astype(np.float32))
    output_directory = f'{work_directory}/energies'
    previous_directory, previous_working_directory = os.getcwd(), GLOBALS['working_directory']
    previous_threads = os.environ.get('OMP_NUM_THREADS')
    logs = io.StringIO()
    # The workflow gets the structure and charges ready, so building them is not part of the analysis' time
    structure = build_structure(data)
    charges = [float(charge) for charge in data['charges']]
    try:
        os.chdir(work_directory)
        os.environ['OMP_NUM_THREADS'] = str(threads)
        GLOBALS['working_directory'] = work_directory
        start = time.perf_counter()
        with contextlib.redirect_stdout(logs):
            energies(
                trajectory_file=File(trajectory_filepath),
                output_directory=output_directory,
                structure=structure,
                interactions=[interaction],
                charges=charges,
                cg_selection=Selection([]),
                dummy_selection=Selection([]),
                snapshots=1,
                **({'energies_engine': engine} if engine else {}),
            )
        seconds = time.perf_counter() - start
    finally:
        os.chdir(previous_directory)
        GLOBALS['working_directory'] = previous_working_directory
        if previous_threads is None:
            os.environ.pop('OMP_NUM_THREADS', None)
        else:
            os.environ['OMP_NUM_THREADS'] = previous_threads
    output = json.load(open(f'{output_directory}/mda.energies_00.json'))
    result = {'seconds': seconds, 'logs': logs.getvalue()}
    # Atoms set to 0 for being above 100 kcal/mol (clashes), as counted by the analysis since version 1.2.0
    result['clashes'] = (sum(output[agent]['zeroed'][0] for agent in ('agent1', 'agent2'))
                         if 'zeroed' in output['agent1'] else None)
    for agent in ('agent1', 'agent2'):
        result[agent] = np.array([[np.nan if value is None else value for value in output[agent][key]]
                                  for key in ENERGY_KEYS]).T
    # Atoms outside the CMIP grid of the last interaction (B2), from the files the analysis leaves behind
    if not os.path.exists(f'{output_directory}/cmip.in'):
        result['outside_grid'] = 0
    else:
        low, high = read_grid_box(f'{output_directory}/cmip.in')
        result['outside_grid'] = sum(
            count_atoms_outside(f'{output_directory}/{name.replace(" ", "_").replace("/", "_")}_guest.cmip.pdb', low, high)
            for name in (interaction['agent_1'], interaction['agent_2']))
    if keep:
        print(f'  {label} work directory kept: {work_directory}')
    else:
        shutil.rmtree(work_directory)
    return result


def measures(result: dict) -> dict:
    """Totals per agent, the two-direction electrostatic mismatch, the atoms with all energies 0 and the clashes."""
    totals = {agent: np.nansum(result[agent], axis=0) for agent in ('agent1', 'agent2')}
    es_1, es_2 = totals['agent1'][1], totals['agent2'][1]
    zeroed = sum(int(np.all(result[agent] == 0, axis=1).sum()) for agent in ('agent1', 'agent2'))
    return dict(es_1=es_1, es_2=es_2, vdw_1=totals['agent1'][0], vdw_2=totals['agent2'][0],
                mismatch=abs(es_1 - es_2) / (abs(es_1 + es_2) / 2), zeroed=zeroed, clashes=result['clashes'],
                outside_grid=result['outside_grid'],
                seconds=result['seconds'])


def repeat(systems: list[str], threads: int, runs: int, engine: str) -> int:
    """Run each system several times and report the largest differences between runs (determinism, bug B7)."""
    print('| system | runs | threads | max per-atom difference between runs (kcal/mol) | '
          'max agent total difference (kcal/mol) | time (s) |')
    print('|---|---|---|---|---|---|')
    for system in systems:
        results = [run_system(system, threads, engine=engine) for _ in range(runs)]
        per_atom = total = 0.0
        for agent in ('agent1', 'agent2'):
            values = np.array([result[agent] for result in results])
            per_atom = max(per_atom, float(np.nanmax(np.nanmax(values, 0) - np.nanmin(values, 0))))
            totals = np.nansum(values, axis=1)
            total = max(total, float(np.max(totals.max(0) - totals.min(0))))
        seconds = np.mean([result['seconds'] for result in results])
        print(f'| {system} | {runs} | {threads} | {per_atom:.2g} | {total:.2g} | {seconds:.1f} |', flush=True)
    return 0


def compare(result: dict, reference) -> tuple[bool, float]:
    """Check whether per-atom results equal an earlier run (atol 1e-6, None in the same places) and get the largest difference."""
    differences = [np.nanmax(np.abs(result[agent] - reference[agent]), initial=0) for agent in ('agent1', 'agent2')]
    same_nones = all(np.array_equal(np.isnan(result[agent]), np.isnan(reference[agent])) for agent in ('agent1', 'agent2'))
    same_shape = all(result[agent].shape == reference[agent].shape for agent in ('agent1', 'agent2'))
    return same_shape and same_nones and max(differences) <= ATOL, float(max(differences))


def main(arguments: list[str]) -> int:
    """Run T1, and store the per-atom results under a label or compare them with a labelled run."""
    parser = argparse.ArgumentParser(prog='python -m test.benchmarks.energies t1', description=__doc__.split('\n')[0])
    parser.add_argument('--systems', nargs='+', default=SYSTEMS, choices=SYSTEMS)
    parser.add_argument('--threads', type=int, default=1, help='OpenMP threads for CMIP (default 1)')
    parser.add_argument('--keep', action='store_true', help='keep the work directories, for debugging')
    parser.add_argument('--engine', default=DEFAULT_ENERGIES_ENGINE, choices=ENERGIES_ENGINES,
                        help=f'energies engine (default {DEFAULT_ENERGIES_ENGINE})')
    parser.add_argument('--repeat', type=int, default=1,
                        help='run each system this many times and report how much the results vary between runs')
    parser.add_argument('--label', help='store the per-atom results under this label (downloads/t1/results/<label>/)')
    parser.add_argument('--baseline', help='compare with the per-atom results of this label (atol 1e-6)')
    args = parser.parse_args(arguments)
    if args.repeat > 1:
        return repeat(args.systems, args.threads, args.repeat, args.engine)
    failures = []
    rows = []
    for system in args.systems:
        result = run_system(system, args.threads, args.keep, args.engine)
        if args.label:
            (LABELS_DIRECTORY / args.label).mkdir(parents=True, exist_ok=True)
            np.savez_compressed(LABELS_DIRECTORY / args.label / f'{system}.npz', agent1=result['agent1'],
                                agent2=result['agent2'], seconds=result['seconds'], threads=args.threads,
                                engine=args.engine)
        status = ''
        if args.baseline:
            identical, difference = compare(result, np.load(LABELS_DIRECTORY / args.baseline / f'{system}.npz'))
            status = (f'identical to {args.baseline} (max diff {difference:.1e})' if identical
                      else f'DIFFERENT from {args.baseline} (max diff {difference:.3g})')
            if not identical:
                failures.append(system)
        elif args.label:
            status = f'stored as {args.label}'
        rows.append((system, measures(result), status))
        sys.stdout.flush()
    print(f'\nT1, {args.engine}, {args.threads} OpenMP thread(s)')
    print('| system | es agent 1 | es agent 2 | mismatch | vdw agent 1 | vdw agent 2 | atoms with all energies 0 | '
          'clash atoms set to 0 | atoms outside grid | time (s) | comparison |')
    print('|---|---|---|---|---|---|---|---|---|---|---|')
    for system, m, status in rows:
        print(f'| {system} | {m["es_1"]:.2f} | {m["es_2"]:.2f} | {100 * m["mismatch"]:.1f}% | {m["vdw_1"]:.2f} | '
              f'{m["vdw_2"]:.2f} | {m["zeroed"]} | {"-" if m["clashes"] is None else m["clashes"]} | {m["outside_grid"]} | '
              f'{m["seconds"]:.3g} | {status} |')
    return 1 if failures else 0
