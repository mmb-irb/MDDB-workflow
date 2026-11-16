"""T2 and T3 accuracy: the workflow's energies analysis on experimental binding benchmarks.

- T2: PL-REX (10 targets, 164 complexes) and SKEMPI (7 complexes, 376 alanine mutations).
- T3: T2 plus the FEP+ benchmark (37 series, 725 restrained-minimized complexes) and 33 more SKEMPI complexes (40 in all,
  1138 mutations). A T3 run computes only the extra systems; its scores also use the T2 results of the same label.
- T4 (nucleic): protein-nucleic acid alanine scanning (mmCSM-NA, 39 complexes with 6 or more mutations), scored like SKEMPI.
- T5 (sealed hold-out): the SKEMPI alanine-scan complexes not in T2 or T3 (scripts/t5_select.py), scored like SKEMPI.
  It is kept for final decisions: `t5 score` needs `--decision` and records each use in downloads/t5/ledger.jsonl.

Every complex is run through `mddb_workflow.analyses.energies.energies` as one frame (see t1.py), so a change is
measured on the code users run. Inputs are in downloads/t2 and downloads/t3 (ignored by git; see README.md).

Run in mwf_env, from the repository root (hours with CMIP; run it in its own memory-capped unit):
    python -m test.benchmarks.energies t2 run <label> [--jobs 3] [--threads 2] [--engine cmip|gridfree] [--only text]
    python -m test.benchmarks.energies t2 score <label> [--baseline <label>]
and the same with t3, t4 and t5 (`t5 score <label> --decision '<decision and rule>'`). `run` writes per-atom results
to downloads/<tier>/results/<label>/ and skips the systems already done. `score` prints the metrics (and the paired
differences with a baseline) and writes them to downloads/<tier>/scores/<label>.json.
"""
import os
import sys
import json
import glob
import argparse
from datetime import datetime
import multiprocessing
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np

DIRECTORY = Path(__file__).parent
# Side chains are the atoms of a residue other than these (SKEMPI alanine scanning)
BACKBONE = {'N', 'H', 'CA', 'HA', 'HA2', 'HA3', 'C', 'O', 'OXT', 'H1', 'H2', 'H3'}
HOT_SPOT_DDG = 2.0
# T2 guard of the 1.2.0 bug fixes: no metric may drop by more than this
GUARD = 0.10


def inputs_directory(tier: str) -> Path:
    """Get the inputs directory of a tier."""
    return DIRECTORY / 'downloads' / tier


def list_systems(tier: str = 't2') -> list[str]:
    """List the systems a tier runs, as '<benchmark>__<id>'."""
    directory = inputs_directory(tier)
    if tier == 't4':
        return sorted(f'nucleic__{Path(p).stem}' for p in glob.glob(str(directory / 'nucleic' / '*_*_*.npz')))
    if tier == 't5':
        return sorted(f'skempi__{Path(p).stem}' for p in glob.glob(str(directory / 'skempi' / '*_*_*.npz')))
    ligand_benchmark = 'plrex' if tier == 't2' else 'fep'
    complexes = sorted(f'{ligand_benchmark}__{Path(p).parent.parent.name}__{Path(p).parent.name}'
                       for p in glob.glob(str(directory / ligand_benchmark / '*' / '*' / 'atoms.npz')))
    skempi = sorted(f'skempi__{Path(p).stem}' for p in glob.glob(str(directory / 'skempi' / '*_*_*.npz')))
    return complexes + skempi


def residue_indices(chains, resnums, resnames) -> np.ndarray:
    """Get residue indices from the changes of chain, residue number or residue name along the atoms."""
    keys = list(zip(chains, resnums, resnames))
    changes = [0] + [int(keys[i] != keys[i - 1]) for i in range(1, len(keys))]
    return np.cumsum(changes)


def bond_pairs(bonds_json) -> np.ndarray:
    """Get bonded atom pairs from a JSON list of bonded atoms per atom."""
    bonds = json.loads(str(bonds_json))
    return np.array(sorted({(min(a, b), max(a, b)) for a, partners in enumerate(bonds) for b in partners}), dtype=np.int32)


def load_system(system: str, tier: str = 't2') -> tuple[dict, dict]:
    """Load the arrays of a system in the T1 input format, and its interaction."""
    benchmark, *ids = system.split('__')
    directory = inputs_directory(tier)
    if benchmark in ('plrex', 'fep'):
        # Series names may contain '__' too
        group, complex_id = '__'.join(ids[:-1]), ids[-1]
        source = np.load(directory / benchmark / group / complex_id / 'atoms.npz')
        n_receptor, atom_count = int(source['n_rec']), len(source['charges'])
        names = np.array([str(name)[:4] for name in source['names']])
        data = dict(coords=source['coords'], elements=source['elements'], names=names, resnames=source['resnames'],
                    resnums=source['resnums'], chains=source['chains'], icodes=np.full(atom_count, ''),
                    resindex=residue_indices(source['chains'], source['resnums'], source['resnames']),
                    charges=source['charges'], bonds=bond_pairs(source['bonds']),
                    agent1=np.arange(n_receptor), agent2=np.arange(n_receptor, atom_count))
        interaction = dict(name='Receptor-ligand interaction', agent_1='Receptor', agent_2='Ligand', strong_bonds=[])
    else:
        source = np.load(directory / benchmark / f'{ids[0]}.npz')
        data = {key: source[key] for key in ('elements', 'names', 'resnames', 'resnums', 'chains', 'icodes', 'resindex',
                                             'charges', 'agent1', 'agent2')}
        data['coords'] = source['coords'][0]
        data['bonds'] = bond_pairs(source['bonds'])
        chains_1, chains_2 = ids[0].split('_')[1:]
        interaction = dict(name=f'{ids[0]} interaction', agent_1=f'Chains {chains_1}', agent_2=f'Chains {chains_2}',
                           strong_bonds=[])
    return data, interaction


def save_result(output_filepath: Path, **arrays):
    """Save a result in a temporary file and then rename it, so an interrupted run leaves no truncated result."""
    temporary_filepath = output_filepath.with_name(output_filepath.name + '.part')
    with open(temporary_filepath, 'wb') as file:
        np.savez_compressed(file, **arrays)
    os.replace(temporary_filepath, output_filepath)


def run_one(system: str, label: str, threads: int, tier: str, engine: str) -> str:
    """Run one system and save its per-atom energies."""
    from test.benchmarks.energies.t1 import run_workflow_energies
    data, interaction = load_system(system, tier)
    result = run_workflow_energies(data, interaction, threads, label=system, engine=engine)
    warnings = [line for line in result['logs'].splitlines() if 'WARNING' in line]
    output_filepath = inputs_directory(tier) / 'results' / label / f'{system}.npz'
    save_result(output_filepath, agent1=result['agent1'], agent2=result['agent2'], seconds=result['seconds'],
                outside_grid=result['outside_grid'], warnings=np.array(json.dumps(warnings)))
    return f'{system} {result["seconds"]:.1f} s, {result["outside_grid"]} atoms outside the grid'


def run(label: str, jobs: int, threads: int, only: str = '', tier: str = 't2', engine: str = None):
    """Run every system not done yet for this label, or only those whose name contains a text."""
    output_directory = inputs_directory(tier) / 'results' / label
    output_directory.mkdir(parents=True, exist_ok=True)
    todo = [system for system in list_systems(tier)
            if only in system and not (output_directory / f'{system}.npz').exists()]
    print(f'{len(todo)} systems to run with {jobs} jobs of {threads} threads', flush=True)
    with ProcessPoolExecutor(max_workers=jobs, mp_context=multiprocessing.get_context('spawn')) as pool:
        futures = {pool.submit(run_one, system, label, threads, tier, engine): system for system in todo}
        for future in as_completed(futures):
            try:
                print(future.result(), flush=True)
            except Exception as error:
                print(f'{futures[future]} FAILED: {error!r}', flush=True)


def spearman(x, y) -> float:
    """Get the Spearman correlation."""
    from scipy.stats import spearmanr
    return float(spearmanr(x, y)[0])


def auc(labels, scores) -> float:
    """AUC as the Mann-Whitney statistic: the probability that a hot spot outranks a non hot spot."""
    from scipy.stats import rankdata
    labels = np.asarray(labels, bool)
    ranks = rankdata(scores)
    positives, negatives = labels.sum(), (~labels).sum()
    return float((ranks[labels].sum() - positives * (positives + 1) / 2) / (positives * negatives))


def read_experimental_dg(filepath: Path) -> dict:
    """Read a PL-REX experimental dG table."""
    return {line.split()[0]: float(line.split()[1]) for line in open(filepath) if line.strip() and not line.startswith('#')}


def load_results(label: str, tier: str) -> dict:
    """Load the per-atom results of a label in a tier."""
    return {Path(p).stem: np.load(p) for p in glob.glob(str(inputs_directory(tier) / 'results' / label / '*.npz'))}


def score_ligand_series(results: dict, benchmark: str, tier: str) -> dict:
    """Spearman between the ligand total energy (the ligand as the guest) and the experimental dG, per target or series."""
    scores = {}
    for group_directory in sorted((inputs_directory(tier) / benchmark).iterdir()):
        if benchmark == 'plrex':
            experimental = read_experimental_dg(group_directory / 'experimental_dG.txt')
        else:
            ligands = json.load(open(group_directory / 'series.json'))['ligands']
            experimental = {ligand: values['exp'] for ligand, values in ligands.items()}
        pairs = [(float(np.nansum(results[key]['agent2'][:, 2])), dg) for complex_id, dg in experimental.items()
                 if (key := f'{benchmark}__{group_directory.name}__{complex_id}') in results]
        if len(pairs) >= 5:
            x, y = np.array(pairs).T
            scores[group_directory.name] = dict(n=len(pairs), spearman=spearman(x, y))
    return scores


def score_skempi(results: dict, tier: str, benchmark: str = 'skempi') -> dict:
    """Spearman and hot-spot AUC of the side-chain energy of each mutated residue (both agents) against the ddG.

    The AUC is NaN in complexes with no hot spot or only hot spots.
    """
    scores = {}
    for complex_filepath in sorted((inputs_directory(tier) / benchmark).glob('*_*_*.npz')):
        key = f'{benchmark}__{complex_filepath.stem}'
        if key not in results:
            continue
        # Read every array once: an NpzFile reads the whole array again at each access
        with np.load(complex_filepath) as npz:
            source = {key: npz[key] for key in npz.files}
        ddg = json.load(open(str(complex_filepath)[:-4] + '_ddg.json'))
        order = np.concatenate([source['agent1'], source['agent2']])
        residue_keys = np.array([f'{source["chains"][i]}{source["resnums"][i]}{source["icodes"][i]}' for i in order])
        side_chain = ~np.isin(source['names'][order], list(BACKBONE))
        both = np.concatenate([results[key]['agent1'][:, 2], results[key]['agent2'][:, 2]])
        values, experimental = [], []
        for residue, value in ddg.items():
            mask = (residue_keys == residue) & side_chain
            if mask.any():
                values.append(float(np.nansum(both[mask])))
                experimental.append(value)
        values, experimental = -np.array(values), np.array(experimental)
        hot_spots = experimental >= HOT_SPOT_DDG
        scores[complex_filepath.stem] = dict(n=len(values), spearman=spearman(values, experimental),
                                             auc=auc(hot_spots, values) if 0 < hot_spots.sum() < len(values) else np.nan)
    return scores


def other_measures(results: dict) -> dict:
    """Two-direction mismatch, zeroed atoms, atoms outside the grid and time."""
    mismatches, zeroed, outside, seconds = [], 0, 0, []
    for result in results.values():
        es_1, es_2 = np.nansum(result['agent1'][:, 1]), np.nansum(result['agent2'][:, 1])
        mismatches.append(abs(es_1 - es_2) / max(abs(es_1 + es_2) / 2, 1.0))
        zeroed += sum(int(np.all(result[agent] == 0, axis=1).sum()) for agent in ('agent1', 'agent2'))
        outside += int(result['outside_grid'])
        seconds.append(float(result['seconds']))
    return dict(median_mismatch=float(np.median(mismatches)), mismatch_90th=float(np.percentile(mismatches, 90)),
                zeroed_atoms=zeroed, atoms_outside_grid=outside, hours=sum(seconds) / 3600,
                median_seconds=float(np.median(seconds)))


def mean_of(scores: dict, key: str) -> float:
    """Mean of a score over targets, series or complexes, leaving out those where it is not defined."""
    return float(np.nanmean([value[key] for value in scores.values()])) if scores else float('nan')


def score(label: str, tier: str = 't2') -> dict:
    """Score the per-atom results of a label."""
    if tier == 't4':
        results = load_results(label, 't4')
        report = {'label': label, 'tier': tier, 'nucleic': score_skempi(results, 't4', 'nucleic'), 'systems': len(results)}
        report['summary'] = dict(nucleic_spearman=mean_of(report['nucleic'], 'spearman'),
                                 nucleic_auc=mean_of(report['nucleic'], 'auc'), **other_measures(results))
        return report
    if tier == 't5':
        results = load_results(label, 't5')
        report = {'label': label, 'tier': tier, 'skempi': score_skempi(results, 't5'), 'systems': len(results)}
        report['summary'] = dict(skempi_spearman=mean_of(report['skempi'], 'spearman'),
                                 skempi_auc=mean_of(report['skempi'], 'auc'), **other_measures(results))
        return report
    t2_results = load_results(label, 't2')
    report = {'label': label, 'tier': tier, 'plrex': score_ligand_series(t2_results, 'plrex', 't2'),
              'skempi': score_skempi(t2_results, 't2')}
    results = dict(t2_results)
    if tier == 't3':
        t3_results = load_results(label, 't3')
        results.update(t3_results)
        report['fep'] = score_ligand_series(t3_results, 'fep', 't3')
        report['skempi'] = {**report['skempi'], **score_skempi(t3_results, 't3')}
    report['systems'] = len(results)
    report['summary'] = dict(plrex_spearman=mean_of(report['plrex'], 'spearman'),
                             skempi_spearman=mean_of(report['skempi'], 'spearman'),
                             skempi_auc=mean_of(report['skempi'], 'auc'), **other_measures(results))
    if tier == 't3':
        report['summary']['fep_spearman'] = mean_of(report['fep'], 'spearman')
    return report


def bootstrap_mean(values, resamples: int = 10000) -> tuple[float, float, float]:
    """Mean and 95% bootstrap interval."""
    values = np.asarray(values, float)
    rng = np.random.default_rng(0)
    means = values[rng.integers(0, len(values), (resamples, len(values)))].mean(1)
    return float(values.mean()), float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def paired(report: dict, baseline: dict) -> dict:
    """Paired differences (report minus baseline) per metric, with bootstrap intervals and Wilcoxon p-values."""
    from scipy.stats import wilcoxon
    differences = {}
    metrics = [('plrex_spearman', 'plrex', 'spearman'), ('skempi_spearman', 'skempi', 'spearman'),
               ('skempi_auc', 'skempi', 'auc'), ('nucleic_spearman', 'nucleic', 'spearman'),
               ('nucleic_auc', 'nucleic', 'auc')]
    if 'fep' in report and 'fep' in baseline:
        metrics.insert(0, ('fep_spearman', 'fep', 'spearman'))
    for metric, benchmark, key in metrics:
        if benchmark not in report or benchmark not in baseline:
            continue
        common = sorted(k for k in set(report[benchmark]) & set(baseline[benchmark])
                        if np.isfinite(report[benchmark][k][key]) and np.isfinite(baseline[benchmark][k][key]))
        if not common:
            continue
        values = [report[benchmark][k][key] - baseline[benchmark][k][key] for k in common]
        mean, low, high = bootstrap_mean(values)
        p_value = float(wilcoxon(values).pvalue) if np.any(np.array(values) != 0) else 1.0
        differences[metric] = dict(difference=mean, low=low, high=high, p=p_value, n=len(values),
                                   better=int(sum(v > 0 for v in values)), worse=int(sum(v < 0 for v in values)),
                                   guard_passed=bool(mean >= -GUARD))
    return differences


def record_t5_use(decision: str, label: str, baseline: str | None):
    """Append a use of T5 to its ledger, before scoring, so every look at the hold-out is on record."""
    ledger_filepath = inputs_directory('t5') / 'ledger.jsonl'
    ledger_filepath.parent.mkdir(parents=True, exist_ok=True)
    entry = dict(time=datetime.now().isoformat(timespec='seconds'), decision=decision, label=label, baseline=baseline)
    with open(ledger_filepath, 'a') as file:
        file.write(json.dumps(entry) + '\n')
    print(f'T5 use recorded in {ledger_filepath}')


def main(arguments: list[str], tier: str = 't2') -> int:
    """Run or score a tier."""
    from mddb_workflow.analyses.energies import ENERGIES_ENGINES
    parser = argparse.ArgumentParser(prog=f'python -m test.benchmarks.energies {tier}', description=__doc__.split('\n')[0])
    subparsers = parser.add_subparsers(dest='command', required=True)
    run_parser = subparsers.add_parser('run')
    run_parser.add_argument('label')
    run_parser.add_argument('--jobs', type=int, default=3)
    run_parser.add_argument('--threads', type=int, default=2)
    run_parser.add_argument('--only', default='', help='run only the systems whose name contains this text')
    run_parser.add_argument('--engine', default=None, choices=ENERGIES_ENGINES,
                            help='energies engine (default: the analysis default)')
    score_parser = subparsers.add_parser('score')
    score_parser.add_argument('label')
    score_parser.add_argument('--baseline')
    if tier == 't5':
        score_parser.add_argument('--decision', required=True,
                                  help='the decision this score settles, with its rule (T5 is a sealed hold-out)')
    args = parser.parse_args(arguments)
    if args.command == 'run':
        run(args.label, args.jobs, args.threads, args.only, tier, args.engine)
        return 0
    if tier == 't5':
        record_t5_use(args.decision, args.label, args.baseline)
    report = score(args.label, tier)
    if args.baseline:
        report['baseline'] = args.baseline
        report['paired'] = paired(report, score(args.baseline, tier))
    scores_directory = inputs_directory(tier) / 'scores'
    scores_directory.mkdir(parents=True, exist_ok=True)
    with open(scores_directory / f'{args.label}.json', 'w') as file:
        json.dump(report, file, indent=1)
    s = report['summary']
    print(f'{tier.upper()} {args.label}: {report["systems"]} systems, {s["hours"]:.2f} h, '
          f'median {s["median_seconds"]:.2f} s per system')
    columns = (['nucleic_spearman', 'nucleic_auc'] if tier == 't4' else ['skempi_spearman', 'skempi_auc'] if tier == 't5'
               else (['fep_spearman'] if tier == 't3' else []) + ['plrex_spearman', 'skempi_spearman', 'skempi_auc'])
    print('| ' + ' | '.join(columns) + ' | median mismatch (90th) | zeroed atoms | atoms outside grid |')
    print('|' + '---|' * (len(columns) + 3))
    print('| ' + ' | '.join(f'{s[c]:.3f}' for c in columns) + f' | {100 * s["median_mismatch"]:.1f}% '
          f'({100 * s["mismatch_90th"]:.1f}%) | {s["zeroed_atoms"]} | {s["atoms_outside_grid"]} |')
    if args.baseline:
        print(f'\nPaired difference with {args.baseline} [95% CI], Wilcoxon p:')
        for metric, d in report['paired'].items():
            print(f'  {metric:16s} {d["difference"]:+.3f} [{d["low"]:+.3f}, {d["high"]:+.3f}]  p={d["p"]:.3f}  '
                  f'better {d["better"]}/{d["n"]}, worse {d["worse"]}')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
