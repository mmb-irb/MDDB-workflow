"""Benchmark of the recognition of coarse grained lipids by bead names (see README.md).

The corpus (corpus.json, made by build_corpus.py) holds the residues of the parameter sets of Martini and other coarse
grained force fields, labelled as lipids or not. The residues which the rule gets wrong must be the recorded ones, so
both regressions and improvements make the test fail and the record is kept up to date.

Run this file to print the results by source as a markdown table:
    python test/benchmarks/cg_lipid_naming/test_cg_lipid_naming.py
"""
import json
import os
from collections import OrderedDict
from mddb_workflow.tools.cg_lipids import coarse_grain_lipid_tails

CORPUS = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'corpus.json')
# Lipids which are not recognised by their bead names, by class, so they must be found by residue name or forced
MISSED_LIPIDS = {
    # SPICA and SIRAH name their beads in other ways
    'SPICA lipid': 'ASM BMGP BSM CER1 CER2 CHL1 CLOL CLR DAPC DLIP DLIPC DLPC DMPC DMPG DOPC DOPE DOPG DOPS DPPC DSPC LSM '
        'POPC POPE POPG POPI POPS PSM SPM2 SSM',
    'SIRAH lipid': 'CMM CPO CPP EPO SPO',
    # No ROH head bead, and many small molecules have R1, R2... beads too
    'hopanoid': 'HBHT HDPT HHOP HOPR',
    # Only one or two tail beads
    'lyso-PC': 'CPC TPC',
    'detergent': 'BOG',
    # Martini 2 single residue LPS with V11, V12... beads
    'LPS': 'LPS',
}
# Molecules which are not lipids but are recognised as lipids by their bead names, by class
FALSE_POSITIVES = {
    # Two tailed voltage-sensitive dyes, which are membrane-embedded amphiphiles anyway
    'membrane dye': 'ALT5 DYE4 DYE5 DYE7 DYE8',
    # Parameterisation fragments of ionizable lipids, never in production systems
    'IL fragment': 'ALC Ester',
    # One acyl residue of a multi-residue lipid A molecule
    'lipid A / LPS residue': 'LPB',
}


def load_corpus() -> list[dict]:
    with open(CORPUS) as handle:
        return json.load(handle)


def is_recognised(case: dict) -> bool:
    return coarse_grain_lipid_tails(case['beads'].split()) is not None


def wrong_cases(cases: list[dict], lipids: bool) -> dict[str, str]:
    """Get the residue names of the lipids which are not recognised (or the non-lipids which are) by class."""
    wrong = {}
    for case in cases:
        if (case['expected'] is not None) != lipids or is_recognised(case) == lipids: continue
        wrong.setdefault(case['class'], set()).add(case['resname'])
    return {label: ' '.join(sorted(resnames)) for label, resnames in wrong.items()}


def results_by_group(cases: list[dict]) -> OrderedDict:
    """Get the number of lipids, recognised lipids, non-lipids and false positives of every source group and in total."""
    results = {}
    for case in cases:
        lipid, recognised = case['expected'] is not None, is_recognised(case)
        for group in case['groups'] + ['Total']:
            counts = results.setdefault(group, [0, 0, 0, 0])
            counts[0 if lipid else 2] += 1
            counts[1 if lipid else 3] += recognised
    return OrderedDict(sorted(results.items(), key=lambda item: (item[0] == 'Total', item[0].lower())))


def test_corpus_size():
    cases = load_corpus()
    lipids = [case for case in cases if case['expected'] is not None]
    assert (len(lipids), len(cases) - len(lipids)) == (1305, 665)


def test_missed_lipids():
    assert wrong_cases(load_corpus(), lipids=True) == MISSED_LIPIDS


def test_false_positives():
    assert wrong_cases(load_corpus(), lipids=False) == FALSE_POSITIVES


if __name__ == '__main__':
    print('| source | lipids | recognised | non-lipids | false positives |')
    print('|---|---|---|---|---|')
    for group, (n_lipids, recognised, n_others, false_positives) in results_by_group(load_corpus()).items():
        print(f'| {group} | {n_lipids} | {recognised} | {n_others} | {false_positives} |')
