"""Select and check the SKEMPI complexes of T5, the sealed hold-out tier (2026-11-13).

Rule, fixed before any structure was prepared or scored: every SKEMPI 2.0 complex not used in T2 or T3 with at least
MIN_MUTATIONS single X -> Ala mutations (X not Ala or Gly, parsed as in skempi_prepare.py), and, once prepared with
skempi_prepare.py, at least MIN_MUTATIONS of those residues present with side-chain atoms in the structure. Each complex
counts once in the mean per-complex Spearman, whatever its number of mutations.

Usage, from the repository root:
    python test/benchmarks/energies/scripts/t5_select.py select <skempi_v2.csv> <output dir>
    (then skempi_prepare.py on each selected complex, into the same output dir)
    python test/benchmarks/energies/scripts/t5_select.py check <output dir>
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

MIN_MUTATIONS = 6
R = 1.9872e-3
# The SKEMPI complexes of T2 (7) and T3 (33 more)
USED = set('''1A22_A_B 1EMV_A_B 1JRH_LH_I 1JTD_A_B 1JTG_A_B 1VFB_AB_C 3S9D_A_B 1A4Y_A_B 1AO7_ABC_DE 1BP3_A_B 1BRS_A_D
1DAN_HL_UT 1DQJ_AB_C 1DVF_AB_CD 1F47_A_B 1K8R_A_B 1KTZ_A_B 1LFD_A_B 1OGA_ABC_DE 1Z7X_W_X 2G2U_A_B 2NYY_DC_A 2WPT_A_B
3BN9_B_CD 3C60_CD_AB 3EQS_A_B 3EQY_A_C 3HFM_HL_Y 3M62_A_B 3M63_A_B 3NGB_HL_G 3QHY_A_B 3QIB_ABP_CD 3SZK_AB_C 4G0N_A_B
4OFY_A_D 4P23_CD_AB 4RS1_A_B 5C6T_HL_A 5M2O_A_B'''.split())
BACKBONE = {'N', 'CA', 'C', 'O', 'H', 'HA', 'OXT', 'H1', 'H2', 'H3'}


def alanine_mutations(skempi_filepath: str) -> dict[str, set[str]]:
    """Get the residues with a single X -> Ala mutation and a usable affinity, per SKEMPI complex."""
    rows = list(csv.reader(open(skempi_filepath), delimiter=';'))[1:]
    mutations = {}
    for row in rows:
        if ',' in row[1] or not row[1].endswith('A') or row[1][0] in 'AG':
            continue
        try:
            mutant, wild_type = float(row[7]), float(row[9])
            float(row[13].replace('(assumed)', '').strip() or 298)
        except ValueError:
            continue
        if mutant > 0 and wild_type > 0 and math.isfinite(R * math.log(mutant / wild_type)):
            mutations.setdefault(row[0], set()).add(row[1][1:-1])
    return mutations


def select(skempi_filepath: str, output_directory: str):
    """Write the candidate complexes, before any preparation."""
    mutations = alanine_mutations(skempi_filepath)
    candidates = {complex_id: len(residues) for complex_id, residues in sorted(mutations.items())
                  if complex_id not in USED and len(residues) >= MIN_MUTATIONS}
    Path(output_directory).mkdir(parents=True, exist_ok=True)
    with open(Path(output_directory) / 'candidates.json', 'w') as file:
        json.dump(dict(rule=__doc__.split('\n\n')[1], candidates=candidates), file, indent=1)
    print(f'{len(candidates)} candidates, {sum(candidates.values())} mutations')
    print(' '.join(candidates))


def check(output_directory: str):
    """Keep the prepared complexes with enough mutated residues in their structure, and write the selection."""
    directory = Path(output_directory)
    candidates = json.load(open(directory / 'candidates.json'))['candidates']
    selected, excluded = {}, {}
    for complex_id in candidates:
        structure_filepath = directory / f'{complex_id}.npz'
        if not structure_filepath.exists():
            excluded[complex_id] = 'preparation failed'
            continue
        with np.load(structure_filepath) as source:
            names, chains, resnums, icodes = (source[key] for key in ('names', 'chains', 'resnums', 'icodes'))
        side_chain = ~np.isin(names, list(BACKBONE))
        present = {f'{c}{r}{i}' for c, r, i, s in zip(chains, resnums, icodes, side_chain) if s}
        ddg = json.load(open(directory / f'{complex_id}_ddg.json'))
        matched = sum(1 for residue in ddg if residue in present)
        if matched >= MIN_MUTATIONS:
            selected[complex_id] = matched
        else:
            excluded[complex_id] = f'{matched} mutated residues in the structure'
    with open(directory / 'selection.json', 'w') as file:
        json.dump(dict(selected=selected, excluded=excluded), file, indent=1)
    for complex_id in excluded:
        for suffix in ('.npz', '.pdb', '_ddg.json'):
            (directory / f'{complex_id}{suffix}').unlink(missing_ok=True)
    print(f'{len(selected)} complexes selected ({sum(selected.values())} mutations), {len(excluded)} excluded')
    for complex_id, reason in excluded.items():
        print(f'  {complex_id}: {reason}')


if __name__ == '__main__':
    command, *arguments = sys.argv[1:]
    {'select': select, 'check': check}[command](*arguments)
