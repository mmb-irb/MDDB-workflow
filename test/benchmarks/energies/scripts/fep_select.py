"""Select the FEP+ benchmark series and ligands to score, before any energy is computed.

Rules (fixed up front): a series needs >= 10 ligands with an experimental dG and a structure, and an
experimental range >= 2 kcal/mol (below that, ranks are dominated by assay noise). Series with more than
25 ligands are randomly subsampled to 25 (seed 0) to bound the CMIP cost. A ligand's structure is the SDF
entry whose title equals the ligand name (the first alternate rotamer/protomer if there is no exact one).

Run: python fep_select.py <fepplus repo> <output dir>
"""
import os
import sys
import csv
import json
import glob
import random

repo, out = sys.argv[1:3]
inputs = f'{repo}/fep_benchmark_inputs/structure_inputs'
preds = f'{repo}/21_4_results/ligand_predictions'
MIN_LIGANDS, MIN_RANGE, CAP = 10, 2.0, 25
# Prediction files whose names do not follow the input naming scheme
RENAMED = {'p2y1_meta_sub': 'py21_merged.csv', 'p2y1_ortho_sub': 'py21_merged.csv',
           'bace_ciordia_retro': 'retrospective_custom_core_extra_pkacorr_out.csv',
           'syk_4puz_fullmap': 'syk_pkacorr_out.csv', 'frag_mcl1_noweak': 'frag_mcl1_manual_uvt_out.csv'}


def split_sdf(path):
    entries, block = {}, []
    for line in open(path):
        block.append(line)
        if line.startswith('$$$$'):
            title = block[0].strip()
            entries.setdefault(title, ''.join(block))
            block = []
    return entries


selected, skipped = [], []
for row in csv.DictReader(open(f'{inputs}/../benchmark_metadata.csv')):
    group, scheme = row['Group abbreviation'], row['Input file naming scheme']
    name = f'{group}__{scheme}'
    sdf, pdb = f'{inputs}/{group}/{scheme}_ligands.sdf', f'{inputs}/{group}/{scheme}_protein.pdb'
    files = sorted(glob.glob(f'{preds}/{group}/*.csv'))
    exact = [f for f in files if os.path.basename(f) == f'{scheme}_out.csv']
    prefix = [f for f in files if os.path.basename(f).startswith(scheme + '_')]
    renamed = [f for f in files if os.path.basename(f) == RENAMED.get(scheme)]
    pred_file = (exact or prefix or renamed or [None])[0]
    if pred_file is None or not os.path.exists(sdf) or not os.path.exists(pdb):
        skipped.append((name, 'no prediction file' if pred_file is None else 'no structure'))
        continue
    entries = split_sdf(sdf)
    ligands = {}
    for r in csv.DictReader(open(pred_file)):
        lig = (r.get('Ligand name') or r.get('Ligand') or '').strip()
        try:
            dg, fep = float(r['Exp. dG (kcal/mol)']), float(r['Pred. dG (kcal/mol)'])
        except (ValueError, KeyError):
            continue
        block = entries.get(lig) or next((v for k, v in entries.items() if k.split()[0] == lig), None)
        if block:
            ligands[lig] = dict(exp=dg, fep=fep, sdf=block)
    span = max(v['exp'] for v in ligands.values()) - min(v['exp'] for v in ligands.values()) if ligands else 0
    if len(ligands) < MIN_LIGANDS or span < MIN_RANGE:
        skipped.append((name, f'{len(ligands)} ligands, range {span:.1f}'))
        continue
    keys = sorted(ligands)
    if len(keys) > CAP:
        keys = sorted(random.Random(0).sample(keys, CAP))
    sdir = f'{out}/{name}'
    os.makedirs(sdir, exist_ok=True)
    table = {}
    for i, lig in enumerate(keys):
        ldir = f'{sdir}/lig{i:02d}'
        os.makedirs(ldir, exist_ok=True)
        open(f'{ldir}/lig.sdf', 'w').write(ligands[lig]['sdf'])
        table[f'lig{i:02d}'] = dict(name=lig, exp=ligands[lig]['exp'], fep=ligands[lig]['fep'])
    json.dump(dict(protein=os.path.abspath(pdb), prediction_file=os.path.abspath(pred_file),
                   n_available=len(ligands), range=span, ligands=table), open(f'{sdir}/series.json', 'w'), indent=1)
    selected.append((name, len(ligands), len(keys), round(span, 1)))

print(f'selected {len(selected)} series, {sum(s[2] for s in selected)} complexes')
for s in selected:
    print('  ', *s)
print(f'skipped {len(skipped)}')
for s in skipped:
    print('  ', *s)
