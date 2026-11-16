"""Prepare SKEMPI complexes (crystal structure -> Amber ff14SB with hydrogens) and the experimental
alanine-scanning ddG table.

Run in openmm_env: python skempi_prepare.py <skempi_v2.csv> <output dir> <PDB_chains1_chains2> ...
Writes <id>.npz (coords, charges, LJ, OBC2 radii/scales, elements, residue info, bonds, agents),
<id>.pdb and <id>_ddg.json (WT residue -> mean experimental ddG of the mutation to Ala, kcal/mol).
"""
import sys
import os
import io
import csv
import json
import math
import urllib.request
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as u
from pdbfixer import PDBFixer

R = 1.9872e-3
csv_file, outdir = sys.argv[1:3]
os.makedirs(outdir, exist_ok=True)
rows = list(csv.reader(open(csv_file), delimiter=';'))[1:]

for complex_id in sys.argv[3:]:
    pdb, c1, c2 = complex_id.split('_')
    # Experimental ddG for single X -> Ala mutations (X not Ala/Gly), averaged over repeated entries
    ddg = {}
    for r in rows:
        if r[0] != complex_id or ',' in r[1] or not r[1].endswith('A') or r[1][0] in 'AG':
            continue
        try:
            km, kw = float(r[7]), float(r[9])
            T = float(r[13].replace('(assumed)', '').strip() or 298)
        except ValueError:
            continue
        if km > 0 and kw > 0:
            key = r[1][1:-1]  # chain + residue number (+ insertion code)
            ddg.setdefault(key, []).append(R * T * math.log(km / kw))
    ddg = {k: float(np.mean(v)) for k, v in ddg.items()}
    json.dump(ddg, open(f'{outdir}/{complex_id}_ddg.json', 'w'), indent=1)

    text = urllib.request.urlopen(f'https://files.rcsb.org/download/{pdb}.pdb', timeout=60).read().decode()
    fixer = PDBFixer(pdbfile=io.StringIO(text))
    keep = set(c1) | set(c2)
    fixer.removeChains([i for i, c in enumerate(fixer.topology.chains()) if c.id not in keep])
    fixer.findMissingResidues()
    fixer.missingResidues = {}  # do not model missing loops
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    ff = app.ForceField('amber14-all.xml', 'implicit/obc2.xml')
    topology, positions = fixer.topology, fixer.positions
    system = ff.createSystem(topology, nonbondedMethod=app.NoCutoff)
    nb = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
    gb = next(f for f in system.getForces() if isinstance(f, (mm.GBSAOBCForce, mm.CustomGBForce)))
    atoms = list(topology.atoms())
    nbp = [nb.getParticleParameters(i) for i in range(len(atoms))]
    if isinstance(gb, mm.GBSAOBCForce):
        gbp = [gb.getParticleParameters(i) for i in range(len(atoms))]
        radii = [p[1].value_in_unit(u.angstrom) for p in gbp]
        scales = [p[2] for p in gbp]
    else:
        # amber14 implicit/obc2.xml CustomGBForce: (charge, offset radius 'or' nm, scaled offset radius 'sr' nm)
        gbp = [gb.getParticleParameters(i) for i in range(len(atoms))]
        radii = [p[1] * 10 + 0.09 for p in gbp]
        scales = [p[2] / p[1] for p in gbp]
    bonds = [[] for _ in atoms]
    for b in topology.bonds():
        bonds[b.atom1.index].append(b.atom2.index)
        bonds[b.atom2.index].append(b.atom1.index)
    chain_ids = np.array([a.residue.chain.id for a in atoms])
    agent1 = np.where(np.isin(chain_ids, list(c1)))[0]
    agent2 = np.where(np.isin(chain_ids, list(c2)))[0]
    np.savez(f'{outdir}/{complex_id}.npz',
             coords=np.array(positions.value_in_unit(u.angstrom))[None],
             charges=np.array([p[0].value_in_unit(u.elementary_charge) for p in nbp]),
             sigma=np.array([p[1].value_in_unit(u.angstrom) for p in nbp]),
             epsilon=np.array([p[2].value_in_unit(u.kilocalorie_per_mole) for p in nbp]),
             gb_radius=np.array(radii), gb_scale=np.array(scales),
             elements=np.array([a.element.symbol for a in atoms]), names=np.array([a.name for a in atoms]),
             resnames=np.array([a.residue.name for a in atoms]),
             # Residue numbers keep their sign, since SKEMPI names mutations at negative numbers too (e.g. D-13)
             resnums=np.array([int(''.join(ch for ch in a.residue.id if ch.isdigit() or ch == '-')) for a in atoms]),
             icodes=np.array([a.residue.insertionCode.strip() for a in atoms]),
             chains=chain_ids, resindex=np.array([a.residue.index for a in atoms]),
             bonds=np.array(json.dumps(bonds)), agent1=agent1, agent2=agent2,
             interaction=np.array(json.dumps({'name': complex_id})))
    with open(f'{outdir}/{complex_id}.pdb', 'w') as f:
        app.PDBFile.writeFile(topology, positions, f, keepIds=True)
    print(complex_id, 'atoms', len(atoms), 'agents', len(agent1), len(agent2),
          'charges', round(sum(nbp[i][0].value_in_unit(u.elementary_charge) for i in agent1), 2),
          round(sum(nbp[i][0].value_in_unit(u.elementary_charge) for i in agent2), 2), 'ala ddG', len(ddg), flush=True)
