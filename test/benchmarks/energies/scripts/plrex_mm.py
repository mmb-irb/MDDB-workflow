"""MM interaction energies of PL-REX complexes with OpenMM: gas phase, GBn2 (with and without the
surface-area term), and an export of charges / LJ / coordinates for the other methods.

Run in openmm_env: python plrex_mm.py <plrex work dir> <output json> [platform CPU|CUDA]
Each complex dir needs com/rec/lig prmtop+inpcrd (plrex_build.sh).
"""
import os
import sys
import json
import glob
import time
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as u

work, out_file = sys.argv[1:3]
platform = mm.Platform.getPlatformByName(sys.argv[3] if len(sys.argv) > 3 else 'CPU')
results = json.load(open(out_file)) if os.path.exists(out_file) else {}
KCAL = u.kilocalorie_per_mole


def energy(prmtop, inpcrd, **kwargs):
    system = prmtop.createSystem(nonbondedMethod=app.NoCutoff, constraints=None, **kwargs)
    context = mm.Context(system, mm.VerletIntegrator(0.001), platform)
    context.setPositions(inpcrd.positions)
    start = time.perf_counter()
    e = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KCAL)
    return e, time.perf_counter() - start


for cdir in sorted(glob.glob(f'{work}/0*/*/')):
    name = '/'.join(cdir.rstrip('/').split('/')[-2:])
    if name in results or not os.path.exists(f'{cdir}/com.prmtop'):
        continue
    parts = {k: (app.AmberPrmtopFile(f'{cdir}/{k}.prmtop'), app.AmberInpcrdFile(f'{cdir}/{k}.inpcrd')) for k in ('com', 'rec', 'lig')}
    row = {}
    for label, kwargs in (('vacuum', {}), ('gbn2', dict(implicitSolvent=app.GBn2, gbsaModel=None)),
                          ('gbn2_sa', dict(implicitSolvent=app.GBn2, gbsaModel='ACE'))):
        e = {k: energy(p, c, **kwargs) for k, (p, c) in parts.items()}
        row[label] = e['com'][0] - e['rec'][0] - e['lig'][0]
        row[f'{label}_seconds'] = sum(v[1] for v in e.values())
    # Export per-atom data of the complex (receptor atoms first, then ligand) for the other methods
    prmtop, inpcrd = parts['com']
    system = prmtop.createSystem(nonbondedMethod=app.NoCutoff)
    nb = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
    params = [nb.getParticleParameters(i) for i in range(system.getNumParticles())]
    atoms = list(prmtop.topology.atoms())
    bonds = [[] for _ in atoms]
    for b in prmtop.topology.bonds():
        bonds[b.atom1.index].append(b.atom2.index)
        bonds[b.atom2.index].append(b.atom1.index)
    n_rec = parts['rec'][0].topology.getNumAtoms()
    np.savez(f'{cdir}/atoms.npz',
             coords=np.array(inpcrd.positions.value_in_unit(u.angstrom)),
             charges=np.array([p[0].value_in_unit(u.elementary_charge) for p in params]),
             sigma=np.array([p[1].value_in_unit(u.angstrom) for p in params]),
             epsilon=np.array([p[2].value_in_unit(KCAL) for p in params]),
             elements=np.array([a.element.symbol for a in atoms]),
             names=np.array([a.name for a in atoms]), resnames=np.array([a.residue.name for a in atoms]),
             resnums=np.array([int(a.residue.id) for a in atoms]), chains=np.array(['A' if a.index < n_rec else 'L' for a in atoms]),
             bonds=np.array(json.dumps(bonds)), n_rec=n_rec)
    row['atoms'] = len(atoms)
    results[name] = row
    print(name, {k: round(v, 2) for k, v in row.items()}, flush=True)
    json.dump(results, open(out_file, 'w'), indent=1)
