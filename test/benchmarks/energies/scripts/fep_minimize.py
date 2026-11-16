"""Restrained minimization of each FEP+ complex before scoring.

The FEP+ inputs are unrelaxed starting poses (ligands placed by core alignment), and many substituents clash
with the protein. MD frames never clash, so every complex is relaxed the same way before any method scores it:
ff14SB + GAFF2/AM1-BCC, reaction-field electrostatics (1.2 nm, solvent eps 78.5), all heavy atoms restrained
to the input positions (5 kcal/mol/A^2), hydrogens free. atoms.npz gets the minimized coordinates; the input
ones are kept in atoms_raw.npz.

Run in openmm_env: python fep_minimize.py <fep work dir> [series ...]
"""
import os
import sys
import json
import glob
import shutil
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as u
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.argv, args = sys.argv[:1], sys.argv[1:]
from fep_receptor import build_receptor

work = args[0]
series = args[1:] or sorted(os.path.basename(os.path.dirname(p)) for p in glob.glob(f'{work}/*/series.json'))
K = 5.0 * 418.4  # 5 kcal/mol/A^2 in kJ/mol/nm^2
platform = mm.Platform.getPlatformByName('CUDA')


def merge(rec_system, lig_system):
    """Append the ligand's particles and forces to a copy of the receptor system."""
    system = mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(rec_system))
    n = system.getNumParticles()
    for i in range(lig_system.getNumParticles()):
        system.addParticle(lig_system.getParticleMass(i))
    nb = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
    lnb = next(f for f in lig_system.getForces() if isinstance(f, mm.NonbondedForce))
    for i in range(lnb.getNumParticles()):
        nb.addParticle(*lnb.getParticleParameters(i))
    for i in range(lnb.getNumExceptions()):
        a, b, q, s, e = lnb.getExceptionParameters(i)
        nb.addException(a + n, b + n, q, s, e)
    for f in lig_system.getForces():
        if isinstance(f, mm.HarmonicBondForce):
            g = mm.HarmonicBondForce()
            for i in range(f.getNumBonds()):
                a, b, *p = f.getBondParameters(i)
                g.addBond(a + n, b + n, *p)
        elif isinstance(f, mm.HarmonicAngleForce):
            g = mm.HarmonicAngleForce()
            for i in range(f.getNumAngles()):
                a, b, c, *p = f.getAngleParameters(i)
                g.addAngle(a + n, b + n, c + n, *p)
        elif isinstance(f, mm.PeriodicTorsionForce):
            g = mm.PeriodicTorsionForce()
            for i in range(f.getNumTorsions()):
                a, b, c, d, *p = f.getTorsionParameters(i)
                g.addTorsion(a + n, b + n, c + n, d + n, *p)
        else:
            continue
        system.addForce(g)
    nb.setNonbondedMethod(mm.NonbondedForce.CutoffNonPeriodic)
    nb.setCutoffDistance(1.2 * u.nanometer)
    return system


for name in series:
    sdir = f'{work}/{name}'
    out_file = f'{sdir}/minimize.json'
    done = json.load(open(out_file)) if os.path.exists(out_file) else {}
    info = json.load(open(f'{sdir}/series.json'))
    todo = [l for l in sorted(info['ligands']) if l not in done and os.path.exists(f'{sdir}/{l}/atoms.npz')]
    if not todo:
        continue
    topology, positions, rec_system, meta = build_receptor(info['protein'], name)
    rec = np.load(f'{sdir}/rec.npz')
    # Same atoms in the same order; added hydrogens may sit slightly differently (Modeller relaxes them),
    # so the system is only used for parameters and the coordinates come from atoms.npz
    elements = np.array([a.element.symbol for a in topology.atoms()])
    assert (elements == rec['elements']).all(), name
    heavy_rec = elements != 'H'
    assert np.abs(np.array(positions.value_in_unit(u.angstrom))[heavy_rec] - rec['coords'][heavy_rec]).max() < 1e-3, name
    for lig in todo:
        ldir = f'{sdir}/{lig}'
        if not os.path.exists(f'{ldir}/atoms_raw.npz'):
            shutil.copy(f'{ldir}/atoms.npz', f'{ldir}/atoms_raw.npz')
        raw = dict(np.load(f'{ldir}/atoms_raw.npz'))
        lig_system = app.AmberPrmtopFile(f'{ldir}/lig.prmtop').createSystem(nonbondedMethod=app.NoCutoff)
        system = merge(rec_system, lig_system)
        xyz0 = raw['coords'] / 10.0  # nm
        heavy = np.where(raw['elements'] != 'H')[0]
        restraint = mm.CustomExternalForce('0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
        restraint.addGlobalParameter('k', K)
        for p in ('x0', 'y0', 'z0'):
            restraint.addPerParticleParameter(p)
        for i in heavy:
            restraint.addParticle(int(i), xyz0[i].tolist())
        system.addForce(restraint)
        context = mm.Context(system, mm.VerletIntegrator(0.001), platform)
        context.setPositions(xyz0 * u.nanometer)
        e0 = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole)
        mm.LocalEnergyMinimizer.minimize(context, 10.0, 5000)
        state = context.getState(getEnergy=True, getPositions=True)
        e1 = state.getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole)
        xyz = state.getPositions(asNumpy=True).value_in_unit(u.angstrom)
        n_rec = int(raw['n_rec'])
        lig_heavy = heavy[heavy >= n_rec]
        rmsd = lambda idx: float(np.sqrt(((xyz[idx] - xyz0[idx] * 10) ** 2).sum(1).mean()))
        done[lig] = dict(e0=e0, e1=e1, heavy_rmsd=rmsd(heavy), ligand_heavy_rmsd=rmsd(lig_heavy),
                         max_heavy_shift=float(np.linalg.norm(xyz[heavy] - xyz0[heavy] * 10, axis=1).max()))
        raw['coords'] = np.array(xyz)
        np.savez(f'{ldir}/atoms.npz', **raw)
        json.dump(done, open(out_file, 'w'), indent=1)
        print(name, lig, {k: round(v, 2) for k, v in done[lig].items()}, flush=True)
