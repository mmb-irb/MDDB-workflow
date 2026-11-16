"""Receptor (ff14SB) for each selected FEP+ series, then per-complex atoms.npz files (receptor first, then ligand).

Run in openmm_env: python fep_receptor.py <fep work dir> [series ...]
Receptor: the prepared FEP+ protein with its hydrogens kept (Modeller only adds missing ones), modified
residues mapped to their standard parents, and waters, ions, glycans and cofactors removed, as the energies
analysis keeps only the two interacting agents. Ligand: lig.prmtop from fep_ligand.sh (GAFF2, AM1-BCC).
"""
import os
import sys
import io
import json
import glob
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as u
from pdbfixer import PDBFixer

work = sys.argv[1] if len(sys.argv) > 1 else ''
series = sys.argv[2:] or sorted(os.path.basename(os.path.dirname(p)) for p in glob.glob(f'{work}/*/series.json'))
KCAL = u.kilocalorie_per_mole
ff = app.ForceField('amber14-all.xml')


def clean_caps(path):
    """Maestro writes N-methyl caps as NMA (CA, 1HA..) and some caps as UNK; rename both to Amber ACE/NME."""
    # HXT on a C-terminus without OXT (Maestro H-capped carbonyl) has no Amber template: drop it, OXT is added later
    lines = [l for l in open(path) if l.startswith(('ATOM', 'HETATM', 'TER', 'END')) and not (l.startswith(('ATOM', 'HETATM')) and l[12:16].strip() == 'HXT')]
    groups = {}
    for i, l in enumerate(lines):
        if l.startswith(('ATOM', 'HETATM')) and l[17:20] in ('NMA', 'UNK'):
            groups.setdefault(l[17:27], []).append(i)
    for key, idx in groups.items():
        atoms = {i: (lines[i][12:16].strip(), lines[i][76:78].strip(), np.array([float(lines[i][30 + 8 * k:38 + 8 * k]) for k in range(3)])) for i in idx}
        heavy = {i: a for i, a in atoms.items() if a[1] != 'H'}
        has_o = any(a[1] == 'O' for a in heavy.values())
        rename = {}
        if has_o:  # ACE: carbonyl C (bonded to O), methyl C, O
            o = next(i for i, a in heavy.items() if a[1] == 'O')
            c = min((i for i, a in heavy.items() if a[1] == 'C'), key=lambda i: np.linalg.norm(heavy[i][2] - heavy[o][2]))
            rename[o], rename[c] = 'O', 'C'
            parent_names = {i: 'CH3' for i, a in heavy.items() if a[1] == 'C' and i != c}
            rename.update(parent_names)
            hnames = iter(['H1', 'H2', 'H3'])
            for i, a in atoms.items():
                if a[1] == 'H':
                    rename[i] = next(hnames)
            res = 'ACE'
        else:  # NME: N, methyl C; H on N is H, H on C are H1-H3
            n = next(i for i, a in heavy.items() if a[1] == 'N')
            for i, a in heavy.items():
                rename[i] = 'N' if i == n else 'C'
            hnames = iter(['H1', 'H2', 'H3'])
            for i, a in atoms.items():
                if a[1] == 'H':
                    near = min(heavy, key=lambda j: np.linalg.norm(heavy[j][2] - a[2]))
                    rename[i] = 'H' if near == n else next(hnames)
            res = 'NME'
        for i, new in rename.items():
            l = lines[i]
            lines[i] = 'ATOM  ' + l[6:12] + (' ' + new.ljust(3) if len(new) < 4 else new) + l[16] + res + l[20:]
    return io.StringIO(''.join(lines))


def params(system, topology, positions):
    nb = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
    p = [nb.getParticleParameters(i) for i in range(system.getNumParticles())]
    atoms = list(topology.atoms())
    bonds = [[] for _ in atoms]
    for b in topology.bonds():
        bonds[b.atom1.index].append(b.atom2.index)
        bonds[b.atom2.index].append(b.atom1.index)
    return dict(coords=np.array(positions.value_in_unit(u.angstrom)),
                charges=np.array([x[0].value_in_unit(u.elementary_charge) for x in p]),
                sigma=np.array([x[1].value_in_unit(u.angstrom) for x in p]),
                epsilon=np.array([x[2].value_in_unit(KCAL) for x in p]),
                elements=np.array([a.element.symbol for a in atoms]), names=np.array([a.name for a in atoms]),
                resnames=np.array([a.residue.name for a in atoms]),
                resnums=np.array([int(''.join(c for c in a.residue.id if c.isdigit() or c == '-') or 0) for a in atoms]),
                bonds=bonds)


def build_receptor(protein, name=''):
    """Receptor topology, positions and ff14SB system; rebuilds all hydrogens if the kept ones match no template."""
    for strip_h in (False, True):
        fixer = PDBFixer(pdbfile=clean_caps(protein))
        if strip_h:  # fallback: rebuild all hydrogens at pH 7
            modeller = app.Modeller(fixer.topology, fixer.positions)
            modeller.delete([a for a in modeller.topology.atoms() if a.element is not None and a.element.symbol == 'H'])
            fixer.topology, fixer.positions = modeller.topology, modeller.positions
        n_h0 = sum(a.element.symbol == 'H' for a in fixer.topology.atoms())
        fixer.findMissingResidues()
        fixer.missingResidues = {}  # do not model missing loops
        fixer.findNonstandardResidues()
        replaced = sorted({r.name for r, _ in fixer.nonstandardResidues})
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(keepWater=False)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)
        n_h1 = sum(a.element.symbol == 'H' for a in fixer.topology.atoms())
        try:
            system = ff.createSystem(fixer.topology, nonbondedMethod=app.NoCutoff)
        except ValueError as error:
            print(name, 'template error' + (' after rebuilding H' if strip_h else ', rebuilding H'), str(error)[:150], flush=True)
            continue
        return fixer.topology, fixer.positions, system, dict(strip_h=strip_h, n_h0=n_h0, n_h1=n_h1, replaced=replaced)
    return None


def main():
  for name in series:
    sdir = f'{work}/{name}'
    info = json.load(open(f'{sdir}/series.json'))
    if not os.path.exists(f'{sdir}/rec.npz'):
        built = build_receptor(info['protein'], name)
        if built is None:
            continue
        topology, positions, system, meta = built
        rec = params(system, topology, positions)
        # Duplicate hydrogens (a name mismatch would make Modeller add a second copy) show up as close H pairs
        xyz = rec['coords']
        h = np.where(rec['elements'] == 'H')[0]
        dh = np.linalg.norm(xyz[h][:, None] - xyz[h][None], axis=-1)
        close = int((np.triu(dh < 0.9, 1)).sum())
        np.savez(f'{sdir}/rec.npz', **{k: v for k, v in rec.items() if k != 'bonds'}, bonds=np.array(json.dumps(rec['bonds'])),
                 chains=np.array([a.residue.chain.id for a in topology.atoms()]))
        print(name, 'atoms', len(xyz), 'H kept', meta['n_h0'], 'H after', meta['n_h1'], 'close H pairs', close, 'charge',
              round(rec['charges'].sum(), 2), 'replaced', meta['replaced'], 'rebuilt H' if meta['strip_h'] else '', flush=True)
    if not os.path.exists(f'{sdir}/rec.npz'):
        continue
    rec = np.load(f'{sdir}/rec.npz')
    rec_bonds = json.loads(str(rec['bonds']))
    n_rec = len(rec['charges'])
    for lig in sorted(info['ligands']):
        ldir = f'{sdir}/{lig}'
        if os.path.exists(f'{ldir}/atoms.npz') or not os.path.exists(f'{ldir}/lig.prmtop'):
            continue
        prmtop, inpcrd = app.AmberPrmtopFile(f'{ldir}/lig.prmtop'), app.AmberInpcrdFile(f'{ldir}/lig.inpcrd')
        system = prmtop.createSystem(nonbondedMethod=app.NoCutoff)
        l = params(system, prmtop.topology, inpcrd.positions)
        bonds = rec_bonds + [[j + n_rec for j in b] for b in l['bonds']]
        np.savez(f'{ldir}/atoms.npz',
                 **{k: np.concatenate([rec[k], l[k]]) for k in ('coords', 'charges', 'sigma', 'epsilon', 'elements', 'names', 'resnames', 'resnums')},
                 chains=np.concatenate([np.where(rec['chains'] == 'L', 'A', rec['chains']), ['L'] * len(l['charges'])]),
                 bonds=np.array(json.dumps(bonds)), n_rec=n_rec)


if __name__ == '__main__':
    main()
