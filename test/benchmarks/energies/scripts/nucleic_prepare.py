"""Prepare the protein-nucleic acid alanine scanning set of T4 (mmCSM-NA single mutations) like the SKEMPI complexes.

Run in an environment with openmm and pdbfixer:
    python nucleic_prepare.py <single_mutations.csv> <output dir> [--minimum 6]
The mutations are those of mmCSM-NA (Nguyen et al., Bioinformatics 2021; data from ProNIT, dbAMEPNI and the
literature), downloaded from https://biosig.lab.uq.edu.au/mmcsm_na/static/datasets/single_mutations.csv. Their ddG is
negative for destabilizing mutations, so its sign is changed to the SKEMPI convention (positive = weaker binding).

For each PDB entry with at least `minimum` X -> Ala mutations (X not Ala or Gly), the agents are the mutated protein
chains and the nucleic acid chains with a heavy atom within 5 Å of them; other chains, waters and heterogens are
removed, missing atoms and hydrogens are added (pH 7) and charges come from Amber ff14SB, OL15 (DNA) and OL3 (RNA).
Writes <PDB>_<protein chains>_<nucleic chains>.npz (the SKEMPI format of skempi_prepare.py), .pdb and _ddg.json.
"""
import io
import os
import csv
import json
import argparse
import urllib.request
from collections import defaultdict
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as u
from pdbfixer import PDBFixer

NUCLEIC_RESIDUES = {'DA', 'DC', 'DG', 'DT', 'DU', 'DI', 'A', 'C', 'G', 'U', 'I'}
CONTACT_DISTANCE = 5.0
ONE_LETTER = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H',
              'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
              'TYR': 'Y', 'VAL': 'V'}


def read_mutations(filepath: str) -> dict:
    """Get the mean ddG (SKEMPI sign) of every X -> Ala mutation, per PDB entry: {pdb: {(chain, wt, residue): ddg}}."""
    values = defaultdict(lambda: defaultdict(list))
    for row in csv.DictReader(open(filepath)):
        chain, mutation = row['MUTATION'].split()
        wild_type, residue, mutant = mutation[0], mutation[1:-1], mutation[-1]
        if mutant != 'A' or wild_type in 'AG':
            continue
        values[row['PDB'].upper()][(chain, wild_type, residue)].append(-float(row['DDG']))
    return {pdb: {key: float(np.mean(v)) for key, v in mutations.items()} for pdb, mutations in values.items()}


def download(pdb: str) -> PDBFixer:
    """Read an entry from the RCSB, as PDB or else as mmCIF."""
    try:
        text = urllib.request.urlopen(f'https://files.rcsb.org/download/{pdb}.pdb', timeout=60).read().decode()
        return PDBFixer(pdbfile=io.StringIO(text))
    except Exception:
        text = urllib.request.urlopen(f'https://files.rcsb.org/download/{pdb}.cif', timeout=60).read().decode()
        return PDBFixer(pdbxfile=io.StringIO(text))


def is_nucleic(chain) -> bool:
    """Tell whether most residues of a chain are nucleotides."""
    residues = list(chain.residues())
    return sum(residue.name.strip() in NUCLEIC_RESIDUES for residue in residues) > len(residues) / 2


def heavy_coords(chain, positions) -> np.ndarray:
    """Get the heavy atom coordinates (Å) of a chain."""
    return np.array([positions[atom.index].value_in_unit(u.angstrom) for atom in chain.atoms()
                     if atom.element is not None and atom.element.symbol != 'H'])


def select_chains(fixer: PDBFixer, protein_chains: set) -> tuple[list, list]:
    """Get the mutated protein chains and the nucleic acid chains in contact with them."""
    from scipy.spatial import cKDTree
    chains = list(fixer.topology.chains())
    protein = [c for c in chains if c.id in protein_chains and not is_nucleic(c)]
    if not protein:
        return [], []
    protein_tree = cKDTree(np.concatenate([heavy_coords(c, fixer.positions) for c in protein]))
    nucleic = [c for c in chains if is_nucleic(c) and len(heavy_coords(c, fixer.positions)) and
               protein_tree.query(heavy_coords(c, fixer.positions), distance_upper_bound=CONTACT_DISTANCE)[0].min()
               <= CONTACT_DISTANCE]
    # A chain id may appear in several OpenMM chains (split by TER records)
    return list(dict.fromkeys(c.id for c in protein)), list(dict.fromkeys(c.id for c in nucleic))


def remove_five_prime_phosphates(topology, positions):
    """Remove the phosphate of the first nucleotide of each nucleic acid chain (Amber 5' templates have none)."""
    modeller = app.Modeller(topology, positions)
    atoms = [atom for chain in modeller.topology.chains() if is_nucleic(chain)
             for atom in next(iter(chain.residues())).atoms() if atom.name in ('P', 'OP1', 'OP2', 'OP3', 'O1P', 'O2P')]
    modeller.delete(atoms)
    return modeller.topology, modeller.positions


def build(fixer: PDBFixer, force_field):
    """Fix the structure and parametrize it, removing 5' phosphates when a template does not match."""
    fixer.findMissingResidues()
    fixer.missingResidues = {}  # do not model missing loops
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    topology, positions = fixer.topology, fixer.positions
    try:
        return topology, positions, force_field.createSystem(topology, nonbondedMethod=app.NoCutoff)
    except ValueError:
        topology, positions = remove_five_prime_phosphates(topology, positions)
        modeller = app.Modeller(topology, positions)
        modeller.addHydrogens(force_field, pH=7.0)
        topology, positions = modeller.topology, modeller.positions
        return topology, positions, force_field.createSystem(topology, nonbondedMethod=app.NoCutoff)


def prepare(pdb: str, mutations: dict, output_directory: str, minimum: int, force_field) -> str:
    """Prepare one entry; return a line of report."""
    fixer = download(pdb)
    protein_chains, nucleic_chains = select_chains(fixer, {chain for chain, _, _ in mutations})
    if not protein_chains or not nucleic_chains:
        return f'{pdb} skipped: no protein chain or no nucleic acid chain in contact'
    keep = set(protein_chains) | set(nucleic_chains)
    fixer.removeChains([i for i, c in enumerate(fixer.topology.chains()) if c.id not in keep])
    topology, positions, system = build(fixer, force_field)
    atoms = list(topology.atoms())
    # Mutations whose residue is in the structure with its wild type
    residue_names = {(a.residue.chain.id, a.residue.id.strip() + a.residue.insertionCode.strip()): a.residue.name
                     for a in atoms}
    ddg = {}
    for (chain, wild_type, residue), value in mutations.items():
        name = residue_names.get((chain, residue))
        if name is not None and ONE_LETTER.get({'HIE': 'HIS', 'HID': 'HIS', 'HIP': 'HIS', 'CYX': 'CYS'}.get(name, name)) \
                == wild_type:
            ddg[f'{chain}{residue}'] = value
    if len(ddg) < minimum:
        return f'{pdb} skipped: {len(ddg)} mutations match the structure'
    complex_id = f'{pdb}_{"".join(protein_chains)}_{"".join(nucleic_chains)}'
    nonbonded = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
    charges = np.array([nonbonded.getParticleParameters(i)[0].value_in_unit(u.elementary_charge)
                        for i in range(len(atoms))])
    bonds = [[] for _ in atoms]
    for bond in topology.bonds():
        bonds[bond.atom1.index].append(bond.atom2.index)
        bonds[bond.atom2.index].append(bond.atom1.index)
    chain_ids = np.array([a.residue.chain.id for a in atoms])
    agent1 = np.where(np.isin(chain_ids, protein_chains))[0]
    agent2 = np.where(np.isin(chain_ids, nucleic_chains))[0]
    np.savez(f'{output_directory}/{complex_id}.npz',
             coords=np.array(positions.value_in_unit(u.angstrom))[None], charges=charges,
             elements=np.array([a.element.symbol for a in atoms]), names=np.array([a.name for a in atoms]),
             resnames=np.array([a.residue.name for a in atoms]),
             resnums=np.array([int(''.join(ch for ch in a.residue.id if ch.isdigit() or ch == '-')) for a in atoms]),
             icodes=np.array([a.residue.insertionCode.strip() for a in atoms]),
             chains=chain_ids, resindex=np.array([a.residue.index for a in atoms]),
             bonds=np.array(json.dumps(bonds)), agent1=agent1, agent2=agent2,
             interaction=np.array(json.dumps({'name': complex_id})))
    with open(f'{output_directory}/{complex_id}_ddg.json', 'w') as file:
        json.dump(ddg, file, indent=1)
    with open(f'{output_directory}/{complex_id}.pdb', 'w') as file:
        app.PDBFile.writeFile(topology, positions, file, keepIds=True)
    return (f'{complex_id} atoms {len(atoms)}, agents {len(agent1)} {len(agent2)}, charges '
            f'{charges[agent1].sum():.1f} {charges[agent2].sum():.1f}, {len(ddg)} mutations')


def main():
    """Prepare every entry with enough mutations."""
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('mutations')
    parser.add_argument('output_directory')
    parser.add_argument('--minimum', type=int, default=6)
    args = parser.parse_args()
    os.makedirs(args.output_directory, exist_ok=True)
    force_field = app.ForceField('amber14-all.xml')
    for pdb, mutations in sorted(read_mutations(args.mutations).items()):
        if len(mutations) < args.minimum:
            continue
        try:
            print(prepare(pdb, mutations, args.output_directory, args.minimum, force_field), flush=True)
        except Exception as error:
            print(f'{pdb} FAILED: {error!r}'[:300], flush=True)


if __name__ == '__main__':
    main()
