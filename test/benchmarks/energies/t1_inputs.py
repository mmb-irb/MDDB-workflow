"""Fetch and build the T1 input frames when they are first needed, in downloads/t1/<system>.npz (ignored by git).

- A01TW (kinase-crizotinib) and A0221 (TNF-antibody): MDDB projects of the irb-dev node, the workflow's default API.
  The workflow downloads a 5-frame sample of each replica (`mwf run -proj <accession> -i download -smp 5`) and runs its
  interactions and charges tasks on the first one (`mwf run -mdir <replica> -i inter charges -t`). The frame is read as
  the energies analysis reads it, through a PDB file, so coordinates have 3 decimals. A01TW uses frame 3 (in frames
  0-2 the crizotinib overlaps the protein, because the deposited trajectory has a transposed box) and A0221 frame 2
  (the largest two-direction mismatch of CMIP 1.1.0).
- 1VFB (D1.3 Fv-lysozyme): the SKEMPI complex 1VFB_AB_C of T2, downloads/t2/skempi/1VFB_AB_C.npz, which
  scripts/skempi_prepare.py writes (crystal structure, PDBFixer at pH 7, Amber ff14SB charges).

Each file has the frame's coordinates, atoms, residues, chains, charges, bonds and agents.
"""
import os
import sys
import json
import subprocess
from pathlib import Path
import numpy as np

DIRECTORY = Path(__file__).parent / 'downloads' / 't1'
REPOSITORY = Path(__file__).parents[3]
MDDB_SYSTEMS = {
    'A01TW': dict(md='replica_1', frame=3),
    'A0221': dict(md='310_K_37C_seed_1', frame=2),
}
SKEMPI_SYSTEMS = {
    '1VFB': dict(complex='1VFB_AB_C', agents=('D1.3 Fv', 'Lysozyme')),
}
SNAPSHOTS = 5
INTERACTION_KEYS = ('name', 'agent_1', 'agent_2', 'strong_bonds')


def get_input_filepath(system: str) -> Path:
    """Get the input file of a T1 system, fetching and building it if it is not there yet."""
    filepath = DIRECTORY / f'{system}.npz'
    if not filepath.exists():
        DIRECTORY.mkdir(parents=True, exist_ok=True)
        data = build_mddb_system(system) if system in MDDB_SYSTEMS else build_skempi_system(system)
        save(filepath, data)
    return filepath


def save(filepath: Path, data: dict):
    """Save an input file through a temporary file, so an interrupted build leaves no truncated file."""
    temporary_filepath = filepath.with_name(filepath.name + '.part')
    with open(temporary_filepath, 'wb') as file:
        np.savez_compressed(file, **data)
    os.replace(temporary_filepath, filepath)


def run_workflow(arguments: list[str], directory: Path, log_filepath: Path):
    """Run the workflow's command line of this repository in a directory."""
    command = [sys.executable, '-c', 'from mddb_workflow.console import main; main()', *arguments]
    environment = dict(os.environ, PYTHONPATH=str(REPOSITORY))
    print(f'  mwf {" ".join(arguments)} (log in {log_filepath})', flush=True)
    with open(log_filepath, 'w') as log:
        process = subprocess.run(command, cwd=directory, env=environment, stdout=log, stderr=subprocess.STDOUT)
    if process.returncode != 0:
        raise RuntimeError(f'mwf {" ".join(arguments)} failed, see {log_filepath}')


def build_mddb_system(system: str) -> dict:
    """Download an MDDB project with the workflow, run the tasks the energies analysis needs and take one frame."""
    from mddb_workflow.utils.structures import Structure
    from mddb_workflow.tools.get_pdb_frames import get_pdb_frames
    md, frame = MDDB_SYSTEMS[system]['md'], MDDB_SYSTEMS[system]['frame']
    print(f'Fetching the T1 system {system} from MDDB', flush=True)
    project_directory = DIRECTORY / 'sources' / system
    project_directory.mkdir(parents=True, exist_ok=True)
    md_directory = project_directory / md
    if not (md_directory / 'trajectory.xtc').exists():
        run_workflow(['run', '-proj', system, '-i', 'download', '-smp', str(SNAPSHOTS)], project_directory,
                     project_directory / 'download.log')
    run_workflow(['run', '-mdir', md, '-i', 'inter', 'charges', '-t'], project_directory, project_directory / 'tasks.log')
    cache = json.load(open(md_directory / '.mwf_cache.json'))
    interaction = cache['inter_task_output'][0]
    structure = Structure.from_mmcif_file(str(md_directory / 'structure.cif'))
    # The frames come as one PDB file that is rewritten for each frame, so the frame is read while iterating
    frames, _, _ = get_pdb_frames(structure, str(md_directory / 'trajectory.xtc'), cache['snapshots'], SNAPSHOTS)
    for index, frame_pdb in enumerate(frames):
        if index == frame:
            coordinates = np.array([atom.coords for atom in Structure.from_pdb_file(frame_pdb).atoms], dtype=float)
    atoms = structure.atoms
    return dict(
        coords=coordinates,
        elements=np.array([atom.element for atom in atoms]),
        names=np.array([atom.name for atom in atoms]),
        resnames=np.array([atom.residue.name for atom in atoms]),
        resnums=np.array([atom.residue.number for atom in atoms]),
        icodes=np.full(len(atoms), ''),
        chains=np.array([atom.chain.name for atom in atoms]),
        resindex=np.array([atom.residue.index for atom in atoms]),
        charges=np.array(cache['charges_task_output'], dtype=float),
        bonds=bond_pairs(cache['refbonds_task_output']),
        agent1=np.array(interaction['atom_indices_1']),
        agent2=np.array(interaction['atom_indices_2']),
        interaction=np.array(json.dumps({key: interaction[key] for key in INTERACTION_KEYS})),
        source=np.array(f'MDDB {system} (irb-dev), {md}, frame {frame} of a {SNAPSHOTS}-frame sample'),
    )


def build_skempi_system(system: str) -> dict:
    """Take the SKEMPI complex from the T2 inputs."""
    complex_id, (agent_1, agent_2) = SKEMPI_SYSTEMS[system]['complex'], SKEMPI_SYSTEMS[system]['agents']
    skempi_filepath = DIRECTORY.parent / 't2' / 'skempi' / f'{complex_id}.npz'
    if not skempi_filepath.exists():
        raise FileNotFoundError(f'{skempi_filepath} is missing: prepare it with scripts/skempi_prepare.py, in an '
                                f'environment with OpenMM and PDBFixer (see README.md, T2 and T3, SKEMPI)')
    with np.load(skempi_filepath) as source:
        data = {key: source[key] for key in source.files}
    interaction = dict(name=f'{agent_1}-{agent_2} interaction', agent_1=agent_1, agent_2=agent_2, strong_bonds=[])
    return dict(
        coords=data['coords'][0],
        **{key: data[key] for key in ('elements', 'names', 'resnames', 'resnums', 'icodes', 'chains', 'resindex',
                                      'charges', 'agent1', 'agent2')},
        bonds=bond_pairs(json.loads(str(data['bonds']))),
        interaction=np.array(json.dumps(interaction)),
        source=np.array(f'SKEMPI 2.0 {complex_id}: PDB 1VFB chains A, B and C, PDBFixer at pH 7, Amber ff14SB charges'),
    )


def bond_pairs(bonds: list[list[int]]) -> np.ndarray:
    """Get the bonded atom pairs from a list of bonded atoms per atom."""
    return np.array(sorted({(min(a, b), max(a, b)) for a, partners in enumerate(bonds) for b in partners}), dtype=np.int32)
