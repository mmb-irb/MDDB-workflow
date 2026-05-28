from pathlib import Path
import MDAnalysis as mda
import subprocess
from mddb_workflow.tools.get_reduced_trajectory import get_reduced_trajectory
from mddb_workflow.utils.auxiliar import save_json, load_json
from mddb_workflow.utils.file import File
from mddb_workflow.utils.constants import OUTPUT_CHANNELS_FILENAME
from mddb_workflow.utils.type_hints import *
import json


def _parse_chap_residue_fractions(pdb_file: Path) -> dict[str, list[float]]:
    """Parse residue-level pore-lining and pore-facing fractions from CHAP output."""
    # https://www.channotation.org/docs/contents_pdb_file/
    residue_fractions = {}
    if not Path(pdb_file).exists():
        return residue_fractions
    pore_lining = []
    pore_facing = []
    residue_keys = []
    with open(pdb_file, 'r') as pdb_handle:
        for line in pdb_handle:
            if not line.startswith(('ATOM', 'HETATM')):
                continue

            chain_id = line[21].strip() or '_'
            residue_index = line[22:26].strip()
            insertion_code = line[26].strip()
            residue_key = f'{chain_id}:{residue_index}{insertion_code}' if insertion_code else f'{chain_id}:{residue_index}'

            if residue_key in residue_keys:
                continue

            residue_keys.append(residue_key)
            pore_lining.append(float(line[54:60]))
            pore_facing.append(float(line[60:66]))

    return {'pore_lining': pore_lining, 'pore_facing': pore_facing}


def channels(
    membrane_map: dict,
    universe: 'Universe',
    output_directory: str,
    thickness_analysis: dict,
    snapshots: int,
    frames_limit: int
):
    """Analyze channels in a membrane protein using MDAnalysis mda_hole."""
    if membrane_map is None or membrane_map['n_mems'] == 0:
        print('  -> No membranes found, skipping channels analysis')
        return
    if universe.select_atoms('protein').n_atoms == 0:
        print('  -> No protein found, skipping channels analysis')
        return
    print('  -> Running channels analysis')
    frame = 0
    top = thickness_analysis['data']['mean_positive'][frame] + thickness_analysis['data']['midplane_z'][frame]
    bot = thickness_analysis['data']['mean_negative'][frame] + thickness_analysis['data']['midplane_z'][frame]
    # Select protein between membrane planes
    TM_at = universe.select_atoms(f'protein and prop z > {bot} and prop z < {top}').residues.atoms
    with mda.selections.gromacs.SelectionWriter(f'{output_directory}/membrane_atoms.ndx', mode='w') as ndx:
        ndx.write(universe.select_atoms('protein').atoms, name='Protein')
        ndx.write(TM_at, name='membrane_atoms')
    pdb = Path(universe.filename).absolute()
    reduced_trajectory_file, frame_step, n_frames = get_reduced_trajectory(
        input_topology_file=File(universe.filename),
        input_trajectory_file=File(universe.trajectory.filename),
        snapshots=snapshots,
        reduced_trajectory_frames_limit=frames_limit
    )
    # TODO: find better fallbacks values
    cmd = ("chap",
           "-f", str(Path(reduced_trajectory_file).absolute()),
           "-s", str(pdb),
           "-n", "membrane_atoms.ndx",
           "-sel-pathway", "0",
           "-pf-sel-ipp", "1",
           "-hydrophob-fallback", "-0.568",  # hardcoded value for histidine forms
           "-pf-vdwr-fallback", "0.2",       # hardcoded value for biggest radiues
           "-out-num-points", "200",
           "-out-detailed"
    )
    cmd = ' '.join(cmd)
    print('  -> Running command:', cmd)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=output_directory)
    if result.returncode != 0:
        failt_str = 'Pore radius at initial probe position is infinite'
        if failt_str in result.stderr:
            print('  -> No channels found.')
            return
        # Temporal raise to detect and fix possible issues with CHAP execution,
        # but we should handle this more gracefully in the future
        raise RuntimeError(f'CHAP command failed with return code {result.returncode}. Stderr: {result.stderr}')
    output = load_json(f'{output_directory}/output.json')
    stream = f'{output_directory}/stream_output.json'
    stream_data = []
    # Each line in the stream output is a JSON object
    with open(stream, 'r') as f:
        for line in f:
            stream_data.append(json.loads(line))

    data_to_save = {
        'data': {
            'output': output,
            'stream': stream_data,
            'pore_residues': _parse_chap_residue_fractions(
                f"{output_directory}/output.pdb"),
        },
        'version': '0.1.0',
    }
    save_json(data_to_save, output_directory + '/' + OUTPUT_CHANNELS_FILENAME)
