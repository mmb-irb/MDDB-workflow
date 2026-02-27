import re
import numpy as np
from mdahole2.analysis import HoleAnalysis
from mddb_workflow.tools.get_reduced_trajectory import calculate_frame_step
from mddb_workflow.utils.auxiliar import save_json
from mddb_workflow.utils.constants import OUTPUT_CHANNELS_FILENAME
from mddb_workflow.utils.type_hints import *


def channels(
    membrane_map: dict,
    universe: 'Universe',
    output_directory: str,
    snapshots: int,
    frames_limit: int
):
    """Analyze channels in a membrane protein using MDAnalysis mda_hole."""
    if membrane_map is None or membrane_map['n_mems'] == 0:
        print('-> No membranes found, skipping channels analysis')
        return
    if universe.select_atoms('protein').n_atoms == 0:
        print('-> No protein found, skipping channels analysis')
        return
    frame_step, frame_count = calculate_frame_step(snapshots, frames_limit)
    hole = HoleAnalysis(universe)
    hole.run(step=frame_step)
    hole.create_vmd_surface(output_directory + '/hole.vmd', dot_density=13)
    hole.delete_temporary_files()

    with open(output_directory + '/hole.vmd', 'r') as f:
        lines = f.readlines()

    # Find lines with triangle coordinates
    trinorms = []
    for i, line in enumerate(lines):
        if i > 3 and 'set triangle' in line:
            vmd_set = re.sub(r'set triangles\(\d+\)', '', line)  # Remove set triangles(i)
            vmd_set = re.sub(r'\{(\s*-?\d[^\s]*)(\s*-?\d[^\s]*)(\s*-?\d[^}]*)\}', r'[\1,\2,\3]', vmd_set)  # Convert { x y z } to [x,y,z]
            vmd_set = vmd_set.replace('{', '[').replace('}', ']')  # Convert { to [ and } to ]
            vmd_set = re.sub(r'\]\s*\[', '], [', vmd_set)  # Add commas between brackets
            vmd_set = eval(vmd_set.strip())  # Evaluate string as list
            # different hole colors
            trinorms.append(vmd_set)
    # Create a list of positions, colors, and normals
    assert frame_count == len(trinorms), f'Frame count {frame_count} does not match trinorms length {len(trinorms)}'
    data = {}
    for frame in range(frame_count):
        poss = []
        cols = [0, 0, 0]
        z_range = 0
        for i in range(3):  # RGB
            # Get triangle coordinates for this frame and color
            # Red are low radiues, green medium, blue high
            trinorms_cl = np.array(trinorms[frame][i])  # (N, 6, 3) # 6: 3 positions + 3 normals, 3: vertex
            if len(trinorms_cl) == 0:
                continue
            pos = trinorms_cl[:, :3, :]
            z_range = max(pos[..., 2].max() - pos[..., 2].min(), z_range)
            poss.append(pos.flatten())
            cols[i] = pos.size
        poss = np.concatenate(poss)
        if z_range < 40:
            # failed to find channel on this frame
            continue
        # compute cumulative offsets for colors
        data[frame] = {
            'position': poss.tolist(),
            'color_offset': np.cumsum(cols).tolist(),
        }
    if len(data) == 0:
        print('-> No channels found')
        return
    data_to_save = {'data': data, 'n_frames': frame_count}
    save_json(data_to_save, output_directory + '/' + OUTPUT_CHANNELS_FILENAME)
