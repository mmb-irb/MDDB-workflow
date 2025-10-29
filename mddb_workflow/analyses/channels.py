from biobb_mem.mdanalysis_biobb.mda_hole import mda_hole
from mddb_workflow.tools.get_reduced_trajectory import calculate_frame_step
from mddb_workflow.utils.auxiliar import save_json
from mddb_workflow.utils.constants import OUTPUT_CHANNELS_FILENAME
from mddb_workflow.utils.type_hints import *
import re
import numpy as np

def channels (
    structure_file : 'File',
    trajectory_file : 'File',
    output_directory : str,
    membrane_map: dict,
    snapshots: int,
    frames_limit: int):

    if membrane_map is None or membrane_map['n_mems'] == 0:
        print('-> Skipping channels analysis')
        return
    
    frame_step, frame_count = calculate_frame_step(snapshots, frames_limit)
    
    prop = {
        'select': 'protein',
        'steps': frame_step,
        'sample': 0.2,
        'dot_density': 13,
        'disable_logs': True,
        'disable_sandbox': True,
    }

    mda_hole(input_top_path=structure_file.path,
            input_traj_path=trajectory_file.path,
            output_hole_path=output_directory+'/hole.vmd',
            output_csv_path=output_directory+'/mda.hole_profile.csv',
            properties=prop)
    
    with open(output_directory+'/hole.vmd', 'r') as f:
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
        cols = [0,0,0]
        z_range = 0
        for i in range(3): # RGB
            # Get triangle coordinates for this frame and color
            # Red are low radiues, green medium, blue high
            trinorms_cl = np.array(trinorms[frame][i]) # (N, 6, 3) # 6: 3 positions + 3 normals, 3: vertex
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
