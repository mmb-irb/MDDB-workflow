from biobb_mem.lipyphilic_biobb.lpp_zpositions import lpp_zpositions, frame_df
from model_workflow.utils.pyt_spells import get_reduced_pytraj_trajectory
from model_workflow.utils.auxiliar import save_json
from model_workflow.utils.constants import OUTPUT_THICKNESS_FILENAME
from model_workflow.utils.type_hints import *
from contextlib import redirect_stdout
import os


def thickness (
    structure_file : 'File',
    trajectory_file : 'File',
    output_directory : str,
    membrane_map: dict,
    snapshots : int,
    frames_limit: int = 100):
    """Membrane thickness analysis."""

    if membrane_map is None or membrane_map['n_mems'] == 0:
        print('-> Skipping thickness analysis')
        return
    
    # Set the main output filepath
    output_analysis_filepath = f'{output_directory}/{OUTPUT_THICKNESS_FILENAME}'

    tj, frame_step, frames_count = get_reduced_pytraj_trajectory(
        structure_file.path, trajectory_file.path, snapshots, frames_limit)
    head_sel = []
    for n in range(membrane_map['n_mems']):
        head_sel.extend(membrane_map['mems'][str(n)]['polar_atoms']['top'])
        head_sel.extend(membrane_map['mems'][str(n)]['polar_atoms']['bot'])
    head_sel_mda = 'index ' + " ".join(map(str,(head_sel)))
    # Run the analysis on the whole membrane
    prop = {
        'lipid_sel': head_sel_mda,
        'steps': frame_step,
        'height_sel': head_sel_mda,
        'ignore_no_box': True,
        'disable_logs': True
    }
    print(' Running BioBB LiPyphilic ZPositions')
    with redirect_stdout(None):
        lpp_zpositions(input_top_path=structure_file.path,
                    input_traj_path=trajectory_file.path,
                    output_positions_path='.zpositions.csv',
                    properties=prop)
    df = frame_df('.zpositions.csv') # Per frame data
    os.remove('.zpositions.csv')

    # Calculate the mean z position of midplane wrt the box axis using pytraj
    midplane_z = []
    for i in range(0, frames_count):
        frame = tj[i]
        selected_atoms = frame[head_sel]
        mean_z = selected_atoms.mean(axis=0)[2]
        midplane_z.append(float(mean_z))
    # Save the data
    data = { 'data':{
        'frame': df.index.tolist(),
        'mean_positive': df['mean_positive'].tolist(),
        'mean_negative': df['mean_negative'].tolist(),
        'std_positive': df['std_positive'].tolist(),
        'std_negative': df['std_negative'].tolist(),
        'thickness': df['thickness'].tolist(),
        'std_thickness': df['std_thickness'].tolist(),
        'midplane_z': midplane_z,
        'step': frame_step
        }
    }
    save_json(data, output_analysis_filepath)