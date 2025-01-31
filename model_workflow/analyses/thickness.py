from biobb_mem.lipyphilic_biobb.lpp_zpositions import lpp_zpositions, frame_df
from model_workflow.tools.get_pytraj_trajectory import get_reduced_pytraj_trajectory
from model_workflow.utils.auxiliar import save_json
from model_workflow.utils.type_hints import *
import os


def thickness (
    input_structure_filepath : str,
    input_trajectory_filepath : str,
    output_analysis_filepath : str,
    membrane_map: dict,
    snapshots : int,
    frames_limit: int = 1000):
    print('-> Running thickness analysis')

    if membrane_map['n_mems'] == 0:
        # Do something special for density analysis for 
        # membranes like leaflets separation, polargroups, etc.
        print(' No membranes found in the structure. Skipping thickness analysis.')
        return

    # Load. Only used to get the frame step. TODO: extract the frame step on a separate function
    tj, frame_step, frames_count = get_reduced_pytraj_trajectory(input_structure_filepath, input_trajectory_filepath, snapshots, frames_limit)

    head_sel = membrane_map['polar_atoms']
    head_sel_mda = 'index ' + " ".join(map(str,(head_sel)))
    # Run the analysis on the whole membrane
    prop = {
        'lipid_sel': head_sel_mda,
        'steps': frame_step,
        'height_sel': head_sel_mda,
        'ignore_no_box': True,
        'disable_logs': True
    }
    lpp_zpositions(input_top_path=input_structure_filepath,
                   input_traj_path=input_trajectory_filepath,
                   output_positions_path='.zpositions.csv',
                   properties=prop)
    df = frame_df('.zpositions.csv') # Per frame data
    os.remove('.zpositions.csv')

    # Calculate the mean z position of midplane wrt the box axis using pytraj
    midplane_z = []
    for i in range(0, frames_count, frame_step):
        frame = tj[i]
        selected_atoms = frame[head_sel]
        mean_z = selected_atoms.mean(axis=0)[2]
        midplane_z.append(float(mean_z))
    # Save the data
    data = {
        'frame': df.index.tolist(),
        'mean_positive': df['mean_positive'].tolist(),
        'mean_negative': df['mean_negative'].tolist(),
        'std_positive': df['std_positive'].tolist(),
        'std_negative': df['std_negative'].tolist(),
        'thickness': df['thickness'].tolist(),
        'std_thickness': df['std_thickness'].tolist(),
        'midplane_z': midplane_z,

    }
    save_json(data, output_analysis_filepath)