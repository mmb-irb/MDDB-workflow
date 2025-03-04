from biobb_mem.fatslim.fatslim_apl import fatslim_apl
from model_workflow.utils.auxiliar import save_json
from model_workflow.utils.type_hints import *
from scipy.interpolate import griddata
from contextlib import redirect_stdout
import pandas as pd
import numpy as np
import os


def area_per_lipid (
    input_structure_filepath : str,
    input_trajectory_filepath : str,
    output_analysis_filepath : str,
    membrane_map: dict,):

    if membrane_map is None or membrane_map['n_mems'] == 0:
        print('-> Skipping area per lipid analysis')
        return
    print('-> Running area per lipid analysis')

    head_sel = []
    for n in range(membrane_map['n_mems']):
        head_sel.extend(membrane_map['mems'][str(n)]['polar_atoms']['top'])
        head_sel.extend(membrane_map['mems'][str(n)]['polar_atoms']['bot'])
    head_sel_mda = 'index ' + " ".join(map(str,(head_sel)))
    # Run the analysis on the whole membrane
    prop = {
    'lipid_selection': head_sel_mda,
    'ignore_no_box': True,
    'disable_logs': True,
    }
    apl_tmp = '.apl.csv'
    print(' Running BioBB FATSLiM APL')
    with redirect_stdout(None):
        fatslim_apl(input_top_path=input_structure_filepath,
                    input_traj_path=input_trajectory_filepath,
                    output_csv_path=apl_tmp,
                    properties=prop)
    grids, grid_x, grid_y, m, s = process_apl(apl_tmp)
    os.remove(apl_tmp)
    # Replace NaNs with -1 in the grids so the loader don't break
    grids = [np.nan_to_num(grid, nan=-1) for grid in grids]
    # Save the data
    data = { 'data':{
        'lower leaflet': grids[0].tolist(),
        'upper leaflet': grids[1].tolist(),
        'grid_x': grid_x[:,0].tolist(),
        'grid_y': grid_y[0,:].tolist(),
        'median': m,
        'std': s,
        }
    }
    save_json(data, output_analysis_filepath)


def process_apl(output_csv_path, res=100j):
    df = pd.read_csv(output_csv_path)
    grids = []
    # Create separate plots for each leaflet
    df['Area per lipid'] *= 100  # Convert to A^2
    m = df['Area per lipid'].median()
    s = df['Area per lipid'].std()
    
    # Define common grid for both plots
    x_all = df['X coords']
    y_all = df['Y coords']
    grid_x, grid_y = np.mgrid[min(x_all):max(x_all):res,
                             min(y_all):max(y_all):res]
    for leaflet in ['lower leaflet', 'upper leaflet']:
        df_leaflet = df[df['leaflet'] == leaflet]
        points = np.stack((np.array(df_leaflet['X coords']).T, np.array(df_leaflet['Y coords']).T), axis=-1)
        values = np.array(df_leaflet['Area per lipid'])
        grid = griddata(points, values, (grid_x, grid_y), method='cubic')
        grids.append(grid)
    return grids, grid_x, grid_y, m, s
