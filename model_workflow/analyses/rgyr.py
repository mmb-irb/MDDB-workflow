from os.path import exists
from os import remove
from subprocess import run, PIPE, Popen
from numpy import mean, std

from model_workflow.tools.get_reduced_trajectory import get_reduced_trajectory
from model_workflow.tools.xvg_parse import xvg_parse
from model_workflow.utils.auxiliar import save_json
from model_workflow.utils.constants import GROMACS_EXECUTABLE
from model_workflow.utils.type_hints import *

# Set an auxiliar data filename
rgyr_data_filename = '.rgyr_data.xvg'

# Radius of gyration (Rgyr)
# 
# Perform the RMSd analysis 
# Use the first trajectory frame in .pdb format as a reference
def rgyr (
    input_topology_file : 'File',
    input_trajectory_file : 'File',
    output_analysis_filepath : str,
    snapshots : int,
    frames_limit : int,
    structure : 'Structure',
    pbc_residues : List[int]):

    print('-> Running RGYR analysis')

    # Use a reduced trajectory in case the original trajectory has many frames
    reduced_trajectory_filepath, step, frames = get_reduced_trajectory(
        input_topology_file,
        input_trajectory_file,
        snapshots,
        frames_limit,
    )

    # Generate a custom index file to exclude PBC residues from the analysis
    pbc_selection = structure.select_residue_indices(pbc_residues)
    not_pbc_selection = structure.invert_selection(pbc_selection)
    if not not_pbc_selection:
        print(' No selection to run the analysis')
        return
    selection_name = 'pbc_residues'
    ndx_selection = not_pbc_selection.to_ndx(selection_name)
    ndx_filename = '.not_pbc.ndx'
    with open(ndx_filename, 'w') as file:
        file.write(ndx_selection)
    
    # Run Gromacs
    p = Popen([
        "echo",
        selection_name,
    ], stdout=PIPE)
    process = run([
        GROMACS_EXECUTABLE,
        "gyrate",
        "-s",
        input_topology_file.path,
        "-f",
        reduced_trajectory_filepath,
        '-o',
        rgyr_data_filename,
        '-n',
        ndx_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE, stderr=PIPE)
    logs = process.stdout.decode()
    p.stdout.close()

    # If the output does not exist at this point it means something went wrong with gromacs
    if not exists(rgyr_data_filename):
        print(logs)
        error_logs = process.stderr.decode()
        print(error_logs)
        raise SystemExit('Something went wrong with GROMACS')

    # Read the output file and parse it
    raw_rgyr_data = xvg_parse(rgyr_data_filename, ['times', 'rgyr', 'rgyrx', 'rgyry', 'rgyrz'])

    # Format data
    rgyr_data = {
        'start': 0,
        'step': step,
        'y': {
            'rgyr': {
                'average': mean(raw_rgyr_data['rgyr']),
                'stddev': std(raw_rgyr_data['rgyr']),
                'min': min(raw_rgyr_data['rgyr']),
                'max': max(raw_rgyr_data['rgyr']),
                'data': raw_rgyr_data['rgyr']
            },
            'rgyrx': {
                'average': mean(raw_rgyr_data['rgyrx']),
                'stddev': std(raw_rgyr_data['rgyrx']),
                'min': min(raw_rgyr_data['rgyrx']),
                'max': max(raw_rgyr_data['rgyrx']),
                'data': raw_rgyr_data['rgyrx']
            },
            'rgyry': {
                'average': mean(raw_rgyr_data['rgyry']),
                'stddev': std(raw_rgyr_data['rgyry']),
                'min': min(raw_rgyr_data['rgyry']),
                'max': max(raw_rgyr_data['rgyry']),
                'data': raw_rgyr_data['rgyry']
            },
            'rgyrz': {
                'average': mean(raw_rgyr_data['rgyrz']),
                'stddev': std(raw_rgyr_data['rgyrz']),
                'min': min(raw_rgyr_data['rgyrz']),
                'max': max(raw_rgyr_data['rgyrz']),
                'data': raw_rgyr_data['rgyrz']
            }
        }
    }

    # Export formatted data to a json file
    save_json(rgyr_data, output_analysis_filepath)

    # Remove residual files
    remove(rgyr_data_filename)
    remove(ndx_filename)