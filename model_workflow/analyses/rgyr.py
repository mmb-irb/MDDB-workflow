from os.path import exists
from os import remove
from subprocess import run, PIPE, Popen
from numpy import mean, std
from json import dump
from model_workflow.tools.get_reduced_trajectory import get_reduced_trajectory
from model_workflow.tools.xvg_parse import xvg_parse

# Set an auxiliar data filename
rgyr_data_filename = '.rgyr_data.xvg'

# Radius of gyration (Rgyr)
# 
# Perform the RMSd analysis 
# Use the first trajectory frame in .pdb format as a reference
def rgyr (
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_analysis_filename : str,
    snapshots : int,
    frames_limit : int):

    # Use a reduced trajectory in case the original trajectory has many frames
    reduced_trajectory_filename, step, frames = get_reduced_trajectory(
        input_topology_filename,
        input_trajectory_filename,
        snapshots,
        frames_limit,
    )
    
    # Run Gromacs
    p = Popen([
        "echo",
        "System",
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "gyrate",
        "-s",
        input_topology_filename,
        "-f",
        reduced_trajectory_filename,
        '-o',
        rgyr_data_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE, stderr=PIPE).stdout.decode()
    p.stdout.close()

    # If the output does not exist at this point it means something went wrong with gromacs
    if not exists(rgyr_data_filename):
        print(logs)
        raise SystemExit('Something went wrong with GROMACS')

    # Read the output file and parse it
    raw_rgyr_data = xvg_parse(rgyr_data_filename, ['times', 'rgyr', 'rgyrx', 'rgyry', 'rgyrz'])

    # Cleanup the auxiliar file
    remove(rgyr_data_filename)

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
    with open(output_analysis_filename, 'w') as file:
        dump(rgyr_data, file)