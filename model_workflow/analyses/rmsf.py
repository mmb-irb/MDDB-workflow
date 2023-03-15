# Generic analyses
# Easy and fast trajectory analyses carried by Gromacs

from subprocess import run, PIPE, Popen
from os.path import exists
from os import remove
from numpy import mean, std
from json import dump

from model_workflow.tools.xvg_parse import xvg_parse

# Set an auxiliar data filename
rmsf_data_filename = '.rmsf_data.xvg'

# Set a residual data filename
# This analysis produces a 'noelem' file which is never used and thus removed
output_noelem_filename = '.noelem.pdb'

# Fluctuation
# 
# Perform the fluctuation analysis
def rmsf (
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_analysis_filename : str,
):
    
    # Run Gromacs
    p = Popen([
        "echo",
        "System",
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "rmsf",
        "-s",
        input_topology_filename,
        "-f",
        input_trajectory_filename,
        '-o',
        rmsf_data_filename,
        '-oq',
        output_noelem_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE, stderr=PIPE).stdout.decode()
    p.stdout.close()

    # If the output does not exist at this point it means something went wrong with gromacs
    if not exists(rmsf_data_filename):
        print(logs)
        raise SystemExit('Something went wrong with GROMACS')

    # Read the output file and parse it
    raw_rmsf_data = xvg_parse(rmsf_data_filename, ['atom', 'rmsf'])

    # Cleanup both the auxiliar and the residual files
    remove(rmsf_data_filename)
    remove(output_noelem_filename)

    # Format data
    rmsf_data = {
        'y': {
            'rmsf': {
                'average': mean(raw_rmsf_data['rmsf']),
                'stddev': std(raw_rmsf_data['rmsf']),
                'min': min(raw_rmsf_data['rmsf']),
                'max': max(raw_rmsf_data['rmsf']),
                'data': raw_rmsf_data['rmsf']
            }
        }
    }

    # Export formatted data to a json file
    with open(output_analysis_filename, 'w') as file:
        dump(rmsf_data, file)