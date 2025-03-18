# Generic analyses
# Easy and fast trajectory analyses carried by Gromacs

from subprocess import run, PIPE, Popen
from os.path import exists
from os import remove
from numpy import mean, std

from model_workflow.tools.xvg_parse import xvg_parse
from model_workflow.utils.auxiliar import save_json
from model_workflow.utils.constants import GROMACS_EXECUTABLE
from model_workflow.utils.type_hints import *

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
    pbc_selection : 'Selection'):

    print('-> Running RMSF analysis')
    
    # Run Gromacs
    p = Popen([
        "echo",
        "System",
    ], stdout=PIPE)
    process = run([
        GROMACS_EXECUTABLE,
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
    ], stdin=p.stdout, stdout=PIPE, stderr=PIPE)
    logs = process.stdout.decode()
    p.stdout.close()

    # If the output does not exist at this point it means something went wrong with gromacs
    if not exists(rmsf_data_filename):
        print(logs)
        error_logs = process.stderr.decode()
        print(error_logs)
        raise SystemExit('Something went wrong with GROMACS')

    # Read the output file and parse it
    raw_rmsf_data = xvg_parse(rmsf_data_filename, ['atom', 'rmsf'])
    rmsf_values = raw_rmsf_data['rmsf']

    # Filter out values from PBC residue atoms since they may have not sense
    for index in pbc_selection.atom_indices:
        rmsf_values[index] = None

    # Get all rmsf values which are not None
    actual_rmsf_values = [ v for v in rmsf_values if v != None ]
    # There may be none if the whole system is in PBC
    if len(actual_rmsf_values) == 0:
        print(' No actual values to do RMSF')
        return

    # Format data
    rmsf_data = {
        'y': {
            'rmsf': {
                'average': mean(actual_rmsf_values),
                'stddev': std(actual_rmsf_values),
                'min': min(actual_rmsf_values),
                'max': max(actual_rmsf_values),
                'data': rmsf_values # Keep all values here to make the list length match the number of atoms
            }
        }
    }

    # Export formatted data to a json file
    save_json(rmsf_data, output_analysis_filename)

    # Cleanup both the auxiliar and the residual files
    remove(rmsf_data_filename)
    remove(output_noelem_filename)