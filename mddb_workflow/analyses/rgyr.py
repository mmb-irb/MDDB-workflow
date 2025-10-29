from os import remove
from numpy import mean, std

from mddb_workflow.tools.get_reduced_trajectory import get_reduced_trajectory
from mddb_workflow.tools.xvg_parse import xvg_parse
from mddb_workflow.utils.auxiliar import save_json
from mddb_workflow.utils.constants import OUTPUT_RGYR_FILENAME
from mddb_workflow.utils.gmx_spells import run_gromacs
from mddb_workflow.utils.type_hints import *

# Set an auxiliar data filename
rgyr_data_filename = '.rgyr_data.xvg'


def rgyr (
    structure_file : 'File',
    trajectory_file : 'File',
    output_directory : str,
    snapshots : int,
    structure : 'Structure',
    pbc_selection : 'Selection',
    frames_limit : int = 5000
):
    """Perform the RMSd analysis. Use the first trajectory frame in .pdb format as a reference."""

    # Set the main output filepath
    output_analysis_filepath = f'{output_directory}/{OUTPUT_RGYR_FILENAME}'

    # Use a reduced trajectory in case the original trajectory has many frames
    reduced_trajectory_filepath, step, frames = get_reduced_trajectory(
        structure_file,
        trajectory_file,
        snapshots,
        frames_limit,
    )

    # Generate a custom index file to exclude PBC residues from the analysis
    not_pbc_selection = structure.invert_selection(pbc_selection)
    if not not_pbc_selection:
        print(' No selection to run the analysis')
        return
    selection_name = 'pbc_atoms'
    ndx_selection = not_pbc_selection.to_ndx(selection_name)
    ndx_filename = '.not_pbc.ndx'
    with open(ndx_filename, 'w') as file:
        file.write(ndx_selection)
    
    # Run Gromacs
    run_gromacs(f'gyrate -s {structure_file.path} -f {reduced_trajectory_filepath} \
        -o {rgyr_data_filename} -n {ndx_filename}', user_input = selection_name)

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