# Generic analyses
# Easy and fast trajectory analyses carried by Gromacs

from numpy import mean, std

from mddb_workflow.utils.auxiliar import save_json
from mddb_workflow.utils.constants import OUTPUT_RMSF_FILENAME
from mddb_workflow.utils.type_hints import *

from MDAnalysis import Universe
from MDAnalysis.analysis import rms


# LORE: Fluctuation analysis was done with Gromacs before
# LORE: However Gromacs required masses and atom radii, which was a problem in coarse grain
def rmsf (
    universe : 'Universe',
    output_directory : str,
    pbc_selection : 'Selection'):
    """Perform the fluctuation analysis."""
    # Set the main output filepath
    output_analysis_filepath = f'{output_directory}/{OUTPUT_RMSF_FILENAME}'

    # Run the analysis in all atoms
    # Even atoms to be excluded (e.g. PBC) are excluded later to don't mess atom results order
    all_atoms = universe.select_atoms('all')
    rmsf_analysis = rms.RMSF(all_atoms).run()
    # Make a list of it because it is a ndarray, which is not JSON serializable
    # Access to the values thorugh 'results.rmsf' instead of a direct 'rmsf' to avoid a deprecation warning
    rmsf_values = list(rmsf_analysis.results.rmsf)

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
        },
        'version': '0.1.0'
    }

    # Export formatted data to a json file
    save_json(rmsf_data, output_analysis_filepath)