from model_workflow.analyses.generic_analyses import rmsd
from model_workflow.tools.xvg_parse import xvg_parse

from subprocess import run, PIPE
import json

# Run multiple RMSD analyses
# A RMSD analysis is run with each reference:
# - First frame
# - Average structure
# A RMSD analysis is run over each rmsd target:
# - Whole protein
# - Heavy atoms
# - Backbone
# - Alpha carbons
def rmsds(
    input_trajectory_filename : str,
    output_analysis_filename : str,
    first_frame_filename : str,
    average_structure_filename : str,
    rmsd_groups : list = ['Protein', 'Protein-H', 'Backbone', 'C-alpha'],
    ):

    rmsd_references = [first_frame_filename, average_structure_filename]

    output_analysis = []
    start = 0
    # The step is always 1 since we are not using a reduced trajectory
    # WARNING: This is the frames step, no the time step. The time step is calculated in the client
    # WARNING: Do never mine the time step from a xvg file, since gromacs may have wrong times
    step = 1

    # Iterate over each reference and group
    for reference in rmsd_references:
        # Get a standarized reference name
        reference_name = reference[0:-4].lower()
        for group in rmsd_groups:
            # Get a standarized group name
            group_name = group.lower()
            # Set the analysis filename
            rmsd_analysis = 'rmsd.' + reference_name + '.' + group_name + '.xvg'
            # Run the rmsd
            rmsd(reference, input_trajectory_filename, group, rmsd_analysis)
            # Read and parse the output file
            rmsd_data = xvg_parse(rmsd_analysis)
            # Format the mined data and append it to the overall output
            # Multiply by 10 since rmsd comes in nanometers (nm) and we want it in Ångstroms (Å)
            rmsd_values = [ v*10 for v in rmsd_data['values'] ]
            data = {
                'values': rmsd_values,
                'reference': reference_name,
                'group': group_name
            }
            output_analysis.append(data)
            # Remove the analysis xvg file since it is not required anymore
            run([
                "rm",
                rmsd_analysis,
            ], stdout=PIPE).stdout.decode()

    # Export the analysis in json format
    with open(output_analysis_filename, 'w') as file:
        json.dump({ 'start': start, 'step': step, 'data': output_analysis }, file)