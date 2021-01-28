from model_workflow.analyses.generic_analyses import rmsd

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
    trajectory_filename : str,
    output_analysis_filename : str,
    rmsd_references : list,
    rmsd_groups : list = ['Protein', 'Protein-H', 'Backbone', 'C-alpha'],
    ):

    output_analysis = []

    # Iterate over each reference and group
    for reference in rmsd_references:
        for group in rmsd_groups:
            # Get standarized reference and group names
            reference_name = reference[0:-4].lower()
            group_name = group.lower()
            # Set the analysis filename
            rmsd_analysis = 'rmsd.' + reference_name + '.' + group_name + '.xvg'
            # Run the rmsd
            rmsd(reference, trajectory_filename, group, rmsd_analysis)
            # Read and parse the output file
            rmsd_data = xvg_parse(rmsd_analysis)
            # Name the mined data and append it to the overall output
            rmsd_data['reference'] = reference_name
            rmsd_data['group'] = group_name
            output_analysis.append(rmsd_data)

    # Export the analysis in json format
    with open(output_analysis_filename, 'w') as file:
        json.dump({ 'data': output_analysis }, file)

# Read and parse a xvg file
# xvg file example:
# # This is a comment line
# @ This is a comment line
#    0.0000000    0.3644627
#   10.0000000    0.3536768
#   20.0000000    0.3509805
def xvg_parse (filename : str):
    times = []
    values = []
    # Read the specified file line per line
    with open(filename, 'r') as file:
        lines = list(file)
        for line in lines:
            # Skip comment lines
            first_character = line[0]
            if first_character in ['#','@']:
                continue
            # Useful lines are splitted by spaces
            # Splits are saved in 2 columns: times and values
            [time, value] = line.split()
            times.append(float(time))
            values.append(float(value))
    # Finally set the output object
    start = times[0]
    step = times[1] - start
    return { 'start': start, 'step': step, 'values': values }