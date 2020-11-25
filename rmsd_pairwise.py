# RMSD pairwise analysis
# 
# Perform the RMSD analysis for pair of frames in the trajectory
# The analysis is carried by pytraj

import pytraj as pt
import re

import json

# The pytraj trajectory may be reduced
# DANI: Todavía no está acabado
def rmsd_pairwise (
    pytraj_trajectory,
    output_analysis_filename : str,
    interfaces ):

    # Run the analysis
    data = pt.pairwise_rmsd(pytraj_trajectory, '@CA')
    # Convert data to a normal list, since numpy ndarrays are not json serializable
    data = data.tolist()

    # Set the final structure data 
    output_analysis = [
        {
            'name': 'Overall',
            'data': data,
        }
    ]

    # Repeat the analysis for each seletion of interface residues
    for interface in interfaces:
        # Hay que ver en que formato vienen las interfaces
        print(interface)


    return
    # Export the analysis in json format
    with open(output_analysis_filename, 'w') as file:
        json.dump(output_analysis, file)
    