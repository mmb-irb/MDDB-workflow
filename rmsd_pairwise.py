# RMSD pairwise analysis
# 
# Perform the RMSD analysis for pair of frames in the trajectory
# The analysis is carried by pytraj

import pytraj as pt
import re

import json

# Perform an analysis for the overall structure and then one more analysis for each interface
# The 'interfaces' input is mandatory but it may be an empty list (i.e. there are no interfaces)
# The pytraj trajectory ('pt_trajectory') may be reduced
def rmsd_pairwise (
    pt_trajectory,
    output_analysis_filename : str,
    interfaces : list ):

    # Run the analysis
    data = pt.pairwise_rmsd(pt_trajectory, '@CA')
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
        
        # Select all interface residues in pytraj notation
        pt_interface = interface['pt_interface_1'] + interface['pt_interface_2']
        pt_selection = ':' + ','.join(map(str, pt_interface)) + ' @CA'
        
        # Run the analysis
        data = pt.pairwise_rmsd(pt_trajectory, pt_selection)
        # Convert data to a normal list, since numpy ndarrays are not json serializable
        data = data.tolist()

        output_analysis.append(
            {
                'name': interface['name'],
                'data': data,
            }
        )

    # Export the analysis in json format
    with open(output_analysis_filename, 'w') as file:
        json.dump({ 'data': output_analysis }, file)
    