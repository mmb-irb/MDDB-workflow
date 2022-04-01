# RMSD pairwise analysis
#
# Perform the RMSD analysis for pair of frames in the trajectory
# The analysis is carried by pytraj

import pytraj as pt
import re

import json

from model_workflow.tools.get_pytraj_trajectory import get_reduced_pytraj_trajectory

# Perform an analysis for the overall structure and then one more analysis for each interaction
# The 'interactions' input is mandatory but it may be an empty list (i.e. there are no interactions)
# The pytraj trajectory ('pt_trajectory') may be reduced


def rmsd_pairwise(
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_analysis_filename : str,
    interactions : list,
    frames_limit : int):

    # Parse the trajectory intro ptraj
    # Reduce it in case it exceeds the frames limit
    pt_trajectory = get_reduced_pytraj_trajectory(input_topology_filename, input_trajectory_filename, frames_limit)

    # Run the analysis
    overall_selection = '@CA'
    data = pt.pairwise_rmsd(pt_trajectory, overall_selection)
    # Convert data to a normal list, since numpy ndarrays are not json serializable
    data = data.tolist()
    # In case there is no data we stop here
    if len(data) == 0:
        raise SystemExit('The rmsd-pairwise overall selection (' + overall_selection + ') is empty')

    # Set the final structure data
    output_analysis = [
        {
            'name': 'Overall',
            'rmsds': data,
        }
    ]

    # Repeat the analysis with the interface residues of each interaction
    for interaction in interactions:

        # Select all interface residues in pytraj notation
        pt_interface = interaction['pt_interface_1'] + \
            interaction['pt_interface_2']
        pt_selection = ':' + ','.join(map(str, pt_interface))

        # Run the analysis
        data = pt.pairwise_rmsd(pt_trajectory, pt_selection)
        # Convert data to a normal list, since numpy ndarrays are not json serializable
        data = data.tolist()

        output_analysis.append(
            {
                'name': interaction['name'],
                'rmsds': data,
            }
        )

    # Write in the analysis the starting frame and the step between frames
    # By default the first frame in the reduced trajectory is the first frame (0)
    start = 0
    step = pt_trajectory.step

    # Export the analysis in json format
    with open(output_analysis_filename, 'w') as file:
        json.dump({'data': output_analysis, 'start': start, 'step': step}, file)
