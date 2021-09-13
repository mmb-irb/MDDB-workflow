# Distance per residue analysis
# 
# Perform the distance per residue analysis between each pair of interacting agents
# The analysis is carried by pytraj

import pytraj as pt
import numpy

import json

from model_workflow.tools.get_pytraj_trajectory import get_reduced_pytraj_trajectory

# Calculate the distance mean and standard deviation of each pair of residues*
# * Where each residue is from a different agent
# Note that the distances are calculated for all residues in the agent, not only the interface residues
def distance_per_residue (
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_analysis_filename : str,
    interactions : list,
    frames_limit : int):

    # Return before doing anything if there are no interactions
    if not interactions or len(interactions) == 0:
        print('No interactions were specified')
        return

    # Parse the trajectory intro ptraj
    # Reduce it in case it exceeds the frames limit
    pt_trajectory = get_reduced_pytraj_trajectory(input_topology_filename, input_trajectory_filename, frames_limit)

    output_analysis = []
    for interaction in interactions:

        # Contact Matrix -- Initialization
        # Create 2 lists filled with 0s with the length of the residue number arrays respectively
        h,w = len(interaction['residues_2']), len(interaction['residues_1'])
        means_matrix = [[0 for x in range(w)] for y in range(h)]
        stdvs_matrix = [[0 for x in range(w)] for y in range(h)]

        # Contact Matrix -- Calculation
        for i, r2 in enumerate(interaction['pt_residues_2']):
            for j, r1 in enumerate(interaction['pt_residues_1']):
                txt = ":" + str(r2) + " " + ":" + str(r1)
                ptdist = pt.distance(pt_trajectory, txt)
                mean = numpy.mean(ptdist)
                stdv = numpy.std(ptdist)
                means_matrix[i][j] = mean
                stdvs_matrix[i][j] = stdv

        output = {
            'name': interaction['name'],
            'means': means_matrix,
            'stdvs': stdvs_matrix,
        }
        output_analysis.append(output)

    # Export the analysis in json format
    with open(output_analysis_filename, 'w') as file:
        json.dump({ 'data': output_analysis }, file)