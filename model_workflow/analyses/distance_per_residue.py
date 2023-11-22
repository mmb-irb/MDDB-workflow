# Distance per residue analysis
# 
# Perform the distance per residue analysis between each pair of interacting agents
# The analysis is carried by pytraj

import pytraj as pt
import numpy

from model_workflow.utils.auxiliar import save_json

from model_workflow.tools.get_pytraj_trajectory import get_reduced_pytraj_trajectory

# Set a limit of values that can be stored in the analysis without exceeding the mongo limit of 16Mb
# Also, this limit is a good reference to avoid loading huge analyses which would lead to long response time in the client
# This is an aproximation below the limit which has been observed experimentally
n_values_limit = 400000

# Calculate the distance mean and standard deviation of each pair of residues*
# * Where each residue is from a different agent
# Note that the distances are calculated for all residues in the agent, not only the interface residues
def distance_per_residue (
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_analysis_filename : str,
    interactions : list,
    snapshots : int,
    frames_limit : int
):
    # Return before doing anything if there are no interactions
    if not interactions or len(interactions) == 0:
        print('No interactions were specified')
        return
    # First of all, calculate the number of values for each interaction
    # If we have more values than we can store then flag some interactions as reduced
    # Flag interactions one by one starting by the biggest one until we are below the limit of values
    # Those interactions will store only residues in the interface in order to reduce the number of values
    interaction_n_values = []
    for interaction in interactions:
        n_values = len(interaction['residues_1']) * len(interaction['residues_2'])
        interaction_n_values.append(n_values)
    reduced_analyses = [ False ] * len(interactions)
    # We must be below the limit in order to proceed
    while sum(interaction_n_values) > n_values_limit:
        # Get the number of values for each interaction which has not been yet reduced
        non_reduced_n_values = [ n_values for i, n_values in enumerate(interaction_n_values) if not reduced_analyses[i] ]
        # If there are not more analysis to reduce then we can not respect the limit
        if len(non_reduced_n_values) == 0:
            # Prepare a small report
            for i, interaction in enumerate(interactions):
                n_values = interaction_n_values[i]
                print('     ' + interaction['name'] + ' -> ' + str(n_values) + ' values')
            raise ValueError('We are still over the limit of values after reducing all analyses')
        # Find the non-reduced interaction with the greatest number of values
        maximum_non_reduced_n_values = max(non_reduced_n_values)
        biggest_interaction_index = interaction_n_values.index(maximum_non_reduced_n_values)
        interaction = interactions[biggest_interaction_index]
        # Calculate the reduced number of values and set the interaction as reduced
        reduced_n_values = len(interaction['interface_1']) * len(interaction['interface_2'])
        interaction_n_values[biggest_interaction_index] = reduced_n_values
        reduced_analyses[biggest_interaction_index] = True
    # Parse the trajectory intro ptraj
    # Reduce it in case it exceeds the frames limit
    pt_trajectory, frame_step, frames_count = get_reduced_pytraj_trajectory(input_topology_filename, input_trajectory_filename, snapshots, frames_limit)
    # Run the analysis for each interaction
    output_analysis = []
    for n, interaction in enumerate(interactions):
        interaction_name = interaction['name']
        print('Running distance per residue analysis for ' + interaction_name)
        # Check if the analysis has been reduced for this interaction
        reduced = reduced_analyses[n]
        # Contact Matrix -- Initialization
        # Create 2 lists filled with 0s with the length of the residue number arrays respectively
        if reduced:
            print('     The analysis has been reduced to interface residues only for this interaction')
            residues_1, residues_2 = interaction['interface_1'], interaction['interface_2']
            pt_residues_1, pt_residues_2 = interaction['pt_interface_1'], interaction['pt_interface_2']
        else:
            residues_1, residues_2 = interaction['residues_1'], interaction['residues_2']
            pt_residues_1, pt_residues_2 = interaction['pt_residues_1'], interaction['pt_residues_2']
        h,w = len(residues_2), len(residues_1)
        print('     ' + str(h) + 'x' + str(w) + ' residues')
        means_matrix = [[0 for x in range(w)] for y in range(h)]
        stdvs_matrix = [[0 for x in range(w)] for y in range(h)]
        # Contact Matrix -- Calculation
        for i, r2 in enumerate(pt_residues_2):
            for j, r1 in enumerate(pt_residues_1):
                txt = ":" + str(r2) + " " + ":" + str(r1)
                ptdist = pt.distance(pt_trajectory, txt)
                mean = numpy.mean(ptdist)
                stdv = numpy.std(ptdist)
                means_matrix[i][j] = mean
                stdvs_matrix[i][j] = stdv
        # Set the output data for this interaction
        output = {
            'name': interaction_name,
            'means': means_matrix,
            'stdvs': stdvs_matrix,
        }
        output_analysis.append(output)
    # Export the analysis in json format
    save_json({ 'data': output_analysis }, output_analysis_filename)