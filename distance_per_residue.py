# Distance per residue analysis
# 
# Perform the distance per residue analysis between each pair of interface residues
# The analysis is carried by pytraj

import pytraj as pt
import numpy

import json

# Perform a distance per residue over all residues from different interface agents
# Note that the distances are calculated for all residues in the agent, not only the interface residues
def distance_per_residue (
    pt_trajectory,
    output_analysis_filename : str,
    interfaces : list ):

    # Return before doing anything if there are no interfaces
    if len(interfaces) == 0:
        return

    output_analysis = []
    for interface in interfaces:

        # Contact Matrix -- Initialization
        # Create 2 lists filled with 0s with the length of the residue number arrays respectively
        h,w = len(interface['residues_2']), len(interface['residues_1'])
        means_matrix = [[0 for x in range(w)] for y in range(h)]
        stdvs_matrix = [[0 for x in range(w)] for y in range(h)]

        # Contact Matrix -- Calculation
        for i, r2 in enumerate(interface['pt_residues_2']):
            for j, r1 in enumerate(interface['pt_residues_1']):
                txt = ":" + str(r2) + " " + ":" + str(r1)
                ptdist = pt.distance(pt_trajectory, txt)
                mean = numpy.mean(ptdist)
                stdv = numpy.std(ptdist)
                means_matrix[i][j] = mean
                stdvs_matrix[i][j] = stdv

        output = {
            'name': interface['name'],
            'means': means_matrix,
            'stdvs': stdvs_matrix,
        }
        output_analysis.append(output)

    # Export the analysis in json format
    with open(output_analysis_filename, 'w') as file:
        json.dump({ 'data': output_analysis }, file)