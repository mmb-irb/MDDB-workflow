# Distance per residue analysis
# 
# Perform the distance per residue analysis between each pair of interface residues
# The analysis is carried by pytraj

import pytraj as pt
import numpy

import json

# Perform an hydrogen bonds analysis for each interface
# The 'interfaces' input may be an empty list (i.e. there are no interfaces)
# In case there are no interfaces the analysis stops
def distance_per_residue (
    pt_trajectory,
    output_analysis_filename : str,
    interfaces : list ):

    output_analysis = []
    for interface in interfaces:

        # Contact Matrix -- Initialization
        # Create 2 lists filled with 0s with the length of the residue number arrays respectively
        h,w = len(interface['residues_2']), len(interface['residues_1'])
        mat_mean = [[0 for x in range(w)] for y in range(h)]
        mat_stdv = [[0 for x in range(w)] for y in range(h)]

        # Contact Matrix -- Calculation
        for i, r2 in enumerate(interface['pt_residues_2']):
            for j, r1 in enumerate(interface['pt_residues_1']):
                txt = ":" + str(r2) + " " + ":" + str(r1)
                ptdist = pt.distance(pt_trajectory, txt)
                mean = numpy.mean(ptdist)
                stdv = numpy.std(ptdist)
                mat_mean[i][j] = mean
                mat_stdv[i][j] = stdv

        output = {
            'name': interface['name'],
            'xLabels': list(map(str, interface['residues_1'])),
            'yLabels': list(map(str, interface['residues_2'])),
            'mean': mat_mean,
            'stdv': mat_stdv,
        }
        output_analysis.append(output)

    # Export the analysis in json format
    with open(output_analysis_filename, 'w') as file:
        json.dump({ 'data': output_analysis }, file)