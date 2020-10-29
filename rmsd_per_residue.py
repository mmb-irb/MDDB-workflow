# RMSD per resiude analysis
# 
# Perform the RMSD analysis for each residue, which is carried by pytraj

import pytraj as pt
import re

def rmsd_per_residue (
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_analysis_filename : str,
    topology_reference ):

    # Load the trajectory to pytraj and set a reduced trajectory
    pytrajectory = pt.iterload(input_trajectory_filename, input_topology_filename)
    reduced_pytrajectory = pytrajectory[0:2000:10]

    # Run the analysis in pytraj
    data = pt.rmsd_perres(reduced_pytrajectory)

    # Write the output to a new filename in a standarized format
    with open(output_analysis_filename,'w') as file:

        for d in data:
            # Key format: SER:1, TYR:2, ...
            key = d.key
            #print(key)
            tag = key
            match = re.match('(.*):(.*)', key)
            if match != None:
                #print(key)
                id = match.groups(0)[1]
                tag = str(topology_reference.pytraj2source(int(id)))
                
            file.write("@ key " + tag + "\n")

        for i in range(reduced_pytrajectory.n_frames):
            line = str(reduced_pytrajectory[i].time) + '    '
            for d in data:
                line = line + str(d[i]) + '    '
            file.write(line + "\n")

    # It is not possible to represent the rmsd per residue with a classical graph