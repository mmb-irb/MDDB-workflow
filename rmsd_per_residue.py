# RMSD per resiude analysis
# 
# Perform the RMSD analysis for each residue
# The analysis is carried by pytraj

import pytraj as pt
import re

# The pytraj trajectory may be reduced
def rmsd_per_residue (
    pt_trajectory,
    output_analysis_filename : str,
    topology_reference ):
    
    # Run the analysis in pytraj
    # DANI: Esto devuelve "Error: Range::SetRange(None): Range is -1 for None"
    # DANI: No se por que pasa pero aparentemente funciona bien
    data = pt.rmsd_perres(pt_trajectory)
    
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

        for i in range(pt_trajectory.n_frames):
            line = str(pt_trajectory[i].time) + '    '
            for d in data:
                line = line + str(d[i]) + '    '
            file.write(line + "\n")

    # It is not possible to represent the rmsd per residue with a classical graph