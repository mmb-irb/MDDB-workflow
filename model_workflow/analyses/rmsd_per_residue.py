# RMSD per resiude analysis
# 
# Perform the RMSD analysis for each residue
# The analysis is carried by pytraj

import pytraj as pt
import re

import json

from model_workflow.tools.get_pytraj_trajectory import get_reduced_pytraj_trajectory

# The pytraj trajectory may be reduced
def rmsd_per_residue (
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_analysis_filename : str,
    topology_reference,
    frames_limit : int):

    # Parse the trajectory intro ptraj
    # Reduce it in case it exceeds the frames limit
    pt_trajectory = get_reduced_pytraj_trajectory(input_topology_filename, input_trajectory_filename, frames_limit)
    
    # Run the analysis in pytraj
    # The result data is a custom pytraj class: pytraj.datasets.datasetlist.DatasetList
    # This class has keys but its attributes can not be accessed through the key
    # They must be accessed thorugh the index
    # DANI: Esto devuelve "Error: Range::SetRange(None): Range is -1 for None"
    # DANI: No se por que pasa pero aparentemente funciona bien
    data = pt.rmsd_perres(pt_trajectory)

    # We remove the first result, which is meant to be the whole rmsd and whose key is 'RMSD_00001'
    del data[0]

    # Mine the analysis data
    output_analysis = []
    for residue in data:
        # Convert pytraj residue keys to source notation
        # Key format: SER:1, TYR:2, ...
        match = re.match('(.*):(.*)', residue.key)
        id = match.groups(0)[1]
        source_residue = topology_reference.pytraj2source(int(id))
        # Residue tag may be not found in excepcional cases which Prody fails to handle:
        # e.g. Residues repeated like: 1, 2, 2, 2, 3, ... (e.g. NMA and ACE residues)
        # In this cases we just ignore the residue
        # If this residue MUST be recorded then change the topology manually adding icodes
        # e.g. Residues repeated like: 1, 2, 2A, 2B, 3, ...
        if not source_residue:
            continue
        residue_tag = str(source_residue)
        # Write data to the output file
        # The 'residue' DataArray contains numeric values (rmsds) and it is not JSON serializable
        output_analysis.append(
            {
                'name': residue_tag,
                'rmsds': list(residue),
            }
        )

    # Export the analysis in json format
    with open(output_analysis_filename, 'w') as file:
        json.dump({ 'data': output_analysis }, file)
    
    # It is not possible to represent the whole rmsd per residue with a classical graph