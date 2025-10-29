# Distance per residue analysis
# 
# Perform the distance per residue analysis between each pair of interacting agents
# The analysis is carried by pytraj

from os.path import exists

import pytraj as pt
import numpy as np

from mddb_workflow.utils.auxiliar import save_json, get_analysis_name, numerate_filename, reprint, warn
from mddb_workflow.utils.constants import OUTPUT_DIST_PERRES_FILENAME
from mddb_workflow.utils.pyt_spells import get_reduced_pytraj_trajectory
from mddb_workflow.utils.type_hints import *

# Set a limit of values that can be stored in the analysis without exceeding the mongo limit of 16Mb
# Also, this limit is a good reference to avoid loading huge analyses which would lead to long response time in the client
# This is an aproximation below the limit which has been observed experimentally
N_VALUES_LIMIT = 400000

def distance_per_residue (
    structure_file : 'File',
    trajectory_file : 'File',
    output_directory : str,
    structure : 'Structure',
    interactions : list,
    snapshots : int,
    frames_limit : int
):
    """Calculate the distance mean and standard deviation of each pair of residues of different agents.
    Note that the distances are calculated for all residues in the agent, not only the interface residues."""
    
    # Return before doing anything if there are no interactions
    if not interactions or len(interactions) == 0:
        print('No interactions were specified')
        return
    
    # Set the main output filepath
    output_analysis_filepath = f'{output_directory}/{OUTPUT_DIST_PERRES_FILENAME}'
    
    # Set a reference system to handle conversions to pytraj residue numeration
    # First set the pytraj topology
    pytraj_topology = structure.get_pytraj_topology()
    pytraj_residues = list(pytraj_topology.residues)
    # Transform a structure residue to the pytraj residue numeration (1, 2, ... n)
    def residue_2_pytraj_residue_index (residue_index : int) -> int:
        residue = structure.residues[residue_index]
        residue_number = residue.number
        residue_name = residue.name[0:3]
        # And check that this residue data matches the pytraj residues data
        pytraj_residue = pytraj_residues[residue_index]
        if (residue_number == pytraj_residue.original_resid and residue_name == pytraj_residue.name):
            return residue_index + 1
        # If not, we must iterate over all pytraj residues to find a match
        for index, pytraj_residue in enumerate(pytraj_residues):
            if (residue_number == pytraj_residue.original_resid and residue_name == pytraj_residue.name):
                return index + 1
        # Return None if there is no match
        return None

    # Parse the trajectory intro ptraj
    # Reduce it in case it exceeds the frames limit
    pt_trajectory, frame_step, frames_count = get_reduced_pytraj_trajectory(
        structure_file.path, trajectory_file.path, snapshots, frames_limit)
    # Save each analysis to a dict which will be parsed to json
    # Here we keep track of the summary, which will also be parsed to json at the end
    output_summary = []
    # Run the analysis for each interaction
    print()
    for n, interaction in enumerate(interactions):
        # Get the interaction name
        name = interaction['name']
        # Set a filename for the current interaction data
        numbered_output_analysis_filepath = numerate_filename(output_analysis_filepath, n)
        # Add the root of the output analysis filename to the run data
        analysis_name = get_analysis_name(numbered_output_analysis_filepath)
        # Append current interaction to the summary
        output_summary.append({
            'name': name,
            'analysis': analysis_name
        })
        # If the file already exists then skip to the next
        if exists(numbered_output_analysis_filepath): continue
        reprint(f' Analyzing {name} ({n+1}/{len(interactions)})')
        # Get the residues to be used for this interaction
        residues_1 = interaction['residue_indices_1']
        residues_2 = interaction['residue_indices_2']
        # First of all, calculate the number of values for this interaction
        # If the interaction has more values than we can store then it will be reduced
        # Reduced interactions will store only residues in the interface in order to fit
        n_values = len(residues_1) * len(residues_2)
        reduced = n_values > N_VALUES_LIMIT        
        # Contact Matrix -- Initialization
        # Create 2 lists filled with 0s with the length of the residue number arrays respectively
        if reduced:
            warn('The analysis has been reduced to interface residues only for this interaction')
            residues_1 = interaction['interface_residue_indices_1']
            residues_2 = interaction['interface_residue_indices_2']
            n_values = len(residues_1) * len(residues_2)
            if n_values > N_VALUES_LIMIT: raise ValueError('Too many values, even after reducing')
        # Show the size of the matrix
        h,w = len(residues_2), len(residues_1)
        print(f'     {h}x{w} residues')
        # Convert residues to pytraj residue numbers
        pt_residues_1 = list(map(residue_2_pytraj_residue_index, residues_1))
        pt_residues_2 = list(map(residue_2_pytraj_residue_index, residues_2))
        
        # Contact Matrix -- Calculation (Vectorized)
        # Create all mask pairs at once
        mask_pairs = []
        for r2 in pt_residues_2:
            for r1 in pt_residues_1:
                mask_pairs.append(f":{r2} :{r1}")
        
        # Calculate all distances in one call
        all_distances = pt.distance(pt_trajectory, mask_pairs)
        
        # Reshape the results into matrices
        # all_distances shape: (n_pairs, n_frames)
        means_matrix = np.zeros((h, w))
        stdvs_matrix = np.zeros((h, w))

        pair_idx = 0
        for i in range(h):
            for j in range(w):
                distances = all_distances[pair_idx]
                means_matrix[i][j] = np.mean(distances)
                stdvs_matrix[i][j] = np.std(distances)
                pair_idx += 1
        # Set the output data for this interaction
        save_json({
            'name': name,
            'means': means_matrix.tolist(),
            'stdvs': stdvs_matrix.tolist(),
        }, numbered_output_analysis_filepath)
    # Export the analysis in json format
    save_json(output_summary, output_analysis_filepath)