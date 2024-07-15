import numpy as np

import mdtraj as mdt

from model_workflow.utils.auxiliar import round_to_thousandths, save_json, otherwise
from model_workflow.tools.get_screenshot import get_screenshot

AUXILIAR_PDB_FILENAME = '.model.pdb'

# Run the cluster analysis
def clusters_analysis (
    input_structure_filename : str,
    input_trajectory_filename : str,
    interactions : list,
    structure : 'Structure',
    output_analysis_filename : str,
    output_screenshots_filename : str,
    # Set the number of steps between the maximum and minimum RMSD so set how many cutoff are tried and how far they are
    n_steps : int = 100,
    # Set the final amount of desired clusters
    desired_n_clusters : int = 20,
    # Set the atom selection for the overall clustering
    overall_selection : str = "name CA or name C5'",
):

    # Load the whole trajectory
    traj = mdt.load(input_trajectory_filename, top=input_structure_filename)

    # The cluster analysis is run for the overall structure and then once more for every interaction
    # We must set the atom selection of every run in atom indices, for MDtraj
    runs = []

    # Start with the overall selection
    runs.append({
        'name': 'Overall',
        'selection': structure.select(overall_selection)
    })

    # Now setup the interaction runs
    for interaction in interactions:
        interface_residue_indices = interaction['interface_indices_1'] + interaction['interface_indices_2']
        runs.append({
            'name': interaction['name'],
            'selection': structure.select_residue_indices(interface_residue_indices)
        })

    # Now iterate over the different runs
    for r, run in enumerate(runs):
        name = run['name']
        print(f'Calculating distances for {name}')
        # Get the run selection atom indices
        atom_indices = run['selection'].atom_indices
        
        # Calculate the RMSD matrix
        distance_matrix = np.empty((traj.n_frames, traj.n_frames))
        for i in range(traj.n_frames):
            print(f' Frame {i+1} out of {traj.n_frames}', end='\r')
            # Calculate the RMSD between every frame in the trajectory and the frame 'i'
            distance_matrix[i] = mdt.rmsd(traj, traj, i, atom_indices=atom_indices)

        # Get the maximum RMSD value in the whole matrix
        maximum_rmsd = np.max(distance_matrix)
        # Get the minimum RMSD value in the whole matrix
        # Discard 0s from frames against themselves
        minimum_rmsd = np.min(distance_matrix[distance_matrix != 0])

        # Set the difference between the minimum and maximum to determine the cutoffs step
        rmsd_difference = maximum_rmsd - minimum_rmsd
        rmsd_step = rmsd_difference / n_steps

        # Set the initial RMSD cutoff
        cutoff = round_to_thousandths(minimum_rmsd + rmsd_difference / 2)
        # Keep a register of already tried cutoffs so we do not repeat
        already_tried_cutoffs = set()

        # Adjust the RMSD cutoff until we get the desired amount of clusters
        # Note that final clusters will be ordered by the time they appear
        clusters = None
        n_clusters = 0
        while n_clusters != desired_n_clusters:
            # Find clusters
            print(f'Trying with cutoff {cutoff}', end='')
            clusters = clustering(distance_matrix, cutoff)
            n_clusters = len(clusters)
            print(f' -> Found {n_clusters} clusters')
            # Update the cutoff
            already_tried_cutoffs.add(cutoff)
            if n_clusters > desired_n_clusters:
                cutoff = round_to_thousandths(cutoff + rmsd_step)
            if n_clusters < desired_n_clusters:
                cutoff = round_to_thousandths(cutoff - rmsd_step)
            # If we already tried the updated cutoff then we are close enough to the desired number of clusters
            if cutoff in already_tried_cutoffs:
                break

        # Count the number of frames per cluster
        cluster_lengths = [ len(cluster) for cluster in clusters ]

        # Resort clusters in a "cluster per frame" structure
        frame_clusters = np.empty(traj.n_frames, dtype=int)
        for c, cluster in enumerate(clusters):
            for frame in cluster:
                frame_clusters[frame] = c

        # Count the transitions between clusters
        transitions = []

        # Iterate over the different frames
        previous_cluster = frame_clusters[0]
        for cluster in frame_clusters[1:]:
            # If this is the same cluster then there is no transition here
            if previous_cluster == cluster:
                continue
            # Otherwise save the transition
            transition = previous_cluster, cluster
            transitions.append(transition)
            previous_cluster = cluster

        print(f'Found {len(transitions)} transitions')

        # Count every different transition
        transition_counts = {}
        for transition in transitions:
            current_count = transition_counts.get(transition, 0)
            transition_counts[transition] = current_count + 1

        # Now for every cluster find the most representative frame (i.e. the one with less RMSD distance to its neighbours)
        # Then make a screenshot for this specific frame
        representative_frames = []
        # Save the screenshot parameters so we can keep images coherent between clusters
        screenshot_parameters = None
        for c, cluster in enumerate(clusters):
            most_representative_frame = None
            min_distance = float('inf') # Positive infinity
            for frame, neighbour_frames in otherwise(cluster):
                # Calculate the sum of all rmsd distances
                total_distance = 0
                for neighbour_frame in neighbour_frames:
                    total_distance += distance_matrix[frame][neighbour_frame]
                # If the distance is inferior the current minimum then set this frame as the most representative
                if total_distance < min_distance:
                    most_representative_frame = frame
                    min_distance = total_distance
            # Save the most representative frame in the list
            representative_frames.append(most_representative_frame)
            # Once we have the most representative frame we take a screenshot
            # This screenshots will be then uploaded to the database as well
            # Generate a pdb with coordinates from the most representative frame
            coordinates = traj[most_representative_frame]
            coordinates.save(AUXILIAR_PDB_FILENAME)
            # Set the screenshot filename from the input template
            screenshot_filename = output_screenshots_filename.replace('*', str(r).zfill(2)).replace('?', str(c).zfill(2))
            # Generate the screenshot
            screenshot_parameters = get_screenshot(AUXILIAR_PDB_FILENAME, screenshot_filename, parameters=screenshot_parameters)

        # Set the output clusters which include all frames in the cluster and the main or more representative frame
        output_clusters = []
        for frames, most_representative_frame in zip(clusters, representative_frames):
            output_clusters.append({ 'frames': frames, 'main': most_representative_frame })

        # Set the output transitions in a hashable and json parseable way
        output_transitions = []
        for transition, count in transition_counts.items():
            # Set frames as regular ints to make them json serializable
            output_transitions.append({ 'from': int(transition[0]), 'to': int(transition[1]), 'count': count })

        # Set the output analysis
        output_analysis = {
            'subname': name,
            'cutoff': cutoff,
            'clusters': output_clusters,
            'transitions': output_transitions
        }

        # Set the analysis filename from the input template
        analysis_filename = output_analysis_filename.replace('*', str(r).zfill(2))

        # The output filename must be different for every run to avoid overwritting previous results
        # However the filename is not important regarding the database since this analysis is found by its 'run'
        save_json(output_analysis, analysis_filename)

# Set a function to cluster frames in a RMSD matrix given a RMSD cutoff
# https://github.com/boneta/RMSD-Clustering/blob/master/rmsd_clustering/clustering.py
def clustering (rmsd_matrix : np.ndarray, cutoff : float) -> list:
    clusters = []
    for i in range(rmsd_matrix.shape[0]):
        for cluster in clusters:
            if all(rmsd_matrix[i,j] < cutoff for j in cluster):
                cluster.append(i)
                break
        else:
            clusters.append([i])
    return clusters