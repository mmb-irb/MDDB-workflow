from itertools import product

import mdtraj as mdt
import pytraj as pt
import numpy as np
# from scipy.special import expit
from sklearn.decomposition import PCA

from mddb_workflow.utils.auxiliar import save_json


def pytraj_residue_pairs(resid_1, resid_2):
    # get list of residue pairs
    # residue_pairs = list(product(*residue_lists))
    residue_pairs = [f":{i} :{j}" for i, j in product(resid_1, resid_2)]
    return residue_pairs


def mdtraj_residue_pairs(resid_1, resid_2):
    residue_pairs = list(product(resid_1, resid_2))
    return residue_pairs


def pytraj_distances(traj, residue_pairs):
    # get pairwise minimal distances
    dists = pt.distance(traj, residue_pairs)
    dists = dists.reshape(-1, dists.shape[0])
    return dists


def mdtraj_distances(traj, residue_pairs):
    dists, _ = mdt.compute_contacts(traj, contacts=residue_pairs)
    return dists


def get_most_frequent_pairs(
        traj,
        dists,
        residue_pairs,
        frequency_threshold,
        distance_threshold):
    # how frequent does the contact need to be: frequency_threshold
    threshold = frequency_threshold * len(traj)
    # set contact definitions: distance_threshold
    close_contacts = np.sum(
        dists < distance_threshold, axis=0)
    # finding all pairs where contact frequency exceeds threshold
    contact_is_frequent = close_contacts > threshold

    # selecting indices for which the above is true
    selection = np.where(contact_is_frequent)[0]

    # select the residue indices pairs corresponding to the selection
    frequent_pairs = np.array(residue_pairs)[selection]

    return frequent_pairs


def sigmoid(x):
    return 1/(1+np.exp(-x))


def pca_contacts(
        trajectory: str,
        topology: str,
        interactions: list,
        output_analysis_filename: str,
        use_pytraj=True,
        n_components=2,
        distance_threshold=15.0,
        frequency_threshold=0.05,
        smooth=5.0):

    print('-> Running PCA contacts analysis')

    # Return before doing anything if there are no interactions
    if len(interactions) == 0:
        return

    output_analysis = []
    for interaction in interactions:
        # DANI: Estos campos ya no están en interactions
        # DANI: Se pueden recuperar tal y como se hace en distance_per_residue
        # DANI: No lo hice en su día porque este análisis nunca se ha llegado a usar
        residue_lists = (interaction["pt_residues_1"],
                         interaction["pt_residues_2"])
        # get list of residue pairs and pairwise minimal distances
        if use_pytraj:
            traj = pt.load(trajectory, top=topology)
            residue_pairs = pytraj_residue_pairs(*residue_lists)
            dists = pytraj_distances(traj, residue_pairs)
        else:
            traj = mdt.load(trajectory, top=topology)
            residue_pairs = mdtraj_residue_pairs(*residue_lists)
            dists = mdtraj_distances(traj, residue_pairs)

        frequent_pairs = get_most_frequent_pairs(
            traj,
            dists,
            residue_pairs,
            frequency_threshold,
            distance_threshold)

        # compute new only with frequent pairs distances
        if use_pytraj:
            dists = pytraj_distances(traj, frequent_pairs)
        else:
            dists = mdtraj_distances(traj, frequent_pairs)

        # smooth distances
        smooth_distances = 1 - sigmoid(smooth * (dists - distance_threshold))

        # compute PCA
        pca = PCA(n_components=n_components)
        transformed = pca.fit_transform(smooth_distances)

        # most important pairs indices ordered in descending order
        pca_components_argsort = pca.components_.argsort()
        most_important_contacts = [frequent_pairs[pca_components_argsort[i][::-1]]
                                   for i in range(n_components)]
        most_important_contacts = [tuple(map(int, s.replace(":", "").split()))
                                   for arr in most_important_contacts
                                   for s in arr]

        # writing data
        # the transformed distances can be plotted individually
        # as a function of the trajectory time (frames)
        transformed_dists = {
            f"transformed_dist_{i+1}": list(transformed.T[i])
            for i in range(n_components)}

        # each of the components (axes in the feature space)
        # sorted in descending order
        components_values = {
            f"component_{i+1}": list(pca.components_[i][pca_components_argsort[i]])
            for i in range(n_components)}

        # residue pairs sorted according to the corresponding components
        ordered_residues = {
            f"ordered_residues_{i+1}": list(most_important_contacts[i])
            for i in range(n_components)}

        # residues that are part of both interaction groups used for the analysis
        interaction_residues = {
            f"interaction_residues": residue_lists}

        output_analysis.append(transformed_dists)
        output_analysis.append(components_values)
        output_analysis.append(ordered_residues)
        output_analysis.append(interaction_residues)

    if output_analysis:
        # Export the analysis in json format
        save_json(output_analysis, output_analysis_filename)
