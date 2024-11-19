import numpy as np
import MDAnalysis
from MDAnalysis.transformations.boxdimensions import set_dimensions
from lipyphilic.lib.neighbours import Neighbours
from lipyphilic.lib.assign_leaflets import AssignLeaflets
from scipy.sparse.csgraph import connected_components
from model_workflow.utils.topology_converter import to_MDAnalysis_topology
from model_workflow.tools.get_inchi_keys import get_inchi_keys, is_in_LIPID_MAPS
from model_workflow.utils.type_hints import *
from typing import List


def generate_membrane_mapping(structure : 'Structure',
                              topology_file : 'File',
                              structure_file : 'File',
                              debug=False
                              ) -> List[dict]:
    """
    Generates a list of residue numbers of membrane components from a given structure and topology file.
    Args:
        structure (Structure): The molecular structure to analyze.
        topology_file (File): The topology file in JSON format.
        structure_file (File): The structure file.
        debug (bool, optional): If True, additional debug information is returned. Defaults to False.
    
    Returns:
        List[dict]: A list containing the membrane mapping. If debug is True, additional information is returned.
    
    Raises:
        AssertionError: If the topology file is not in JSON format.
    
    Notes:
        - The function identifies lipid and non-lipid residues based on InChI keys.
        - It classifies residues and checks for potential misclassifications.
        - Lipid residues are selected and neighboring lipids are found.
        - Clusters of lipids are identified, and clusters with more than 30 lipids are considered as membranes.
        - If debug is enabled, the function returns additional information including lipid residues, neighbors, counts, and clusters.
    """
    print('Calculando la membrana...')
    assert topology_file.extension == 'json', 'Input topology file must be in json format'
    mda_top = to_MDAnalysis_topology(topology_file.absolute_path)
    u = MDAnalysis.Universe(mda_top, structure_file.absolute_path)
    # Get InChI keys of non-proteic/non-nucleic residues
    inchi_keys = get_inchi_keys(u, structure)
    # TO-DO: identify substructures like carbohydrates by neutralizing bonds:
    # ich = inchi_keys['ZNJXAXZHPQQZRY-LXGUWJNJSA-N']['inchi'] # A01IP
    # for args in [[ich], [ich, 1], [ich, 0, 1], [ich, 1, 1]]:
    #    display(inchi_2_mol(*args))
    # Classsify the residues as lipid or not
    lipid, not_lipid = [], []
    for inchikey, data in inchi_keys.items():
        lipid_data = is_in_LIPID_MAPS(inchikey)
        if lipid_data:
            lipid.extend(data['residues'])
            if data['residues'][0].classification != 'fatty':
                print('WARNING: The residue', data['residues'][0].name, 'is not classified as fatty')
        else:
            not_lipid.extend(data['residues'])
            if data['residues'][0].classification == 'fatty':
                print('WARNING: The residue', data['residues'][0].name, 'is classified as fatty')
    """ esto será para los glucolipidos
    # Propierty to make easier to see if the residue is lipid
    # without having to iterate over all the residues in mem_dict['lipid]
    for res in mem_dict['lipid']:
        res.is_lipid = True

     # Check if non lipid are bonded to lipid like in    
    for res in not_lipid:
        for bonded_res in res.get_bonded_residues():
            if getattr(res, 'is_lipid', False):
                mem_dict['other'].append(res)
                break """
    
    # Select only the lipids and potential membrane members
    res_idx = [res.index for res in lipid]
    all_res_sele = '(resindex ' + ' '.join(map(str,(res_idx)))+')'

    # Find neighbouring lipids
    neighbours = Neighbours(
    universe=u,
    lipid_sel=all_res_sele + ' and not element H' # remove hydrogens to make it faster
    )
    neighbours.run(verbose=False)

    # Find clusters
    frame_idx = 0 # TO-DO: clustering in more than one frame
    n_clusters, labels, clusters = find_clusters(neighbours.neighbours[frame_idx])

    # cluster con más de 30 lipidos agrupados es la membrana:
    # https://pythonhosted.org/fatslim/documentation/leaflets.html#membrane-identification
    membranes_map = {}
    for n, c in clusters.items():
        if len(c) > 30:
            membranes_map[str(n)] = neighbours.membrane.residues.resindices[clusters[n]] 
    
    charges = abs(np.array([atom.charge for atom in u.atoms]))
    for key, mem in membranes_map.items():
        # for better leaflet assignation we only use polar atoms
        polar_atoms = []
        for res_idx in mem:
            res = u.residues[res_idx]
            res_ch = charges[[at.index for at in res.atoms]]
            max_ch_idx = np.argmax(res_ch)
            polar_atoms.append(res.atoms[max_ch_idx].index)

        # Add dimensions to the universe so AssignLeaflets does not crash
        if u.dimensions is None:
            print('WARNING trajectory probably has no box variable. Setting dimensions for lipyphilic')
            u.trajectory.add_transformations(set_dimensions([ 100,100,100, 90, 90, 90]))

        sele = 'index ' + ' '.join(map(str,polar_atoms))

        leaflets = AssignLeaflets(
            universe=u,
            lipid_sel=sele,
            #n_bins=10
        )
        leaflets.run(verbose=False)

        top_idx = leaflets.membrane.resindices[leaflets.leaflets[:,0] == 1]
        bot_idx = leaflets.membrane.resindices[leaflets.leaflets[:,0] == -1]
        membranes_map[key] = {
            'res_idx': mem,
            'top': top_idx,
            'bot': bot_idx
        }

    
    if debug:
        counts = neighbours.count_neighbours()
        return membranes_map, lipid, neighbours, counts, clusters
    else:
        return membranes_map


def find_clusters(contact_matrix):
    """
    Find clusters of molecules in contact using a sparse matrix.
    
    Parameters:
    contact_matrix (csr_matrix): A sparse matrix where each entry (i, j) indicates contact 
                                 between molecule i and molecule j.
    
    Returns:
    tuple: 
        - n_components (int): Number of clusters found.
        - labels (ndarray): An array of shape (n_samples,) where labels[i] is the cluster 
                            index for the ith molecule.
        - clusters (dict): A dictionary where keys are cluster labels and values are lists of 
                           molecule indices that belong to that cluster.
    """
    # Find the connected components (clusters) in the graph
    n_components, labels = connected_components(contact_matrix)
    
    # Create a dictionary to store molecules in each cluster
    clusters = {}
    for idx, label in enumerate(labels):
        if label not in clusters:
            clusters[label] = []
        clusters[label].append(idx)
    
    return n_components, labels, clusters

