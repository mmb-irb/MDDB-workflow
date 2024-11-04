import MDAnalysis
from lipyphilic.lib.neighbours import Neighbours
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
    print('Calculando la membrana...')
    assert topology_file.extension == 'json', 'Input topology file must be in json format'
    mda_top = to_MDAnalysis_topology(topology_file.absolute_path)
    u = MDAnalysis.Universe(mda_top, structure_file.absolute_path)
    # Get InChI keys of non-proteic/non-nucleic residues
    print("dadaa")
    inchi_keys = get_inchi_keys(u, structure)
    # TO-DO: identify substructures like carbohydrates by neutralizing bonds:
    # ich = inchi_keys['ZNJXAXZHPQQZRY-LXGUWJNJSA-N']['inchi'] # A01IP
    # for args in [[ich], [ich, 1], [ich, 0, 1], [ich, 1, 1]]:
    #    display(inchi_2_mol(*args))
    print(inchi_keys)
    # Classsify the residues as lipid or not
    lipid, not_lipid = [], []
    for inchikey, data in inchi_keys.items():
        lipid_data = is_in_LIPID_MAPS(inchikey)
        if lipid_data:
            lipid.extend(data['residues'])
            # print(f'Residue {data["resname"]} is a lipid in LIPID MAPS', lipid_data)
        else:
            not_lipid.extend(data['residues'])
            # print(f'Residue {data["resname"]} not found in LIPID MAPS')
    print(lipid)
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
    res_num = [res.number for res in lipid]
    all_res_sele = '(resnum ' + ' '.join(map(str,(res_num)))+')'

    # Find neighbouring lipids
    neighbours = Neighbours(
    universe=u,
    lipid_sel=all_res_sele + ' and not element H' # remove hydrogens to make it faster
    )
    neighbours.run()

    # Find clusters
    frame_idx = 0 # TO-DO: clustering in more than one frame
    n_clusters, labels, clusters = find_clusters(neighbours.neighbours[frame_idx])

    # cluster con más de 10 lipidos agrupados es la membrana
    membranes_map = {}
    for n, c in clusters.items():
        if len(c) > 10:
            membranes_map[str(n)] = neighbours.membrane.residues.resnums[clusters[n]] 
    print(membranes_map)
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

