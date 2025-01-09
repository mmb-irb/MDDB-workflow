import os
import numpy as np
import MDAnalysis
from biobb_mem.fatslim.fatslim_membranes import fatslim_membranes, parse_index
from model_workflow.utils.constants import MEMBRANE_MAPPING_FILENAME
from model_workflow.utils.auxiliar import save_json
from model_workflow.utils.topology_converter import to_MDAnalysis_topology
from model_workflow.tools.get_inchi_keys import get_inchi_keys, is_in_LIPID_MAPS
from model_workflow.utils.type_hints import *
from typing import List


def generate_membrane_mapping(structure : 'Structure',
                              topology_file : 'File',
                              structure_file : 'File',
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
    # RUBEN: esto será para los glucolipidos
    """ 
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
    mem_candidates = f'(resindex {" ".join(map(str,(res_idx)))})'

    # for better leaflet assignation we only use polar atoms
    charges = abs(np.array([atom.charge for atom in u.atoms]))
    polar_atoms = []
    for res_idx in u.select_atoms(mem_candidates).resindices:
        res = u.residues[res_idx]
        res_ch = charges[[at.index for at in res.atoms]]
        max_ch_idx = np.argmax(res_ch)
        polar_atoms.append(res.atoms[max_ch_idx].index)
    headgroup_sel = f'(index {" ".join(map(str,(polar_atoms)))})'

    # Run FATSLiM to find the membranes
    prop = {
        'selection': headgroup_sel,
        'cutoff': 2.2,
        'ignore_no_box': True,
        'disable_logs': True,
        'return_hydrogen':True
    }
    output_ndx_path = "tmp_mem_map.ndx"
    fatslim_membranes(input_top_path=structure_file.absolute_path,
                    output_ndx_path=output_ndx_path,
                    properties=prop)
    # Parse the output to get the membrane mapping
    mem_map = parse_index(output_ndx_path)
    os.remove(output_ndx_path)
    # Save the membrane mapping as a JSON file
    n_mems = len(mem_map)//2
    mem_map_js = {'n_mems': n_mems, 'mems': {}}
    for i in range(n_mems):
        mem_map_js['mems'][(str(i))] = {
            'leaflets': {
                'top': mem_map[f'membrane_{2*i+1}_leaflet_1'], 
                'bot': mem_map[f'membrane_{2*i+1}_leaflet_2']
                }
            }
    save_json(mem_map_js, MEMBRANE_MAPPING_FILENAME)
    return mem_map_js

