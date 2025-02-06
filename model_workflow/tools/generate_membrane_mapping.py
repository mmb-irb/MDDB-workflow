import os, warnings
import numpy as np
import MDAnalysis
from biobb_mem.fatslim.fatslim_membranes import fatslim_membranes, parse_index
from model_workflow.utils.constants import MEMBRANE_MAPPING_FILENAME
from model_workflow.utils.auxiliar import save_json
from model_workflow.utils.topology_converter import to_MDAnalysis_topology
from model_workflow.tools.get_inchi_keys import get_inchi_keys, is_in_swiss_lipids
from model_workflow.utils.type_hints import *
from model_workflow.utils.warnings import warn


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
     # Patch case where there no internet
    try:
        # This would return a ConnectionError
        is_in_swiss_lipids('test')
    except:
        # Then we map the lipids/membrane
        warn('There is a problem connecting to the SwissLipids database.\n'
             'Skipping membrane mapping.')
        return None
    print('-> Generating membrane mapping')
    assert topology_file.extension == 'json', 'Input topology file must be in json format: '+ topology_file.extension
    mda_top = to_MDAnalysis_topology(topology_file.absolute_path)
    u = MDAnalysis.Universe(mda_top, structure_file.absolute_path)
    # Get InChI keys of non-proteic/non-nucleic residues
    inchi_keys = get_inchi_keys(u, structure)
    # Classsify the residues as lipid or not
    lipid_ridx, glclipid_ridx = [], []
    for inchikey, res_data in inchi_keys.items():
        lipid_data = is_in_swiss_lipids(inchikey)
        # We don't use lipid data for now, if we have it it is present in LIPID MAPS
        if lipid_data:
            lipid_ridx.extend(res_data['resindices'])
            if all('fatty' not in classes for classes in res_data['classification']):
                warn(f'The InChIKey {inchikey} of {str(res_data["resname"])} is not '
                     f'classified as fatty {res_data["classification"]} but it is a lipid')
            # Glucolipids are saved separately to solve a problem with FATSLiM later
            if all(len(classes)>1 for classes in res_data['classification']):
                for grp in res_data['resgroups']:
                    if len(grp)>1:
                        glclipid_ridx.append(grp)
        else:
            if any('fatty' in classes for classes in res_data['classification']):
                warn(f'The InChIKey {inchikey} of {str(res_data["resname"])} is '
                     f'classified as fatty but is not a lipid.\n'
                     f'Resindices: {str(res_data["resindices"])}\n'
                     'In case it is a lipid, please add it to the LIPID MAPS database: https://www.lipidmaps.org/new/reg/')
    # Prepare the membrane mapping OBJ/JSON
    mem_map_js = {'n_mems': 0, 'mems': {}, 'no_mem_lipid': {}}
    # if no lipids are found, we save the empty mapping and return
    if len(lipid_ridx) == 0:
        print('No lipids found in the structure.')
        save_json(mem_map_js, MEMBRANE_MAPPING_FILENAME)
        return mem_map_js
    
    # Select only the lipids and potential membrane members
    mem_candidates = u.select_atoms(f'(resindex {" ".join(map(str,(lipid_ridx)))})')
    # for better leaflet assignation we only use polar atoms
    charges = abs(np.array([atom.charge for atom in u.atoms]))
    polar_atoms = []
    for ridx in mem_candidates.residues.resindices:
        res = u.residues[ridx]
        res_ch = charges[res.atoms.ix]
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
    print('Running BioBB FATSLiM Membranes:')
    fatslim_membranes(input_top_path=structure_file.absolute_path,
                    output_ndx_path=output_ndx_path,
                    properties=prop)
    # Parse the output to get the membrane mapping
    mem_map = parse_index(output_ndx_path)
    os.remove(output_ndx_path)
    # Save the membrane mapping as a JSON file
    n_mems = len(mem_map)//2
    mem_map_js['n_mems'] = n_mems
    no_mem_lipids = set(mem_candidates.atoms.indices)
    for i in range(n_mems):
        # and they are not assigned to any membrane. FATSLiM indexes start at 1
        top = (np.array(mem_map[f'membrane_{i+1}_leaflet_1'])-1).tolist()  # JSON does not support numpy arrays
        bot = (np.array(mem_map[f'membrane_{i+1}_leaflet_2'])-1).tolist()
        top_ridx = u.atoms[top].residues.resindices
        bot_rdix = u.atoms[bot].residues.resindices
        # Some polar atoms from the glucids are to far from the polar atoms of the lipids
        top = set(top)
        bot = set(bot)
        remove_ridx = []
        for grp in glclipid_ridx:
            if np.in1d(grp, top_ridx).any():
                top.update(u.residues[grp].atoms.indices)
                remove_ridx.append(grp)
                continue
            if np.in1d(grp,bot_rdix).any():
                bot.update(u.residues[grp].atoms.indices)
                remove_ridx.append(grp)
        for grp in remove_ridx:
            glclipid_ridx.remove(grp)
        mem_map_js['mems'][(str(i))] = {
            'leaflets': {
                'top': list(map(int, top)),
                'bot': list(map(int, bot))
                }
            }
        # Remove lipids not assigned to any membrane
        no_mem_lipids -= top
        no_mem_lipids -= bot
    mem_map_js['no_mem_lipid'] = list(map(int, no_mem_lipids))
    # Print the results and save the membrane mapping
    no_mem_lipids_str = f'{len(glclipid_ridx)} lipid/s not assigned to any membrane.' if len(glclipid_ridx)>0 else ''
    print(f'{n_mems} membrane/s found. ' + no_mem_lipids_str)
    save_json(mem_map_js, MEMBRANE_MAPPING_FILENAME)
    return mem_map_js

