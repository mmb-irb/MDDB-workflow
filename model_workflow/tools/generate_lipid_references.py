from model_workflow.utils.auxiliar import save_json

from model_workflow.tools.get_inchi_keys import get_inchi_keys, is_in_swiss_lipids, is_in_LIPID_MAPS
from model_workflow.utils.type_hints import *
from model_workflow.utils.warnings import warn


def generate_lipid_references(structure: 'Structure',
                              universe: 'Universe',
                              lipid_map_filepath: str
                              ) -> List[dict]:
    # Patch case where there no internet
    try:
        # This would return a ConnectionError
        is_in_swiss_lipids('test')
    except:
        # Then we map the lipids/membrane
        warn('There was a problem connecting to the SwissLipids database.')
        return None

    if universe.atoms.charges is None:
        print('Topology file does not have charges, cannot generate lipid references.')
        return save_json([], lipid_map_filepath)
    
    print('-> Getting lipid references')

    # Get InChI keys of non-proteic/non-nucleic residues
    inchi_keys = get_inchi_keys(universe, structure)

    lipid_references = []
    for inchikey, res_data in inchi_keys.items():
        SL_data = is_in_swiss_lipids(inchikey)
        LM_data = is_in_LIPID_MAPS(inchikey)

        # We don't use lipid data for now, if we have it it is present in LIPID MAPS
        if SL_data or LM_data:
            lipid_references.append({'inchikey': inchikey,
                                     'resname': list(res_data['resname'])[0],
                                     'resindices': list(map(int, res_data['resindices'])),
                                     'indices': structure.select_residue_indices(res_data['resindices']).to_ngl(),
                                     'swisslipids': SL_data,
                                     'lipidmaps': LM_data,
                                     })
            cls = res_data['classification']
            # If the residue is a lipid, we check if it is classified as fatty
            if all('fatty' not in classes for classes in cls):
                warn(f'The InChIKey {inchikey} of {str(res_data["resname"])} is not '
                     f'classified as fatty {cls} but it is a lipid')

            # Glucolipids are saved separately to solve a problem with FATSLiM later (ex: A01IR)
            if type(cls) == list and all(len(classes) > 1 for classes in cls):
                lipid_references[-1]['resgroups'] = list(map(int, res_data['resgroups']))

        else:
            # If the InChIKey is not in SwissLipids or LIPID MAPS, we check if it is classified as fatty
            if any('fatty' in classes for classes in res_data['classification']):
                warn(f'The InChIKey {inchikey} of {str(res_data["resname"])} is '
                     f'classified as fatty but is not a lipid.\n'
                     f'Resindices: {str(res_data["resindices"])}')

    save_json(lipid_references, lipid_map_filepath)
    return lipid_references
