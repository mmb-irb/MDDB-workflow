from model_workflow.utils.auxiliar import save_json

from model_workflow.tools.get_inchi_keys import get_inchi_keys, is_in_swiss_lipids, is_in_LIPID_MAPS
from model_workflow.utils.auxiliar import warn
from model_workflow.utils.type_hints import *

def generate_lipid_references(structure: 'Structure',
                              universe: 'Universe',
                              output_filepath : str,
                              ) -> List[dict]:
    """Generate the lipid references."""
    # Patch case where there no internet
    try:
        # This would return a ConnectionError
        is_in_swiss_lipids('test')
    except Exception as e:
        # Then we map the lipids/membrane
        warn(f'There was a problem connecting to the SwissLipids database: {e}')
        return None

    if not universe.universe.atoms.charges.any(): # AGUS: he añadido .any() porque el error me lo indicaba, pero también me sugería .all() , no sé cuál encaja mejor
        print('Topology file does not have charges, cannot generate lipid references.')
        return save_json([], output_filepath)

    # Get InChI keys of non-proteic/non-nucleic residues
    inchi_keys = get_inchi_keys(universe, structure)

    lipid_references = []
    lipid_map = []
    for inchikey, res_data in inchi_keys.items():
        SL_data = is_in_swiss_lipids(inchikey)
        LM_data = is_in_LIPID_MAPS(inchikey)
        # If we dont find it, we try without stereochemistry
        if not SL_data:
            SL_data = is_in_swiss_lipids(inchikey, only_first_layer=True)
        if not LM_data:
            LM_data = is_in_LIPID_MAPS(inchikey, only_first_layer=True)

        # We don't use lipid data for now, if we have it it is present in LIPID MAPS
        if SL_data or LM_data:
            lipid_references.append({'inchikey': inchikey,
                                     'swisslipids': SL_data,
                                     'lipidmaps': LM_data,
                                     })
            # Format needed for generate_residue_mapping
            lipid_map.append({ 
                'name': list(res_data['resname'])[0], 
                'residue_indices': list(map(int, res_data['resindices'])),
                'resgroups': res_data['resgroups'], 
                'match': { 
                    'ref': { 'inchikey': inchikey } } 
                    })
 
            # QUALITY CHECKS
            cls = res_data['classification']
            # If the residue is a lipid, we check if it is classified as fatty
            if all('fatty' not in classes for classes in cls):
                warn(f'The InChIKey {inchikey} of {str(res_data["resname"])} is not '
                     f'classified as fatty {cls} but it is a lipid')

        else:
            # If the InChIKey is not in SwissLipids or LIPID MAPS, we check if it is classified as fatty
            if any('fatty' in classes for classes in res_data['classification']):
                warn(f'The InChIKey {inchikey} of {str(res_data["resname"])} is '
                     f'classified as fatty but is not a lipid.\n'
                     f'Resindices: {str(res_data["resindices"])}')

    save_json(lipid_references, output_filepath)
    return lipid_map
