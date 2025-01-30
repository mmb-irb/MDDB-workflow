from model_workflow.tools.get_pytraj_trajectory import get_reduced_pytraj_trajectory
from model_workflow.utils.auxiliar import save_json
from model_workflow.utils.type_hints import *

import pytraj as pt

def density (
    input_structure_filepath : str,
    input_trajectory_filepath : str,
    output_analysis_filepath : str,
    membrane_map: dict,
    structure : 'Structure',
    snapshots : int,
    density_type = 'mass',
    frame_limit = 1000):
    print('-> Running density analysis')

    if membrane_map['n_mems'] == 0:
        # Do something special for density analysis for 
        # membranes like leaflets separation, polargroups, etc.
        print(' No membranes found in the structure. Skipping density analysis.')
        return

    # Load
    tj, frame_step, frames_count = get_reduced_pytraj_trajectory(input_structure_filepath, input_trajectory_filepath, snapshots, frame_limit)
    
    # Set every selections to be analyzed separately
    components = []
    for chain in structure.chains:
        components.append({
            'name': chain.name,
            'selection': chain.get_selection()
        })
    
    # Parse selections to pytraj masks
    pytraj_masks = [ component['selection'].to_pytraj() for component in components ]
    # Run pytraj
    out = pt.density(tj, pytraj_masks, density_type)
    # Iterate pytraj results
    results = iter(out.values())
    for component in components:
        # Mine pytraj results
        component['dens'] = list(next(results))
        component['stdv'] = list(next(results))
        # Parse the selection to atom indices
        component['selection'] = component['selection'].atom_indices
    # Export results
    data = {'data': { 'comps': components, 'z': list(out['z']) } }
    save_json(data, output_analysis_filepath)