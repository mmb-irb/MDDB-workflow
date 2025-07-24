from model_workflow.utils.pyt_spells import get_reduced_pytraj_trajectory
from model_workflow.utils.auxiliar import save_json
from model_workflow.utils.constants import OUTPUT_DENSITY_FILENAME
from model_workflow.utils.type_hints import *

import pytraj as pt

def density (
    structure_file : 'File',
    trajectory_file : 'File',
    output_directory : str,
    membrane_map: dict,
    structure : 'Structure',
    snapshots : int,
    density_types = ['number', 'mass', 'charge', 'electron'],
    frames_limit = 1000):
    """Membrane density analysis."""
    
    if membrane_map is None or membrane_map['n_mems'] == 0:
        print('-> Skipping density analysis')
        return
    
    # Set the main output filepath
    output_analysis_filepath = f'{output_directory}/{OUTPUT_DENSITY_FILENAME}'

    # Load
    tj, frame_step, frames_count = get_reduced_pytraj_trajectory(
        structure_file.path, trajectory_file.path, snapshots, frames_limit)
    
    # Set every selections to be analyzed separately
    components = []
    for chain in structure.chains:
        components.append({
            'name': chain.name,
            'selection': chain.get_selection(),
            'number': {},
            'mass': {},
            'charge': {}, # charge will be all 0 because we cannot add charges to pytraj topology
            'electron': {}
        })
    # Parse selections to pytraj masks
    pytraj_masks = [ component['selection'].to_pytraj() for component in components ]
    # Add polar atoms selection
    polar_atoms = []
    for n in range(membrane_map['n_mems']):
        polar_atoms.extend(membrane_map['mems'][str(n)]['polar_atoms']['top'])
        polar_atoms.extend(membrane_map['mems'][str(n)]['polar_atoms']['bot'])
    components.append({
        'name': 'polar',
        'selection': polar_atoms,
        'number': {},'mass': {},'charge': {},'electron': {}
    }) 
    pytraj_masks.append('@'+', '.join(map(str,polar_atoms)))

    # Run pytraj
    for density_type in density_types:
        out = pt.density(tj, pytraj_masks, density_type)
        # Iterate pytraj results
        results = iter(out.values())
        for component in components:
            # Mine pytraj results
            component[density_type]['dens'] = list(next(results))
            component[density_type]['stdv'] = list(next(results))

    # Parse the selection to atom indices
    # Selections could be removed to make the file smaller
    for component in components:
        if component['name'] == 'polar': continue
        component['selection'] = component['selection'].atom_indices
    # Export results
    data = {'data': { 'comps': components, 'z': list(out['z']) } }
    save_json(data, output_analysis_filepath)