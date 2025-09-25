from model_workflow.tools.get_reduced_trajectory import calculate_frame_step
from model_workflow.utils.auxiliar import save_json
from model_workflow.utils.constants import OUTPUT_LIPID_INTERACTIONS_FILENAME
from model_workflow.utils.type_hints import *
import numpy as np
import MDAnalysis

def lipid_interactions (
    universe : 'MDAnalysis.Universe',
    output_directory : str,
    lipid_map: List[dict],
    snapshots : int,
    frames_limit: int = 100):
    """
        Lipid-protein interactions analysis.
    """
    if lipid_map is None or len(lipid_map) == 0:
        print('-> Skipping lipid-protein interactions analysis')
        return
    
    # Set the main output filepath
    output_analysis_filepath = f'{output_directory}/{OUTPUT_LIPID_INTERACTIONS_FILENAME}'
    
    frame_step, frame_count = calculate_frame_step(snapshots, frames_limit)
    lipids_residx = [residue_idx for lipid in lipid_map for residue_idx in lipid['residue_indices']]
    # Create a mapping from lipid residue indices to array indices
    residx_2_idx = {residx: i for i, residx in enumerate(lipids_residx)}
    lipid_group = universe.select_atoms(f'resindex {" ".join(map(str, lipids_residx))}')
    # A counter for each pair of protein-lipid residues
    protein_residx = universe.select_atoms('protein').residues.resindices
    ocupancy_arrs = np.zeros((len(protein_residx), len(lipids_residx)))

    # Only iterate through the frames you need. TODO: parallel frames
    for ts in universe.trajectory[0:snapshots:frame_step]:
        # Select lipid atoms near the protein once
        lipid_near_prot = universe.select_atoms(f'(around 6 protein) and group lipid_group and not protein', 
                                                lipid_group=lipid_group)
        for i, residx in enumerate(protein_residx):
            # Find lipid atoms near a specific residue
            residue_atoms = universe.select_atoms(f'resindex {residx}')
            lipid_near_res = lipid_near_prot.select_atoms(f'around 6 global group residuegroup', 
                                                residuegroup=residue_atoms)
            # Add the count of lipids
            for lipid_residx in lipid_near_res.residues.resindices:
                ocupancy_arrs[i, residx_2_idx[lipid_residx]] += 1

    # Normalize the occupancy arrays by dividing by the number of frames
    ocupancy_arrs /= frame_count
    
    # Save the data
    data = { 'residue_indices': protein_residx.tolist()}
    for lipid in lipid_map:
        # Convert lipid residue indices to array indices
        lipid_idx = [residx_2_idx[residx] for residx in lipid['residue_indices']]
        data[lipid['match']['ref']['inchikey']] = ocupancy_arrs[:, lipid_idx].sum(1).tolist()
    # Wrap the data in a dictionary
    data = { 'data': data}
    save_json(data, output_analysis_filepath)