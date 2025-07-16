from model_workflow.tools.get_reduced_trajectory import calculate_frame_step
from model_workflow.utils.mda_spells import to_MDAnalysis_topology
from model_workflow.utils.auxiliar import save_json
from model_workflow.utils.constants import OUTPUT_LIPID_INTERACTIONS_FILENAME
from model_workflow.utils.type_hints import *
from collections import Counter
import numpy as np
import MDAnalysis

def lipid_interactions (
    trajectory_file : 'File',
    standard_topology_file : 'File',
    output_directory : str,
    membrane_map: dict,
    lipid_map: dict,
    snapshots : int,
    frames_limit: int = 100):
    """
        Lipid-protein interactions analysis.
    """
    if membrane_map is None or membrane_map['n_mems'] == 0:
        print('-> Skipping lipid-protein interactions analysis')
        return
    
    # Set the main output filepath
    output_analysis_filepath = f'{output_directory}/{OUTPUT_LIPID_INTERACTIONS_FILENAME}'
    
    mda_top = to_MDAnalysis_topology(standard_topology_file)
    u = MDAnalysis.Universe(mda_top, trajectory_file.path)
    frame_step, frame_count = calculate_frame_step(snapshots, frames_limit)
    lipids = set([ lipid['resname'] for lipid in lipid_map])
    lipids_str = " ".join(lipids)
    resids = np.unique(u.select_atoms('protein').resindices)
    ocupancy = {resid: Counter() for resid in resids}

    # Only iterate through the frames you need
    for ts in u.trajectory[0:snapshots:frame_step]:
        # Select non-protein atoms near the protein once
        non_protein_near = u.select_atoms(f'(around 6 protein) and (resname {lipids_str}) and not protein')
        # Make residue selections only once
        for resid in resids:
            residue_atoms = u.select_atoms(f'resid {resid}')
            # Find non-protein atoms near this specific residue
            nearby = non_protein_near.select_atoms(f'around 6 global group residuegroup', 
                                                residuegroup=residue_atoms)
            # Count residue names
            ocupancy[resid].update(Counter(nearby.residues.resnames))

    ocupancy_arrs = {lipid: np.zeros(len(ocupancy)) for lipid in lipids}
    for resid, counter in ocupancy.items():
        for lipid, count in counter.items():
            ocupancy_arrs[lipid][resid] = count

    # Normalize the occupancy arrays by dividing by the number of frames
    for lipid in lipids:
        ocupancy_arrs[lipid] /= frame_count
        ocupancy_arrs[lipid] = ocupancy_arrs[lipid].tolist()
    # Save the data
    data = { 'data': ocupancy_arrs}
    save_json(data, output_analysis_filepath)