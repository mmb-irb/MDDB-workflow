from model_workflow.tools.get_pytraj_trajectory import calculate_frame_step
from model_workflow.utils.topology_converter import to_MDAnalysis_topology
from model_workflow.utils.auxiliar import save_json
from model_workflow.utils.type_hints import *
from collections import Counter
import numpy as np
import MDAnalysis

def lipid_interactions (
    input_trajectory_filepath : str,
    topology_file : 'File',
    output_analysis_filepath : str,
    membrane_map: dict,
    snapshots : int,
    frames_limit: int = 100):
    """
        Lipid-protein interactions analysis.
    """
    if membrane_map is None or membrane_map['n_mems'] == 0:
        print('-> Skipping lipid-protein interactions analysis')
        return
    print('-> Running lipid-protein interactions analysis')
    
    mda_top = to_MDAnalysis_topology(topology_file.absolute_path)
    u = MDAnalysis.Universe(mda_top, input_trajectory_filepath)
    frame_step, frame_count = calculate_frame_step(snapshots, frames_limit)
    lipids = [ data['resname'] for data in membrane_map['references'].values()]
    lipids_str = " ".join(lipids)
    resids = np.unique(u.select_atoms('protein').resindices)
    ocupancy = {resid: Counter() for resid in resids}

    # Only iterate through the frames you need
    for ts in u.trajectory[0:1:frame_step]:
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
    print()
    for resid, counter in ocupancy.items():
        for lipid, count in counter.items():
            ocupancy_arrs[lipid][resid] = count

    # Normalize the occupancy arrays by dividing by the number of frames
    for lipid in lipids:
        print(lipids,type(ocupancy_arrs[lipid]))
        ocupancy_arrs[lipid] /= frame_count
        ocupancy_arrs[lipid] = ocupancy_arrs[lipid].tolist()
    # Save the data
    data = { 'data': ocupancy_arrs}
    save_json(data, output_analysis_filepath)