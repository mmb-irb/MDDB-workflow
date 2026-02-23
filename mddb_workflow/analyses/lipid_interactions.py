from mddb_workflow.tools.get_reduced_trajectory import calculate_frame_step
from mddb_workflow.utils.auxiliar import save_json, load_json
from mddb_workflow.utils.constants import OUTPUT_LIPID_INTERACTIONS_FILENAME, LIPIDS_RESIDUE_NAMES
from mddb_workflow.utils.type_hints import *
import numpy as np
import MDAnalysis


def lipid_interactions(
    universe: 'MDAnalysis.Universe',
    output_directory: str,
    inchikey_map: list[dict],
    cg_residues: list[int],
    snapshots: int,
    frames_limit: int = 100
):
    """Lipid-protein interactions analysis."""
    if universe.select_atoms('protein').n_atoms == 0:
        print('-> No protein found, skipping channels analysis')
        return
    output_analysis_filepath = f'{output_directory}/{OUTPUT_LIPID_INTERACTIONS_FILENAME}'

    # Check if we're dealing with coarse-grain simulations
    lipid_map = [ref for ref in inchikey_map if ref['is_lipid']]
    if len(lipid_map) and len(lipid_map[0]['residue_indices']) == 1:
        # This is probably a false positive lipid
        return
    if len(cg_residues) > 0:
        data = cg_lipid_interactions(universe, snapshots, frames_limit)
    elif inchikey_map and len(inchikey_map) > 0 and lipid_map:
        data = aa_lipid_interactions(universe, snapshots, frames_limit, lipid_map)
    else:
        print('-> Skipping lipid-protein interactions analysis')
        return

    # Wrap the data in a dictionary
    data = {'data': data}
    save_json(data, output_analysis_filepath)


def aa_lipid_interactions(universe, snapshots, frames_limit, lipid_map):
    """Atomistic lipid-protein interactions analysis."""
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
        lipid_near_prot = universe.select_atoms('(around 6 protein) and group lipid_group and not protein',
                                                lipid_group=lipid_group)
        for i, residx in enumerate(protein_residx):
            # Find lipid atoms near a specific residue
            residue_atoms = universe.select_atoms(f'resindex {residx}')
            lipid_near_res = lipid_near_prot.select_atoms('around 6 global group residuegroup',
                                                residuegroup=residue_atoms)
            # Add the count of lipids
            for lipid_residx in lipid_near_res.residues.resindices:
                ocupancy_arrs[i, residx_2_idx[lipid_residx]] += 1

    # Normalize the occupancy arrays by dividing by the number of frames
    ocupancy_arrs /= frame_count

    # Save the data
    data = {'residue_indices': protein_residx.tolist()}
    for lipid in lipid_map:
        # Convert lipid residue indices to array indices
        lipid_idx = [residx_2_idx[residx] for residx in lipid['residue_indices']]
        data[lipid['inchikey']] = ocupancy_arrs[:, lipid_idx].sum(1).tolist()

    return data


def cg_lipid_interactions(universe, snapshots, frames_limit):
    """Coarse-grain lipid-protein interactions analysis."""
    print('-> Performing coarse-grain lipid-protein interactions analysis')
    # Identify lipid residues by name for coarse-grain simulations
    lipid_residues = []
    for resname in LIPIDS_RESIDUE_NAMES:
        lipids = universe.select_atoms(f'resname {resname}')
        if len(lipids) > 0:
            lipid_residues.append(resname)

    if not lipid_residues:
        print('-> No lipid residues found for coarse-grain analysis')
        return {}

    frame_step, frame_count = calculate_frame_step(snapshots, frames_limit)
    protein_residx = universe.select_atoms('protein').residues.resindices
    assert len(protein_residx) > 0

    # Create occupancy arrays for each lipid type
    lipid_occupancy = {resname: np.zeros(len(protein_residx)) for resname in lipid_residues}

    # Iterate through frames
    for ts in universe.trajectory[0:snapshots:frame_step]:
        for i, residx in enumerate(protein_residx):
            residue_atoms = universe.select_atoms(f'resindex {residx}')

            # Check interactions with each lipid type
            for resname in lipid_residues:
                lipids_near_res = universe.select_atoms(f'(around 6 global group residuegroup) and resname {resname}',
                                                        residuegroup=residue_atoms)
                if len(lipids_near_res) > 0:
                    lipid_occupancy[resname][i] += len(lipids_near_res.residues)

    # Normalize by frame count
    for resname in lipid_residues:
        lipid_occupancy[resname] /= frame_count

    # Save the data in the same format as atomistic analysis
    data = {'residue_indices': protein_residx.tolist()}
    for resname in lipid_residues:
        data[resname] = lipid_occupancy[resname].tolist()

    return data


def plot_lipid_interactions(output_analysis_filepath: str):
    """Plot lipid-protein interactions analysis results."""
    import matplotlib.pyplot as plt
    data = load_json(output_analysis_filepath)['data']

    protein_residues = data['residue_indices']
    lipid_keys = [key for key in data.keys() if key != 'residue_indices']

    # Create a figure for the plot
    plt.figure(figsize=(10, 6))

    for lipid_key in lipid_keys:
        interactions = data[lipid_key]
        plt.plot(protein_residues, interactions, label=lipid_key)

    plt.xlabel('Protein Residue Indices')
    plt.ylabel('Interaction Count (Normalized)')
    plt.title('Lipid-Protein Interactions')
    plt.legend(); plt.grid(True); plt.tight_layout(); plt.show()
