# RMSD per resiude analysis
# 
# Perform the RMSD analysis for each residue
# The analysis is carried by pytraj

import pytraj as pt
import re
from typing import List

from distutils.version import StrictVersion

from model_workflow.tools.get_pytraj_trajectory import get_reduced_pytraj_trajectory
from model_workflow.utils.auxiliar import delete_previous_log, save_json

# The pytraj trajectory may be reduced
def rmsd_per_residue (
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_analysis_filename : str,
    structure : 'Structure',
    pbc_residues : List[int],
    snapshots : int,
    frames_limit : int):

    # Parse the trajectory intro ptraj
    # Reduce it in case it exceeds the frames limit
    pt_trajectory, frame_step, frames_count = get_reduced_pytraj_trajectory(input_topology_filename, input_trajectory_filename, snapshots, frames_limit)

    # We must filter out residues which only have 1 atom (e.g. ions)
    # This is because sometimes pytraj does not return results for them and then the number of results and residues does not match
    # More info: https://github.com/Amber-MD/pytraj/issues/1580
    ion_atom_indices = []
    for residue in structure.residues:
        if len(residue.atom_indices) == 1:
            ion_atom_indices += residue.atom_indices

    # We must exclude here PBC residues from the analysis
    # PBC residues may jump through boundaries thus having non-sense high RMSD values
    pbc_atom_indices = []
    if pbc_residues and len(pbc_residues) > 0:
        pbc_selection = structure.select_residue_indices(pbc_residues)
        pbc_atom_indices = pbc_selection.atom_indices

    # Filter the trajectory with the specified residue indices
    filter_out_atom_indices = ion_atom_indices + pbc_atom_indices
    filter_out_selection = structure.select_atom_indices(filter_out_atom_indices)
    filter_in_selection = structure.invert_selection(filter_out_selection)

    # Calculate the residue indices of the overall structure remaining in the filtered trajectory
    residue_indices = structure.get_selection_residue_indices(filter_in_selection)
    # Count the number of residues
    residue_count = len(residue_indices)

    # Make sure we have more than 1 residue left at this point
    # Otherwise it makes not sense running this analysis for a single residue
    if residue_count <= 1:
        tail = 'no resiudes at all' if residue_count == 0 else 'one residue only'
        print('  The RMSD per residue analysis will be skipped since it would be applied to ' + tail)
        return

    # Convert the selection to a pytraj mask
    pytraj_selection = filter_in_selection.to_pytraj()
    # Apply the mask to the trajectory
    filtered_pt_trajectory = pt_trajectory[pytraj_selection]
    # WARNING: This extra line prevents the error "Segment violation (core dumped)" in some pdbs
    # This happens with some random pdbs which pytraj considers to have 0 Mols
    # More info: https://github.com/Amber-MD/cpptraj/pull/820
    # IMPORTANT: This is critical for pytraj <= 2.0.5 but it makes the code fail for pytraj 2.0.6
    if StrictVersion(pt.__version__) <= StrictVersion('2.0.5'):
        filtered_pt_trajectory.top.start_new_mol()

    # Run the analysis in pytraj
    # The result data is a custom pytraj class: pytraj.datasets.datasetlist.DatasetList
    # This class has keys but its attributes can not be accessed through the key
    # They must be accessed thorugh the index
    # DANI: When the 'resname' argument is missing it prints "Error: Range::SetRange(None): Range is -1 for None"
    # DANI: However there is no problem and the analysis runs flawlessly
    # DANI: For this reason we call this function with no resname and then we remove the log
    data = pt.rmsd_perres(filtered_pt_trajectory)
    # We remove the previous error log
    delete_previous_log()

    # We remove the first result, which is meant to be the whole rmsd and whose key is 'RMSD_00001'
    del data[0]

    # Check the expected output number of residues to match with the pytraj results
    if residue_count != len(data):
        raise ValueError(f'Number of residues in structure ({residue_count}) does not match number of residues in RMSD-perres results ({len(data)})')

    # Mine the analysis data
    whole_structure_residues_count = len(structure.residues)
    rmsd_per_residue = [ None ] * whole_structure_residues_count
    for index, residue_data in enumerate(data):
        # Get the actual residue index of the current data
        residue_index = residue_indices[index]
        # Add current data to its corresponding position in the 'per residue' list
        rmsd_per_residue[residue_index] = list(residue_data)

    # Export the analysis in json format
    output_analysis = { 'step': frame_step, 'rmsdpr': rmsd_per_residue }
    save_json(output_analysis, output_analysis_filename)