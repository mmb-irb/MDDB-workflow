# RMSD per resiude analysis
# 
# Perform the RMSD analysis for each residue
# The analysis is carried by pytraj

import pytraj as pt
import re

import json

from distutils.version import StrictVersion

from model_workflow.tools.get_pytraj_trajectory import get_reduced_pytraj_trajectory

# The pytraj trajectory may be reduced
def rmsd_per_residue (
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_analysis_filename : str,
    structure : 'Structure',
    membranes,
    frames_limit : int):

    # Parse the trajectory intro ptraj
    # Reduce it in case it exceeds the frames limit
    pt_trajectory = get_reduced_pytraj_trajectory(input_topology_filename, input_trajectory_filename, frames_limit)

    # We must filter out residues which only have 1 atom (e.g. ions)
    # This is because sometimes pytraj does not return results for them and then the number of results and residues does not match
    # More info: https://github.com/Amber-MD/pytraj/issues/1580
    ion_atom_indices = []
    for residue in structure.residues:
        if len(residue.atom_indices) == 1:
            ion_atom_indices += residue.atom_indices

    # We must exclude here membranes from the analysis
    # Membrane lipids close to boundaries are use to jump so the RMSD values of those residues would eclipse the protein
    membrane_atom_indices = []
    if membranes and len(membranes) > 0:
        membrane_selection_string = ' and '.join([ '( ' + membrane['selection'] + ' )' for membrane in membranes ])
        membrane_selection = structure.select(membrane_selection_string, syntax='vmd')
        if not membrane_selection:
            raise SystemExit('ERROR: Membrane selection "' + membrane_selection_string + '" is empty')
        membrane_atom_indices = membrane_selection.atom_indices

    # Filter the trajectory with the specified residue indices
    filter_out_atom_indices = ion_atom_indices + membrane_atom_indices
    filter_out_selection = structure.select_atom_indices(filter_out_atom_indices)
    filter_in_selection = structure.invert_selection(filter_out_selection)
    pytraj_selection = filter_in_selection.to_pytraj()
    filtered_pt_trajectory = pt_trajectory[pytraj_selection]
    # WARNING: This extra line prevents the error "Segment violation (core dumped)" in some pdbs
    # This happens with some random pdbs which pytraj considers to have 0 Mols
    # More info: https://github.com/Amber-MD/cpptraj/pull/820
    # IMPORTANT: This is critical for pytraj <= 2.0.5 but it makes the code fail for pytraj 2.0.6
    if StrictVersion(pt.__version__) <= StrictVersion('2.0.5'):
        filtered_pt_trajectory.top.start_new_mol()

    # Calculate the residue indices of the overall structure remaining in the filtered trajectory
    residue_indices = structure.get_selection_residue_indices(filter_in_selection)

    # Run the analysis in pytraj
    # The result data is a custom pytraj class: pytraj.datasets.datasetlist.DatasetList
    # This class has keys but its attributes can not be accessed through the key
    # They must be accessed thorugh the index
    # DANI: Esto devuelve "Error: Range::SetRange(None): Range is -1 for None"
    # DANI: No se por que pasa pero aparentemente funciona bien
    data = pt.rmsd_perres(filtered_pt_trajectory)

    # We remove the first result, which is meant to be the whole rmsd and whose key is 'RMSD_00001'
    del data[0]

    # Check the expected output number of residues to match with the pytraj results
    if len(residue_indices) != len(data):
        raise ValueError('Number of residues in structure (' + str(len(residue_indices)) + ') does not match number of residues in RMSD-perres results (' + str(len(data)) + ')')

    # Mine the analysis data
    whole_structure_residues_count = len(structure.residues)
    rmsd_per_residue = [ None ] * whole_structure_residues_count
    for index, residue_data in enumerate(data):
        # Get the actual residue index of the current data
        residue_index = residue_indices[index]
        # Add current data to its corresponding position in the 'per residue' list
        rmsd_per_residue[residue_index] = list(residue_data)

    # Export the analysis in json format
    output_analysis = { 'step': pt_trajectory.step, 'rmsdpr': rmsd_per_residue }
    with open(output_analysis_filename, 'w') as file:
        json.dump(output_analysis, file)