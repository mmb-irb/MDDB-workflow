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

    # Set the indices of residues in the analysis
    # This is important since some residues may be excluded (membrane residues)
    residue_indices = list(range(len(structure.residues)))

    # We must exclude here membranes from the analysis
    # Membrane lipids close to boundaries are use to jump so the RMSD values of those residues would eclipse the protein
    if membranes and len(membranes) > 0:
        selection = ' and '.join([ '( not ' + membrane['selection'] + ' )' for membrane in membranes ])
        selection = structure.select(selection, syntax='vmd')
        pytraj_selection = selection.to_pytraj()
        filtered_pt_trajectory = pt_trajectory[pytraj_selection]
        # WARNING: This extra line prevents the error "Segment violation (core dumped)" in some pdbs
        # This happens with some random pdbs which pytraj considers to have 0 Mols
        # More info: https://github.com/Amber-MD/cpptraj/pull/820
        # IMPORTANT: This is critical for pytraj <= 2.0.5 but it makes the code fail for pytraj 2.0.6
        if StrictVersion(pt.__version__) <= StrictVersion('2.0.5'):
            filtered_pt_trajectory.top.start_new_mol()
        # Update the residue indices list with the selected atom residues
        residue_indices = list(set([ structure.atoms[atom_index].residue_index for atom_index in selection.atom_indices ]))
    else:
        filtered_pt_trajectory = pt_trajectory

    # Run the analysis in pytraj
    # The result data is a custom pytraj class: pytraj.datasets.datasetlist.DatasetList
    # This class has keys but its attributes can not be accessed through the key
    # They must be accessed thorugh the index
    # DANI: Esto devuelve "Error: Range::SetRange(None): Range is -1 for None"
    # DANI: No se por que pasa pero aparentemente funciona bien
    data = pt.rmsd_perres(filtered_pt_trajectory)

    # We remove the first result, which is meant to be the whole rmsd and whose key is 'RMSD_00001'
    del data[0]

    # Check the structure to match in number of residues with the pytraj results
    if len(residue_indices) != len(data):
        raise ValueError('Number of residues in structure does not match number of residues in RMSD-perres results')

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