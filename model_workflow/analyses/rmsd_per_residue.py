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
        # DANI: Esto es útil en pytraj <= 2.0.5 pero hace fallar el código a partir de pytraj 2.0.6
        if StrictVersion(pt.__version__) <= StrictVersion('2.0.5'):
            filtered_pt_trajectory.top.start_new_mol()
        # Create a filtered topology with the same selection than pytraj
        control_structure = structure.filter(selection)
    else:
        filtered_pt_trajectory = pt_trajectory
        control_structure = structure

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
    if len(control_structure.residues) != len(data):
        raise ValueError('Number of residues in structure does not match number of residues in RMSD-perres results')

    # Mine the analysis data
    output_analysis = []
    for residue_index, residue_data in enumerate(data):
        # DEPRECATED: Now residues are found using the residue index and the control structure
        # Convert pytraj residue keys to source notation
        # Key format: SER:1, TYR:2, ...
        # match = re.match('(.*):(.*)', residue.key)
        # num = match.groups(0)[1]
        residue = control_structure.residues[residue_index]
        residue_tag = residue.chain.name + ':' + str(residue.number) + residue.icode
        # Write data to the output file
        # The 'residue' DataArray contains numeric values (rmsds) and it is not JSON serializable
        output_analysis.append(
            {
                'name': residue_tag,
                'rmsds': list(residue_data),
            }
        )

    # Export the analysis in json format
    with open(output_analysis_filename, 'w') as file:
        json.dump({ 'data': output_analysis }, file)
    
    # It is not possible to represent the whole rmsd per residue with a classical graph