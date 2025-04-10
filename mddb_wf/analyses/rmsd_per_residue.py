# RMSD per resiude analysis
# 
# Perform the RMSD analysis for each residue
# The analysis is carried by pytraj

import pytraj as pt

from distutils.version import StrictVersion

from mddb_wf.utils.pyt_spells import get_reduced_pytraj_trajectory
from mddb_wf.utils.auxiliar import save_json
from mddb_wf.utils.type_hints import *

# The pytraj trajectory may be reduced
def rmsd_per_residue (
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_analysis_filename : str,
    structure : 'Structure',
    pbc_selection : 'Selection',
    snapshots : int,
    frames_limit : int):

    print('-> Running RMSD per residue analysis')

    # Parse the trajectory intro pytraj
    # Reduce it in case it exceeds the frames limit
    pt_trajectory, frame_step, frames_count = get_reduced_pytraj_trajectory(input_topology_filename, input_trajectory_filename, snapshots, frames_limit)

    # We must exclude here PBC residues from the analysis
    # PBC residues may jump through boundaries thus having non-sense high RMSD values
    # Filter the trajectory with atoms to be included in the analysis
    filter_in_selection = structure.invert_selection(pbc_selection)

    # Calculate the residue indices of the overall structure remaining in the filtered trajectory
    residue_indices = structure.get_selection_residue_indices(filter_in_selection)
    # Count the number of residues
    residue_count = len(residue_indices)

    # Make sure we have more than 1 residue left at this point
    # Otherwise it makes not sense running this analysis for a single residue
    if residue_count <= 1:
        tail = 'no residues at all' if residue_count == 0 else 'one residue only'
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
    # They must be accessed through the index
    # IMPORTANT
    # Note that the argument 'resrange' may be redundant, since we always target for all residues
    # If this argument is not passed then it prints "Error: Range::SetRange(None): Range is -1 for None"
    # However, most of the times there is no problem and the analysis runs flawlessly anyway
    # We started to have problems with some ions but then we had no results at all with coarse grain
    # More info: https://github.com/Amber-MD/pytraj/issues/1580
    resrange = f'1-{structure.residue_count + 1}'
    data = pt.rmsd_perres(filtered_pt_trajectory, resrange=resrange)

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