# RMSD pairwise analysis
#
# Perform the RMSD analysis for pair of frames in the trajectory
# The analysis is carried by pytraj

import pytraj as pt

from typing import List

from model_workflow.tools.get_pytraj_trajectory import get_reduced_pytraj_trajectory
from model_workflow.utils.auxiliar import save_json

# Perform an analysis for the overall structure and then one more analysis for each interaction
# The 'interactions' input is mandatory but it may be an empty list (i.e. there are no interactions)
# The pytraj trajectory ('pt_trajectory') may be reduced
# Take a minimal subset of atoms representing both proteins and nucleic acids

def rmsd_pairwise(
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_analysis_filename : str,
    interactions : list,
    snapshots : int,
    frames_limit : int,
    structure : 'Structure',
    pbc_residues : List[int],
    overall_selection : str = "name CA or name C5'", # equivalent to "@CA,C5'" in pytraj
    ):

    print('-> Running RMSD pairwise analysis')

    # Parse the trajectory intro ptraj
    # Reduce it in case it exceeds the frames limit
    pt_trajectory, frame_step, frames_count = get_reduced_pytraj_trajectory(input_topology_filename, input_trajectory_filename, snapshots, frames_limit)

    # Parse the overall selection
    selection = structure.select(overall_selection, syntax='vmd')
    if not selection:
        raise SystemExit('Selection "' + overall_selection + '" is empty')

    # Remove PBC residues from the selection
    pbc_selection = structure.select_residue_indices(pbc_residues)
    selection -= pbc_selection
    if not selection:
        raise SystemExit('Selection "' + overall_selection + '" is empty after substracting PBC residues')
    print(' Analyzing ' + str(len(selection)) + ' atoms')

    # Run the analysis
    data = pt.pairwise_rmsd(pt_trajectory, selection.to_pytraj())
    # Convert data to a normal list, since numpy ndarrays are not json serializable
    data = data.tolist()
    # In case there is no data we stop here
    if len(data) == 0:
        raise SystemExit('Something went wrong with pytraj')

    # Set the final structure data
    output_analysis = [
        {
            'name': 'Overall',
            'rmsds': data,
        }
    ]

    # Repeat the analysis with the interface residues of each interaction
    for interaction in interactions:

        # Select all interface residues in pytraj notation
        pt_interface = interaction['pt_interface_1'] + \
            interaction['pt_interface_2']
        pt_selection = ':' + ','.join(map(str, pt_interface))

        # Run the analysis
        data = pt.pairwise_rmsd(pt_trajectory, pt_selection)
        # Convert data to a normal list, since numpy ndarrays are not json serializable
        data = data.tolist()

        output_analysis.append(
            {
                'name': interaction['name'],
                'rmsds': data,
            }
        )

    # Write in the analysis the starting frame and the step between frames
    # By default the first frame in the reduced trajectory is the first frame (0)

    # Export the analysis in json format
    save_json({'data': output_analysis, 'start': 0, 'step': frame_step}, output_analysis_filename)
