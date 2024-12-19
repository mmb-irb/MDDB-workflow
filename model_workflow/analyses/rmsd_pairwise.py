# RMSD pairwise analysis
#
# Perform the RMSD analysis for pair of frames in the trajectory
# The analysis is carried by pytraj

import pytraj as pt

from model_workflow.tools.get_pytraj_trajectory import get_reduced_pytraj_trajectory
from model_workflow.utils.auxiliar import save_json
from model_workflow.utils.type_hints import *

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
    # If the default selection is empty then use all atoms instead
    if not selection:
        print(f'Default selection "{overall_selection}" is empty -> All atoms will be used instead')
        selection = structure.select_all()
    # If the default selection has one atom only then use all atoms instead
    # Note that a single atom selection will lead to all RMSD values beeing 0 since there is a previous alignment
    if len(selection) == 1:
        print(f'Default selection "{overall_selection}" has 1 atom only -> All atoms will be used instead')
        selection = structure.select_all()

    # Remove PBC residues from the selection
    pbc_selection = structure.select_residue_indices(pbc_residues)
    selection -= pbc_selection
    if not selection:
        raise SystemExit(f'Empty selection after substracting PBC residues')
    print(f' Analyzing {len(selection)} atoms')

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

    # Get all not failed interactions
    valid_interactions = [ interaction for interaction in interactions if not interaction.get('failed', False) ]
    
    # Make sure we have valid interactions
    # DANI: Esto es temporal, lo suyo sería que las interacciones válidas si sean analizadas
    # DANI: Lo que pasa es que pronto cambiaré los análisis de interacciones para que se haga 1 por interacción
    # DANI: De manera que no merece la pena invertir tiempo en dar soporte a esto ahora
    if len(valid_interactions) != len(interactions):
        print('There are no valid interactions -> This analysis will be skipped')
        save_json({'data': output_analysis, 'start': 0, 'step': frame_step}, output_analysis_filename)
        return

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
