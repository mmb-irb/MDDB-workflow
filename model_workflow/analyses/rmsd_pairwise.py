# RMSD pairwise analysis
#
# Perform the RMSD analysis for pair of frames in the trajectory
# The analysis is carried by pytraj

from os.path import exists

import pytraj as pt

from model_workflow.utils.pyt_spells import get_reduced_pytraj_trajectory
from model_workflow.utils.auxiliar import save_json, numerate_filename, get_analysis_name, warn
from model_workflow.utils.constants import OUTPUT_RMSD_PAIRWISE_FILENAME
from model_workflow.utils.type_hints import *

# The 'interactions' input is mandatory but it may be an empty list (i.e. there are no interactions)
# Take a minimal subset of atoms representing both proteins and nucleic acids
def rmsd_pairwise(
    structure_file : str,
    trajectory_file : str,
    output_directory : str,
    interactions : list,
    snapshots : int,
    frames_limit : int,
    structure : 'Structure',
    pbc_selection : 'Selection',
    overall_selection : str = "name CA or name C5'", # equivalent to "@CA,C5'" in pytraj
    ):
    """Perform an analysis for the overall structure and then one more analysis for each interaction."""

    # Set the main output filepath
    output_analysis_filepath = f'{output_directory}/{OUTPUT_RMSD_PAIRWISE_FILENAME}'

    # Parse the trajectory intro ptraj
    # Reduce it in case it exceeds the frames limit
    pt_trajectory, frame_step, frames_count = get_reduced_pytraj_trajectory(
        structure_file.path, trajectory_file.path, snapshots, frames_limit)

    # Parse the overall selection
    parsed_overall_selection = structure.select(overall_selection, syntax='vmd')
    # If the default selection is empty then use all atoms instead
    if not parsed_overall_selection:
        print(f' Default selection "{overall_selection}" is empty -> All atoms will be used instead')
        parsed_overall_selection = structure.select_all()
    # If the default selection has one atom only then use all atoms instead
    # Note that a single atom selection will lead to all RMSD values beeing 0 since there is a previous alignment
    if len(parsed_overall_selection) == 1:
        print(f' Default selection "{overall_selection}" has 1 atom only -> All atoms will be used instead')
        parsed_overall_selection = structure.select_all()

    # Remove PBC residues from the selection
    parsed_overall_selection -= pbc_selection
    if not parsed_overall_selection:
        print(f' Empty selection after substracting PBC atoms')
        return

    # Save each analysis to a dict which will be parsed to json
    output_summary = []

    # Set a filename for the current interaction data
    overall_output_analysis_filepath = numerate_filename(output_analysis_filepath, 0)
    # Get the root name out of the output analysis filepath
    overall_analysis_name = get_analysis_name(overall_output_analysis_filepath)
    # Add it to the summary
    output_summary.append({
        'name': 'Overall',
        'analysis': overall_analysis_name,
    })

    # If the overall output file does not exist then run the overall analysis
    if not exists(overall_output_analysis_filepath):
        print(f' Analyzing overall structure ({len(parsed_overall_selection)} atoms)')
        # Run the analysis
        data = pt.pairwise_rmsd(pt_trajectory, parsed_overall_selection.to_pytraj())
        # Convert data to a normal list, since numpy ndarrays are not json serializable
        data = data.tolist()
        # In case there is no data we stop here
        if len(data) == 0: raise SystemExit('Something went wrong with pytraj')
        # Write the overall analysis output data to disk
        save_json({
            'name': 'Overall',
            'rmsds': data,
            'start': 0,
            'step': frame_step
        }, overall_output_analysis_filepath)

    # Repeat the analysis with the interface residues of each interaction
    for i, interaction in enumerate(interactions, 1):
        # Get the interaction name
        name = interaction['name']
        # Parse the interaction selection
        interaction_atom_indices = interaction['interface_atom_indices_1'] + interaction['interface_atom_indices_2']
        interaction_selection = structure.select_atom_indices(interaction_atom_indices)
        # If the interaction selection matches the overall selection then reuse its data
        # This was very unlinkly until we started to support coarse grain systems
        # Now a CG DNA-DNA hybridization will likely fall here
        if interaction_selection == parsed_overall_selection:
            # Warn the user
            warning = 'Interface selection is identical to overall selection. Data will be reused.'
            warn(f'{name}: {warning}')
            output_summary.append({
                'name': name,
                # Add a note so nobody thinks this is an error
                'note': warning,
                'analysis': overall_analysis_name
            })
            # At this point overall analysis outout should already exists so we skip the analysis
            continue

        # Set a filename for the current interaction data
        numbered_output_analysis_filepath = numerate_filename(output_analysis_filepath, i)
        # Add the root of the output analysis filename to the run data
        analysis_name = get_analysis_name(numbered_output_analysis_filepath)
        # Append current interaction to the summary
        output_summary.append({
            'name': name,
            'analysis': analysis_name
        })

        # If the analysis already exists then proceed to the next interaction
        if exists(numbered_output_analysis_filepath): continue
        print(f' Analyzing "{name}" interface ({len(interaction_selection)} atoms)')
        # Run the analysis
        data = pt.pairwise_rmsd(pt_trajectory, interaction_selection.to_pytraj())
        # Convert data to a normal list, since numpy ndarrays are not json serializable
        data = data.tolist()
        # Write the interaction analysis output data to disk
        save_json({
            'name': name,
            'rmsds': data,
            'start': 0,
            'step': frame_step
        }, numbered_output_analysis_filepath)

    # Write in the analysis the starting frame and the step between frames
    # By default the first frame in the reduced trajectory is the first frame (0)

    # Export the analysis in json format
    save_json(output_summary, output_analysis_filepath)
