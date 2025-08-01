import tmscoring

from os import remove

from model_workflow.tools.get_pdb_frames import get_pdb_frames
from model_workflow.utils.auxiliar import save_json
from model_workflow.utils.constants import REFERENCE_LABELS, OUTPUT_TMSCORES_FILENAME
from model_workflow.utils.gmx_spells import pdb_filter
from model_workflow.utils.type_hints import *

def tmscores (
    trajectory_file : 'File',
    output_directory : str,
    first_frame_file : 'File',
    average_structure_file : 'File',
    structure : 'Structure',
    pbc_selection : 'Selection',
    snapshots : int,
    frames_limit : int = 200):
    """Perform the tm score using the tmscoring package."""

    # Set the main output filepath
    output_analysis_filepath = f'{output_directory}/{OUTPUT_TMSCORES_FILENAME}'

    tmscore_references  = [first_frame_file, average_structure_file]

    # Set the alpha carbon selection
    selection = structure.select('name CA', syntax='vmd')
    if not selection:
        print('WARNING: There are not atoms to be analyzed for the TM score analysis')
        return

    # Remove PBC residues from the selection
    selection -= pbc_selection
    if not selection:
        print('WARNING: There are not atoms to be analyzed for the TM score analysis after PBC substraction')
        return

    # Convert the selection to ndx so we can use it in gromacs
    selection_name = 'alpha-carbons'
    ndx_selection = selection.to_ndx(selection_name)
    ndx_filename = '.tmscore.ndx'
    with open(ndx_filename, 'w') as file:
        file.write(ndx_selection)

    # Set the frame start and step
    start = 0
    step = None
    
    output_analysis = []

    # Iterate over each reference and group
    for reference in tmscore_references:
        print(f' Running TM score using {reference.filename} as reference')
        # Create a reference topology with only the group atoms
        # WARNING: Yes, TM score would work also with the whole reference, but it takes more time!!
        # This has been experimentally tested and it may take more than the double of time
        grouped_reference = 'gref.pdb'
        pdb_filter(reference.path, grouped_reference, ndx_filename, selection_name)
        # Get the TM score of each frame
        # It must be done this way since tmscoring does not support trajectories
        tmscores = []
        frames, step, count = get_pdb_frames(reference.path, trajectory_file.path, snapshots, frames_limit)
        for current_frame in frames:

            # Filter atoms in the current frame
            filtered_frame = 'filtered_frame.pdb'
            pdb_filter(current_frame, filtered_frame, ndx_filename, selection_name)

            # Run the tmscoring over the current frame against the current reference
            # Append the result data for each ligand
            tmscore = tmscoring.get_tm(grouped_reference, filtered_frame)
            tmscores.append(tmscore)

            remove(filtered_frame)

        # Get a standarized reference name
        reference_name = REFERENCE_LABELS[reference.filename]
        # Save the tmscores in the output object
        data = {
            'values': tmscores,
            'reference': reference_name,
            'group': 'c-alpha'
        }
        output_analysis.append(data)

        remove(grouped_reference)
    remove(ndx_filename)
    # Export the analysis in json format
    save_json({ 'start': start, 'step': step, 'data': output_analysis }, output_analysis_filepath)