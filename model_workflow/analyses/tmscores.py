from os import remove
import biotite.structure as struc
import biotite.structure.io as strucio
import biotite.structure.io.xtc as xtc

from model_workflow.utils.auxiliar import save_json
from model_workflow.utils.constants import REFERENCE_LABELS, OUTPUT_TMSCORES_FILENAME
from model_workflow.utils.gmx_spells import pdb_filter
from model_workflow.utils.type_hints import *
from model_workflow.tools.get_reduced_trajectory import calculate_frame_step

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

    selection.to_list()

    # Set the frame start and step
    start = 0
    step = None
    
    output_analysis = []

    template = strucio.load_structure(first_frame_file.path)
    template = template[selection.to_list()]
    frame_step, _ = calculate_frame_step(snapshots, frames_limit)
    xtc_file = xtc.XTCFile.read(trajectory_file.path, atom_i=selection.to_list(), step=frame_step)
    trajectory = xtc_file.get_structure(template)
    
    # Iterate over each reference and group
    for reference in tmscore_references:
        print(f' Running TM score using {reference.filename} as reference')

        # Get the TM score of each frame
        tmscores = []
        reference_frame = trajectory[0]
        for i in range(1, len(trajectory)):
            subject = trajectory[i]
            superimposed, _, ref_indices, sub_indices = struc.superimpose_structural_homologs(
                reference_frame, subject)
            # Run the tmscoring over the current frame against the current reference
            tm = struc.tm_score(reference, superimposed, ref_indices, sub_indices)
            tmscores.append(tm)

        # Get a standarized reference name
        reference_name = REFERENCE_LABELS[reference.filename]
        # Save the tmscores in the output object
        data = {
            'values': tmscores,
            'reference': reference_name,
            'group': 'c-alpha'
        }
        output_analysis.append(data)

    # Export the analysis in json format
    save_json({ 'start': start, 'step': step, 'data': output_analysis }, output_analysis_filepath)