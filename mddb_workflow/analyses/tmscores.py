"""Module to perform the TM score analysis."""

import numpy as np
import biotite.structure as struc
import biotite.structure.io as strucio
import biotite.structure.io.xtc as xtc
from mddb_workflow.utils.auxiliar import save_json
from mddb_workflow.utils.constants import REFERENCE_LABELS, OUTPUT_TMSCORES_FILENAME
from mddb_workflow.utils.type_hints import *
from mddb_workflow.tools.get_reduced_trajectory import calculate_frame_step


def _tm_score_no_filter(reference, subject, reference_indices, subject_indices, reference_length="shorter"):
    # Bypass the filter_amino_acids check by directly computing the score
    ref_length = struc.tm._get_reference_length(
        reference_length,
        struc.get_residue_count(reference),
        struc.get_residue_count(subject)
    )
    distances = struc.geometry.distance(reference[reference_indices], subject[subject_indices])
    return np.sum(struc.tm._tm_score(distances, ref_length)).item() / ref_length


def tmscores(
    trajectory_file: 'File',
    output_directory: str,
    first_frame_file: 'File',
    average_structure_file: 'File',
    structure: 'Structure',
    pbc_selection: 'Selection',
    snapshots: int,
    frames_limit: int = 200):
    """Perform the tm score using the tmscoring package."""
    # Set the main output filepath
    output_analysis_filepath = f'{output_directory}/{OUTPUT_TMSCORES_FILENAME}'
    struc.filter_amino_acids = lambda array: np.ones(array.array_length(), dtype=bool)
    tmscore_references = [first_frame_file, average_structure_file]

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

    output_analysis = []
    template = strucio.load_structure(first_frame_file.path)
    ca_list = selection.to_list()

    ca_mask = np.isin(np.arange(template.array_length()), ca_list)
    template = template[ca_mask]
    frame_step, _ = calculate_frame_step(snapshots, frames_limit)
    xtc_file = xtc.XTCFile.read(trajectory_file.path, atom_i=ca_list, step=frame_step)
    reduced_trajectory = xtc_file.get_structure(template)

    # Iterate over each reference and group
    for reference in tmscore_references:
        print(f' Running TM score using {reference.filename} as reference')

        # Get the TM score of each frame
        tmscores = []
        reference_frame = reduced_trajectory[0]
        for i in range(len(reduced_trajectory)):
            subject = reduced_trajectory[i]
            superimposed, _, ref_indices, sub_indices = struc.superimpose_homologs(
                reference_frame, subject)
            # Run the tmscoring over the current frame against the current reference
            tm = _tm_score_no_filter(reference_frame, superimposed, ref_indices, sub_indices)
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
    save_json({'start': 0, 'step': frame_step, 'data': output_analysis}, output_analysis_filepath)
