from model_workflow.tools.xvg_parse import xvg_parse
from model_workflow.tools.get_reduced_trajectory import get_reduced_trajectory
from model_workflow.utils.auxiliar import save_json
from model_workflow.utils.constants import GROMACS_EXECUTABLE, REFERENCE_LABELS

import os
from subprocess import run, PIPE, Popen

from typing import List

# Run multiple RMSD analyses
# A RMSD analysis is run with each reference:
# - First frame
# - Average structure
# A RMSD analysis is run over each rmsd target:
# - Protein
# - Nucleic acid
def rmsds(
    trajectory_file : 'File',
    first_frame_file : 'File',
    average_structure_file : 'File',
    output_analysis_filepath : str,
    snapshots : int,
    frames_limit : int,
    structure : 'Structure',
    pbc_residues: List[int],
    ligand_map : List[dict],
    selections : List[str] = ['protein', 'nucleic'],
    ):

    # Find PBC residues, which are to be removed from parsed selections
    pbc_selection = structure.select_residue_indices(pbc_residues)

    # Parse the selections to meaningfull atom indices
    parsed_selections = { selection: structure.select(selection, syntax='vmd') for selection in selections }

    # If there is a ligand map then parse them to selections as well
    if ligand_map:
        for ligand in ligand_map:
            selection_name = 'ligand ' + ligand['name']
            parsed_selection = structure.select_residue_indices(ligand['residue_indices'])
            parsed_selections[selection_name] = parsed_selection

    # Remove PBC residues from parsed selections
    pbc_selection = structure.select_residue_indices(pbc_residues)
    non_pbc_selections = {}
    for selection_name, selection in parsed_selections.items():
        # If selection was empty from the begining then discard it
        if not selection:
            continue
        # Substract PBC atoms
        non_pbc_selection = selection - pbc_selection
        # If selection after substracting pbc atoms becomes empty then discard it
        if not non_pbc_selection:
            continue
        # Add the the filtered selection to the dict
        non_pbc_selections[selection_name] = non_pbc_selection

    # The start will be always 0 since we start with the first frame
    start = 0

    # Reduce the trajectory according to the frames limit
    # Use a reduced trajectory in case the original trajectory has many frames
    # Note that it makes no difference which reference is used here
    reduced_trajectory_filepath, step, frames = get_reduced_trajectory(
        first_frame_file,
        trajectory_file,
        snapshots,
        frames_limit,
    )

    # Save results in this array
    output_analysis = []

    # Set the reference structures to run the RMSD against
    rmsd_references = [first_frame_file, average_structure_file]

    # Iterate over each reference and group
    for reference in rmsd_references:
        # Get a standarized reference name
        reference_name = REFERENCE_LABELS[reference.filename]
        for group_name, group_selection in non_pbc_selections.items():
            # Set the analysis filename
            rmsd_analysis = 'rmsd.' + reference_name + '.' + group_name.lower() + '.xvg'
            # Run the rmsd
            rmsd(reference.path, reduced_trajectory_filepath, group_selection, rmsd_analysis)
            # Read and parse the output file
            rmsd_data = xvg_parse(rmsd_analysis, ['times', 'values'])
            # Format the mined data and append it to the overall output
            # Multiply by 10 since rmsd comes in nanometers (nm) and we want it in Ångstroms (Å)
            rmsd_values = [ v*10 for v in rmsd_data['values'] ]
            data = {
                'values': rmsd_values,
                'reference': reference_name,
                'group': group_name
            }
            output_analysis.append(data)
            # Remove the analysis xvg file since it is not required anymore
            os.remove(rmsd_analysis)

    # Export the analysis in json format
    save_json({ 'start': start, 'step': step, 'data': output_analysis }, output_analysis_filepath)

# RMSD
# 
# Perform the RMSd analysis 
def rmsd (
    reference_filepath : str,
    trajectory_filepath : str,
    selection : 'Selection', # This selection will never be empty, since this is checked previously
    output_analysis_filepath : str):

    # Convert the selection to a ndx file gromacs can read
    selection_name = 'rmsd_selection'
    ndx_selection = selection.to_ndx(selection_name)
    ndx_filename = '.rmsd.ndx'
    with open(ndx_filename, 'w') as file:
        file.write(ndx_selection)
    
    # Run Gromacs
    p = Popen([
        "echo",
        selection_name, # Select group for least squares fit
        selection_name, # Select group for RMSD calculation
    ], stdout=PIPE)
    process = run([
        GROMACS_EXECUTABLE,
        "rms",
        "-s",
        reference_filepath,
        "-f",
        trajectory_filepath,
        '-o',
        output_analysis_filepath,
        '-n',
        ndx_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE, stderr=PIPE)
    # Consuming the output makes the process run
    output_logs = process.stdout.decode()
    # Close the input
    p.stdout.close()

    # If the output does not exist at this point it means something went wrong with gromacs
    if not os.path.exists(output_analysis_filepath):
        print(output_logs)
        error_logs = process.stderr.decode()
        print(error_logs)
        raise Exception('Something went wrong with GROMACS')

    # Remove the ndx file
    os.remove(ndx_filename)