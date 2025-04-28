from model_workflow.tools.xvg_parse import xvg_parse
from model_workflow.tools.get_reduced_trajectory import get_reduced_trajectory
from model_workflow.utils.auxiliar import save_json
from model_workflow.utils.constants import GROMACS_EXECUTABLE, REFERENCE_LABELS
from model_workflow.utils.type_hints import *

import os
from subprocess import run, PIPE, Popen

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
    pbc_selection : 'Selection',
    ligand_map : List[dict],
    ):

    print('-> Running RMSDs analysis')

    # Set the default selections to be analyzed
    default_selections = {
        'protein': structure.select_protein(),
        'nucleic': structure.select_nucleic()
    }

    # Set selections to be analyzed
    selections = { **default_selections }

    # If there is a ligand map then parse them to selections as well
    for ligand in ligand_map:
        selection_name = 'ligand ' + ligand['name']
        selection = structure.select_residue_indices(ligand['residue_indices'])
        # If the ligand has less than 3 atoms then gromacs can not fit it so it will fail
        if len(selection) < 3: continue
        # Add current ligand selection to be analyzed
        selections[selection_name] = selection
        # If the ligand selection totally overlaps with a default selection then remove the default
        for default_selection_name, default_selection in default_selections.items():
            if default_selection == selection: del selections[default_selection_name]

    # Remove PBC residues from parsed selections
    non_pbc_selections = {}
    for selection_name, selection in selections.items():
        # Substract PBC atoms
        non_pbc_selection = selection - pbc_selection
        # If selection after substracting pbc atoms becomes empty then discard it
        if not non_pbc_selection:
            continue
        # Add the the filtered selection to the dict
        non_pbc_selections[selection_name] = non_pbc_selection

    # If there is nothing lo analyze at this point then skip the analysis
    if len(non_pbc_selections) == 0:
        print('  The RMSDs analysis will be skipped since there is nothing to analyze')
        return
    
    # Get the coarse grain selection
    cg_selection = structure.select_cg()

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
            rmsd_analysis = f'rmsd.{reference_name}.{group_name.lower()}.xvg'
            # If part of the selection has coarse grain atoms then skip mass weighting
            # Otherwise Gromacs may fail since the atom name may not be in the atommass library
            has_cg = group_selection & cg_selection
            # Run the rmsd
            print(f' Reference: {reference_name}, Selection: {group_name},{" NOT" if has_cg else ""} mass weighted')
            rmsd(reference.path, reduced_trajectory_filepath, group_selection, rmsd_analysis, skip_mass_weighting=has_cg)
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
    output_analysis_filepath : str,
    skip_mass_weighting : bool = False):

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
        # Supress mass weighting if requested
        *([ '-mw', 'no' ] if skip_mass_weighting else []),
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