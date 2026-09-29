from mddb_workflow.tools.xvg_parse import xvg_parse
from mddb_workflow.utils.auxiliar import save_json, get_auxiliar_filepath
from mddb_workflow.utils.constants import REFERENCE_LABELS, OUTPUT_RMSDS_FILENAME
from mddb_workflow.utils.gmx_spells import run_gromacs
from mddb_workflow.utils.type_hints import *

import numpy as np
from math import ceil
from os import remove


def rmsds(
    trajectory_file : 'File',
    first_frame_file : 'File',
    average_structure_file : 'File',
    output_directory : str,
    snapshots : int,
    structure : 'Structure',
    pbc_selection : 'Selection',
    cg_selection : 'Selection',
    dummy_selection : 'Selection',
    inchikey_map : list[dict],
    # Number of splits along the trajectory
    time_splits : int = 100,
    ) -> Optional[tuple[float, float, float]]:
    """Run multiple RMSD analyses. One with each reference (first frame, average structure)
    and each selection (default: global, protein, nucleic).
    RMSD is calculated for every frame in the trajectory but results are summarized in time splits.
    For each split we store the mean, standard deviation and maximum RMSD.
    Return the global RMSD against the first frame along the whole trajectory as (mean, stdv, max).
    Return None if there is no global selection."""
    # Set the main output filepath
    output_analysis_filepath = f'{output_directory}/{OUTPUT_RMSDS_FILENAME}'

    # Set the default selections to be analyzed
    default_selections = {
        'protein': structure.select_protein(),
        'nucleic': structure.select_nucleic()
    }

    # Set selections to be analyzed
    selections = { **default_selections }

    # If there is a ligand map then parse them to selections as well
    for ligand in inchikey_map:
        if ligand['is_lipid']: continue
        selection_name = 'ligand ' + ligand['name']
        selection = structure.select_residue_indices(ligand['residue_indices'])
        # If the ligand has less than 3 atoms then gromacs can not fit it so it will fail
        if len(selection) < 3: continue
        # Add current ligand selection to be analyzed
        selections[selection_name] = selection
        # If the ligand selection totally overlaps with a default selection then remove the default
        for default_selection_name, default_selection in default_selections.items():
            if default_selection == selection: del selections[default_selection_name]

    # If there is anything left apart from proteins, nucleic acids and ligands then select it
    missing_regions = structure.select_all()
    for selection in selections.values():
        missing_regions -= selection
    # WARNING: Atom selections with less than 3 atoms will raise an error from Gromacs
    if missing_regions and len(missing_regions.atom_indices) >= 3:
        selections['other'] = missing_regions

    # Always analyze the whole system (excluding PBC atoms) as a global selection
    non_pbc_selections: dict[str, 'Selection'] = {}
    global_selection = structure.select_all() - pbc_selection
    # WARNING: Atom selections with less than 3 atoms will raise an error from Gromacs
    if len(global_selection) >= 3:
        non_pbc_selections['global'] = global_selection

    # Remove PBC residues from parsed selections
    for selection_name, selection in selections.items():
        # Substract PBC atoms
        non_pbc_selection = selection - pbc_selection
        # If selection after substracting pbc atoms becomes empty then discard it
        if not non_pbc_selection: continue
        # If the selection is identical to the global selection then discard it, since it would be redundant
        if non_pbc_selection == global_selection: continue
        # Add the the filtered selection to the dict
        non_pbc_selections[selection_name] = non_pbc_selection

    # If there is nothing lo analyze at this point then skip the analysis
    if len(non_pbc_selections) == 0:
        print('  The RMSDs analysis will be skipped since there is nothing to analyze')
        return

    # The start will be always 0 since we start with the first frame
    start = 0

    # Calculate how many frames fall in every time split
    step = ceil(snapshots / time_splits)
    # Calculate how many time splits we will have at the end
    # Note that the last split may have less frames than the rest
    nsteps = ceil(snapshots / step)

    # Save results in this array
    output_analysis = []
    # Global RMSD (mean, stdv, max) along the whole trajectory, to be returned
    global_rmsd = None

    # Set the reference structures to run the RMSD against
    rmsd_references = [first_frame_file, average_structure_file]

    # Iterate over each reference and group
    for reference in rmsd_references:
        # Get a standarized reference name
        reference_name = REFERENCE_LABELS[reference.filename]
        for group_name, group_selection in non_pbc_selections.items():
            # Set the analysis filename
            rmsd_analysis_filename = f'rmsd.{reference_name}.{group_name.lower().replace(" ","_")}.xvg'
            rmsd_analysis_filepath = get_auxiliar_filepath(rmsd_analysis_filename)
            # If part of the selection has coarse grain or dummy atoms atoms then skip mass weighting
            # Note that we may have already automatically added masses to the gromas atommass file
            # Thus we could make it work, but the values may be not real since masses are not real
            has_cg = group_selection & cg_selection
            has_dummy = group_selection & dummy_selection
            mass_weighted = not has_cg and not has_dummy
            # Run the rmsd
            print(f' Reference: {reference_name}, Selection: {group_name},{"" if mass_weighted else " NOT"} mass weighted')
            rmsd(reference.path, trajectory_file.path, group_selection, rmsd_analysis_filepath, skip_mass_weighting=not mass_weighted)
            # Read and parse the output file
            rmsd_data = xvg_parse(rmsd_analysis_filepath, ['times', 'values'])
            # Multiply by 10 since rmsd comes in nanometers (nm) and we want it in Ångstroms (Å)
            rmsd_values = np.array(rmsd_data['values']) * 10
            if len(rmsd_values) != snapshots:
                raise ValueError(f'Number of RMSD values ({len(rmsd_values)}) does not match the number of snapshots ({snapshots})')
            # Split the values in time splits and summarize each split
            splits = [ rmsd_values[i*step:(i+1)*step] for i in range(nsteps) ]
            data = {
                'means': [ float(split.mean()) for split in splits ],
                'stdvs': [ float(split.std()) for split in splits ],
                'maxs': [ float(split.max()) for split in splits ],
                # Frame number (1-based) of the maximum RMSD in each split
                'maxframes': [ i*step + int(split.argmax()) + 1 for i, split in enumerate(splits) ],
                'reference': reference_name,
                'group': group_name,
                'massw': mass_weighted,
            }
            output_analysis.append(data)
            # Save the overall global RMSD against the first frame to be returned
            if reference == first_frame_file and group_name == 'global':
                global_rmsd = (float(rmsd_values.mean()), float(rmsd_values.std()), float(rmsd_values.max()))
            # Remove the analysis xvg file since it is not required anymore
            remove(rmsd_analysis_filepath)

    # Export the analysis in json format
    save_json({
        'start': start,
        'step': step,
        'data': output_analysis,
        'version': '1.0.0',
    }, output_analysis_filepath)

    return global_rmsd

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
    ndx_filepath = get_auxiliar_filepath('.rmsd.ndx')
    with open(ndx_filepath, 'w') as file:
        file.write(ndx_selection)

    # Run Gromacs
    run_gromacs(f'rms -s {reference_filepath} -f {trajectory_filepath} \
        -o {output_analysis_filepath} -n {ndx_filepath} {"-mw no" if skip_mass_weighting else ""}',
        user_input = f'{selection_name} {selection_name}')

    # Remove the ndx file
    remove(ndx_filepath)
