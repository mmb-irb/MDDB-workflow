from mddb_workflow.tools.xvg_parse import xvg_parse
from mddb_workflow.tools.get_pdb_frames import get_pdb_frames
from mddb_workflow.utils.auxiliar import save_json
from mddb_workflow.utils.constants import OUTPUT_SASA_FILENAME
from mddb_workflow.utils.gmx_spells import run_gromacs
from mddb_workflow.utils.type_hints import *

import os
import numpy

# This is a residual file produced by the sasa analysis
# It must be deleted after each
RESIDUAL_AREA_FILENAME = 'area.xvg'

def sasa(
    structure_file: 'File',
    trajectory_file: 'File',
    output_directory: str,
    structure : 'Structure',
    pbc_residues : list[int],
    snapshots : int,
    frames_limit : int = 100,
):
    """Perform the Solvent Accessible Surface Analysis."""

    # If all residues are to be excluded since the whole system is in PCB then stop here
    if len(pbc_residues) == len(structure.residues):
        print(' No residues to run the analysis')
        return

    # Set the output filepath
    output_analysis_filepath = f'{output_directory}/{OUTPUT_SASA_FILENAME}'

    # Get the solvent selection so we can further exclude it
    # Note that the solvent would ironically cover the SAS regions
    solvent_selection = structure.select_water_and_counter_ions()
    # Select heavy atoms only thus skipping hydrogens
    sasa_selection = structure.select_heavy_atoms() - solvent_selection

    # Convert the structure selection to a ndx file
    sasa_selection_name = 'sasa'
    sasa_ndx_selection = sasa_selection.to_ndx(sasa_selection_name)
    ndx_filename = f'{output_directory}/indices.ndx'
    with open(ndx_filename, 'w') as file:
        file.write(sasa_ndx_selection)

    # We must exclude PBC residues (e.g. membranes) from sasa results
    # PBC residues close to the boundary will always have unrealistic sasa values
    # WARNING: These results must be exlucded from the analysis afterwards, but the membrane can not be excluded from the structure
    # Otherwise, those residues which are covered by the membrane would be exposed to the solvent during the analysis
    skipped_residue_indices = pbc_residues

    # Calculate the sasa for each frame
    sasa_per_frame = []
    frames, step, count = get_pdb_frames(structure_file.path, trajectory_file.path, snapshots, frames_limit)
    for f, current_frame in enumerate(frames):

        # Run the sasa analysis over the current frame
        # WARNING: We want the SASA per residue, and this could be obtained replacing '-oa' per -'or'
        # WARNING: However, residues are not enumerated the same in Gromacs and other tools (e.g. pytraj)
        # For this reason we must always rely on atom numeration, which is the same along different tools
        # Note that the '-surface' argument sets what atoms are present during the calculation
        # Thus we must use this arguments to remove solvent, but not PBC regions
        # Set the main output filepath
        current_frame_sasa = f'{output_directory}/sasa{f}.xvg'
        # Set a residual output filepath
        residual_area_filepath = f'{output_directory}{RESIDUAL_AREA_FILENAME}'
        # If not defined, then by default it creates a 'area.xvg' file where we are running the workflow
        run_gromacs(f'sasa -s {current_frame} -oa {current_frame_sasa} -o {residual_area_filepath}\
            -n {ndx_filename} -surface {sasa_selection_name}', expected_output_filepath = current_frame_sasa)

        # Mine the sasa results (.xvg file)
        # Areas from excluded atoms are not recorded in the xvg file
        sasa = xvg_parse(current_frame_sasa, ['n', 'area', 'sd'])
        # Restructure data by adding all atoms sas per residue
        atom_numbers = sasa['n']
        atom_areas = sasa['area']
        sas_per_residues = [0.0] * len(structure.residues)
        for atom_number, atom_area in zip(atom_numbers, atom_areas):
            atom_index = atom_number - 1
            atom = structure.atoms[atom_index]
            residue_index = atom.residue_index
            sas_per_residues[residue_index] += atom_area
        sasa_per_frame.append(sas_per_residues)
        # Delete current frame files before going for the next frame
        os.remove(current_frame_sasa)
        os.remove(residual_area_filepath)

    # Remove the indices file since it is not required anymore
    os.remove(ndx_filename)

    # Format output data
    # Sasa values must be separated by residue and then ordered by frame
    saspf = []
    means = []
    stdvs = []
    for r, residue in enumerate(structure.residues):
        # Skip membrane residues
        if r in skipped_residue_indices:
            saspf.append(None)
            means.append(None)
            stdvs.append(None)
            continue
        # Get the number of atoms in current residue
        atom_count = len(residue.atoms)
        # Harvest its sasa along each frame
        residue_saspf = []
        for frame in sasa_per_frame:
            # IMPORTANT: The original SASA value is modified to be normalized
            # We divide the value by the number of atoms
            frame_sas = frame[r]
            normalized_frame_sas = frame_sas / atom_count
            # To make is standard with the rest of analyses we pass the results from nm² to A²
            standard_frame_sas = normalized_frame_sas * 100
            residue_saspf.append(standard_frame_sas)
        # Add current resiude sas per frame to the overall list
        saspf.append(residue_saspf)
        # Calculate the mean and standard deviation of the residue sasa values
        mean = numpy.mean(residue_saspf)
        means.append(mean)
        stdv = numpy.std(residue_saspf)
        stdvs.append(stdv)
    # Set the content of the analysis output file
    output = {
        'step': step,
        'saspf': saspf,
        'means': means,
        'stdvs': stdvs,
        'version': '0.0.1',
    }

    # Export the analysis in json format
    save_json(output, output_analysis_filepath)
