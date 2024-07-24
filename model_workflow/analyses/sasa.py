from model_workflow.tools.xvg_parse import xvg_parse
from model_workflow.tools.get_pdb_frames import get_pdb_frames
from model_workflow.utils.auxiliar import save_json
from model_workflow.utils.constants import GROMACS_EXECUTABLE

import os
import numpy
from subprocess import run, PIPE, Popen
from typing import List

# This is a residual file produced by the sasa analysis
# It must be deleted after each
area_filename = 'area.xvg'

# Perform the Solvent Accessible Surface Analysis
def sasa(
    input_topology_filename: str,
    input_trajectory_filename: str,
    output_analysis_filename: str,
    structure : 'Structure',
    pbc_residues : List[int],
    snapshots : int,
    frames_limit : int,
):

    print('-> Running SAS analysis')

    # For this analysis me must filter out hydrogens
    heavy_atoms_selection = '( not name "H.*" )'

    # Execute the atom selection over the structure
    structure_selection = structure.select(heavy_atoms_selection, syntax='vmd')

    # At this point the number of residues between the original and the non-hydrogens structure should match
    # We rely on this, so we must check
    filtered_structure = structure.filter(structure_selection)
    filtered_residues_count = len(filtered_structure.residues)
    original_residues_count = len(structure.residues)
    if filtered_residues_count != original_residues_count:
        raise ValueError('The number of residues does not match after filtering out hydrogens')

    # Convert the structure selection to a ndx file
    selection_name = 'sasa'
    ndx_selection = structure_selection.to_ndx(selection_name)
    ndx_filename = 'indices.ndx'
    with open(ndx_filename, 'w') as file:
        file.write(ndx_selection)

    # We must exclude PBC residues (e.g. membranes) from sasa results
    # PBC residues close to the boundary will always have unrealistic sasa values
    # WARNING: These results must be exlucded from the analysis afterwards, but the membrane can not be excluded from the structure
    # Otherwise, those residues which are covered by the membrane would be exposed to the solvent during the analysis
    skipped_residue_indices = pbc_residues

    # Calculate the sasa for each frame
    sasa_per_frame = []
    frames, step, count = get_pdb_frames(input_topology_filename, input_trajectory_filename, snapshots, frames_limit)
    for f, current_frame in enumerate(frames):

        # Run the sasa analysis over the current frame
        # WARNING: We want the SASA per residue, and this could be obtained replacing '-oa' per -'or'
        # WARNING: However, residues are not enumerated the same in Gromacs and other tools (e.g. pytraj)
        # For this reason we must always rely on atom numeration, which is the same along different tools
        current_frame_sasa = 'sasa' + str(f) + '.xvg'
        logs = run([
            GROMACS_EXECUTABLE,
            "sasa",
            "-s",
            current_frame,
            '-oa',
            current_frame_sasa,
            "-n",
            ndx_filename,
            "-surface",
            "0",
            '-quiet'
        ], stdout=PIPE, stderr=PIPE).stdout.decode()

        # In case the output file does not exist at this point it means something went wrong with Gromacs
        if not os.path.exists(current_frame_sasa):
            print(logs)
            raise SystemExit('Something went wrong with Gromacs')

        # Mine the sasa results (.xvg file)
        # Hydrogen areas are not recorded in the xvg file
        sasa = xvg_parse(current_frame_sasa, ['n', 'area', 'sd'])
        # Restructure data by adding all atoms sas per residue
        atom_numbers = sasa['n']
        atom_areas = sasa['area']
        sas_per_residues = [0.0] * len(structure.residues)
        for atom_number, atom_area in zip(atom_numbers, atom_areas):
            atom_index = int(atom_number) - 1
            atom = structure.atoms[atom_index]
            residue_index = atom.residue_index
            sas_per_residues[residue_index] += atom_area
        sasa_per_frame.append(sas_per_residues)
        # Delete current frame files before going for the next frame
        run([
            "rm",
            current_frame_sasa,
            area_filename,
        ], stdout=PIPE).stdout.decode()

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
    }

    # Export the analysis in json format
    save_json(output, output_analysis_filename)