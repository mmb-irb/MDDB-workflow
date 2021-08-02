from model_workflow.tools.xvg_parse import xvg_parse
from model_workflow.tools.get_pdb_frames import get_pdb_frames

import os
import json
import numpy
from subprocess import run, PIPE, Popen

# This is a residual file produced by the sasa analysis
# It must be deleted after each
area_filename = 'area.xvg'

# Perform the Solvent Accessible Surface Analysis
def sasa(
    input_topology_filename: str,
    input_trajectory_filename: str,
    output_analysis_filename: str,
    reference,
):

    # Set indexes to select the system without hydrogens
    indexes = 'indexes.ndx'

    # Remove all groups but 0
    remove_others = 'keep 0' + '\n'
    selection_without_hydrogens = '0 & !a H*' + '\n'
    remove_system = 'del 0' + '\n'
    all_commands = remove_others + selection_without_hydrogens + remove_system + 'q'

    p = Popen([
        "echo",
        all_commands,
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "make_ndx",
        "-f",
        input_topology_filename,
        '-o',
        indexes,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

    # Calculate the sasa for each frame
    sasa_per_frame = []
    frames_limit = 200
    frames, step, count = get_pdb_frames(input_topology_filename, input_trajectory_filename, frames_limit)
    for f, current_frame in enumerate(frames):

        # Run the sasa analysis over the current frame
        # WARNING: We want the SASA per residue, and this could be obtained replacing '-oa' per -'or'
        # WARNING: However, residues are not enumerated the same in Gromacs and other tools (e.g. pytraj)
        # For this reason we must always rely on atom numeration, which is the same along different tools
        current_frame_sasa = 'sasa' + str(f) + '.xvg'
        logs = run([
            "gmx",
            "sasa",
            "-s",
            current_frame,
            '-oa',
            current_frame_sasa,
            "-n",
            indexes,
            "-surface",
            "0",
            '-quiet'
        ], stdout=PIPE).stdout.decode()

        # Mine the sasa results (.xvg file)
        # Hydrogen areas are not recorded in the xvg file
        sasa = xvg_parse(current_frame_sasa, ['n', 'area', 'sd'])
        # Restructure data by adding all atoms sas per residue
        atom_numbers = sasa['n']
        atom_areas = sasa['area']
        sas_per_residues = [0.0] * len(reference.residues)
        for atom_number, atom_area in zip(atom_numbers, atom_areas):
            atom_index = int(atom_number) - 1
            residue_index = reference.get_atom_residue_index(atom_index)
            sas_per_residues[residue_index] += atom_area
        sasa_per_frame.append(sas_per_residues)
        # Delete current frame files before going for the next frame
        run([
            "rm",
            current_frame_sasa,
            area_filename,
        ], stdout=PIPE).stdout.decode()

    # Format output data
    # Sasa values must be separated by residue and then ordered by frame
    data = []
    for r, residue in enumerate(reference.residues):
        # Name the residue in the source format
        name = reference.get_residue_name(residue)
        atom_count = residue.numAtoms()
        # Harvest its sasa along each frame
        saspf = []
        for frame in sasa_per_frame:
            # IMPORTANT: The original SASA value is modified to be normalized
            # We divide the value by the number of atoms
            frame_sas = frame[r]
            normalized_frame_sas = frame_sas / atom_count
            # To make is standard with the rest of analyses we pass the results from nm² to A²
            standard_frame_sas = normalized_frame_sas * 100
            saspf.append(standard_frame_sas)
        # Calculate the mean and standard deviation of the residue sasa values
        mean = numpy.mean(saspf)
        stdv = numpy.std(saspf)
        data.append({
            'name': name,
            'saspf': saspf,
            'mean': mean,
            'stdv': stdv
        })

    # Export the analysis in json format
    with open(output_analysis_filename, 'w') as file:
        json.dump({'data': data}, file)

# Read a pdb file and return the residues which gromacs would consider
# This is used for those exotic topologies where the gromacs amount of residues does not match prody's
# def get_gromacs_residues (pdb_filename : str):
#     residues = []
#     # Read the topology line by line and get all residue data (name, chain, number and icode)
#     with open(pdb_filename, "r") as file:
#         lines = file.readlines()
#         last_residue = None
#         for line in lines:
#             # Skip non atom lines
#             if line[0:4] != 'ATOM':
#                 continue
#             # Get the part of the line where the residue data is found
#             residue = line[17:31]
#             # If it is the same residue than the previous line then skip it
#             if residue == last_residue:
#                 continue
#             # Otherwise, record it
#             residues.append(residue)
#             last_residue = residue
#     return residues