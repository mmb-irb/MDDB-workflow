from model_workflow.tools.xvg_parse import xvg_parse
from model_workflow.tools.get_reduced_trajectory import get_reduced_trajectory

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
        snapshots: int):

    # Set the frames where we calculate the sasa
    # WARNING: The gromacs '-fr' option counts frames starting at 1, not at 0
    frames = range(1, snapshots +1)

    # Set a maximum of frames
    # If trajectory has more frames than the limit create a reduced trajectory
    reduced_trajectory_filename = input_trajectory_filename
    frames_number = 200
    if snapshots > frames_number:
        # WARNING: The gromacs '-fr' option counts frames starting at 1, not at 0
        frames = range(1, frames_number +1)  # if frames_number > 1 else [1]
        reduced_trajectory_filename = 'f' + str(frames_number) + '.trajectory.xtc'
        if not os.path.exists(pockets_trajectory):
            get_reduced_trajectory(
                input_topology_filename,
                input_trajectory_filename,
                reduced_trajectory_filename,
                snapshots,
                frames_number,
            )
    else:
        frames_number = snapshots

    # Set indexes to select the system without hydrogens
    indexes = 'indexes.ndx'
    p = Popen([
        "echo",
        "0 & !a H*\nq",
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
    frames_ndx = 'frames.ndx'
    sasa_per_frame = []
    for f in frames:
        print('Frame ' + str(f) + ' / ' + str(frames_number))
        # Extract the current frame
        current_frame = 'frame' + str(f) + '.pdb'
        # The frame selection input in gromacs works with a 'ndx' file
        with open(frames_ndx, 'w') as file:
            file.write('[frames]\n' + str(f))
        p = Popen([
            "echo",
            "System_&_!H*",
        ], stdout=PIPE)
        logs = run([
            "gmx",
            "trjconv",
            "-s",
            input_topology_filename,
            "-f",
            reduced_trajectory_filename,
            '-o',
            current_frame,
            "-fr",
            frames_ndx,
            "-n",
            indexes,
            '-quiet'
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()
        # Run the sasa analysis over the current frame
        current_frame_sasa = 'sasa' + str(f) + '.xvg'
        logs = run([
            "gmx",
            "sasa",
            "-s",
            current_frame,
            "-surface",
            '0',
            '-or',
            current_frame_sasa,
            '-quiet'
        ], stdout=PIPE).stdout.decode()
        # Mine the sasa results (.xvg file)
        sasa = xvg_parse(current_frame_sasa, ['c1', 'c2', 'c3'])
        sasa_per_frame.append(sasa['c2'])
        # Delete current frame files before going for the next frame
        run([
            "rm",
            frames_ndx,
            current_frame,
            current_frame_sasa,
            area_filename,
        ], stdout=PIPE).stdout.decode()

    # Check that the number of sasa values per frame is the same that the number of residues
    sasa_residues_count = len(sasa_per_frame[0])
    if sasa_residues_count != len(reference.residues):
        # If the number does not match we may have a problem
        print('sasa residues: ' + str(sasa_residues_count))
        print('reference residues: ' + str(len(reference.residues)))
        print('WARNING: The reference number of residues does not match the SASA analysis')
        # Try to reconfigure the sasa residues to match the reference
        print('Trying alternative residues configuration...')
        gromacs_residues = get_gromacs_residues(input_topology_filename)
        unique_gromacs_residues = list(set(gromacs_residues))
        # If the new residues configuration neither matches the sasa analysis then surrender
        if sasa_residues_count != len(gromacs_residues) or len(reference.residues) != len(unique_gromacs_residues):
            print('gromacs residues: ' + str(len(gromacs_residues)))
            raise SystemExit('ERROR: The number of residues does not match in SASA analysis')
        print("It's OK")
        # If the new residues configuration matches the sasa analysis
        # Then add sasa values for identical residues
        reconfigured_sasa_per_frame = []
        for sasa in sasa_per_frame:
            reconfigured_sasa = []
            for unique_residue in unique_gromacs_residues:
                sasa_values = [ value for i, value in enumerate(sasa) if gromacs_residues[i] == unique_residue ]
                new_sasa_value = sum(sasa_values)
                reconfigured_sasa.append(new_sasa_value)
            reconfigured_sasa_per_frame.append(reconfigured_sasa)
        sasa_per_frame = reconfigured_sasa_per_frame

    # Format output data
    # Sasa values must be separated by residue and then ordered by frame
    data = []
    for r, residue in enumerate(reference.residues):
        # Name the residue in the source format
        name = reference.get_residue_name(residue)
        # Harvest its sasa along each frame
        saspf = []
        for frame in sasa_per_frame:
            saspf.append(frame[r])
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
# This is used for those exotic topologies where the gromacs amout of residues does not match prody's
def get_gromacs_residues (pdb_filename : str):
    residues = []
    # Read the topology line by line and get all residue data (name, chain, number and icode)
    with open(pdb_filename, "r") as file:
        lines = file.readlines()
        last_residue = None
        for line in lines:
            # Skip non atom lines
            if line[0:4] != 'ATOM':
                continue
            # Get the part of the line where the residue data is found
            residue = line[17:31]
            # If it is the same residue than the previous line then skip it
            if residue == last_residue:
                continue
            # Otherwise, record it
            residues.append(residue)
            last_residue = residue
    return residues