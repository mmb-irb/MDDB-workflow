# Pockets analysis
# 
# The pockets analysis is carried by the 'Fpocket' software. Fpocket is a fast open source protein 
# pocket (cavity) detection algorithm based on Voronoi tessellation.
#
# The analysis is performed over 100 frames selected along the trajectory. First of all, an 'MDpocket'
# analysis is run over these frames thus generating an output grid file. This grid contains a measure 
# of frequency of how many times the pocket was open during the trajectory. The most stable and wide 
# pockets are selected from this grid and then analyzed again with MDpocket individually. This second 
# run returns dynamic data of the pocket volume and drugability score, among others.
#
# Vincent Le Guilloux, Peter Schmidtke and Pierre Tuffery, “Fpocket: An open source platform for ligand 
# pocket detection”, BMC Bioinformatics 2009, 10:168
#
# Peter Schmidtke & Xavier Barril “Understanding and predicting druggability. A high-throughput method 
# for detection of drug binding sites.”, J Med Chem, 2010, 53(15):5858-67
#
# Peter Schmldtke, Axel Bidon-Chanal, Javier Luque, Xavier Barril, “MDpocket: open-source cavity detection 
# and characterization on molecular dynamics trajectories.”, Bioinformatics. 2011 Dec 1;27(23):3276-85

import re
import collections
import math

from subprocess import run, PIPE, Popen

import json

from model_workflow.tools.get_reduced_trajectory import get_reduced_trajectory

# Perform the pockets analysis
def pockets (
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_analysis_filename : str,
    topology_reference,
    snapshots : int ):

    # This anlaysis produces many useless output files
    # Create a new folder to store all ouput files so they do not overcrowd the main directory
    mdpocket_folder = 'mdpocket'
    logs = run([
        "mkdir",
        mdpocket_folder,
    ], stdout=PIPE).stdout.decode()

    # Set a reduced trajectory with only 100 frames
    # Get the step between frames of the new reduced trajectory, since it will be append to the output
    pockets_trajectory = input_trajectory_filename
    frames_number = 100
    step = 1
    if snapshots > frames_number:
        pockets_trajectory = 'pockets.trajectory.xtc'
        step = get_reduced_trajectory(
            input_topology_filename,
            input_trajectory_filename,
            pockets_trajectory,
            snapshots,
            frames_number,
        )

    # Run the mdpocket analysis focusing in this specific pocket
    mdpocket_output = mdpocket_folder + '/mdpout'
    logs = run([
        "mdpocket",
        "--trajectory_file",
        pockets_trajectory,
        "--trajectory_format",
        "xtc",
        "-f",
        input_topology_filename,
        "-o",
        mdpocket_output,
    ], stdout=PIPE).stdout.decode()

    # Read and harvest the gird file
    grid_filename = mdpocket_output + '_freq.dx'
    with open(grid_filename,'r') as file:

        # First, mine the grid dimensions and origin
        dimensions = None
        dimensions_pattern = "counts ([0-9]+) ([0-9]+) ([0-9]+)"
        origin = None
        origin_pattern = "origin ([.0-9-]+) ([.0-9-]+) ([.0-9-]+)"
        for line in file:
            search = re.search(dimensions_pattern, line)
            if search != None:
                dimensions = (int(search.group(1)),int(search.group(2)),int(search.group(3)))
            search = re.search(origin_pattern, line)
            if search != None:
                origin = (float(search.group(1)),float(search.group(2)),float(search.group(3)))
            if origin and dimensions:
                break

        # Next, mine the grid values
        grid_values = []
        grid_values_pattern = "^([.0-9]+) ([.0-9]+) ([.0-9]+) $"
        # The last line of the grid may have less than 3 values
        last_line_grid_values_pattern_1 = "^([.0-9]+) ([.0-9]+) $"
        last_line_grid_values_pattern_2 = "^([.0-9]+) $" 
        for line in file:
            search = re.search(grid_values_pattern, line)
            if (search != None):
                grid_values.append(float(search.group(1)))
                grid_values.append(float(search.group(2)))
                grid_values.append(float(search.group(3)))
                continue
            # This should only happend at the end
            search = re.search(last_line_grid_values_pattern_1, line)
            if (search != None):
                grid_values.append(float(search.group(1)))
                grid_values.append(float(search.group(2)))
                continue
            search = re.search(last_line_grid_values_pattern_2, line)
            if (search != None):
                grid_values.append(float(search.group(1)))
                continue

    # Set a function to get the value of a given 'x, y, z' position
    # Grid values are disposed first in 'z' order, then 'y', and finally 'x'
    xl, yl, zl = dimensions
    def getIndex (x, y, z):
        index = x*yl*zl + y*zl + z
        return index
    def getValue (x, y, z):
        index = getIndex(x,y,z)
        return grid_values[index]

    # Classify all non-zero values by 'pocket' groups according to if they are connected
    pockets = [0] * xl * yl * zl
    pockets_count = 0
    # Save also each point coordinates
    point_coordinates = [None] * xl * yl * zl
    # Set a cuttoff value to consider a point valid
    cuttoff = 0.5
    for x in range(xl):
        for y in range(yl):
            for z in range(zl):
                index = getIndex(x,y,z)
                value = grid_values[index]
                # If it is lower than the cutoff then pass
                if cuttoff and value < cuttoff:
                    continue
                # If the value exists, save the coordinates
                point_coordinates[index] = (x,y,z)
                # Then look for connected values in previously searched values
                # If any of these connected values is marked as a pocket we mark this value with the same number
                pz = z - 1
                if pz >= 0:
                    pIndex = getIndex(x,y,pz)
                    pPocket = pockets[pIndex]
                    if pPocket:
                        pockets[index] = pPocket
                        continue
                py = y - 1
                if py >= 0:
                    pIndex = getIndex(x,py,z)
                    pPocket = pockets[pIndex]
                    if pPocket:
                        pockets[index] = pPocket
                        continue
                px = x - 1
                if px >= 0:
                    pIndex = getIndex(px,y,z)
                    pPocket = pockets[pIndex]
                    if pPocket:
                        pockets[index] = pPocket
                        continue
                # If none of the previous values was marked as a pocket then we set a new pocket number
                pockets_count += 1
                pockets[index] = pockets_count

    # Get the first 10 pockets with more points
    pockets_number = 10

    # Exclude the first result which will always be 0 and it stands for no-pocket points
    biggest_pockets = collections.Counter(pockets).most_common()
    if len(biggest_pockets) > pockets_number:
        biggest_pockets = collections.Counter(pockets).most_common()[1:11]

    # First of all, get all header lines from the original grid file
    # We need them to write grid files further
    with open(grid_filename,'r') as file:
        header_lines = []
        header_pattern = "^[a-z]"
        for line in file:
            if re.match(header_pattern, line):
                header_lines.append(line)
            elif re.match(grid_values_pattern, line):
                break

    # Set the dict where all output data will be stored
    output_analysis = []

    # Next, we analyze each selected pocket independently. The steps for each pocket are:
    # 1 - Create a new grid file
    # 2 - Conver the grid to pdb
    # 3 - Analyze this pdb with mdpocket
    # 4 - Harvest the volumes over time and write them in the pockets analysis file
    for i, p in enumerate(biggest_pockets):
        pocket_name = "p" + str(i+1)
        print('Analyzing ' + pocket_name + ' (' + str(i+1) + '/' + str(pockets_number) + ')')
        # Create the new grid for this pocket, where all values from other pockets are set to 0
        pocket_value = p[0]
        new_grid_values = [str(v).ljust(5,'0') if pockets[i] == pocket_value else '0.000' for i, v in enumerate(grid_values)]
        new_grid_filename = mdpocket_folder + '/' + pocket_name + '.dx'
        with open(new_grid_filename,'w') as file:
            # Write the header lines
            for line in header_lines:
                file.write(line)
            # Write values in lines of 3 values
            count = 0
            line = ''
            for value in new_grid_values:
                if count == 0:
                    line = value
                    count += 1
                elif count == 1:
                    line = line + ' ' + value
                    count += 1
                elif count == 2:
                    line = line + ' ' + value + ' \n'
                    file.write(line)
                    count = 0

        # Convert the grid coordinates to pdb
        new_pdb_lines = []
        lines_count = 0
        new_pdb_filename = pocket_name + '.pdb'
        for i, v in enumerate(point_coordinates):
            if pockets[i] == pocket_value:
                lines_count += 1
                atom_num = str(lines_count).rjust(6,' ')
                x_coordinates = str(round((origin[0] + v[0]) * 1000) / 1000).rjust(8, ' ')
                y_coordinates = str(round((origin[1] + v[1]) * 1000) / 1000).rjust(8, ' ')
                z_coordinates = str(round((origin[2] + v[2]) * 1000) / 1000).rjust(8, ' ')
                line = "ATOM "+ atom_num +"  C   PTH     1    "+ x_coordinates + y_coordinates + z_coordinates +"  0.00  0.00\n"
                new_pdb_lines.append(line)

        # Write the pdb file
        with open(new_pdb_filename,'w') as file:
            for line in new_pdb_lines:
                file.write(line)

        # Run the mdpocket analysis focusing in this specific pocket
        pocket_output = mdpocket_folder + '/' + pocket_name
        logs = run([
            "mdpocket",
            "--trajectory_file",
            pockets_trajectory,
            "--trajectory_format",
            "xtc",
            "-f",
            input_topology_filename,
            "-o",
            pocket_output,
            "--selected_pocket",
            new_pdb_filename,
        ], stdout=PIPE).stdout.decode()

        # Mine data from the mdpocket 'descriptors' output file
        descriptors_data = {}
        mdpocket_descriptors_filename = pocket_output + '_descriptors.txt'
        with open(mdpocket_descriptors_filename,'r') as file:
            # The '[:-1]' is to remove the break line at the end of each line
            entries = re.split("[ ]+", next(file)[:-1])
            for entry in entries:
                descriptors_data[entry] = []
            for line in file:
                line_data = re.split("[ ]+", line[:-1])
                for i, value in enumerate(line_data):
                    descriptors_data[entries[i]].append(value)

        # Mine the atoms implicated in each pocket each frame
        # In this file atoms are listed for each frame, but they are always the same
        # For this reason, we mine only atoms in the first frame
        atoms = []
        atoms_filename = pocket_output + '_mdpocket_atoms.pdb'
        with open(atoms_filename,'r') as file:
            for line in file:
                line_data = re.split("[ ]+", line[:-1])
                if line_data[0] == 'MODEL':
                    continue
                if line_data[0] == 'ATOM':
                    atoms.append(int(line_data[1]))
                if line_data[0] == 'ENDMDL':
                    break

        # Format the mined data and append it to the output data
        # NEVER FORGET: The 'descriptors_data' object contains a lot of data, not only pocket volumes
        # (e.g. drugability score)
        # DANI: Habría que sentarnos un día a ver que otros valores queremos quedarnos
        output = {
            'name': pocket_name,
            'volumes': list(map(float, descriptors_data['pock_volume'])),
            'atoms': atoms,
        }

        output_analysis.append(output)

    # By default, the starting frame is always 0
    start = 0

    # Export the analysis in json format
    with open(output_analysis_filename, 'w') as file:
        json.dump({ 'data': output_analysis, 'start': start, 'step': step }, file)

    # Finally remove the reduced trajectory since it is not required anymore
    if pockets_trajectory == 'pockets.trajectory.xtc':
        logs = run([
            "rm",
            pockets_trajectory,
        ], stdout=PIPE).stdout.decode()