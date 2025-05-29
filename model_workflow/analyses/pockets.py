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

from os.path import exists, getsize, split
from os import mkdir, remove, chdir, getcwd
import re
import collections

from subprocess import run, PIPE

from model_workflow.tools.get_reduced_trajectory import get_reduced_trajectory
from model_workflow.utils.auxiliar import warn, ToolError, save_json
from model_workflow.utils.file import File
from model_workflow.utils.constants import GREY_HEADER, COLOR_END

CURSOR_UP_ONE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'
ERASE_PREVIOUS_LINE = CURSOR_UP_ONE + ERASE_LINE
ERASE_4_PREVIOUS_LINES = ERASE_PREVIOUS_LINE + ERASE_PREVIOUS_LINE + ERASE_PREVIOUS_LINE + ERASE_PREVIOUS_LINE + CURSOR_UP_ONE

# Set some known error logs to check if we are having them
KNOWN_MDPOCKET_ERRORS = set([
    'Error in creating clustering tree, return NULL pointer...breaking up'
])

# Perform the pockets analysis
def pockets (
    structure_file : 'File',
    trajectory_file : 'File',
    output_analysis_filepath : str,
    pockets_prefix : str,
    mdpocket_folder : str,
    pbc_selection : 'Selection',
    snapshots : int,
    frames_limit : int = 100,
    # Get only the 10 first pockets since the analysis is quite slow by now
    # DANI: Cuando hagamos threading y no haya limite de tamaño para cargar en mongo podremos hacer más pockets
    maximum_pockets_number : int = 10):

    print('-> Running pockets analysis')

    # DANI: De momento, no se hacen pockets para simulaciones con residuos en PBC (e.g. membrana)
    # DANI: Esto es debido a que los átomos donde NO queremos que encuentre pockets no se pueden descartar
    # DANI: Descartarlos significa quitarlos, pero si los quitamos entonces podemos encontrar pockets donde están estos átomos
    # DANI: Estamos a la espera que los de mdpocket incluyan un flag para estos casos
    # https://github.com/Discngine/fpocket/issues/77#issuecomment-974193129
    if pbc_selection:
        print(' Pockets analysis will be skipped since we have PBC atoms')
        return

    # Set a reduced trajectory with only 100 frames
    # Get the step between frames of the new reduced trajectory, since it will be append to the output
    pockets_trajectory, step, frames = get_reduced_trajectory(
        structure_file,
        trajectory_file,
        snapshots,
        frames_limit,
    )
    # Save the pockets trajectory as a file
    pockets_trajectory_file = File(pockets_trajectory)    

    # This anlaysis produces many useless output files
    # Create a new folder to store all ouput files so they do not overcrowd the main directory
    if not exists(mdpocket_folder):
        mkdir(mdpocket_folder)

    # WARNING: There is a silent sharp limit of characters here
    # https://github.com/Discngine/fpocket/issues/130
    # To avoid the problem we must change the directory where we run pockets so all the paths are as short as possible

    # Save the MD path
    md_path, mdpocket_folder_name = split(mdpocket_folder)
    # Count the number of directory deep we are now
    # Thus we can know how many directories back we need to jump to get back to the original directory
    if md_path[0] == '/': raise ValueError('This path should not be absolute, the fix below will not work')
    if md_path[-1] == '/': md_path = md_path[0:-1]
    recovery_path = getcwd()

    # Move to the MD path so all relative paths become shorter
    chdir(md_path)

    # Now set the structure and trajectory paths relative to the current new location
    # HARDCODE: We asume the trajectory and structure files are in the current location
    # This is now true, but if this changes in the future this will fail, although the error will be obvious
    auxiliar_trajectory_filepath = f'./{pockets_trajectory_file.filename}'
    auxiliar_structure_filepath = f'./{structure_file.filename}'

    # Run the mdpocket analysis to find new pockets
    mdpocket_output = mdpocket_folder_name + '/mdpout'
    # Set the filename of the fpocket output we are intereseted in
    grid_filename = mdpocket_output + '_freq.dx'
    # Skip this step if the output file already exists and is not empty
    # WARNING: The file is created as soon as mdpocket starts to run
    # WARNING: However the file remains empty until the end of mdpocket
    if not exists(grid_filename) or getsize(grid_filename) == 0:

        print('Searching new pockets')
        print(GREY_HEADER)
        process = run([
            "mdpocket",
            "--trajectory_file",
            # WARNING: There is a silent sharp limit of characters here
            # To avoid the problem we must use the relative path instead of the absolute path
            auxiliar_trajectory_filepath,
            "--trajectory_format",
            "xtc",
            "-f",
            # WARNING: There is a silent sharp limit of characters here
            # To avoid the problem we must use the relative path instead of the absolute path
            auxiliar_structure_filepath,
            "-o",
            mdpocket_output,
        ], stderr=PIPE)
        error_logs = process.stderr.decode()
        print(COLOR_END)

        # If file does not exist or is still empty at this point then somethin went wrong
        if not exists(grid_filename) or getsize(grid_filename) == 0:
            print(error_logs)
            raise ToolError('Something went wrong with mdpocket while searching pockets')

        # Check if we are having concerning error logs
        errored = False
        for error in KNOWN_MDPOCKET_ERRORS:
            matches = re.findall(error, error_logs)
            count = len(matches)
            if count > 0:
                warn(f'We are having "{error}" ({count})')
                errored = True

        # If so we can stop here since we will not be able to mine the grid filename
        if errored:
            if exists(grid_filename):
                remove(grid_filename)
            raise ToolError('We had errors with mdpocket while searching pockets')

    # Read and harvest the gird file
    with open(grid_filename, 'r') as file:

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
        if not dimensions:
            # This may happend when one of the dimensions is negative, which an mdpocket error
            # DANI: El origen de esto no está claro pero cuando me pasó vino acompañado de muchos errores:
            # DANI: '! No Pockets Found while refining' y '! No pocket to reindex.'
            # DANI: Además de algún 'Error in creating clustering tree, return NULL pointer...breaking up!'
            raise Exception('Failed to mine dimensions')

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

    # Set a function to get the value of a given 'x, y, z' grid point
    # Grid points are disposed first in 'z' order, then 'y', and finally 'x'
    xl, yl, zl = dimensions
    def get_index (x : int, y : int, z : int) -> int:
        index = x*yl*zl + y*zl + z
        return index

    # Get the coordinates of the colliding points
    # Set the x, y and z index limit (i.e. the limit - 1)
    xil, yil, zil = (xl-1, yl-1, zl-1)
    def get_colliding_points (point : tuple) -> list:
        x, y, z = point[0], point[1], point[2]
        colliding_points = [
            (x-1, y-1, z-1),
            (x-1, y-1, z  ),
            (x-1, y-1, z+1),
            (x-1, y  , z-1),
            (x-1, y  , z  ),
            (x-1, y  , z+1),
            (x-1, y+1, z-1),
            (x-1, y+1, z  ),
            (x-1, y+1, z+1),
            (x  , y-1, z-1),
            (x  , y-1, z  ),
            (x  , y-1, z+1),
            (x  , y  , z-1),
            (x  , y  , z+1),
            (x  , y+1, z-1),
            (x  , y+1, z  ),
            (x  , y+1, z+1),
            (x+1, y-1, z-1),
            (x+1, y-1, z  ),
            (x+1, y-1, z+1),
            (x+1, y  , z-1),
            (x+1, y  , z  ),
            (x+1, y  , z+1),
            (x+1, y+1, z-1),
            (x+1, y+1, z  ),
            (x+1, y+1, z+1),
        ]
        if x <= 0:
            colliding_points = [ point for point in colliding_points if point[0] != x-1 ]
        if x >= xil:
            colliding_points = [ point for point in colliding_points if point[0] != x+1 ]
        if y <= 0:
            colliding_points = [ point for point in colliding_points if point[1] != y-1 ]
        if y >= yil:
            colliding_points = [ point for point in colliding_points if point[1] != y+1 ]
        if z <= 0:
            colliding_points = [ point for point in colliding_points if point[2] != z-1 ]
        if z >= zil:
            colliding_points = [ point for point in colliding_points if point[2] != z+1 ]
        return colliding_points

    # Set a cuttoff value to consider a point valid
    cuttoff = 0.5

    # Classify all non-zero values by 'pocket' groups according to if they are connected
    pockets = [0] * xl * yl * zl
    pockets_count = 0

    # Given a point which belongs to a pocket, find all pocket points connected to it
    # This connection must be wide enought for a water molecule to fit in
    # Grid points are 1 Ångstrom wide and a water molecule has a diameter of 2.75 Ångstroms aprox
    # So, pocket points must be connected by a net of points at least 3 points wide in all dimensions (x,y,z)
    # Tag all points with the same pocket number
    def set_pocket (start_point : tuple, pocket_number : int):
        points = [start_point]
        for point in points:
            # Get current point colliders
            colliding_points = get_colliding_points(point)
            # Filter only the pocket points
            pocket_points = []
            for point in colliding_points:
                index = get_index(*point)
                value = grid_values[index]
                # If it is a pocket then set it as part of the current pocket
                if value >= cuttoff:
                    # Tag it as the current pocket
                    pockets[index] = pocket_number
                    pocket_points.append(point)
            # In case all its surrounding points are pockets add the points to the list so they can keep expanding the pocket
            # This is done because if all surrounding points are pockets it means we have a wide enought region of the pocket
            # It is a 3 x 3 x 3 region.
            # Otherwise stop here since the surrounding points may not we wide enough to expand the pocket
            if len(pocket_points) < 26:
                continue
            # Filter points which are not in the pockets list already
            new_points = [ point for point in pocket_points if point not in points ]
            # Update the points list
            points += new_points

    # Find out if a point is surrounded by all pocket points so it could be considered a pocket alone
    def is_pocket_base (start_point : tuple) -> bool:
        # Get current point colliders
        colliding_points = get_colliding_points(start_point)
        # Search if these points match the cutoff
        # If only 1 point does not match then this is not a pocket base
        for point in colliding_points:
            index = get_index(*point)
            value = grid_values[index]
            if value < cuttoff:
                return False
        return True
            
    # Save also each point coordinates
    for x in range(xl):
        for y in range(yl):
            for z in range(zl):
                index = get_index(x,y,z)
                # If it is not a pocket then pass
                value = grid_values[index]
                if value < cuttoff:
                    continue
                # If it is a pocket but it has been already tagged then pass
                pocket = pockets[index]
                if pocket:
                    continue
                # If it is not wide enought to start a pocket then pass
                point = (x,y,z)
                if not is_pocket_base(point):
                    continue
                # If none of the previous values was marked as a pocket then we set a new pocket number
                pockets_count += 1
                set_pocket(start_point=point, pocket_number=pockets_count)

    # Exclude the first result which will always be 0 and it stands for no-pocket points
    biggest_pockets = collections.Counter(pockets).most_common()
    if len(biggest_pockets) == 1:
        print(' No pockets were found')
        # Recover the original directory
        chdir(recovery_path)
        return
    biggest_pockets = biggest_pockets[1:]
    pockets_number = len(biggest_pockets)
    # If we have more than the maximum number of pockets then get the first pockets and discard the rest
    if pockets_number > maximum_pockets_number:
        biggest_pockets = biggest_pockets[0:maximum_pockets_number]
        pockets_number = maximum_pockets_number

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
    for p, pock in enumerate(biggest_pockets, 1):
        # WARNING: This name must match the final name of the pocket file once loaded in the database
        pocket_name = 'pocket_' + str(p).zfill(2)
        pocket_output = mdpocket_folder_name + '/' + pocket_name
        # Check if current pocket files already exist and are complete. If so, skip this pocket
        # Output files:
        # - pX.dx: it is created and completed at the begining by this workflow
        # - pX_descriptors.txt: it is created at the begining and completed along the mdpocket progress
        # - pX_mdpocket_atoms.pdb: it is completed at the begining but remains size 0 until the end of mdpocket
        # - pX_mdpocket.pdb: it is completed at the begining but remains size 0 until the end of mdpocket
        # Note that checking pX_mdpocket_atoms.pdb or pX_mdpocket.pdb is enought to know if mdpocket was completed
        checking_filename = pocket_output + '_mdpocket.pdb'
        if not (exists(checking_filename) and getsize(checking_filename) > 0):

            # Update the logs
            print(f' Analyzing pocket {p}/{pockets_number}', end='\r')
            
            # Create the new grid for this pocket, where all values from other pockets are set to 0
            pocket_value = pock[0]
            new_grid_values = [str(value).ljust(5,'0') if pockets[v] == pocket_value else '0.000' for v, value in enumerate(grid_values)]
            new_grid_filename = mdpocket_folder_name + '/' + pocket_name + '.dx'
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

            # Set a function to get the 'x, y, z' coordinates of a given grid point index
            def get_coordinates (index : int) -> tuple:
                x, rem = divmod(index, yl*zl)
                y, z = divmod(rem, zl)
                return (x,y,z)

            # Convert the grid coordinates to pdb
            new_pdb_lines = []
            lines_count = 0
            # HARDCODE: Since we are cding to the current file we must remove the MD path from the prefix
            fixed_pockets_prefix = pockets_prefix.split('/')[-1]
            new_pdb_filename = fixed_pockets_prefix + '_' + str(p).zfill(2) + '.pdb'
            for j, pocket in enumerate(pockets):
                if pocket != pocket_value:
                    continue
                x,y,z = get_coordinates(j)
                lines_count += 1
                atom_num = str(lines_count).rjust(6,' ')
                x_coordinates = str(round((origin[0] + x) * 1000) / 1000).rjust(8, ' ')
                y_coordinates = str(round((origin[1] + y) * 1000) / 1000).rjust(8, ' ')
                z_coordinates = str(round((origin[2] + z) * 1000) / 1000).rjust(8, ' ')
                line = "ATOM "+ atom_num +"  C   PTH     1    "+ x_coordinates + y_coordinates + z_coordinates +"  0.00  0.00\n"
                new_pdb_lines.append(line)

            # Write the pdb file
            with open(new_pdb_filename,'w') as file:
                for line in new_pdb_lines: 
                    file.write(line)

            # Run the mdpocket analysis focusing in this specific pocket
            print(GREY_HEADER)
            error_logs = run([
                "mdpocket",
                "--trajectory_file",
                auxiliar_trajectory_filepath,
                "--trajectory_format",
                "xtc",
                "-f",
                # WARNING: There is a silent sharp limit of characters here
                # To avoid the problem we must use the relative path instead of the absolute path
                auxiliar_structure_filepath,
                "-o",
                pocket_output,
                "--selected_pocket",
                new_pdb_filename,
            ], stderr=PIPE).stderr.decode()
            print(COLOR_END)

            # If file does not exist or is still empty at this point then somethin went wrong
            if not exists(checking_filename) or getsize(checking_filename) == 0:
                print(error_logs)
                raise Exception(f'Something went wrong with mdpocket while analysing pocket {p}')

            # Remove previous lines
            print(ERASE_4_PREVIOUS_LINES)

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
                for v, value in enumerate(line_data):
                    descriptors_data[entries[v]].append(value)

        # Mine the Voronoi vertices which represent the pocket each frame
        # Mine their positions and radius for each frame
        # Example of lines in hte file to be mined:
        # ATOM      0    O STP   105      79.352  75.486  84.079    0.00     4.25
        # ATOM      1    C STP   105      79.735  75.228  84.482    0.00     4.13
        # ATOM      2    C STP   105      79.624  74.251  84.263    0.00     4.52
        # vertices = []
        # vertices_filename = pocket_output + '_mdpocket.pdb'
        # with open(vertices_filename,'r') as file:
        #     frame_vertices = []
        #     for line in file:
        #         line_data = re.split("[ ]+", line[:-1])
        #         if line_data[0] == 'ATOM':
        #             frame_vertices.append({
        #                 'e': line_data[2],
        #                 'n': line_data[4],
        #                 'x': line_data[5],
        #                 'y': line_data[6],
        #                 'z': line_data[7],
        #                 'r': line_data[9],
        #             })
        #         if line_data[0] == 'ENDMDL':
        #             vertices.append(frame_vertices)
        #             frame_vertices = []

        # Mine the atoms implicated in each pocket
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
            #'vertices': vertices,
            'atoms': atoms,
        }

        output_analysis.append(output)

    # Recover the original directory
    chdir(recovery_path)

    # By default, the starting frame is always 0
    start = 0

    # Export the analysis in json format
    save_json({ 'data': output_analysis, 'start': start, 'step': step }, output_analysis_filepath)