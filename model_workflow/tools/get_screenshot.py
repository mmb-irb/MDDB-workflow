from subprocess import run, PIPE
import os
from typing import List, Tuple
from PIL import Image
import math
from model_workflow.utils.structures import Structure
import numpy as np

# The Convex Hull is a polygon which covers all the given point and Convex Hull is the smallest polygon. 
# SciPy provides a method "scipy. spatial. ConvexHull()" to create a Convex Hull
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation
from numpy import cross, eye 
from scipy.linalg import expm, norm

auxiliary_tga_filename = '.transition_screenshot.tga'

# Use this to draw some shapes as references in the screenshot
# This is useful for debugging purposes only
debug = False

# Set a function to write a tuple of numbers (point coorinates) as VMD expect these value: a string separated by space
def tuple_to_vmd (point : tuple) -> str:
    return str(point[0]) + ' ' + str(point[1]) + ' ' + str(point[2])

# Python function to obtain a screenshot from the pdb file using VMD a molecular modelling and visualization computer program
def get_screenshot(
    input_structure_filename : str,
    # input_coordinates_filename : List[float] 
    # Remains like this as we are not able to use a set of values for the rotation, translation and zoom and apply them to the molecule
    output_screenshot_filename : str,
    #input_rotate_matrix : List[float],
    #input_scale_matrix : List[float],
    #input_global_matrix : List[float],
):

    # Check the output screenshot file extension is JPG
    if output_screenshot_filename.split('.')[-1] != 'jpg':
        raise SystemExit('You must provide a .jpg file name!')

    # Number of pixels to scale in x 
    x_number_pixels = 350
    # Number of pixels to scale in y
    y_number_pixels = 350

    # Set an enviornment variable to handle the window size
    # WARNING: Use this system instead of the 'display resize' command
    # This command silently kills VMD when it is run in sbatch (but not in salloc)
    # The cause is not clear but this command may depend on OpenGL rendering features which are not supported in sbatch
    # https://www.ks.uiuc.edu/Training/Tutorials/vmd/tutorial-html/node8.html
    os.environ['VMDSCRSIZE'] = str(x_number_pixels) + " " + str(y_number_pixels)

    # Change the scale factor as wished, since some big molecules may overfill in the VMD window
    scale_factor = 2

    # Prepare a Tcl script in order to execute different commands in the terminal of VMD

    # Generate a file name for the commands file
    commands_filename_1 = '.commands_1.vmd'

    # Set a file to save the result obtained when calculating the center
    center_filename = '.center_point_filename.txt'
    with open(commands_filename_1, "w") as file:
        # Select the whole molecule
        file.write('set sel [atomselect 0 all] \n')
        # Find the center
        file.write('set center [measure center $sel] \n')
        # Export the center to a file
        file.write(f'set center_file [open {center_filename} w] \n')
        file.write('puts $center_file $center \n')
        # Exit VMD 
        file.write('exit \n')
    # Run VMD
    logs = run([
        "vmd",
        input_structure_filename,
        "-e",
        commands_filename_1,
        "-dispdev",
        "none"
    ], stdout=PIPE, stderr=PIPE).stdout.decode()
    # Read the generated file to get the center
    with open(".center_point_filename.txt","r") as file:
        line = file.readline().split()
        line = [float(i) for i in line]
    vmd_center_coordinates = line

    # Set a file name for the commands file
    commands_filename_2 = '.commands_2.vmd'

    # Now write the VMD script for the rotation
    with open(commands_filename_2, "w") as file:
        # Set the Background of the molecule to white color
        file.write('color Display Background white \n')
        # Delete the axes drawing in the VMD window
        file.write('axes location Off \n')
        # Eliminate the molecule to perform changes in the representation, color and material 
        file.write('mol delrep 0 top \n')
        # First add a spcific representation for polymers (protein and nucleic acids)
        # Change the default representation model to Newcartoon
        file.write('mol representation Newcartoon \n')
        # Change the default atom coloring method setting to Chain
        file.write('mol color Chain \n')
        # Set the default atom selection setting to all
        file.write('mol selection "protein or nucleic" \n')
        # Change the current material of the representation of the molecule
        file.write('mol material Opaque \n')
        # Using the new changes performed previously add a new representation to the new molecule
        file.write('mol addrep top \n')
        # Change the default representation model to Newcartoon
        file.write('mol representation cpk \n')
        # Change the default atom coloring method setting to Chain
        file.write('mol color element \n')
        # Set the default atom selection setting to all
        file.write('mol selection "not protein and not nucleic" \n')
        # Change the current material of the representation of the molecule
        file.write('mol material Opaque \n')
        # Using the new changes performed previously add a new representation to the new molecule
        file.write('mol addrep top \n') 
        # Change projection from perspective (used by VMD by default) to orthographic
        file.write('display projection orthographic \n')
        # Select all atoms
        file.write('set sel [atomselect 0 all] \n')

        # Load the whole structure
        # WARNING: Note that custom atoms selection is not yet supported, although we could do it one day
        # WARNING: Not selecting all atoms but a subselection may cause problems when centering the image right now
        structure = Structure.from_pdb_file(input_structure_filename)
        # Obtain all coordinates from each atom of the filtered structure
        coordinates = [list(atom.coords) for atom in structure.atoms]        
        # Convert the list into a Numpy Array, since Scipy library just works with this type of data structure
        coordinates_np = np.array(coordinates)

        # Iterate over all coordinates to store in a dictionary as a key a tuple with the coordinates with the index as its value
        d = {}
        for index,coordinate in enumerate(coordinates):
            d[tuple(coordinate)] = index

        # Compute the hull of the molecule
        hull = ConvexHull(coordinates_np)
        # Extract the points forming the hull
        hullpoints = coordinates_np[hull.vertices,:]
        # Naive way of finding the best pair in O(H^2) time if H is number of points on hull
        hull_best_points = cdist(hullpoints, hullpoints, metric='euclidean')
        # Get the farthest apart points
        bestpair = np.unravel_index(hull_best_points.argmax(), hull_best_points.shape)
        # Convert the first point coordinates from Numpy array to list
        first_point = hullpoints[bestpair[0]].tolist()
        # Now into tuple 
        first_point = tuple(first_point)
        # Convert the second point coordinates from Numpy array to list
        second_point = hullpoints[bestpair[1]].tolist()
        # Now into tuple 
        second_point = tuple(second_point)
        
        ### TRIGONOMETRY TO COMPUTE THE ANGLE WE NEED TO ROTATE

        # FIRST ROTATION
        # Looking at the x-z plane, consider we have a rectangle triangle where the segment between the
        # first and second points is the hypotenuse and the sides of traingle are paralel to x and z axes
        # Calculate the grades we must rotate the molecule to align the two points in the x axis

        # Get the hypotenuse
        # Note that the hypotenuse is calculated as a projection of the segment in the x-z plane
        # This is lower than the actual distance between the points
        hypotenuse = calculate_distance(first_point, second_point, ['x', 'z'])
    
        # Get the z-side which it is just the difference in the z coordinates
        z_side = abs(first_point[2] - second_point[2])

        # Obtain the angle that we are interested in order to rotate the molecule
        # Apply trigonometry, compute the arcsinus of the division between the opposite side and the hypothenuse 
        angle = math.degrees(math.asin(z_side/hypotenuse)) if hypotenuse > 0 else 0

        # Set the right rotation direction
        most_negative_x_point = first_point if first_point[0] <= second_point[0] else second_point
        most_negative_z_point = first_point if first_point[2] <= second_point[2] else second_point
        if most_negative_x_point != most_negative_z_point:
            angle = -angle

        # First rotation of the molecule to set it perpendicular with respect to z axis
        file.write(f'rotate y by {angle} \n')

        # SECOND ROTATION
        # Now looking at the x-y plane, consider we have a rectangle triangle where the segment between the
        # first and second points is the hypotenuse and the sides of traingle are paralel to x and y axes
        # Calculate the grades we must rotate the molecule to align the two points in the x axis

        # Note that once the rotation for this triangle is solved we add an extra 45 grades
        # This is becase we want the molecule to be aligned with the diagonal of the image, not the x axis

        # Get the hypotenuse 
        # Note that here we use all dimensions and not the projection in x-y plane
        # This is because VMD rotates using its current rotation as reference, not the absolute
        # Thus the segment is already in "our" x-y plane after the previous rotation
        # Actually removing the z would be wrong since this z is not "our" z
        # It is hard to understand only with words but trust
        hypotenuse2 = calculate_distance(first_point, second_point, ['x', 'y', 'z'])

        # Get the y-side which it is just the difference in the y coordinates
        y_side = abs(first_point[1] - second_point[1])

        # Obtain the angle that we are interested in order to rotate the molecule
        # Apply trigonometry, compute the arcsinus of the division between the opposite side and the hypothenuse
        angle2 = math.degrees(math.asin(y_side/hypotenuse2)) if hypotenuse2 > 0 else 0

        # Set the right rotation direction
        most_negative_x_point = first_point if first_point[0] <= second_point[0] else second_point
        most_positive_y_point = first_point if first_point[1] >= second_point[1] else second_point
        if most_negative_x_point != most_positive_y_point:
            angle2 = -angle2
        # As we want it diagonal with respect y and x we add 45 degrees
        angle2 += 45

        # Second rotatoin of the molecule to set it diagonal with respect to z axis
        file.write(f'rotate z by {angle2} \n')

        # Vector that we are going to use as x axis 
        absolute_x_axis = (1,0,0)
        # Vector that we are going to use as y axis 
        absolute_y_axis = (0,1,0)

        # Initial normal vector representing the direction of the camera view
        # Note that when we enter VMD the z axis is pointing towards us
        initial_normal_vector = (0,0,-1) 

        # Rotate normal vector in order to have it in the same direction as the camera after the first rotation
        # Note that the angle here and below is negative
        # Rotating in VMD by the y axis is counter intuitive: negative means clockwise
        rotated_normal_vector = rotate_vector(initial_normal_vector, -angle, absolute_y_axis)

        # Rotate x axis around y axis
        first_rotated_x_axis = rotate_vector(absolute_x_axis, -angle, absolute_y_axis)
        # Note that there is no need to apply the first rotation to the y axis since we rotate around the y axis

        # Now, for the second rotation we must use the rotated normal vector as pivot
        # This emulates how VMD rotates around the recently moved camera axis, not the global axis
        # Obtain rotation matrix to perform the second rotation we have done previously 
        rotated_x_axis = rotate_vector(first_rotated_x_axis, angle2, rotated_normal_vector)
        # Perform just the same procedure as before
        rotated_y_axis = rotate_vector(absolute_y_axis, angle2, rotated_normal_vector)

        ##########################################################

        # Project all atom coordinates in each rotated axis to find the range of coordinates and thus its center
        # Formula from https://gamedev.stackexchange.com/questions/72528/how-can-i-project-a-3d-point-onto-a-3d-line
        # A + dot(AP,AB) / dot(AB,AB) * AB

        # The point to set the line does not matter. It could be (0,0,0)
        # However, we use the vmd center so then we can take the projected points as references for debugging
        A_point = vmd_center_coordinates
        
        # As we are going to compute a line in Y axis we just set AB_vector equal to the rotated y axis 
        AB_vector = rotated_y_axis

        # Function to project a point into a line
        # WARNING: It has been observed experimentally that this function is not working as expected
        # WARNING: Points are projected in the line but very shifted
        # WARNING: However the error is the same for the molecule points and the view center
        # WARNING: For this reason the resulting difference is correct and it works
        def get_projected_point (P_point : Tuple[float, float, float]) -> Tuple[float, float, float]:
            AP_vector = [A_point[i] + P_point[i] for i in range(3)]
            AP_AB_dot_product = np.dot(AP_vector, AB_vector)
            AB_AB_dot_product = np.dot(AB_vector, AB_vector)
            scalar = AP_AB_dot_product / AB_AB_dot_product
            projected_vector = [AB_vector[i] * scalar for i in range(3)]
            projected_point = [A_point[i] + projected_vector[i] for i in range(3)]
            return projected_point

        # Project all atom coordinates in the rotated y axis
        projected_points = [ get_projected_point(P_point) for P_point in coordinates ]
        projected_points.sort(key=lambda projected_points: projected_points[0])
        projected_points.sort(key=lambda projected_points: projected_points[1])
        projected_points.sort(key=lambda projected_points: projected_points[2])
        # Get the most distant projected points
        max_ypoint = projected_points[0]
        min_ypoint = projected_points[-1]
        # Find the center between the most distant points
        y_axis_center = [(max_ypoint[i] + min_ypoint[i]) / 2 for i in range(3)]

        # Project the center of the view in the rotated y axis
        y_axis_projected_center = get_projected_point(vmd_center_coordinates)
        # Calculate the difference between the molecule center and the view center in the rotated y axis
        y_axis_difference_vector = [y_axis_projected_center[i] - y_axis_center[i] for i in range(3)]
        # Move to rectify the difference
        file.write('$sel moveby { ' + tuple_to_vmd(y_axis_difference_vector) + ' } \n')
        
        # Repeat the whole process with the rotated x axis

        # Project all atom coordinates in the rotated x axis
        AB_vector = rotated_x_axis
        projected_points2 = [ get_projected_point(P_point) for P_point in coordinates ]
        projected_points2.sort(key=lambda projected_points2: projected_points2[0])
        projected_points2.sort(key=lambda projected_points2: projected_points2[1])
        projected_points2.sort(key=lambda projected_points2: projected_points2[2])
        # Get the most distant projected points
        max_xpoint = projected_points2[0]
        min_xpoint = projected_points2[-1]
        # Find the center between the most distant points
        x_axis_center = [(max_xpoint[i] + min_xpoint[i]) / 2 for i in range(3)]

        # Project the center of the view in the rotated x axis
        x_axis_projected_center = get_projected_point(vmd_center_coordinates)
        # Calculate the difference between the molecule center and the view center in the rotated x axis
        x_axis_difference_vector = [x_axis_projected_center[i] - x_axis_center[i] for i in range(3)]
        
        # Move to rectify the difference
        file.write('$sel moveby { ' + tuple_to_vmd(x_axis_difference_vector) + ' } \n')
        
        # Calculate height and width and get the widest dimension
        width = calculate_distance(max_xpoint, min_xpoint, ['x','y','z'])
        height = calculate_distance(max_ypoint, min_ypoint, ['x','y','z'])
        widest = max(width, height)
        # Set how close to the molecule we want the camera to be
        # This value has been found experimentally and it keeps a small white margin
        zoom = 2.8
        # Set the scale
        scale = zoom / widest
        file.write(f'scale to {scale} \n')

        # Show the theoretical view center
        if debug:
            # Note that all draw commands work with coordinates and thus they are absolute, not relative to camera
            # Draw the center
            file.write('draw sphere {' + tuple_to_vmd(vmd_center_coordinates) + '} radius 5\n')
            # Set a function to display vectors
            def display_vector (vector : tuple, color : str):
                ref_right_vector = multiply_vector(vector, 100)
                ref_right_point = point_adds_vector(vmd_center_coordinates, ref_right_vector)
                ref_left_vector = multiply_vector(vector, -100)
                ref_left_point = point_adds_vector(vmd_center_coordinates, ref_left_vector)
                file.write('draw color ' + color + '\n')
                file.write('draw line {' + tuple_to_vmd(ref_right_point) + '} {' + tuple_to_vmd(ref_left_point) + '}\n')
            # Draw the rotated axes
            #display_vector(first_rotated_x_axis, 'orange')
            display_vector(rotated_x_axis, 'red')
            display_vector(rotated_y_axis, 'green')
            # Draw projected point ranges
            # file.write('draw color red\n')
            # file.write('draw line {' + tuple_to_vmd(min_xpoint) + '} {' + tuple_to_vmd(max_xpoint) + '}\n')
            # file.write('draw color green\n')
            # file.write('draw line {' + tuple_to_vmd(min_ypoint) + '} {' + tuple_to_vmd(max_ypoint) + '}\n')

        # Finally generate the image from the current view
        file.write('render TachyonInternal ' + auxiliary_tga_filename + ' \n')
        # Exit VMD
        file.write('exit\n')

    # Run VMD
    process = run([
        "vmd",
        input_structure_filename,
        "-e",
        commands_filename_2,
        "-dispdev",
        "none"
    ], stdout=PIPE, stderr=PIPE)
    logs = process.stdout.decode()
    # If the output file does not exist at this point then it means something went wrong with VMD
    if not os.path.exists(auxiliary_tga_filename):
        print(logs)
        error_logs = process.stderr.decode()
        print(error_logs)
        raise SystemExit('Something went wrong with VMD')

    im = Image.open(auxiliary_tga_filename)
    # converting to jpg
    rgb_im = im.convert("RGB")
    # exporting the image
    rgb_im.save(output_screenshot_filename)

    # Remove trash files
    trash_files = [ commands_filename_1, commands_filename_2, auxiliary_tga_filename, center_filename ]
    for trash_file in trash_files:
        os.remove(trash_file)

    # Remove the environment variable to not cause problem in possible future uses of VMD
    os.environ['VMDSCRSIZE'] = ""

# Return the rotated vector
# https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
def rotate_vector (vector : tuple, angle : float, axis : tuple) -> tuple:
    radians = math.radians(angle)
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(radians / 2.0)
    b, c, d = -axis * math.sin(radians / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    rotation_matrix = np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
        [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
        [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])
    return np.dot(rotation_matrix, vector)


# Calculate the distance between two atoms
def calculate_distance (first_coordinates : tuple, second_coordinates : tuple, dimensions : List[str]) -> float:
    squared_distances_sum = 0
    for dimension, i in { 'x': 0, 'y': 1, 'z': 2 }.items():
        if dimension not in dimensions:
            continue
        squared_distances_sum += (first_coordinates[i] - second_coordinates[i])**2
    return math.sqrt(squared_distances_sum)

# Set a function to gen a new point from a point and vector
def point_adds_vector (point : tuple, vector : tuple) -> tuple:
    return tuple([ point[d] + vector[d] for d in range(3) ])

def multiply_vector (vector : tuple, multiplier : float) -> tuple:
    return tuple([vector[d] * multiplier for d in range(3)])