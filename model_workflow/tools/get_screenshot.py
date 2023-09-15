from subprocess import run, PIPE
import os
from typing import List, Tuple
from PIL import Image
import math
from mdtoolbelt.structures import Structure
import numpy as np

# The Convex Hull is a polygon which covers all the given point and Convex Hull is the smallest polygon. 
# SciPy provides a method "scipy. spatial. ConvexHull()" to create a Convex Hull
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation
from numpy import cross, eye 
from scipy.linalg import expm, norm

auxiliary_tga_filename = '.transition_screenshot.tga'
# Python function to obtain a screenshot from the pdb file using VMD a molecular modelling and visualization computer program
def get_screenshot(
    input_structure_filename : str,
    # input_coordinates_filename : List[float] 
    # Remains like this as we are not able to use a set of values for the rotation, translation and zoom and apply them to the molecule
    output_screenshot_filename : str,
    #input_rotate_matrix : List[float],
    #input_scale_matrix : List[float],
    #input_global_matrix : List[float],
    #focus_selection : str = 'protein or nucleic',
    focus_selection : str = 'all',
):
    # Number of pixels to scale in x 
    x_number_pixels = 350
    # Number of pixels to scale in y
    y_number_pixels = 350
    # Change the scale factor as wished, since some big molecules may overfill in the VMD window
    scale_factor = 2
    if output_screenshot_filename.split('.')[-1] != 'jpg':
        raise SystemExit('You must provide a .jpg file name!')


    # Prepare a Tcl script in order to execute different commands in the terminal of VMD

    # Generate a file name for the commands file
    commands_filename = 'commands.vmd' 

    # Generate a file to save the result obtained when calculating the center
    center_filename = '.center_point_filename.txt'
    with open(commands_filename, "w") as file:
        
        # Select the whole molecule and save it as variable sel
        file.write('set sel [atomselect 0 all] \n')
        
        # Save in variable called center the value obtained from the command inside the square brackets that computes the center from the variable sel
        file.write('set center [measure center $sel] \n')

        # Open the file using name center_file as variable 
        file.write(f'set center_file [open {center_filename} w] \n')

        # Copy the result into the file
        file.write('puts $center_file $center \n')

        # Exit VMD 
        file.write('exit \n')
        # Run VMD
    logs = run([
        "vmd",
        input_structure_filename,
        "-e",
        commands_filename,
        "-dispdev",
        "none"
    ], stdout=PIPE, stderr=PIPE).stdout.decode()

    with open(".center_point_filename.txt","r") as file:
        line = file.readline()
        line = line.split()
        line = [float(i) for i in line]
    center_point_protein = line

    with open(commands_filename, "w") as file:
        
        # Set the Background of the molecule to white color
        file.write('color Display Background white \n') 
        
        # Delete the axes drawing in the VMD window
        file.write('axes location Off \n') 
        
        # Eliminate the molecule to perform changes in the representation, color and material 
        file.write('mol delrep 0 top \n') 
        
        # Change the default representation model to Newcartoon
        #file.write('mol representation Newcartoon \n') 
        file.write('mol representation CPK \n') 
        
        # Change the default atom coloring method setting to Chain
        #file.write('mol color Chain \n') 
        


        # Set the default atom selection setting to all
        file.write('mol selection {all} \n') 
        
        # Change the current material of the representation of the molecule
        file.write('mol material Opaque \n') 
        
        # Using the new changes performed previously add a new representation to the new molecule
        file.write('mol addrep top \n') 
        #if input_coordinates_filename != []: 
        #    pass

        # Adjust the VMD window to the desired size
        file.write(f'display resize {x_number_pixels} {y_number_pixels} \n') 

        # Change projection from perspective (used by VMD by default) to orthographic

        file.write('display projection orthographic \n')

        file.write('set sel [atomselect 0 all] \n')
        '''
        if input_rotate_matrix != []: # Check if we have a rotation matrix to use it and change the default rotate matrix values in VMD
            file.write(f'molinfo top set rotate_matrix {input_rotate_matrix}')
        if input_scale_matrix != []: # Check if we have a scale matrix to use it and change the default scale matrix values in VMD
            file.write(f'molinfo top set scale_matrix {input_scale_matrix}')
        if input_global_matrix != []: # Check if we have a global matrix to use it and change the default global matrix values in VMD
            file.write(f'molinfo top set global_matrix {input_global_matrix}')
        '''
        structure = Structure.from_pdb_file(input_structure_filename)

        # Parse the selection in VMD selection syntax
        parsed_selection = structure.select(focus_selection, syntax='vmd')

        # If there is nothing to check then warn the user and stop here
        if not parsed_selection:
            raise SystemExit('There are not atoms to be focused for the screenshot')
        
        # Filter focused atoms
        structure = structure.filter(parsed_selection)

        # Obtain all coordinates from each atom of the filtered structure
        coordinates = [list(atom.coords) for atom in structure.atoms]
        d = {}
        
        # Convert the list into a Numpy Array, since Scipy library just works with this type of data structure
        coordinates_np = np.array(coordinates)

        # Iterate over all coordinates to store in a dictionary as a key a tuple with the coordinates with the index as its value
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
        
        #Print them
        #print([hullpoints[bestpair[0]],hullpoints[bestpair[1]]])

        # Convert the first point coordinates from Numpy array to list
        first_point = hullpoints[bestpair[0]].tolist()
        
        # Now into tuple 
        first_point = tuple(first_point)

        # Convert the second point coordinates from Numpy array to list
        second_point = hullpoints[bestpair[1]].tolist()

        # Now into tuple 
        second_point = tuple(second_point)



        #print('index ' + str(d[first_point]) + ' ' + str(d[second_point]))
        
        ### TRIGONOMETRY TO COMPUTE THE ANGLE WE NEED TO ROTATE

        # Consider we have a rectangle triangle if we look at the x-z plane
        # Get the hypotenuse

        
        hypotenuse = calculate_distance(first_point, second_point, ['x', 'z'])
    
        # Get the z-side which it is just the difference in the z coordinates
        z_side = abs(first_point[2] - second_point[2])

        # Obtain the angle that we are interested in order to rotate the protein
        # Apply trigonometry, compute the arcsinus of the division between the opposite side and the hypothenuse 
        angle = math.degrees(math.asin(z_side/hypotenuse)) if hypotenuse > 0 else 0

        # Set the right rotation direction
        most_negative_x_point = first_point if first_point[0] <= second_point[0] else second_point
        most_negative_z_point = first_point if first_point[2] <= second_point[2] else second_point
        if most_negative_x_point != most_negative_z_point:
            angle = -angle


        # First rotation of the molecule to set it perpendicular with respect to z axis
        file.write(f'rotate y by {angle} \n')

        # Consider another triangle but this time if we look at the x-y plane
        # Get the hypotenuse 
        hypotenuse2 = calculate_distance(first_point, second_point, ['x', 'y', 'z'])
        #print(hypotenuse2)

        # Get the y-side which it is just the difference in the y coordinates
        y_side = abs(first_point[1] - second_point[1])

        # Obtain the angle that we are interested in order to rotate the protein
        # Apply trigonometry, compute the arcsinus of the division between the opposite side and the hypothenuse
        angle2 = math.degrees(math.asin(y_side/hypotenuse2)) if hypotenuse2 > 0 else 0

        # Set the right rotation direction
        most_negative_x_point = first_point if first_point[0] <= second_point[0] else second_point
        most_positive_y_point = first_point if first_point[1] >= second_point[1] else second_point
        if most_negative_x_point != most_positive_y_point:
            angle2 = -angle2
        
        # As we want it diagonal with respect y and x we add 45 degrees as the angle computed will rotate the protein to put it perpendicular to y
        angle2 += 45



        # Second rotatoin of the molecule to set it diagonal with respect to z axis
        file.write(f'rotate z by {angle2} \n')

        # Random initial normal vector
        initial_normal_vector = (0,0,-1) 
        
        # Rotate normal vector in order to have it in the same direction as we have rotated the camera with respect y axis
        rotated_normal_vector = rotate(initial_normal_vector, angle, 'y')
        #rotated_normal_vector = rotate(rotated_normal_vector, angle2, 'z')


        # Vector that we are going to use as x axis 
        vector_x_axis = (1,0,0)

        # Vector that we are going to use as y axis 
        vector_y_axis = (0,1,0)

        # Rotate vector with respect y
        rotated_x_axis = rotate(vector_x_axis, angle, 'y')

        # Note that there is no need to apply the first rotation to the y axis

        # Now, for the second rotation we must use the rotated normal vector as pivot
        # This emulates how VMD rotates around the recently moved camera axis, not the global axis

        # Obtain rotation matrix to perform the second rotation we have done previously 
        rotated_x_axis = rotation_vector(rotated_x_axis, angle2, rotated_normal_vector)


        # Perform just the same procedure as before
        rotated_y_axis = rotation_vector(vector_y_axis, angle2, rotated_normal_vector)


        ##########################################################

        # Compute line in Y axis to project all point on them and find the first and last
        # Formula from https://gamedev.stackexchange.com/questions/72528/how-can-i-project-a-3d-point-onto-a-3d-line
        # A + dot(AP,AB) / dot(AB,AB) * AB

        # Point A
        A_point = (0,0,0)
        
        # As we are going to compute a line in Y axis we just set AB_vector equal to the rotated y axis 
        AB_vector = rotated_y_axis

        # Function to project all the set of point into the line in Y axis
        def get_projected_point (P_point : Tuple[float, float, float]) -> Tuple[float, float, float]:
            AP_vector = [A_point[i] + P_point[i] for i in range(3)]
            AP_AB_dot_product = np.dot(AP_vector, AB_vector)
            AB_AB_dot_product = np.dot(AB_vector, AB_vector)
            scalar = AP_AB_dot_product / AB_AB_dot_product
            projected_vector = [AB_vector[i] * scalar for i in range(3)]
            projected_point = [A_point[i] + projected_vector[i] for i in range(3)]
            return projected_point

        projected_points = [ get_projected_point(P_point) for P_point in coordinates ]
        #print(projected_points)
        projected_points.sort(key=lambda projected_points: projected_points[0])
        projected_points.sort(key=lambda projected_points: projected_points[1])
        projected_points.sort(key=lambda projected_points: projected_points[2])
        max_ypoint = projected_points[0]
        min_ypoint = projected_points[-1]


        y_axis_center = [(max_ypoint[i]+min_ypoint[i])/2 for i in range(3)]

        # Project the center of the protein to the Y line  with respect the center Y point computed before
        # After doing these trials we must read the center point protein from a file that we will generate previously using VMD 
        y_axis_projected_center = get_projected_point(center_point_protein)

        y_axis_difference_vector = [y_axis_projected_center[i] - y_axis_center[i] for i in range(3)]

        y_axis_difference = -math.sqrt(sum([y_axis_difference_vector[i]**2 for i in range(3)]))

        # As rotated_y_axis is normalized we can multiplicate each coordinate between both vectors
        final_vector_rotation_y = [rotated_y_axis[i] * y_axis_difference for i in range(3)]

        file.write('$sel moveby { '+str(final_vector_rotation_y[0])+' '+str(final_vector_rotation_y[1])+' '+str(final_vector_rotation_y[2]) +' } \n')
        
        # Compute line in X axis to project all point on them and find the first and last
        # A + dot(AP,AB) / dot(AB,AB) * AB FORMULA THAT WE ARE APPLYING
        AB_vector = rotated_x_axis
        projected_points2 = [ get_projected_point(P_point) for P_point in coordinates ]
        #print(projected_points)
        projected_points2.sort(key=lambda projected_points2: projected_points2[0])
        projected_points2.sort(key=lambda projected_points2: projected_points2[1])
        projected_points2.sort(key=lambda projected_points2: projected_points2[2])

        max_xpoint = projected_points2[0]

        min_xpoint = projected_points2[-1]
        
        x_axis_center = [(max_xpoint[i]+min_xpoint[i])/2 for i in range(3)]
        # Project the center of the protein to the Y line  with respect the center Y point computed before
        # After doing these trials we must read the center point protein from a file that we will generate previously using VMD 
        x_axis_projected_center = get_projected_point(center_point_protein)
        x_axis_difference_vector = [x_axis_projected_center[i] - x_axis_center[i] for i in range(3)]

        x_axis_difference = math.sqrt(sum([x_axis_difference_vector[i]**2 for i in range(3)]))
        # As rotated_y_axis is normalized we can multiplicate each coordinate between both vectors
        final_vector_rotation_x = [rotated_x_axis[i] * x_axis_difference for i in range(3)]

        file.write('$sel moveby { '+str(final_vector_rotation_x[0])+' '+str(final_vector_rotation_x[1])+' '+str(final_vector_rotation_x[2]) +' } \n')
        file.write(f'scale by {scale_factor} \n') 
        # Show the theoretical center
        vmd_center = str(center_point_protein[0]) + ' ' + str(center_point_protein[1]) + ' ' + str(center_point_protein[2])
        file.write('draw sphere {' + vmd_center + '} radius 50\n')

        # Finally generate the image from the current view
        file.write('render TachyonInternal ' + auxiliary_tga_filename + ' \n')

        file.write('exit\n')

    # Run VMD
    logs = run([
        "vmd",
        input_structure_filename,
        "-e",
        commands_filename,
        "-dispdev",
        "none"
    ], stdout=PIPE, stderr=PIPE).stdout.decode()
#    # If any of the output files do not exist at this point then it means something went wrong with vmd
    expected_output_files = [
        auxiliary_tga_filename
    ]
    for output_file in expected_output_files:
        if not os.path.exists(output_file):
            print(logs)
            raise SystemExit('Something went wrong with VMD')

    im = Image.open(auxiliary_tga_filename)
    # converting to jpg
    rgb_im = im.convert("RGB")
    # exporting the image
    rgb_im.save(output_screenshot_filename)

    
    trash_files = [ commands_filename, auxiliary_tga_filename, center_filename ]
    for trash_file in trash_files:
        os.remove(trash_file)

# Return the rotated vector
# https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
def rotation_vector (vector, angle, axis):
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

# Rotate a vector using a given angle 
def rotate (vector,angle,rotation_axis) -> tuple:
    angle_sin = math.sin(math.radians(angle))
    angle_cos = math.cos(math.radians(angle))

    if rotation_axis == 'z':
        x = + angle_cos * vector[0] - angle_sin * vector[1]
        y = + angle_sin * vector[0] + angle_cos * vector[1]
        z = vector[2]

    elif rotation_axis == 'y':
        x = + angle_cos * vector[0] + angle_sin * vector[2]
        y = vector[1]
        z = - angle_sin * vector[0] + angle_cos * vector[2] 

    else:
        x = vector[0]
        y = + angle_cos * vector[1] - angle_sin * vector[2]
        z = + angle_sin * vector[1] + angle_cos * vector[2]

    #print(str(self) + ' (' + str(angle) + ') ' + str(Vector(x,y)))
    return (x,y,z)