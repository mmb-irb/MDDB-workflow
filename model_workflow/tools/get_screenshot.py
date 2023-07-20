from subprocess import run, PIPE
import os
from typing import List
from PIL import Image

auxiliary_tga_filename = '.transition_screenshot.tga'
# Python function to obtain a screenshot from the pdb file using VMD a molecular modelling and visualization computer program
def get_screenshot(
    input_structure_filename : str,
    # input_coordinates_filename : List[float] 
    # Remains like this as we are not able to use a set of values for the rotation, translation and zoom and apply them to the molecule
    output_screenshot_filename : str
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
    commands_filename = 'commands.vmd' # Generate a file name for the commands file
    with open(commands_filename, "w") as file:
        
        # Set the Background of the molecule to white color
        file.write('color Display Background white \n') 
        
        # Delete the axes drawing in the VMD window
        file.write('axes location Off \n') 
        
        # Eliminate the molecule to perform changes in the representation, color and material 
        file.write('mol delrep 0 top \n') 
        
        # Change the default representation model to Newcartoon
        file.write('mol representation Newcartoon \n') 
        
        # Change the default atom coloring method setting to Chain
        file.write('mol color Chain \n') 
        
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
        file.write(f'scale by {scale_factor} \n') 
        file.write('render TachyonInternal ' + auxiliary_tga_filename+ ' \n')
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

    # If any of the output files do not exist at this point then it means something went wrong with vmd
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

    trash_files = [ commands_filename, auxiliary_tga_filename ]
    for trash_file in trash_files:
        os.remove(trash_file)




