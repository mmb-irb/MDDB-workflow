# Markov

# Set the data needed to represent a Markov State Model graph in the client
# This is finding the most populated frames and the transition probabilities matrix between these frames

# DANI: Este análisis, aunque ya funciona, no se usa
# DANI: Los ejemplos que tenemos de MSM incluyen unas transition probabilities plagadas de 0s
# DANI: Esto hace que al hacer un subset con solo las más populadas a penas capturemos ninguna transición
# DANI: Se podrían prerocesar los datos para hacer clusters o algo similar

import mdtraj as mdt

from mddb_workflow.tools.get_screenshot import get_screenshot 
from mddb_workflow.utils.auxiliar import save_json
from mddb_workflow.utils.type_hints import *

def markov (
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_analysis_filename : str,
    structure : 'Structure',
    populations : list[float],
    transitions : list[list[float]],
    nodes_number : int = 20,
):

    print('-> Running Markov analysis')

    # If there is no populations then we stop here
    if populations is None or len(populations) == 0:
        print(' There are no populations')
        return

    # If there is no transitions then we stop here
    if transitions is None or len(transitions) == 0:
        print(' There are no transitions')
        return

    # Get the numbers of frames with highest populations
    population_per_frames = [ (population, frame) for frame, population in enumerate(populations) ]
    highest_populations = []
    highest_population_frames = []
    for population, frame in sorted(population_per_frames, reverse=True)[0:nodes_number]:
        highest_populations.append(population)
        highest_population_frames.append(frame)
    print(' Reading most equilibrium populated frames in trajectory')
    # Read the trajectory frame by frame looking for the specified frames
    trajectory = mdt.iterload(input_trajectory_filename, top=input_topology_filename, chunk=1)
    # Set a generator for the frames to be selected once sorted
    selected_frames = iter(sorted(highest_population_frames))
    next_frame = next(selected_frames)
    # Conserve only the desired frames
    frame_coordinates = {}
    for frame_number, frame in enumerate(trajectory):
        # Update the current frame log
        print(f' Frame {frame_number}', end='\r')
        # Skip the current frame if we do not need it
        if frame_number != next_frame:
            continue
        # Save it otherwise
        frame_coordinates[frame_number] = frame
        # Update the next frame
        next_frame = next(selected_frames, None)
        if next_frame == None:
            break
    print(' Building transition probability matrix')
    # Get a subset of transitions only for the selected frames
    transitions_matrix = []
    for frame in highest_population_frames:
        row = []
        for other_frame in highest_population_frames:
            transition = transitions[frame][other_frame]
            row.append(transition)
        transitions_matrix.append(row)
    # Make a copy of the structure to avoid mutating the original structure
    reference_structure = structure.copy()
    print(' Taking screenshots of selected frames')
    frame_count = len(frame_coordinates)
    # For each frame coordinates, generate PDB file, take a scrrenshot and delete it
    for i, frame in enumerate(frame_coordinates.values(), 1):
        # Update the current frame log
        print(f'  Screenshot {i}/{frame_count}', end='\r')
        # Get the actual coordinates
        coordinates = frame.xyz[0] * 10 # We multiply by to restor Ångstroms
        # Update the reference structure coordinates
        reference_structure.set_new_coordinates(coordinates)
        # Set the screenshot filename
        screenshot_filename = f'markov_screenshot_{str(i).zfill(2)}.jpg'
        # Generate the screenshot
        get_screenshot(reference_structure, screenshot_filename)
    # Export the analysis data to a json file
    data = {
        'frames': highest_population_frames,
        'populations': highest_populations,
        'transitions': transitions_matrix
    }
    save_json(data, output_analysis_filename)