import tmscoring

from subprocess import run, PIPE, Popen
import json

from model_workflow.tools.get_reduced_trajectory import get_reduced_trajectory

# TM scores
# 
# Perform the tm score using the tmscoring package
def tmscores (
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_analysis_filename : str,
    snapshots : int,
    first_frame_filename : str,
    average_structure_filename : str):

    tmscore_references  = [first_frame_filename, average_structure_filename]

    start = 0
    step = None

    # Set the frames where we extract energies to calculate the average
    # WARNING: The gromacs '-fr' option counts frames starting at 1, not at 0
    frames = range(1, snapshots +1)

    # Set a maximum of frames
    # If trajectory has more frames than the limit create a reduced trajectory
    tmscore_trajectory_filename = input_trajectory_filename
    frames_number = 200
    if snapshots > frames_number:
        tmscore_trajectory_filename = 'tmscore.trajectory.xtc'
        step = get_reduced_trajectory(
            input_topology_filename,
            input_trajectory_filename,
            tmscore_trajectory_filename,
            snapshots,
            frames_number,
        )
        # WARNING: The gromacs '-fr' option counts frames starting at 1, not at 0
        frames = range(1, frames_number +1)
    
    output_analysis = []

    # The only possible group in TM score is the alpha carbon
    # They are the only atoms taken in count by the TM score algorithm
    group = 'C-alpha'

    # Get a standarized group name
    group_name = 'c-alpha'

    frames_ndx = 'frames.ndx'
    # Iterate over each reference and group
    for reference in tmscore_references:
        # Get a standarized reference name
        reference_name = reference[0:-4].lower()
        # Create a reference topology with only the group atoms
        # WARNING: Yes, TM score would work also with the whole reference, but it takes more time!!
        # This has been experimentally tested and it may take more than the double of time
        grouped_reference = 'gref.pdb'
        p = Popen([
            "echo",
            group,
        ], stdout=PIPE)
        logs = run([
            "gmx",
            "trjconv",
            "-s",
            reference,
            "-f",
            reference,
            '-o',
            grouped_reference,
            "-dump",
            "0",
            '-quiet'
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()
        # Get the TM score of each frame
        # It must be done this way since tmscoring does not support trajectories
        tmscores = []
        for f in frames:
            # Extract the current frame
            current_frame = 'frame' + str(f) + '.pdb'
            # The frame selection input in gromacs works with a 'ndx' file
            with open(frames_ndx, 'w') as file:
                file.write('[frames]\n' + str(f))
            p = Popen([
                "echo",
                group,
            ], stdout=PIPE)
            logs = run([
                "gmx",
                "trjconv",
                "-s",
                reference,
                "-f",
                tmscore_trajectory_filename,
                '-o',
                current_frame,
                "-fr",
                frames_ndx,
                '-quiet'
            ], stdin=p.stdout, stdout=PIPE).stdout.decode()
            p.stdout.close()

            # Run the tmscoring over the current frame against the current reference
            # Append the result data for each ligand
            tmscore = tmscoring.get_tm(grouped_reference, current_frame)
            tmscores.append(tmscore)

            # Delete current frame files before going for the next frame
            run([
                "rm",
                current_frame,
                frames_ndx,
            ], stdout=PIPE).stdout.decode()

        # Save the tmscores in the output object
        data = {
            'values': tmscores,
            'reference': reference_name,
            'group': group_name
        }
        output_analysis.append(data)

        run([
            "rm",
            grouped_reference,
        ], stdout=PIPE).stdout.decode()

    # Export the analysis in json format
    with open(output_analysis_filename, 'w') as file:
        json.dump({ 'start': start, 'step': step, 'data': output_analysis }, file)

    if tmscore_trajectory_filename == 'tmscore.trajectory.xtc':
        run([
            "rm",
            tmscore_trajectory_filename,
        ], stdout=PIPE).stdout.decode()