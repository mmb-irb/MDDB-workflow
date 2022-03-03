from model_workflow.tools.xvg_parse import xvg_parse
from model_workflow.tools.get_reduced_trajectory import get_reduced_trajectory

from subprocess import run, PIPE, Popen
import json

# Run multiple RMSD analyses
# A RMSD analysis is run with each reference:
# - First frame
# - Average structure
# A RMSD analysis is run over each rmsd target:
# - Whole protein
# - Heavy atoms
# - Backbone
# - Alpha carbons
def rmsds(
    input_trajectory_filename : str,
    output_analysis_filename : str,
    frames_limit : int,
    first_frame_filename : str,
    average_structure_filename : str,
    rmsd_groups : list = ['Protein', 'Protein-H', 'Backbone', 'C-alpha'],
    ):

    rmsd_references = [first_frame_filename, average_structure_filename]

    # First of all, run an RMSd for the whole trajectory and check there are no sudden changes
    # This RMSd is not saved, but this is only a test to check the inegrity of the simulation
    # Note that it makes no difference which reference is used here
    rmsd_check(rmsd_references[0], input_trajectory_filename)

    # The start will be always 0 since we start with the first frame
    start = 0

    # Reduce the trajectory according to the frames limit
    # Use a reduced trajectory in case the original trajectory has many frames
    # Note that it makes no difference which reference is used here
    reduced_trajectory_filename, step, frames = get_reduced_trajectory(
        rmsd_references[0],
        input_trajectory_filename,
        frames_limit,
    )

    # Save results in this array
    output_analysis = []

    # Iterate over each reference and group
    for reference in rmsd_references:
        # Get a standarized reference name
        reference_name = reference[0:-4].lower()
        for group in rmsd_groups:
            # Get a standarized group name
            group_name = group.lower()
            # Set the analysis filename
            rmsd_analysis = 'rmsd.' + reference_name + '.' + group_name + '.xvg'
            # Run the rmsd
            rmsd(reference, reduced_trajectory_filename, group, rmsd_analysis)
            # Read and parse the output file
            rmsd_data = xvg_parse(rmsd_analysis, ['times', 'values'])
            # Format the mined data and append it to the overall output
            # Multiply by 10 since rmsd comes in nanometers (nm) and we want it in Ångstroms (Å)
            rmsd_values = [ v*10 for v in rmsd_data['values'] ]
            data = {
                'values': rmsd_values,
                'reference': reference_name,
                'group': group_name
            }
            output_analysis.append(data)
            # Remove the analysis xvg file since it is not required anymore
            run([
                "rm",
                rmsd_analysis,
            ], stdout=PIPE).stdout.decode()

    # Export the analysis in json format
    with open(output_analysis_filename, 'w') as file:
        json.dump({ 'start': start, 'step': step, 'data': output_analysis }, file)

# RMSD
# 
# Perform the RMSd analysis 
def rmsd (
    input_reference_filename : str,
    input_trajectory_filename : str,
    group_name : str,
    output_analysis_filename : str) -> str:
    
    # Run Gromacs
    p = Popen([
        "echo",
        group_name, # Select group for least squares fit
        group_name, # Select group for RMSD calculation
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "rms",
        "-s",
        input_reference_filename,
        "-f",
        input_trajectory_filename,
        '-o',
        output_analysis_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

    # Return gromacs logs
    return logs

# Look for sudden raises of RMSd values from one frame to another
def rmsd_check (
    input_topology_filename : str,
    input_trajectory_filename : str
    ):

    print('Checking trajectory integrity')

    # Select the whole protein to check the RMSd
    test_group = 'Protein'

    # Set the name for the output of the test rmsd
    test_filename = 'test.rmsd.xvg'

    # Run Gromacs
    p = Popen([
        "echo",
        test_group, # Select group for least squares fit
        test_group, # Select group for RMSD calculation
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "rms",
        "-s",
        input_topology_filename,
        "-f",
        input_trajectory_filename,
        '-o',
        test_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

    # Read the output and do the check
    test = xvg_parse(test_filename, ['Times', 'Values'])
    values = test['Values']
    previous = values[0]
    for i, value in enumerate(values):
        if abs(value - previous) > 1:
            raise ValueError('There is something wrong with RMSd values. Check frame ' + str(i))
        previous = value

    # Remove the test xvg file since it is not required anymore
    run([
        "rm",
        test_filename,
    ], stdout=PIPE).stdout.decode()