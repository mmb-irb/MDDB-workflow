# Principal component analysis (PCA)
import os
import math
from subprocess import run, PIPE, Popen

from model_workflow.tools.get_reduced_trajectory import get_reduced_trajectory

# Set a name for the pca average file
# If not specified, this file is create with the name 'average.pdb'
# This name enters in conflict with the average filename so it must be changed
pca_average_filename = 'pca.average.pdb'

# Perform the PCA of the trajectory
# The PCA is performed with the whole trajectory when the size is reasonable
# When the trajectory is larger than 2000 snapshots we reduce the trajectory
# A new projection trajectory is made for each eigen vector with 1% or greater explained variance
# WARNING: Projection trajectories are backbone only
def pca(
    input_topology_filename: str,
    input_trajectory_filename: str,
    output_eigenvalues_filename: str,
    output_eigenvectors_filename: str,
    frames_limit : int,
    pca_fit_selection : str = 'Protein-H',
    pca_selection : str = 'Backbone'
):

    # By default we set the whole trajectory as PCA trajectory
    pca_trajectory_filename = input_trajectory_filename
    # If trajectory frames number is bigger than the limit we create a reduced trajectory
    pca_trajectory_filename, step, frames = get_reduced_trajectory(
        input_topology_filename,
        input_trajectory_filename,
        frames_limit,
    )

    # Calculate eigen values and eigen vectors with Gromacs
    p = Popen([
        "echo",
        pca_fit_selection,
        pca_selection,
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "covar",
        "-s",
        input_topology_filename,
        "-f",
        pca_trajectory_filename,
        '-o',
        output_eigenvalues_filename,
        '-v',
        output_eigenvectors_filename,
        '-av',
        pca_average_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

    if not os.path.exists(output_eigenvalues_filename):
        print(logs)
        raise SystemExit('Something went wrong with GROMACS')

    # Read the eigen values file and get an array with all eigen values
    values = []
    with open(output_eigenvalues_filename, 'r') as file:
        for line in file:
            if line.startswith(("#", "@")):
                continue
            else:
                values.append(float(line.split()[1]))

    # Get the total eigen value by adding all eigen values
    total = 0
    for value in values:
        total += value

    # Count how many eigen values are greater than 1% of the total eigen value
    # Eigen values are ordered from greater to lower by default, so we stop at the first value lower than 1%
    greater = 0
    cutoff = total / 100
    for value in values:
        if value >= cutoff:
            greater += 1
        else:
            break

    # Now make a projection for each suitable eigen vector
    for ev in range(1, greater+1):
        strev = str(ev)
        # Set the name of the new projection analysis
        projection = 'pca.proj' + strev + '.xvg'
        # Set the name of the new projection trajectory
        projection_trajectory = 'md.pca-' + strev + '.xtc'
        # UNKNOWN USE
        pca_rmsf = 'pca.rmsf' + strev + '.xvg'

        # Perform the projection analysis through the 'anaeig' gromacs command
        p = Popen([
            "echo",
            pca_fit_selection,
            pca_selection,
        ], stdout=PIPE)
        logs = run([
            "gmx",
            "anaeig",
            "-s",
            input_topology_filename,
            "-f",
            pca_trajectory_filename,
            '-eig',
            output_eigenvalues_filename,
            '-v',
            output_eigenvectors_filename,
            '-proj',
            projection,
            '-extr',
            projection_trajectory,
            '-rmsf',
            pca_rmsf,
            '-nframes',
            '20',
            '-first',
            strev,
            '-last',
            strev,
            '-quiet'
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()