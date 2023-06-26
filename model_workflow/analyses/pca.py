# Principal component analysis (PCA)
from sklearn.decomposition import PCA
import numpy as np
import json

from typing import List

import mdtraj as md

from model_workflow.tools.get_reduced_trajectory import get_reduced_trajectory

# Perform the PCA analysis
def pca (
    input_topology_filename: str,
    input_trajectory_filename: str,
    output_analysis_filename: str,
    output_trajectory_projections_prefix : str,
    snapshots : int,
    frames_limit : int,
    structure : 'Structure',
    fit_selection : str,
    analysis_selection : str,
    pbc_residues : List[int],
    projection_frames : int = 20
) -> dict:

    # By default we set the whole trajectory as PCA trajectory
    pca_trajectory_filename = input_trajectory_filename
    # If trajectory frames number is bigger than the limit we create a reduced trajectory
    pca_trajectory_filename, step, frames = get_reduced_trajectory(
        input_topology_filename,
        input_trajectory_filename,
        snapshots,
        frames_limit,
    )

    # Parse PBC residues into a selection to be removed from further PCA selections
    pbc_selection = structure.select_residue_indices(pbc_residues)

    # Parse the string selections
    # VMD selection syntax
    parsed_fit_selection = structure.select(fit_selection, syntax='vmd') - pbc_selection
    parsed_analysis_selection = structure.select(analysis_selection, syntax='vmd') - pbc_selection

    # Load the trajectory
    mdtraj_trajectory = md.load(pca_trajectory_filename, top=input_topology_filename)
    # Fit the trajectory according to the specified fit selection
    mdtraj_trajectory.superpose(mdtraj_trajectory, frame=0, atom_indices=parsed_fit_selection.atom_indices)
    # Filter the atoms to be analized
    atom_indices = parsed_analysis_selection.atom_indices
    mdtraj_trajectory = mdtraj_trajectory.atom_slice(atom_indices=atom_indices)
    # Reshape data to a sklearn-friendly format
    frames_number = mdtraj_trajectory.n_frames
    atoms_number = mdtraj_trajectory.n_atoms
    reshape = mdtraj_trajectory.xyz.reshape(frames_number, atoms_number * 3)

    # Perform the PCA using sklearn
    pca = PCA()
    transformed = pca.fit_transform(reshape)
    # DANI: Al transponer esto las proyecciones parece que tienen más sentido
    transformed = transformed.transpose()

    # Get eigenvalues
    # Multiply values by 100 since values are in namometers (squared) and we want Ångstroms
    eigenvalues = [ ev * 100 for ev in pca.explained_variance_ ]
    # Get eigenvectors
    eigenvectors = [ [ float(v) for v in eigenvector ] for eigenvector in pca.components_ ]

    # Get the mean structure coordinates
    mean = pca.mean_.flatten()

    # Get the total explained variance by adding all eigenvalues
    total = sum(eigenvalues)

    # Now get the projection for those eigenvectors whose eigenvalue is greater than 1% of the total explained variance
    # Eigenvalues are ordered from greater to lower by default, so we stop at the first value lower than 1%
    # Save projections from those eigenvectors
    projections = []
    cutoff = total / 100
    for i, value in enumerate(eigenvalues):
        if value < cutoff:
            break
        # This logic was copied from here:
        # https://userguide.mdanalysis.org/stable/examples/analysis/reduced_dimensions/pca.html
        eigenvector = eigenvectors[i]
        frame_projections = transformed[i]
        offset = np.outer(frame_projections, eigenvector)
        trajectory_projection = mean + offset
        coordinates = trajectory_projection.reshape(len(frame_projections), -1, 3)
        # Now we have the time dependent projection of the principal component
        # However, we will sort frames according to the projection value
        # In addition we will take only a few frames
        max_projection = max(frame_projections)
        min_projection = min(frame_projections)
        projection_step = (max_projection - min_projection) / (projection_frames - 1)
        selected_projections = [ min_projection + projection_step * s for s in range(projection_frames) ]
        selected_coordinates = [ coordinates[get_closer_value_index(frame_projections, p)] for p in selected_projections ]
        # Load coordinates in mdtraj and export the trajectory to xtc
        trajectory_projection = md.Trajectory(selected_coordinates, mdtraj_trajectory.topology)
        trajectory_projection_filename = output_trajectory_projections_prefix + '_' + str(i+1).zfill(2) + '.xtc'
        trajectory_projection.save_xtc(trajectory_projection_filename)
        # Save projections to be further exported to json
        projections.append([ float(v) for v in frame_projections ])

    # DANI: Usa esto para generar una estructura (pdb) que te permita visualizar las trajectory projections
    #trajectory_projection[0].save_pdb('pca.trajectory_projection_structure.pdb')

    # Set the output dict
    data = {
        'framestep': step,
        'atoms': atom_indices,
        'eigenvalues': eigenvalues,
        'projections': projections
    }
    # Finally, export the analysis in json format
    with open(output_analysis_filename, 'w') as file:
        json.dump({'data': data}, file)

# Given an array with numbers, get the index of the value which is closer
def get_closer_value_index (list : list, value : float) -> int:
    distances = [ abs(value-v) for v in list ]
    return min(range(len(distances)), key=distances.__getitem__)
