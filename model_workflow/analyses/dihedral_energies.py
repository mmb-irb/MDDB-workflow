import mdtraj as mdt
import numpy as np

from model_workflow.tools.get_reduced_trajectory import get_reduced_trajectory
from model_workflow.utils.type_hints import *

# Calculate torsions and then dihedral energies for every dihedral along the trajectory
def dihedral_energies (
    input_structure_file : 'File',
    input_trajectory_file : 'File',
    output_analysis_filepath : str,
    dihedrals_data : List[dict],
    snapshots : int,
    frames_limit : int,
):
    print('-> Running dihedral energies analysis')

    # If trajectory frames number is bigger than the limit we create a reduced trajectory
    reduced_trajectory_filepath, step, frames = get_reduced_trajectory(
        input_structure_file,
        input_trajectory_file,
        snapshots,
        frames_limit,
    )

    # Load the reduced trajectory
    traj = mdt.load(reduced_trajectory_filepath, top=input_structure_file.path)

    # MDtraj is prepared to analyze all dihedrals at once
    # To do so, prepare a list with every group of atoms conforming a dihedral
    dihedral_atom_indices = [ data['atom_indices'] for data in dihedrals_data ]

    # Compute dihedral torsions
    trajectory_torsions = mdt.compute_dihedrals(traj, dihedral_atom_indices)

    # Results are returned by frame
    # DANI: Creo, porque he hecho las pruebas con una trayectoria de 1 frame xD
    # Iterate frames and calculate dihedran energies for each frame
    for frame_torsions in trajectory_torsions:

        # Iterate every dihedral
        for dihedral_data, torsion in zip(dihedrals_data, frame_torsions):

            # Calculate the dihedral energy according to the formula
            k = dihedral_data['force_constant']
            n = dihedral_data['periodicity']
            δ = dihedral_data['phase']
            φ = torsion

            print(f'Dihedral number (1-based): {dihedral_data["index"] + 1}')
            print(f'Dihedral atom indices (0-based): {dihedral_data["atom_indices"]}')
            print(f'k: {k}')
            print(f'n: {n}')
            print(f'δ: {δ}')
            print(f'φ: {φ}')
            print(f'n*φ - δ: {n*φ - δ}')
            print(f'cos(n*φ - δ): {np.cos(n*φ - δ)}')

            dihedral_energy = (k/2) * (1 + np.cos(n*φ - δ))

            print(f'Dihedral energy: {dihedral_energy} kcal/mol?')


            raise SystemExit('Hold up')