import mdtraj as mdt
import numpy as np

from model_workflow.tools.get_reduced_trajectory import get_reduced_trajectory
from model_workflow.utils.type_hints import *
from model_workflow.utils.auxiliar import save_json, mean
from model_workflow.utils.constants import OUTPUT_DIHEDRAL_ENERGIES_FILENAME


def compute_dihedral_energies (
    structure_file : 'File',
    trajectory_file : 'File',
    output_directory : str,
    dihedrals : list[dict],
    snapshots : int,
    frames_limit : int,
):
    """Calculate torsions and then dihedral energies for every dihedral along the trajectory."""
    # Set the main output filepath
    output_analysis_filepath = f'{output_directory}/{OUTPUT_DIHEDRAL_ENERGIES_FILENAME}'

    # Set a dict with the results of energies calculations
    # Atom indices are used as keys
    dihedral_energies = {}
    for dihedral_data in dihedrals:
        atom_indices = dihedral_data['atom_indices']
        dihedral_energies[atom_indices] = {
            'atom_indices': list(atom_indices),
            'torsion': [],
            'ee': [],
            'vdw': [],
        }

    # If trajectory frames number is bigger than the limit we create a reduced trajectory
    reduced_trajectory_filepath, step, frames = get_reduced_trajectory(
        structure_file, trajectory_file, snapshots, frames_limit)

    # Load the reduced trajectory
    traj = mdt.load(reduced_trajectory_filepath, top=structure_file.path)

    # MDtraj is prepared to analyze all dihedrals at once
    # To do so, prepare a list with every group of atoms conforming a dihedral
    dihedral_atom_indices = [ data['atom_indices'] for data in dihedrals ]

    # Compute dihedral torsions using MDtraj
    trajectory_torsions = mdt.compute_dihedrals(traj, dihedral_atom_indices)

    # Results are returned by frame
    # DANI: Creo, porque he hecho las pruebas con una trayectoria de 1 frame xD
    # Iterate frames and calculate dihedral energies for each frame
    for frame_torsions in trajectory_torsions:

        # Iterate every dihedral with its torsion
        for dihedral_data, torsion in zip(dihedrals, frame_torsions):
            atom_indices = dihedral_data['atom_indices']
            # Dihedral torsion energy is the sum of the torsion energy of its terms
            torsion_energy = 0
            # Iterate dihedral terms
            for term in dihedral_data['terms']:
                # Get dihedral term parameters
                k = term['force']
                if k == 0: continue
                n = term['period']
                δ = term['phase']
                φ = torsion
                # Calculate the dihedral torsion energy according to the formula
                term_energy = (k/2) * (1 + np.cos(n*φ - δ))
                torsion_energy += term_energy
            # Add torison energies to the dihedral energies object
            dihedral_energies[atom_indices]['torsion'].append(torsion_energy)

    # Now calculate non covalent energies
    # These energies are always computed between atoms 1-4 in the dihedral
    # i.e. atoms in both ends
    dihedral_1_4_indices = [ [ data['atom_indices'][0], data['atom_indices'][3] ] for data in dihedrals ]

    # Compute atom distance along the trajectory using MDtraj
    trajectory_distances = mdt.compute_distances(traj, dihedral_1_4_indices)
            
    # Results are returned by frame
    # Iterate frames and calculate dihedral energies for each frame
    for frame_distances in trajectory_distances:

        # Iterate every dihedral with its 1-4 distance
        for dihedral_data, distance in zip(dihedrals, frame_distances):
            # Get scaling factors for both electrostatic and Van Der Waals parameters
            # Note that these may be missing for some dihedrals since all its terms are to be ignored
            # SANTI: Cuando esto pasa hay que caluclar la energía igualmente, pero sin escalarla
            scee = dihedral_data.get('ee_scaling', 1)
            scnb = dihedral_data.get('vdw_scaling', 1)
            # get Leonard-Johns constants to calculate Van Der Waals energies
            acoef = dihedral_data['lj_acoef']
            bcoef = dihedral_data['lj_bcoef']
            # Note that charges are not actual atom charges, but they are already scaled
            q1 = dihedral_data['atom1_charge']
            q4 = dihedral_data['atom4_charge']
            # MDtraj outputs in nm and we want distance in Angstroms
            r14 = distance * 10
            # Calculate the dihedral electrostatic energy according to the formula
            ee_energy = (1/scee) * (q1 * q4) / r14
            # Calculate the dihedral Var Der Waals energy according to the formula
            vdw_energy = (1/scnb) * (( acoef / (r14 ** 12) ) - ( bcoef / (r14 ** 6) ))
            # Add non covalent energies to the dihedral energies object
            atom_indices = dihedral_data['atom_indices']
            dihedral_energies[atom_indices]['ee'].append(ee_energy)
            dihedral_energies[atom_indices]['vdw'].append(vdw_energy)

    # Reformat output data by calculating average values instead of per-frame values
    output_data = []
    for energies in dihedral_energies.values():
        output_data.append({
            'indices': energies['atom_indices'],
            'torsion': mean(energies['torsion']),
            'ee': mean(energies['ee']) if len(energies['ee']) > 0 else None,
            'vdw': mean(energies['vdw']) if len(energies['vdw']) > 0 else None,
        })

    # Write results to disk
    save_json(output_data, output_analysis_filepath)