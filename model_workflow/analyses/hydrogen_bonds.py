# Hydrogen bonds analysis

# The pytraj 'hbond' analysis finds possible hydrogen bonds interactions between closer atoms and 
# calculates the percent of frames where the hydrogen bond would happen according to the distance 
# between two atoms and the bond angle.
#
# This analysis will return the residue number, residue type and atom name of each pair of bonded 
# atoms. WARNING: In irregular pdb topologies, there may be atoms with identical residue number, 
# residue type and atom name. This makes impossible to know the real resulted atom. For this reason 
# topology has been corrected previously and duplicated atoms have been renamed

import pytraj as pt
import numpy
import re

from model_workflow.utils.pyt_spells import get_reduced_pytraj_trajectory
from model_workflow.utils.auxiliar import save_json
from model_workflow.utils.type_hints import *

# Perform an hydrogen bonds analysis for each interaction interface
# The 'interactions' input may be an empty list (i.e. there are no interactions)
# In case there are no interactions the analysis stops
def hydrogen_bonds (
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_analysis_filename : str,
    structure : 'Structure',
    interactions : list,
    snapshots : int,
    frames_limit : int):

    print('-> Running hydrogen bonds analysis')

    # Return before doing anything if there are no interactions
    if not interactions or len(interactions) == 0:
        print('No interactions were specified')
        return

    # Get all not failed interactions
    valid_interactions = [ interaction for interaction in interactions if not interaction.get('failed', False) ]
    
    # Make sure we have valid interactions
    # DANI: Esto es temporal, lo suyo sería que las interacciones válidas si sean analizadas
    # DANI: Lo que pasa es que pronto cambiaré los análisis de interacciones para que se haga 1 por interacción
    # DANI: De manera que no merece la pena invertir tiempo en dar soporte a esto ahora
    if len(valid_interactions) != len(interactions):
        print('There are no valid interactions -> This analysis will be skipped')
        return

    # Parse the trajectory intro ptraj
    # Reduce it in case it exceeds the frames limit
    pt_trajectory, frames_step, frames_count = get_reduced_pytraj_trajectory(input_topology_filename, input_trajectory_filename, snapshots, frames_limit)

    # Save the reference function to get an absolue atom index from a pytraj atom index
    get_atom_index = structure.get_atom_index
    # Save the reference function to get a source residue from a pytraj residue
    pytraj_residue_index_2_residue = structure.pytraj_residue_index_2_residue

    # Save each analysis to a dict which will be parsed to json
    output_analysis = []

    for interaction in interactions:
        
        # Select all interface residues in pytraj notation
        pt_interface = interaction['pt_interface_1'] + interaction['pt_interface_2']
        if len(pt_interface) == 0:
            raise ValueError('There are no interface residues for interaction "' + interaction['name'] + '"')
        pt_selection = ':' + ','.join(map(str, pt_interface))

        # Run the analysis
        hbonds = pt.hbond(pt_trajectory, mask=pt_selection)

        # Save appart hbond 'old keys' which are not accessible for each hb individually
        # WARNING: These keys are internal pytraj keys but they are crucial
        # e.g. current key: 'CYS44_O-THR25_OG1-HG1', old key: 'CYS_44@O-THR_25@OG1-HG1'
        # WARNING: The fact that residue name and residue number is separated by '_' is critical
        # Some topologies have residue names with different length (2 or 3 characters)
        # Some topologies have residue names which end in numeric characters (e.g. 'P1', 'P2', etc.)
        # For this reason we can not use current keys where residue name and number are merged
        # This makes impossible to know what is the name and what is the number of the residue sometimes
        hbond_keys = hbonds._old_keys

        # Get residues in each interaction agent interface
        interface_1 = interaction['interface_1']
        interface_2 = interaction['interface_2']

        acceptor_atom_index_list = []
        donor_atom_index_list = []
        hydrogen_atom_index_list = []
        hbond_values = []
        
        # Search all predicted hydrogen bonds
        for i, hb in enumerate(hbonds):
            key = hbond_keys[i]
            # hb.key example: "ASN15_OD1-LEU373_N-H"
            # matchObj = re.match( r'\w{3}(\d*)_(.*)-\w{3}(\d*)_(.*)-(.*)', hb.key, re.M|re.I)
            # hb.key example: "ASN_15@OD1-LEU_373@N-H"
            matchObj = re.match( r'\w*_(\d*)@(.*)-\w*_(\d*)@(.*)-(.*)', key, re.M|re.I)
            # Mine all data from the parsed results
            if matchObj is not None:
                acceptor_resnum = matchObj.group(1)
                acceptor_atom = matchObj.group(2)
                donor_resnum = matchObj.group(3)
                donor_atom = matchObj.group(4)
                hydrogen_atom = matchObj.group(5)

                # Get the acceptor and donor residues in source notation
                acceptor = pytraj_residue_index_2_residue(int(acceptor_resnum))
                donor = pytraj_residue_index_2_residue(int(donor_resnum))
                
                # WARNING: The analysis may return hydrogen bonds between residues from the same agent
                # Accept the hydrogen bond only if its residues belong to different interaction agents
                if ((acceptor in interface_1 and donor in interface_2)
                or (acceptor in interface_2 and donor in interface_1)):
                    
                    # Get the absolute index of each atom
                    acceptor_atom_index_list.append(get_atom_index(acceptor, acceptor_atom))
                    donor_atom_index_list.append(get_atom_index(donor, donor_atom))
                    hydrogen_atom_index_list.append(get_atom_index(donor, hydrogen_atom))
                    # List 'hb.values' because it is a ndarray, which is not JSON serializable
                    # Then convert the 0s and 1s into trues and falses
                    # This is because booleans are lighter once parsed in hte database
                    values = list( map(bool, (map(int, hb.values))) )
                    #hbond_values.append(' '.join(map(str, hb.values)))
                    hbond_values.append(values) 

        # Write 
        output_analysis.append(
            {
                'name': interaction['name'],
                'acceptors': acceptor_atom_index_list,
                'donors': donor_atom_index_list,
                'hydrogens': hydrogen_atom_index_list,
                'hbonds': hbond_values,
            }
        )

    # Export the analysis in json format
    save_json({ 'data': output_analysis }, output_analysis_filename)
