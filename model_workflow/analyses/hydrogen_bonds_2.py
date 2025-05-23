# Hydrogen bonds analysis

# The pytraj 'hbond' analysis finds possible hydrogen bonds interactions between closer atoms and 
# calculates the percent of frames where the hydrogen bond would happen according to the distance 
# between two atoms and the bond angle.
#
# Pytraj analysis will return the residue number, residue type and atom name of each pair of bonded 
# atoms. WARNING: In irregular pdb topologies, there may be atoms with identical residue number, 
# residue type and atom name. This makes impossible to know the real resulted atom. For this reason 
# topology has been corrected previously and duplicated atoms have been renamed

import pytraj as pt
import re
from math import ceil
from os.path import exists

from model_workflow.utils.pyt_spells import get_pytraj_trajectory
from model_workflow.utils.auxiliar import save_json, numerate_filename, get_analysis_name
from model_workflow.utils.type_hints import *

# Perform an hydrogen bonds analysis for each interaction interface
# The 'interactions' input may be an empty list (i.e. there are no interactions)
# In case there are no interactions the analysis stops
# Note that this analysis find hydrogen bonds in a subset of frames along the trajectory
# Storing the results for the whole trajectory is not possible due to storage limits
def hydrogen_bonds (
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_analysis_filepath : str,
    structure : 'Structure',
    interactions : list,
    #snapshots : int,
    is_time_dependend : bool,
    time_splits : int,
    populations : Optional[List[float]],
    # Explicit values for the most populated frames is saved apart
    # Set how many frames (in order of population) are saved
    most_populated_frames_number : int = 5,
    ):

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
    pt_trajectory = get_pytraj_trajectory(input_topology_filename, input_trajectory_filename)

    # Save the reference function to get an absolue atom index from a pytraj atom index
    get_atom_index = structure.get_atom_index
    # Save the reference function to get a source residue from a pytraj residue
    pytraj_residue_index_2_residue = structure.pytraj_residue_index_2_residue

    # Save each analysis to a dict which will be parsed to json
    output_summary = []

    # Iterate interactions
    for i, interaction in enumerate(interactions):

        # If the interaction has coarse grain atoms then just skip it
        # Note that this analysis makes not sense if the interaction is not atomistic
        # FUN FACT: If tried, pytraj fails so spectacularly it is is not even able to report the error
        # RuntimeError: Failed to setup action. Use pytraj._verbose() to turn on the error report.
        if interaction.get('has_cg', False): continue

        # Get the interaction name
        name = interaction['name']
        # Set a filename for the current interaction data
        numbered_output_analysis_filepath = numerate_filename(output_analysis_filepath, i)

        # Add the root of the output analysis filename to the run data
        analysis_name = get_analysis_name(numbered_output_analysis_filepath)
        # Append current interaction to the summary
        output_summary.append({
            'name': name,
            'analysis': analysis_name
        })

        # If the analysis already exists then proceed to the next interaction
        if exists(numbered_output_analysis_filepath):
            continue
        
        # Select all interface residues in pytraj notation
        pt_interface = interaction['pt_interface_1'] + interaction['pt_interface_2']
        if len(pt_interface) == 0:
            raise ValueError(f'There are no interface residues for interaction "{name}"')
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

        # Get the number of snapshots from one analysis result
        snapshots = len(hbonds[0].values)
        # Calculate in how many frames we must split every time slot
        step = ceil(snapshots / time_splits)
        # Calculate how many time slots we will have at the end
        # Note that there is a residual number of frames at the end which is discarded
        nsteps = ceil(snapshots / step)
         # In case we have populations
        if populations:
            # Set a function to calculate the population-weighted value of a given frame
            def get_populated_value (hbond_value : float, frame : int) -> float:
                population = populations[frame]
                multiplier = population * snapshots
                return hbond_value * multiplier
            # Get the most populated frames numbers
            population_per_frames = [ (population, frame) for frame, population in enumerate(populations) ]
            population_sorted_frames = [frame for population, frame in sorted(population_per_frames, reverse=True)]
            most_populated_frames = population_sorted_frames[0:most_populated_frames_number]

        # Get residues in each interaction agent interface
        interface_1 = interaction['interface_1']
        interface_2 = interaction['interface_2']

        acceptor_atom_index_list = []
        donor_atom_index_list = []
        hydrogen_atom_index_list = []
        hbond_overall = []
        hbond_timed = []
        hbond_framed = []
        
        # Search all predicted hydrogen bonds
        for i, hb in enumerate(hbonds):
            key = hbond_keys[i]
            # hb.key example: "ASN15_OD1-LEU373_N-H"
            # matchObj = re.match( r'\w{3}(\d*)_(.*)-\w{3}(\d*)_(.*)-(.*)', hb.key, re.M|re.I)
            # hb.key example: "ASN_15@OD1-LEU_373@N-H"
            matchObj = re.match( r'\w*_(\d*)@(.*)-\w*_(\d*)@(.*)-(.*)', key, re.M|re.I)
            if matchObj == None:
                continue
            # Mine all data from the parsed results
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
            if (acceptor in interface_1 and donor in interface_1) or (acceptor in interface_2 and donor in interface_2):
                continue
            # Get the absolute index of each atom
            acceptor_atom_index_list.append(get_atom_index(acceptor, acceptor_atom))
            donor_atom_index_list.append(get_atom_index(donor, donor_atom))
            hydrogen_atom_index_list.append(get_atom_index(donor, hydrogen_atom))
            # List 'hb.values' because it is a ndarray, which is not JSON serializable
            hbond_values = list(map(int, hb.values))
            # In case we have populations
            if populations:
                # Save if there is or not hydrogen bond for the most populated frames
                framed_hbonds = { frame: bool(hbond_values[frame]) for frame in most_populated_frames }
                hbond_framed.append(framed_hbonds)
                # Weigth hbond values according to populations
                hbond_values = [ get_populated_value(hbond_value, frame) for frame, hbond_value in enumerate(hbond_values) ]
            # Get the overall percent of frames where the hydrogen bond exists
            overall_percent = sum(hbond_values) / snapshots
            hbond_overall.append(overall_percent)
            # The last part is done only when the simulation is time dependent
            if not is_time_dependend:
                continue
            # Now split the whole trajectory in time slots and calculate hbonds percents for every slot
            temporal_percents = []
            for i in range(nsteps):
                slot_hbond_values = hbond_values[i*step:(i+1)*step]
                temporal_percent = sum(slot_hbond_values) / snapshots
                temporal_percents.append(temporal_percent)
            hbond_timed.append(temporal_percents)

        # Set the interaction output
        interaction_data = {
            'name': interaction['name'],
            'acceptors': acceptor_atom_index_list,
            'donors': donor_atom_index_list,
            'hydrogens': hydrogen_atom_index_list,
            'hbonds_overall': hbond_overall
        }
        if is_time_dependend:
            interaction_data['hbonds_timed'] = hbond_timed
        if populations:
            interaction_data['hbonds_framed'] = hbond_framed

        # Write the interaction analysis output to a file
        save_json(interaction_data, numbered_output_analysis_filepath)

    # Export the analysis in json format
    if len(output_summary) > 0:
        save_json(output_summary, output_analysis_filepath)