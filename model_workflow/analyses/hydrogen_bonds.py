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
from model_workflow.utils.auxiliar import save_json, numerate_filename, get_analysis_name, reprint
from model_workflow.utils.constants import OUTPUT_HBONDS_FILENAME
from model_workflow.utils.type_hints import *

# WARNING: the output file size depends on the number of hydrogen bonds
# WARNING: analyses must be no heavier than 16Mb in BSON format
# WARNING: In case of large surface interaction the output analysis may be larger than the limit
def hydrogen_bonds (
    structure_file : 'File',
    trajectory_file : 'File',
    output_directory : str,
    populations : Optional[List[float]],
    structure : 'Structure',
    interactions : list,
    is_time_dependent : bool,
    # Number of splits along the trajectory
    time_splits : int = 100,
    # Explicit values for the most populated frames is saved apart
    # Set how many frames (in order of population) are saved
    most_populated_frames_number : int = 5,
    ):
    """Perform an hydrogen bonds analysis for each interaction interface.
    The 'interactions' input may be an empty list (i.e. there are no interactions).
    In case there are no interactions the analysis stops.
    Note that this analysis find hydrogen bonds in a subset of frames along the trajectory.
    Storing the results for the whole trajectory is not possible due to storage limits."""
    # Return before doing anything if there are no interactions
    if not interactions or len(interactions) == 0:
        print('No interactions were specified')
        return
    
    # Set the main output filepath
    output_analysis_filepath = f'{output_directory}/{OUTPUT_HBONDS_FILENAME}'
    
    # Set a reference system to handle conversions to pytraj topology
    # First set the pytraj topology
    pytraj_topology = structure.get_pytraj_topology()
    pytraj_residues = list(pytraj_topology.residues)

    # Transform a pytraj residue numeration (1-based) and an atom name to a 0-based atom index
    # Double check the residue is the one we expect to find and, if not, use another strategy
    # Remember that sometimes pytraj may parse residues differently for chaotic topologies
    def pytraj_residue_number_2_residue (pytraj_residue_index : int) -> 'Residue':
        residue_index = pytraj_residue_index - 1
        pytraj_residue = pytraj_residues[residue_index]
        expected_number = pytraj_residue.original_resid
        expected_name = pytraj_residue.name
        # In the canonical way this index is equivalent to the structure resiude index
        if residue_index < len(structure.residues):
            residue = structure.residues[residue_index]
            if residue.number == expected_number and residue.name[0:3] == expected_name:
                return residue
        # Pytraj index may not match the structure index in caotic topologies
        # (i.e. when heavy atoms and hydrogen are listed independently)
        # When this happens we can try to find the residue by comparing resnum and resname
        # WARNING: Note that this alternative method is nos sensitive to chains or icodes
        # WARNING: This is because pytraj does not deal with chains or icodes
        for residue in structure.residues:
            if residue.number == expected_number and residue.name[0:3] == expected_name:
                return residue
        # Return None if there are no results    
        return None

    # Parse the trajectory intro ptraj
    # Reduce it in case it exceeds the frames limit
    pt_trajectory = get_pytraj_trajectory(structure_file.path, trajectory_file.path)

    # Save each analysis to a dict which will be parsed to json
    output_summary = []

    # Iterate interactions
    print()
    for i, interaction in enumerate(interactions):

        # If the interaction has coarse grain atoms then just skip it
        # Note that this analysis makes not sense if the interaction is not atomistic
        # FUN FACT: If tried, pytraj fails so spectacularly it is is not even able to report the error
        # RuntimeError: Failed to setup action. Use pytraj._verbose() to turn on the error report.
        if interaction.get('has_cg', False): continue

        # Get the interaction name
        name = interaction['name']
        reprint(f' Processing {name} ({i+1}/{len(interactions)})')
        # Set a filename for the current interaction data
        numbered_output_analysis_filepath = numerate_filename(output_analysis_filepath, i)

        # Add the root of the output analysis filename to the run data
        analysis_name = get_analysis_name(numbered_output_analysis_filepath)
        # Append current interaction to the summary
        analysis_entry = { 'name': name, 'analysis': analysis_name }
        output_summary.append(analysis_entry)

        # If the analysis already exists then proceed to the next interaction
        if exists(numbered_output_analysis_filepath): continue
        
        # Get interface atom indices
        interface_atom_indices_1 = interaction['interface_atom_indices_1']
        interface_atom_indices_2 = interaction['interface_atom_indices_2']
        # Select all interface residues in pytraj notation
        interface_atom_indices = interface_atom_indices_1 + interface_atom_indices_2
        if len(interface_atom_indices) == 0:
            raise ValueError(f'There are no interface atoms for interaction "{name}"')
        interface_selection = structure.select_atom_indices(interface_atom_indices)
        interface_pt_mask = interface_selection.to_pytraj()

        # Run the analysis
        hbonds = pt.hbond(pt_trajectory, mask=interface_pt_mask)

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

        # Save all atom indices for the final output
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
            if matchObj == None: continue
            # Mine all data from the parsed results
            acceptor_resnum = int(matchObj.group(1))
            acceptor_atom_name = matchObj.group(2)
            donor_resnum = int(matchObj.group(3))
            donor_atom_name = matchObj.group(4)
            hydrogen_atom_name = matchObj.group(5)
            # Get the acceptor and donor parsed atoms
            acceptor_residue = pytraj_residue_number_2_residue(acceptor_resnum)
            acceptor_atom = acceptor_residue.get_atom_by_name(acceptor_atom_name)
            donor_residue = pytraj_residue_number_2_residue(donor_resnum)
            donor_atom = donor_residue.get_atom_by_name(donor_atom_name)
            hydrogen_atom = donor_residue.get_atom_by_name(hydrogen_atom_name)
            # WARNING: The analysis may return hydrogen bonds between residues from the same agent
            # Accept the hydrogen bond only if its residues belong to different interaction agents
            if acceptor_atom.index in interface_atom_indices_1 and donor_atom.index in interface_atom_indices_1: continue
            if acceptor_atom.index in interface_atom_indices_2 and donor_atom.index in interface_atom_indices_2: continue
            # Get the absolute index of each atom
            acceptor_atom_index_list.append(acceptor_atom.index)
            donor_atom_index_list.append(donor_atom.index)
            hydrogen_atom_index_list.append(hydrogen_atom.index)
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
            if not is_time_dependent: continue
            # Now split the whole trajectory in time slots and calculate hbonds percents for every slot
            temporal_percents = []
            for i in range(nsteps):
                slot_hbond_values = hbond_values[i*step:(i+1)*step]
                temporal_percent = sum(slot_hbond_values) / step
                temporal_percents.append(temporal_percent)
            hbond_timed.append(temporal_percents)

        # If no hydrogen bonds were found then do not append the analysis
        if len(acceptor_atom_index_list) == 0:
            output_summary.pop()
            continue

        # Set the interaction output
        interaction_data = {
            'name': interaction['name'],
            'acceptors': acceptor_atom_index_list,
            'donors': donor_atom_index_list,
            'hydrogens': hydrogen_atom_index_list,
            'hbonds': hbond_overall,
            'version': '1.0.0'
        }
        if is_time_dependent:
            interaction_data['hbonds_timed'] = hbond_timed
        # DANI: Esto el cliente aún no lo soporta y me he quedado sin tiempo
        # DANI: De momento todo lo que sea ensemble que tire de mostar el overall y listo
        # if populations:
        #    interaction_data['hbonds_framed'] = hbond_framed

        # Write the interaction analysis output to a file
        save_json(interaction_data, numbered_output_analysis_filepath)

    # Export the analysis in json format
    if len(output_summary) > 0:
        save_json(output_summary, output_analysis_filepath)