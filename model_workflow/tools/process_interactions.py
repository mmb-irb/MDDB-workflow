# Processing the interactions means finding the residues of each interacting agent
# In addition, interface residues are listed appart

from subprocess import run, PIPE, Popen
import os
import json

from model_workflow.tools.topology_manager import sourceResidue
from model_workflow.tools.get_reduced_trajectory import get_reduced_trajectory
from model_workflow.tools.xvg_parse import xvg_parse

# Since interactions are heavy to calculate they are stored in a json file
# This file is only a backup, which is not uploaded to the database
# In addition, this file may be used to force interactions with custom interface residues manually
interactions_backup = 'interactions.json'

# The processing uses gromacs to find interface residues
# Residues are filtered by minimum distance along the trajecotry
def process_interactions (
    topology_filename : str,
    trajectory_filename : str,
    interactions : list,
    topology_reference,
    snapshots : int) -> list:

    # If there is a backup then use it
    # Load the backup and return its content as it is
    if os.path.exists(interactions_backup):
        with open(interactions_backup, 'r') as file:
            interactions = json.load(file)
        # Parse the residues in string format to 'sourceResidue' format (e.g. 'A:1')
        for interaction in interactions:
            for key in ['residues_1','residues_2','interface_1','interface_2']:
                interaction[key] = [ sourceResidue.from_tag(residue) for residue in interaction[key] ]
        return interactions

    if not interactions or len(interactions) == 0:
        return []
    
    for interaction in interactions:
        # residues_1 is the list of all residues in the first agent
        interaction['residues_1'] = topology_reference.topology_selection(
            interaction['selection_1']
        )
        # residues_2 is the list of all residues in the second agent
        interaction['residues_2'] = topology_reference.topology_selection(
            interaction['selection_2']
        )
        # Translate the residues selections to pytraj notation
        # These values are used along the workflow but not added to metadata
        interaction.update(
            {
                'pt_residues_1': list(map(topology_reference.source2pytraj, interaction['residues_1'])),
                'pt_residues_2': list(map(topology_reference.source2pytraj, interaction['residues_2'])),
            }
        )
        # Find out the interface resdiues for each agent
        # Update the interaction dict by adding both interfaces
        get_interface_residues(topology_filename, trajectory_filename, interaction, snapshots)
        
        # Translate the residues selections back from pytraj notation
        interaction['interface_1'] = list(map(topology_reference.pytraj2source, interaction['pt_interface_1']))
        interaction['interface_2'] = list(map(topology_reference.pytraj2source, interaction['pt_interface_2']))

        print(interaction['name'] +
            ' -> ' + str(interaction['interface_1'] + interaction['interface_2']))

    # Save a backup for interactions
    with open(interactions_backup, 'w') as file:
        # Create a new interactions object with all 'sourceResidue' values as string
        # e.g. 'A:1'
        serializable_interactions = []
        for interaction in interactions:
            serializable_interaction = {}
            for key, value in interaction.items():
                if type(value) == str:
                    serializable_interaction[key] = value
                else:
                    serializable_interaction[key] = [ str(element) for element in value ]
            serializable_interactions.append(serializable_interaction)
        json.dump(serializable_interactions, file, indent=4)

    return interactions

# There is a file which is always produced by mindist
# We remove it each time for the files to do not disturb us
residual_mindist_file = 'mindist.xvg'

# Find out the interface residues from an interaction and add them to the input dict
# Obtain minimum distances between each agents pair of residues and get the closer ones
# This is powered by the gromacs command 'mindist' and the distance cutoff is 5 Angstroms
# WARNING: The cutoff distance is in nanometers (nm)
def get_interface_residues(
    topology_filename : str,
    trajectory_filename : str,
    interaction : dict,
    snapshots : int,
    cutoff_distance : float = 0.5) -> dict:

    # Use a reduced trajectory since this process is slow
    frames_limit = 100
    reduced_trajectory_filename, step, frames = get_reduced_trajectory(
        topology_filename,
        trajectory_filename,
        snapshots,
        frames_limit,
    )

    # First of all we must set up gromacs index to select each interacting agent separately
    # This is performed using the make_ndx gromacs command with specific commands

    # In order to remove all default groups we run first 'keep 0' and then 'del 0'
    # i.e. first remove all groups but the group 0 and then remove the group 0
    remove_default_groups = 'keep 0' + '\n' + 'del 0' + '\n'

    # Create a new group with the agent 1 residues and name
    residues_1 = interaction['pt_residues_1']
    name_1 = interaction['agent_1'].replace(' ','_') # Avoid white spaces
    string_residues_1 = [ str(residue) for residue in residues_1 ]
    select_1 = 'ri ' + ' '.join(string_residues_1) + '\n'
    rename_1 = 'name 0 ' + name_1 + '\n'

    # Create a new group with the agent 2 residues and name
    residues_2 = interaction['pt_residues_2']
    name_2 = interaction['agent_2'].replace(' ','_') # Avoid white spaces
    string_residues_2 = [ str(residue) for residue in residues_2 ]
    select_2 = 'ri ' + ' '.join(string_residues_2) + '\n'
    rename_2 = 'name 1 ' + name_2 + '\n'

    # Join all commands with the 'save and exit' command at the end
    all_commands = remove_default_groups + select_1 + rename_1 + select_2 + rename_2 + 'q'

    # Run gromacs make_ndx to create the index file
    index = 'index.ndx'
    p = Popen([
        "echo",
        all_commands,
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "make_ndx",
        "-f",
        topology_filename,
        '-o',
        index,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

    # Now calculate the minimum distance for each residue in each agent using the mindist gromacs command

    # Run gromacs mindist to create the minimum distances per residue file
    mindistres_filename = 'mindistres.xvg'
    p = Popen([
        "echo",
        "0",
        "1",
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "mindist",
        "-s",
        topology_filename,
        "-f",
        reduced_trajectory_filename,
        "-n",
        index,
        '-or',
        mindistres_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

    # Mine the distances from the gromacs output file
    data_1 = xvg_parse(mindistres_filename, ['residue', 'distance'])

    # Set which residues are considered as interface residues according to distances
    interface_1 = []
    distances = data_1['distance']
    for i, distance in enumerate(distances):
        if distance <= cutoff_distance:
            interface_1.append(residues_1[i])

    # Then delete the output files
    run([
        "rm",
        mindistres_filename,
        residual_mindist_file,
    ], stdout=PIPE).stdout.decode()

    # Repeat the process with the second agent

    # Run gromacs mindist to create the minimum distances per residue file
    mindistres_filename = 'mindistres.xvg'
    p = Popen([
        "echo",
        "1",
        "0",
    ], stdout=PIPE)
    logs = run([
        "gmx",
        "mindist",
        "-s",
        topology_filename,
        "-f",
        reduced_trajectory_filename,
        "-n",
        index,
        '-or',
        mindistres_filename,
        '-quiet'
    ], stdin=p.stdout, stdout=PIPE).stdout.decode()
    p.stdout.close()

    # Mine the distances from the gromacs output file
    data_2 = xvg_parse(mindistres_filename, ['residue', 'distance'])

    # Set which residues are considered as interface residues according to distances
    interface_2 = []
    distances = data_2['distance']
    for i, distance in enumerate(distances):
        if distance <= cutoff_distance:
            interface_2.append(residues_2[i])

    # Then delete the output files
    run([
        "rm",
        mindistres_filename,
        residual_mindist_file,
    ], stdout=PIPE).stdout.decode()

    # Update the interaction dict with both interfaces
    interaction['pt_interface_1'] = interface_1
    interaction['pt_interface_2'] = interface_2

    # Remove the index file
    run([
        "rm",
        index,
    ], stdout=PIPE).stdout.decode()

    return interaction

# The easy processing uses the cutoff in the topology / structure to set the interface residues
# Interface residues are defined as by a cuttoff distance in Angstroms
# DEPRECATED
def easy_process_interactions (
    interactions : list,
    topology_reference,
    cutoff_distance : float = 5) -> list:

    if not interactions or len(interactions) == 0:
        return []
    
    for interaction in interactions:
        # residues_1 is the list of all residues in the first agent
        interaction['residues_1'] = topology_reference.topology_selection(
            interaction['selection_1']
        )
        # residues_2 is the list of all residues in the second agent
        interaction['residues_2'] = topology_reference.topology_selection(
            interaction['selection_2']
        )
        # interface_1 is the list of residues from the agent 1 which are close to the agent 2
        interaction['interface_1'] = topology_reference.topology_selection(
            '(' + interaction['selection_1'] +
            ') and same residue as exwithin ' +
            str(cutoff_distance) +
            ' of (' +
            interaction['selection_2'] + ')')
        # interface_2 is the list of residues from agent 2 which are close to the agent 1
        interaction['interface_2'] = topology_reference.topology_selection(
            '(' + interaction['selection_2'] +
            ') and same residue as exwithin ' +
            str(cutoff_distance) +
            ' of (' +
            interaction['selection_1'] + ')')

        # Translate all residues selections to pytraj notation
        # These values are used along the workflow but not added to metadata
        interaction.update(
            {
                'pt_residues_1': list(map(topology_reference.source2pytraj, interaction['residues_1'])),
                'pt_residues_2': list(map(topology_reference.source2pytraj, interaction['residues_2'])),
                'pt_interface_1': list(map(topology_reference.source2pytraj, interaction['interface_1'])),
                'pt_interface_2': list(map(topology_reference.source2pytraj, interaction['interface_2'])),
            }
        )

        print(interaction['name'] +
            ' -> ' + str(interaction['interface_1'] + interaction['interface_2']))

    return interactions