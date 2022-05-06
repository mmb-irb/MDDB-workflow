# Processing the interactions means finding the residues of each interacting agent
# In addition, interface residues are listed appart

from subprocess import run, PIPE, Popen
import os
import json

from model_workflow.tools.topology_manager import TopologyReference, sourceResidue
from model_workflow.tools.get_reduced_trajectory import get_reduced_trajectory
from model_workflow.tools.get_pdb_frames import get_pdb_frames
from model_workflow.tools.get_last_frame import get_last_frame
from model_workflow.tools.xvg_parse import xvg_parse

# The cutoff distance is in Ångstroms (Å)
cutoff_distance : float = 5

# Find interfaces by computing a minimum distance between residues along the trajectory
# Residues are filtered by minimum distance along the trajectory
# The heavy results of interactions are stored in a json file which is uploaded to the database independently
# This file is also used as a backup here, since calculating interactions is a heavy calculation
# In addition, this file may be used to force interactions with custom interface residues manually
def process_interactions (
    interactions : list,
    topology_filename : str,
    trajectory_filename : str,
    topology_reference,
    interactions_file : str) -> list:

    # If there is a backup then use it
    # Load the backup and return its content as it is
    if os.path.exists(interactions_file):
        loaded_interactions = load_interactions(interactions_file, topology_reference)
        # Merge the loaded interactions with the input interactions to cover all fields
        complete_interactions = []
        for i, input_interaction in enumerate(interactions):
            loaded_interaction = loaded_interactions[i]
            complete_interaction = { **input_interaction, **loaded_interaction }
            complete_interactions.append(complete_interaction)
        return complete_interactions

    # If there are no interactions return an empty list
    if not interactions or len(interactions) == 0:
        return []
    
    for interaction in interactions:
        selection_1 = interaction['selection_1']
        selection_2 = interaction['selection_2']
        # Get residues and residue indices
        # residues_1 is the list of all residues in the first agent
        residues_1, residue_indices_1 = topology_reference.prody_selection_2_residues(selection_1)
        interaction['residues_1'] = residues_1
        residue_indices_1 = [ int(i) for i in residue_indices_1 ]
        interaction['residue_indices_1'] = residue_indices_1
        # residues_2 is the list of all residues in the second agent
        residues_2, residue_indices_2 = topology_reference.prody_selection_2_residues(selection_2)
        interaction['residues_2'] = residues_2
        residue_indices_2 = [ int(i) for i in residue_indices_2 ]
        interaction['residue_indices_2'] = residue_indices_2
        # Get 100 frames (or less) along the trajectory
        frames_limit = 100
        frames, step, count = get_pdb_frames(topology_filename, trajectory_filename, frames_limit)
        # Find out the interaction residues for each frame and save all residues as the overall interface
        interface_1_residues = []
        interface_2_residues = []
        interface_1_indices = []
        interface_2_indices = []
        for current_frame in frames:
            frame_structure = TopologyReference(current_frame)
            # interface_1 is the list of residues from the agent 1 which are close to the agent 2
            interface_selection_1 = ('(' + selection_1 + ') and same residue as exwithin ' +
                str(cutoff_distance) + ' of (' + selection_2 + ')')
            frame_interface_residues_1, frame_interface_residue_indices_1 = frame_structure.prody_selection_2_residues(
                interface_selection_1
            )
            interface_1_residues += frame_interface_residues_1
            interface_1_indices += frame_interface_residue_indices_1
            # interface_2 is the list of residues from agent 2 which are close to the agent 1
            interface_selection_2 = ('(' + selection_2 + ') and same residue as exwithin ' +
                str(cutoff_distance) + ' of (' + selection_1 + ')')
            frame_interface_residues_2, frame_interface_residue_indices_2 = frame_structure.prody_selection_2_residues(
                interface_selection_2
            )
            interface_2_residues += frame_interface_residues_2
            interface_2_indices += frame_interface_residue_indices_2
        
        # Remove duplicates and sort residues
        interaction['interface_1'] = sorted(list(set(interface_1_residues)))
        interaction['interface_2'] = sorted(list(set(interface_2_residues)))
        # Conver indices in normal integers
        interface_1_indices = [ int(i) for i in interface_1_indices ]
        interface_2_indices = [ int(i) for i in interface_2_indices ]
        interaction['interface_indices_1'] = sorted(list(set(interface_1_indices)))
        interaction['interface_indices_2'] = sorted(list(set(interface_2_indices)))

        # Find strong bonds between residues in different interfaces
        # Use the last trajectory frame to find them
        last_frame_filename = 'last_frame.pdb'
        get_last_frame(topology_filename, trajectory_filename, last_frame_filename)
        vmd_atom_selection_1 = topology_reference.get_prody_selection(interaction['selection_1']).to_vmd()
        vmd_atom_selection_2 = topology_reference.get_prody_selection(interaction['selection_2']).to_vmd()
        strong_bonds = get_strong_bonds(last_frame_filename, vmd_atom_selection_1, vmd_atom_selection_2)

        # Translate all residues selections to pytraj notation
        # These values are used along the workflow but not added to metadata
        interaction.update(
            {
                'pt_residues_1': list(map(topology_reference.source2pytraj, interaction['residues_1'])),
                'pt_residues_2': list(map(topology_reference.source2pytraj, interaction['residues_2'])),
                'pt_interface_1': list(map(topology_reference.source2pytraj, interaction['interface_1'])),
                'pt_interface_2': list(map(topology_reference.source2pytraj, interaction['interface_2'])),
                'strong_bonds': strong_bonds
            }
        )

        print(interaction['name'] +
            ' -> ' + str(interaction['interface_1'] + interaction['interface_2']))

    # Write the interactions file with the fields to be uploaded to the database only
    # i.e. strong bonds and residue indices
    file_keys = [
        'name',
        'agent_1',
        'agent_2',
        'residue_indices_1',
        'residue_indices_2',
        'interface_indices_1',
        'interface_indices_2',
        'strong_bonds'
    ]
    with open(interactions_file, 'w') as file:
        # Create a new interactions object with all 'sourceResidue' values as string
        # e.g. 'A:1'
        file_interactions = []
        for interaction in interactions:
            file_interaction = { key: value for key, value in interaction.items() if key in file_keys }
            file_interactions.append(file_interaction)
        json.dump(file_interactions, file, indent=4)

    return interactions

# Load interactions from an already existinf interactions file
def load_interactions (interactions_file : str, topology_reference) -> list:
    with open(interactions_file, 'r') as file:
        # The stored interactions should carry only residue indices and strong bonds
        interactions = json.load(file)
        # Now we must complete every interactions dict by adding residues in source format and pytraj format
        for interaction in interactions:
            # Get residues from their indices
            converter = topology_reference.index_2_prody
            interaction['residues_1'] = [ converter(index) for index in interaction['residue_indices_1'] ]
            interaction['residues_2'] = [ converter(index) for index in interaction['residue_indices_2'] ]
            interaction['interface_1'] = [ converter(index) for index in interaction['interface_indices_1'] ]
            interaction['interface_2'] = [ converter(index) for index in interaction['interface_indices_2'] ]
            # Transform to pytraj
            converter = topology_reference.source2pytraj
            interaction['pt_residues_1'] = list(map(converter, interaction['residues_1']))
            interaction['pt_residues_2'] = list(map(converter, interaction['residues_2']))
            interaction['pt_interface_1'] = list(map(converter, interaction['interface_1']))
            interaction['pt_interface_2'] = list(map(converter, interaction['interface_2']))
        return interactions

# The easy processing uses the cutoff in the topology / structure to set the interface residues
# Interface residues are defined as by a cuttoff distance in Angstroms
# DEPRECATED
# DANI: El nuevo sistema hace lo mismo, pero en lugar de solo en la primera frame lo hace en varias
# DANI: De manera que al final, cada residuo que ha salido en una frame como mínimo se considera apto
def easy_process_interactions (
    interactions : list,
    topology_reference,
    cutoff_distance : float = 5) -> list:

    if not interactions or len(interactions) == 0:
        return []
    
    for interaction in interactions:
        # residues_1 is the list of all residues in the first agent
        interaction['residues_1'] = topology_reference.residues_selection(
            interaction['selection_1']
        )
        # residues_2 is the list of all residues in the second agent
        interaction['residues_2'] = topology_reference.residues_selection(
            interaction['selection_2']
        )
        # interface_1 is the list of residues from the agent 1 which are close to the agent 2
        interaction['interface_1'] = topology_reference.residues_selection(
            '(' + interaction['selection_1'] +
            ') and same residue as exwithin ' +
            str(cutoff_distance) +
            ' of (' +
            interaction['selection_2'] + ')')
        # interface_2 is the list of residues from agent 2 which are close to the agent 1
        interaction['interface_2'] = topology_reference.residues_selection(
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

# The processing uses gromacs to find interface residues
# Residues are filtered by minimum distance along the trajecotry
# DEPRECATED
# DANI: El problema de este sistema es que puede haber diferencias en la numeración de residuos
# DANI: Gromacs puede contar residuos diferente que pytraj y prody (e.g. hidrógenos al final)
# DANI: Pytraj y prody también pueden contar diferente, pero ya tengo funciones para parsear
# DANI: Hacer estas funciones para parsear de gromacs a prody o pytraj no es banal
# DANI: Esto se debe a que gromacs es mucho menos accesible que prody o pytraj
def process_interactions_gromacs (
    topology_filename : str,
    trajectory_filename : str,
    interactions : list,
    topology_reference) -> list:

    # If there is a backup then use it
    # Load the backup and return its content as it is
    if os.path.exists(interactions_file):
        with open(interactions_file, 'r') as file:
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
        interaction['residues_1'] = topology_reference.residues_selection(
            interaction['selection_1']
        )
        # residues_2 is the list of all residues in the second agent
        interaction['residues_2'] = topology_reference.residues_selection(
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
        get_interface_residues_gromacs(topology_filename, trajectory_filename, interaction)
        
        # Translate the residues selections back from pytraj notation
        interaction['interface_1'] = list(map(topology_reference.pytraj2source, interaction['pt_interface_1']))
        interaction['interface_2'] = list(map(topology_reference.pytraj2source, interaction['pt_interface_2']))

        print(interaction['name'] +
            ' -> ' + str(interaction['interface_1'] + interaction['interface_2']))

    # Save a backup for interactions
    with open(interactions_file, 'w') as file:
        # Create a new interactions object with all 'sourceResidue' values as string
        # e.g. 'A:1'
        serializable_interactions = []
        for interaction in interactions:
            serializable_interaction = {}
            for key, value in interaction.items():
                if type(value) == sourceResidue:
                    serializable_interaction[key] = [ str(element) for element in value ]
                else:
                    serializable_interaction[key] = value
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
# DEPRECATED
def get_interface_residues_gromacs(
    topology_filename : str,
    trajectory_filename : str,
    interaction : dict,
    cutoff_distance : float = 0.5) -> dict:

    # Use a reduced trajectory since this process is slow
    frames_limit = 100
    reduced_trajectory_filename, step, frames = get_reduced_trajectory(
        topology_filename,
        trajectory_filename,
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


# Set a function to retrieve strong bonds between 2 atom selections
# Atom selections must be in VMD selection syntax
def get_strong_bonds (structure_filename : str, atom_selection_1 : str, atom_selection_2 : str) -> list:

    # Prepare a script for the VMD to automate the commands. This is Tcl lenguage
    commands_filename = 'commands.vmd'
    output_index_1_file = 'index1.text'
    output_index_2_file = 'index2.text'
    output_bonds_file = 'bonds.text'
    with open(commands_filename, "w") as file:
        # Select the specified atoms in selection 1
        file.write('set sel1 [atomselect top "' + atom_selection_1 + '"]\n')
        # Save all atom index in the selection
        file.write('set index1 [$sel1 list]\n')
        # Write those index to a file
        file.write('set indexfile1 [open ' + output_index_1_file + ' w]\n')
        file.write('puts $indexfile1 $index1\n')
        # Save all strong atoms in the selection
        file.write('set bonds [$sel1 getbonds]\n')
        # Write those bonds to a file
        file.write('set bondsfile [open ' + output_bonds_file + ' w]\n')
        file.write('puts $bondsfile $bonds\n')
        # Select the specified atoms in selection 2
        file.write('set sel2 [atomselect top "' + atom_selection_2 + '"]\n')
        # Save all atom index in the selection
        file.write('set index2 [$sel2 list]\n')
        # Write those index to a file
        file.write('set indexfile2 [open ' + output_index_2_file + ' w]\n')
        file.write('puts $indexfile2 $index2\n')
        file.write('exit\n')
        
    # Run VMD
    logs = run([
        "vmd",
        structure_filename,
        "-e",
        commands_filename,
        "-dispdev",
        "none"
    ], stdout=PIPE).stdout.decode()
    
    # Read the VMD output
    with open(output_index_1_file, 'r') as file:
        raw_index_1 = file.read()
    with open(output_bonds_file, 'r') as file:
        raw_bonds = file.read()
    with open(output_index_2_file, 'r') as file:
        raw_index_2 = file.read()

    # Remove vmd files since they are no longer usefull
    for f in [ commands_filename, output_index_1_file, output_index_2_file, output_bonds_file ]:
        os.remove(f)

    # Sometimes there is a breakline at the end of the raw bonds string and it must be removed
    # Add a space at the end of the string to make the parser get the last character
    raw_bonds = raw_bonds.replace('\n', '') + ' '
    
    # Raw indexes is a string with all indexes separated by spaces
    index_1 = [ int(i) for i in raw_index_1.split() ]
    index_2 = [ int(i) for i in raw_index_2.split() ]
    
    # Parse the raw bonds string to a list of atom bonds (i.e. a list of lists of integers)
    # Raw bonds format is (for each atom in the selection):
    # '{index1, index2, index3 ...}' with the index of each connected atom
    # 'index' if there is only one connected atom
    # '{}' if there are no connected atoms
    bonds_per_atom = []
    last_atom_index = ''
    last_atom_bonds = []
    in_brackets = False
    for character in raw_bonds:
        if character == ' ':
            if len(last_atom_index) > 0:
                if in_brackets:
                    last_atom_bonds.append(int(last_atom_index))
                else:
                    bonds_per_atom.append([int(last_atom_index)])
                last_atom_index = ''
            continue
        if character == '{':
            in_brackets = True
            continue
        if character == '}':
            if last_atom_index == '':
                bonds_per_atom.append([])
                in_brackets = False
                continue
            last_atom_bonds.append(int(last_atom_index))
            last_atom_index = ''
            bonds_per_atom.append(last_atom_bonds)
            last_atom_bonds = []
            in_brackets = False
            continue
        last_atom_index += character
        
    # At this point indexes and bonds from the first selection should match in number
    if len(index_1) != len(bonds_per_atom):
        raise ValueError('Indexes (' + str(len(index_1)) + ') and atom bonds (' +  str(len(bonds_per_atom)) + ') do not match in number')
        
    # Now get all strong bonds which include an index from the atom selection 2
    crossed_bonds = []
    for i, index in enumerate(index_1):
        bonds = bonds_per_atom[i]
        for bond in bonds:
            if bond in index_2:
                crossed_bond = (index, bond)
                crossed_bonds.append(crossed_bond)
                
    return crossed_bonds