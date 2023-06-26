# Processing the interactions means finding the residues of each interacting agent
# In addition, interface residues are listed appart

from subprocess import run, PIPE, Popen
import os
import json
from typing import List

from model_workflow.tools.get_reduced_trajectory import get_reduced_trajectory

# The cutoff distance is in Ångstroms (Å)
distance_cutoff : float = 5

# Find interfaces by computing a minimum distance between residues along the trajectory
# Residues are filtered by minimum distance along the trajectory
# The heavy results of interactions are stored in a json file which is uploaded to the database independently
# This file is also used as a backup here, since calculating interactions is a heavy calculation
# In addition, this file may be used to force interactions with custom interface residues manually
def process_interactions (
    input_interactions : list,
    topology_filename : str,
    trajectory_filename : str,
    structure : 'Structure',
    snapshots : int,
    interactions_file : str,
    frames_limit : int) -> list:

    # If there is a backup then use it
    # Load the backup and return its content as it is
    if os.path.exists(interactions_file):
        loaded_interactions = load_interactions(interactions_file, structure)
        # Merge the loaded interactions with the input interactions to cover all fields
        complete_interactions = []
        for input_interaction, loaded_interaction in zip(input_interactions, loaded_interactions):
            complete_interaction = { **input_interaction, **loaded_interaction }
            complete_interactions.append(complete_interaction)
        return complete_interactions

    # If there are no interactions return an empty list
    if not input_interactions or len(input_interactions) == 0:
        return []

    # If trajectory frames number is bigger than the limit we create a reduced trajectory
    reduced_trajectory, step, frames = get_reduced_trajectory(
        topology_filename,
        trajectory_filename,
        snapshots,
        frames_limit,
    )
    
    # Duplicate the input interactions to avoid modifying the originals
    interactions = [ { k:v for k,v in interaction.items() } for interaction in input_interactions ]

    # Iterate over each defined interaction
    for interaction in interactions:
        # Find out the interaction residues for each frame and save all residues as the overall interface
        interface_results = get_interface_atom_indices_vmd(
            topology_filename,
            reduced_trajectory,
            interaction['selection_1'],
            interaction['selection_2'],
            distance_cutoff
        )
        # For each agent in the interaction, get the residues in the interface from the previously calculated atom indices
        for agent in ['1','2']:
            # First with all atoms/residues
            atom_indices = interface_results['selection_' + agent + '_atom_indices']
            residue_indices = sorted(list(set([ structure.atoms[atom_index].residue_index for atom_index in atom_indices ])))
            # Check residue lists to not be empty, which should never happen
            if len(residue_indices) == 0:
                agent_name = interaction['agent_' + agent]
                raise ValueError('Empty selection for agent "' + agent_name + '" in interaction "' + interaction['name'] + '": ' + interaction['selection_' + agent])
            interaction['residue_indices_' + agent] = residue_indices
            interaction['residues_' + agent] = [ structure.residues[residue_index] for residue_index in residue_indices ]
            # Then with interface atoms/residues
            interface_atom_indices = interface_results['selection_' + agent + '_interface_atom_indices']
            interface_residue_indices = list(set([ structure.atoms[atom_index].residue_index for atom_index in interface_atom_indices ]))
            interaction['interface_indices_' + agent] = interface_residue_indices
            interaction['interface_' + agent] = [ structure.residues[residue_index] for residue_index in interface_residue_indices ]

        # Find strong bonds between residues in different interfaces
        # Use the main topology, which is corrected and thus will retrieve the right bonds
        strong_bonds = get_strong_bonds(topology_filename, interaction['selection_1'], interaction['selection_2'])

        # Translate all residues selections to pytraj notation
        # These values are used along the workflow but not added to metadata
        converter = structure.residue_2_pytraj_residue_index
        interaction.update(
            {
                'pt_residues_1': list(map(converter, interaction['residues_1'])),
                'pt_residues_2': list(map(converter, interaction['residues_2'])),
                'pt_interface_1': list(map(converter, interaction['interface_1'])),
                'pt_interface_2': list(map(converter, interaction['interface_2'])),
                'strong_bonds': strong_bonds
            }
        )

        print(interaction['name'] +
            ' -> ' + str(interaction['interface_indices_1'] + interaction['interface_indices_2']))

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
        file_interactions = []
        for interaction in interactions:
            file_interaction = { key: value for key, value in interaction.items() if key in file_keys }
            file_interactions.append(file_interaction)
        json.dump(file_interactions, file, indent=4)

    return interactions

# Load interactions from an already existinf interactions file
def load_interactions (interactions_file : str, structure : 'Structure') -> list:
    with open(interactions_file, 'r') as file:
        # The stored interactions should carry only residue indices and strong bonds
        interactions = json.load(file)
        # Now we must complete every interactions dict by adding residues in source format and pytraj format
        for interaction in interactions:
            # Get residues from their indices
            residues = structure.residues
            interaction['residues_1'] = [ residues[index] for index in interaction['residue_indices_1'] ]
            interaction['residues_2'] = [ residues[index] for index in interaction['residue_indices_2'] ]
            # Check residue lists to not be empty, which should never happen
            if len(interaction['residues_1']) == 0:
                raise ValueError('Empty selection for agent "' + interaction['agent_1'] + '" in interaction "' + interaction['name'] + '"')
            if len(interaction['residues_2']) == 0:
                raise ValueError('Empty selection for agent "' + interaction['agent_2'] + '" in interaction "' + interaction['name'] + '"')
            interaction['interface_1'] = [ residues[index] for index in interaction['interface_indices_1'] ]
            interaction['interface_2'] = [ residues[index] for index in interaction['interface_indices_2'] ]
            # Transform to pytraj
            converter = structure.residue_2_pytraj_residue_index
            interaction['pt_residues_1'] = list(map(converter, interaction['residues_1']))
            interaction['pt_residues_2'] = list(map(converter, interaction['residues_2']))
            interaction['pt_interface_1'] = list(map(converter, interaction['interface_1']))
            interaction['pt_interface_2'] = list(map(converter, interaction['interface_2']))
        return interactions

# Given two atom selections, find interface atoms and return their indices
# Interface atoms are those atoms closer than the cutoff in at least 1 frame along a trajectory
# Return also atom indices for the whole selections
def get_interface_atom_indices_vmd (
    input_structure_filename : str,
    input_trajectory_filename : str,
    selection_1 : str,
    selection_2 : str,
    distance_cutoff : float,
) -> List[int]:

    # Set the interface selections
    interface_selection_1 = ('(' + selection_1 + ') and within ' + str(distance_cutoff) + ' of (' + selection_2 + ')')
    interface_selection_2 = ('(' + selection_2 + ') and within ' + str(distance_cutoff) + ' of (' + selection_1 + ')')
    
    # Set the output txt files for vmd to write the atom indices
    selection_1_filename = '.selection_1.txt'
    selection_2_filename = '.selection_2.txt'
    interface_selection_1_filename = '.interface_selection_1.txt'
    interface_selection_2_filename = '.interface_selection_2.txt'

    # Prepare a script for VMD to run. This is Tcl language
    commands_filename = 'commands.vmd'
    with open(commands_filename, "w") as file:
        # -------------------------------------------
        # First get the whole selection atom indices
        # -------------------------------------------
        # Select the specified atoms
        file.write('set selection [atomselect top "' + selection_1 + '"]\n')
        # Save atom indices from the selection
        file.write('set indices [$selection list]\n')
        # Write atom indices to a file
        file.write('set indices_file [open ' + selection_1_filename + ' w]\n')
        file.write('puts $indices_file $indices\n')
        # Select the specified atoms
        file.write('set selection [atomselect top "' + selection_2 + '"]\n')
        # Save atom indices from the selection
        file.write('set indices [$selection list]\n')
        # Write atom indices to a file
        file.write('set indices_file [open ' + selection_2_filename + ' w]\n')
        file.write('puts $indices_file $indices\n')
        # -------------------------------------------
        # Now get the interface selection atom indices
        # -------------------------------------------
        # Capture indices for each frame in the trajectory
        file.write('set accumulated_interface1_atom_indices []\n')
        file.write('set accumulated_interface2_atom_indices []\n')
        file.write('set interface1 [atomselect top "' + interface_selection_1 + '"]\n')
        file.write('set interface2 [atomselect top "' + interface_selection_2 + '"]\n')
        # Get the number of frames in the trajectory
        file.write('set nframes [molinfo top get numframes]\n')
        # Iterate over each frame
        file.write('for { set i 1 } { $i < $nframes } { incr i } {\n')
        # Update the selection in the current frame
        file.write('    $interface1 frame $i\n')
        file.write('    $interface1 update\n')
        # Add its atom indices to the acumulated atom indices
        file.write('    set interface1_atom_indices [$interface1 list]\n')
        file.write('    set accumulated_interface1_atom_indices [concat $accumulated_interface1_atom_indices $interface1_atom_indices ]\n')
        # Repeat with the selection 2
        file.write('    $interface2 frame $i\n')
        file.write('    $interface2 update\n')
        file.write('    set interface2_atom_indices [$interface2 list]\n')
        file.write('    set accumulated_interface2_atom_indices [concat $accumulated_interface2_atom_indices $interface2_atom_indices ]\n')
        file.write('}\n')
        # Remove duplicated indices
        file.write('lsort -unique $accumulated_interface1_atom_indices\n')
        file.write('lsort -unique $accumulated_interface2_atom_indices\n')
        # Write indices to files
        file.write('set indices_file [open ' + interface_selection_1_filename + ' w]\n')
        file.write('puts $indices_file $accumulated_interface1_atom_indices\n')
        file.write('set indices_file [open ' + interface_selection_2_filename + ' w]\n')
        file.write('puts $indices_file $accumulated_interface2_atom_indices\n')
        file.write('exit\n')

    # Run VMD
    logs = run([
        "vmd",
        input_structure_filename,
        input_trajectory_filename,
        "-e",
        commands_filename,
        "-dispdev",
        "none"
    ], stdout=PIPE, stderr=PIPE).stdout.decode()

    # If any of the output files do not exist at this point then it means something went wrong with vmd
    expected_output_files = [
        selection_1_filename,
        selection_2_filename,
        interface_selection_1_filename,
        interface_selection_2_filename
    ]
    for output_file in expected_output_files:
        if not os.path.exists(output_file):
            print(logs)
            raise SystemExit('Something went wrong with VMD')
    
    # Set a function to read the VMD output and parse the atom indices string to an array of integers
    def process_vmd_output (output_filename : str) -> List[int]:
        with open(output_filename, 'r') as file:
            raw_atom_indices = file.read()
        return [ int(i) for i in raw_atom_indices.split() ]

    # Read the VMD output
    selection_1_atom_indices = process_vmd_output(selection_1_filename)
    selection_2_atom_indices = process_vmd_output(selection_2_filename)
    selection_1_interface_atom_indices = process_vmd_output(interface_selection_1_filename)
    selection_2_interface_atom_indices = process_vmd_output(interface_selection_2_filename)
    
    # Remove trash files
    trash_files = [ commands_filename ] + expected_output_files
    for trash_file in trash_files:
        os.remove(trash_file)

    # Return the results
    return {
        'selection_1_atom_indices': selection_1_atom_indices,
        'selection_2_atom_indices': selection_2_atom_indices,
        'selection_1_interface_atom_indices': selection_1_interface_atom_indices,
        'selection_2_interface_atom_indices': selection_2_interface_atom_indices
    }

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
    ], stdout=PIPE, stderr=PIPE).stdout.decode()
    
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