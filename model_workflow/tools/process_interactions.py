# Processing the interactions means finding the residues of each interacting agent
# In addition, interface residues are listed appart

# Interface residues are defined as by a cuttoff distance in Angstroms
cutoff_distance = 5

# The easy processing uses the cutoff in the topology / structure to set the interface residues
def process_interactions (
    interactions : list,
    topology_reference) -> list:

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