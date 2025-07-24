from model_workflow.utils.auxiliar import warn, reprint, InputError
from model_workflow.utils.type_hints import *

def find_interaction_types (
    input_interactions : Optional[list],
    structure : 'Structure',) -> dict:
    """Given a list of interactions, find their types."""
    # If there are no interactions return an empty list
    # The interactions could be a null value
    if input_interactions is None: return {}
    interaction_count = len(input_interactions)
    if interaction_count == 0: return {}
    # Since this is the first time we read the input interactions, we must make sure they are correct
    # Make sure there are no interactions with the same name
    interaction_names = [ interaction['name'] for interaction in input_interactions ]
    if len(set(interaction_names)) < len(interaction_names):
        raise InputError('Interactions must have unique names')
    # Check input interactions to have every expected field and make sure they are coherent
    for i, interaction in enumerate(input_interactions, 1):
        name = interaction["name"]
        # Check agents have different names
        if interaction['agent_1'] == interaction['agent_2']:
            raise InputError(f'Interaction agents must have different names at {name}')
        # Check agents have different selections
        if interaction['selection_1'] == interaction['selection_2']:
            raise InputError(f'Interaction agents must have different selections at {name}')
        # Make sure both agents have valid selections
        agent_1_selection = structure.select(interaction['selection_1'])
        if not agent_1_selection:
            raise InputError(f'Interaction "{name}" has a non valid (or empty) selection for agent 1 ({interaction["agent_1"]}): {interaction["selection_1"]}')
        agent_2_selection = structure.select(interaction['selection_2'])
        if not agent_2_selection:
            raise InputError(f'Interaction "{name}" has a non valid (or empty) selection for agent 2 ({interaction["agent_2"]}): {interaction["selection_2"]}')
        # Make sure selections do not overlap at all
        # This makes not sense as interactions are implemented in this workflow
        overlap = agent_1_selection & agent_2_selection
        if overlap:
            raise InputError(f'Agents in interaction "{name}" have {len(overlap)} overlapping atoms')
    # Save in a dict the results
    interaction_types = {}
    # Print an empty line for the next reprint
    print()
    # Check input interactions to be correct
    for i, interaction in enumerate(input_interactions, 1):
        name = interaction["name"]
        agent_1_selection = structure.select(interaction['selection_1'])
        agent_2_selection = structure.select(interaction['selection_2'])
        reprint(f' Finding interaction type in {name} ({i}/{interaction_count})')
        # Check if there was a type already assigned to the interaction
        # This is not supported anymore since the interaction type is set automatically
        if 'type' in interaction:
            warn(f'Interaction type "{interaction["type"]}" is set for interaction "{name}".\n'
                 'Interaction type is now calculated and the input interaction type is no longer supported.\n'
                 'Note that the input value will be ignored')
        # Set the interaction type
        # LORE: The type was a user input back in time but now we find it automatically
        # WARNING: Do not calculate the type from the interface residue instead of the whole agent
        # WARNING: This seems more coherent BUT the type will be written in the PROJECT metadata
        # WARNING: Interaction type is a valuable search parameter and thus it must remain in project metadata
        # WARNING: However we could have different types in different MDs, if the interaction is different
        agent_1_classification = structure.get_selection_classification(agent_1_selection)
        agent_2_classification = structure.get_selection_classification(agent_2_selection)
        alphabetically_sorted = sorted([agent_1_classification, agent_2_classification])
        interaction_types[name] = f'{alphabetically_sorted[0]}-{alphabetically_sorted[1]}'
    # Return the interaction types
    return interaction_types