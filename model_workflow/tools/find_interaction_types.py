from model_workflow.utils.auxiliar import warn, reprint
from model_workflow.utils.type_hints import *

# Given a list of interactions, find their types
def find_interaction_types (
    input_interactions : Optional[list],
    structure : 'Structure',) -> dict:
    # If there are no interactions return an empty list
    interaction_count = len(input_interactions)
    if interaction_count == 0: return {}
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