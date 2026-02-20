import itertools
import shlex

from mddb_workflow.tools.get_reduced_trajectory import get_reduced_trajectory
from mddb_workflow.utils.auxiliar import InputError, TestFailure, save_json, warn, reprint
from mddb_workflow.utils.constants import STABLE_INTERACTIONS_FLAG, OUTPUT_INTERACTIONS_FILENAME
from mddb_workflow.utils.type_hints import *
from mddb_workflow.utils.vmd_spells import get_covalent_bonds_between, get_interface_atom_indices

# Set the default distance cutoff in Ångstroms (Å)
# This is useful for atomistic simulations
DEFAULT_DISTANCE_CUTOFF = 5

# Set which fields are to be uploaded to the database only
# i.e. not vmd selections or residue numbers
UPLOAD_FIELDS = {
    'name',
    'agent_1',
    'agent_2',
    'atom_indices_1',
    'atom_indices_2',
    'interface_atom_indices_1',
    'interface_atom_indices_2',
    'version',
    'strong_bonds',
    'has_cg'
}

# Set the flag used to label failed interactions
FAILED_INTERACTION_FLAG = 'failed'


# Find interfaces by computing a minimum distance between residues along the trajectory
# Residues are filtered by minimum distance along the trajectory
# Interface residues are listed apart
# The heavy results of interactions are stored in a json file which is uploaded to the database independently
# In addition, this file may be used to force interactions with custom interface residues manually
def process_interactions(
    input_interactions: Optional[list],
    structure_file: 'File',
    trajectory_file: 'File',
    structure: 'Structure',
    inchikey_map: list[dict],
    snapshots: int,
    output_directory: str,
    mercy: list[str],
    interactions_auto: str,
    pbc_selection: 'Selection',
    cg_selection: 'Selection',
    frames_limit: int = 1000,
    # Percent of frames where an interaction must have place (from 0 to 1)
    # If the interactions fails to pass the cutoff then the workflow is killed and the user is warned
    interaction_cutoff: float = 0.1,
    ) -> list:
    """Find the residues of each interacting agent.
    It can automatically detect interactions based on chain names or ligand information, or use a predefined list of interactions.
    """
    # Set our own internal interactions
    interactions = []

    # If there are no interactions then stop here
    if input_interactions is None: return []
    # If interactions is not a list then make it a list of it
    # if it is already a list then copy it to avoid mutating the original
    if type(input_interactions) == list:
        interactions = [*input_interactions]
    else:
        interactions = [input_interactions]
    # If interactions is an empty list then stop here
    if len(input_interactions) == 0: return []

    # Get explicit interactions (not keywords)
    # Duplicate input interactions to avoid modifying the originals
    explicit_interactions = [
        {k: v for k, v in inter.items()} for inter in interactions if type(inter) == dict]

    # Since this is the first time we read the input interactions, we must make sure they are correct
    # Make sure there are no interactions with the same name
    interaction_names = [interaction['name'] for interaction in explicit_interactions]
    if len(set(interaction_names)) < len(interaction_names):
        raise InputError('Interactions must have unique names')
    # Check input interactions to have every expected field and make sure they are coherent
    for interaction in explicit_interactions:
        name = interaction["name"]
        # Check agents have different names
        if interaction['agent_1'] == interaction['agent_2']:
            raise InputError(f'Interaction agents must have different names at {name}')
        # Check agents have different selections
        input_agent_1_selection = interaction['selection_1']
        input_agent_2_selection = interaction['selection_2']
        if input_agent_1_selection == input_agent_2_selection:
            raise InputError(f'Interaction agents must have different selections at {name}')
        # Make sure both agents have valid selections
        agent_1_selection = structure.select(input_agent_1_selection)
        if not agent_1_selection:
            raise InputError(f'Interaction "{name}" has a non valid (or empty) selection for agent 1 ({interaction["agent_1"]}): {interaction["selection_1"]}')
        agent_2_selection = structure.select(input_agent_2_selection)
        if not agent_2_selection:
            raise InputError(f'Interaction "{name}" has a non valid (or empty) selection for agent 2 ({interaction["agent_2"]}): {interaction["selection_2"]}')
        # Make sure selections do not overlap at all
        # This makes not sense as interactions are implemented in this workflow
        overlap = agent_1_selection & agent_2_selection
        if overlap: raise InputError(
            f'Agents in interaction "{name}" have {len(overlap)} overlapping atoms.\n' +
            f'This means that the atom selection in agent 1 ({input_agent_1_selection})' +
            f' and the atom selection in agent 2 ({input_agent_2_selection}) aim for' +
            f' {len(overlap)} atoms in common.' +
             ' This is not supported since it makes no sense in this context ' +
             'to check the interaction of a group of atoms against themselves.\n'
             'Please change the atom selections for this interaction.')

    # Get keywords from the input interactions
    # Right now the only supported instruction is 'auto'
    keyword_interactions = set([inter for inter in interactions if type(inter) == str])

    # If the iauto argument is used from the console then add it here
    if interactions_auto:
        keyword_interactions.add(f'auto "{interactions_auto}"')

    # Now process keyword interactions
    for keyword in keyword_interactions:
        # Parse the keyword
        arguments = shlex.split(keyword)
        # Depending on the header we set the behaviour
        header = arguments[0]
        options = arguments[1:]
        # Find all possible interactions in a specific region
        # All the structure is used by default if no selection is passed
        # Note that an interaction is set by every possible combination of pairs of chains
        if header == 'auto':
            print(' Processing interactions automatically')
            if len(options) > 1: raise InputError('Automatic interactions support one selection only.\n'+
                f' Your instruction was: auto {" ".join(options)}\n'+
                f' Did you forget to add the quotes maybe? Try this: auto "{" ".join(options)}"')
            # Get the structure regions which are not in PBC
            reference_structure = structure.filter(pbc_selection) if pbc_selection else structure
            # Set the selection of atoms where this logic is to be applied, if any
            # If a selection is passed then filter the structure before finding its chains
            if len(options) == 1:
                selection = options[0]
                parsed_selection = structure.select(selection)
                if not parsed_selection:
                    raise InputError(f'Selection for automatic interaction "{selection}" is empty')
                reference_structure = structure.filter(parsed_selection)
            # Iterate possible combinations of chains
            for chain1, chain2 in itertools.combinations(reference_structure.chains, 2):
                interaction = {
                    "name": f"chain {chain1.name}-chain {chain2.name} interaction",
                    "agent_1": f"chain {chain1.name}",
                    "agent_2": f"chain {chain2.name}",
                    "selection_1": f"chain {chain1.name}",
                    "selection_2": f"chain {chain2.name}",
                    "auto": True,
                }
                # Now add the explicit interaction to the list
                explicit_interactions.append(interaction)
        # If we do not recognize the keyword header then we complain
        else:
            raise InputError(f'Instruction "{header}" not supported as interactions input')

    # If there are no interactions at this point then return an empty list
    interaction_count = len(explicit_interactions)
    if interaction_count == 0: return []

    # Set the output filepath
    output_analysis_filepath = f'{output_directory}/{OUTPUT_INTERACTIONS_FILENAME}'

    # If trajectory frames number is bigger than the limit we create a reduced trajectory
    reduced_trajectory_filepath, step, frames = get_reduced_trajectory(
        structure_file,
        trajectory_file,
        snapshots,
        frames_limit,
    )

    # Set if we are to have mercy when an interaction fails
    have_mercy = STABLE_INTERACTIONS_FLAG in mercy

    # Iterate over each defined interaction
    for interaction in explicit_interactions:
        interaction_name = interaction["name"]
        # Set the distance cutoff
        distance_cutoff = interaction.get('distance_cutoff', DEFAULT_DISTANCE_CUTOFF)
        # Find if this interaction has coarse grain atoms involved
        agent_1_selection = structure.select(interaction['selection_1'])
        agent_2_selection = structure.select(interaction['selection_2'])
        has_agent_1_cg = bool(agent_1_selection & cg_selection)
        has_agent_2_cg = bool(agent_2_selection & cg_selection)
        interaction['has_cg'] = has_agent_1_cg or has_agent_2_cg
        # Check if we are using the defualt atomistic distance while selections are coarse grain
        # If this is the case then warn the user
        if interaction['has_cg'] and distance_cutoff == DEFAULT_DISTANCE_CUTOFF:
            warn(f'Using atomistic default distance cutoff ({distance_cutoff}Å) with coarse grain agent(s)\n'
            f'  You may need to manually specify the distance cutoff in the inputs file for interaction "{interaction_name}"')
        # Find out the interaction residues for each frame and save all residues as the overall interface
        interface_results = get_interface_atom_indices(
            structure_file.path,
            reduced_trajectory_filepath,
            interaction['selection_1'],
            interaction['selection_2'],
            distance_cutoff
        )
        # Check if the interaction is respecting the frames percent cutoff and if it fails then kill it
        frames_percent = interface_results['interacting_frames'] / interface_results['total_frames']
        pretty_frames_percent = str(round(frames_percent * 10000) / 100)
        if frames_percent < interaction_cutoff:
            # If this interaction was set automatically then simply skip it silently
            if interaction.get('auto', False):
                interaction[FAILED_INTERACTION_FLAG] = True
                continue
            meaning_log = 'is not happening at all' if frames_percent == 0 else 'is happening only in a small percent of the trajectory'
            print(f'Interaction "{interaction_name}" is not reaching the frames percent cutoff of {interaction_cutoff} ({pretty_frames_percent}).\n'
                f'This means the interaction {meaning_log}.\n'
                'Check agent selections are correct or consider removing this interaction from the inputs.\n'
                f'   - Agent 1 selection: {interaction["selection_1"]}\n'
                f'   - Agent 2 selection: {interaction["selection_2"]}')
            # If we are not to have mercy in case of interaction failure then stop here
            if not have_mercy: raise TestFailure('An interaction failed to be set.\n'
                'Use the "--mercy interact" flag for the workflow to continue.\n'
                'Failed interactions will be ignored and will not appear in further analyses.')
            # If the workflow is not to be killed then just remove this interaction from the interactions list
            # Thus it will not be considered in interaction analyses and it will not appear in the metadata
            # To not mess the interactions iteration we simply flag the interaction
            # It will be removed further
            warn(f'Interaction "{interaction_name}" will be ignored and will not appear in further analyses')
            interaction[FAILED_INTERACTION_FLAG] = True
            continue

        # Iterate interaction agents
        for agent in ['1', '2']:
            # Get agent name and selection for logging purposes
            agent_name = interaction['agent_' + agent]
            agent_selection = interaction['selection_' + agent]
            # Save atom indices in the interaction object
            atom_indices = interface_results[f'selection_{agent}_atom_indices']
            # This should never happen, but make sure they are not empty
            if len(atom_indices) == 0:
                raise ValueError(f'Empty agent "{agent_name}" in interaction "{interaction_name}": {agent_selection}')
            interaction[f'atom_indices_{agent}'] = atom_indices
            # Save interface atom indices in the interaction object
            interface_atom_indices = interface_results[f'selection_{agent}_interface_atom_indices']
            # This should never happen, but make sure they are not empty
            if len(interface_atom_indices) == 0:
                raise ValueError(f'Empty interface for agent "{agent_name}" in interaction "{interaction_name}": {agent_selection}')
            interaction[f'interface_atom_indices_{agent}'] = interface_atom_indices

        # Add residue notations
        add_residues_indices(interaction, structure)

        # Find strong bonds between residues in different interfaces
        # Use the main topology, which is corrected and thus will retrieve the right bonds
        strong_bonds = get_covalent_bonds_between(structure_file.path, interaction['selection_1'], interaction['selection_2'])
        interaction['strong_bonds'] = strong_bonds

        # If one of the agents is fully made of strong bonds then the interaction is not valid
        strong_bonded_atoms = set(sum(strong_bonds, []))
        if all(index in strong_bonded_atoms for index in agent_1_selection.atom_indices) or \
           all(index in strong_bonded_atoms for index in agent_2_selection.atom_indices):
            warn(f'Interaction "{interaction_name}" is not valid since one of the agents is fully bonded to the other agent.\n'
                 'This may be due to wrong selections or a wrong interaction to be considered.\n'
                 'This interaction will be ignored and will not appear in further analyses')
            interaction[FAILED_INTERACTION_FLAG] = True
            continue

        # Save the interactions version
        interaction['version'] = '2.0.0'

        # Log the final results
        interface_residue_indices = sorted(interaction["interface_residue_indices_1"]
            + interaction["interface_residue_indices_2"])
        print(f'{interaction_name} (time: {pretty_frames_percent} %) -> {interface_residue_indices}')

    # Filter away interactions which have failed
    valid_interactions = [inte for inte in explicit_interactions if not inte.get(FAILED_INTERACTION_FLAG, False)]
    interaction_count = len(valid_interactions)
    print(f'There is a total of {interaction_count} valid interactions')

    # If there are no valid interactions then stop here
    if interaction_count == 0: return []

    # Print an empty line for the next reprint
    print()
    ligands_selection = structure.select_ligands(inchikey_map)
    # Check input interactions to be correct
    for i, interaction in enumerate(valid_interactions, 1):
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
        if ligands_selection and (agent_1_selection & ligands_selection):
            agent_1_classification = 'ligand'
        if ligands_selection and (agent_2_selection & ligands_selection):
            agent_2_classification = 'ligand'
        alphabetically_sorted = sorted([agent_1_classification, agent_2_classification])
        interaction['type'] = f'{alphabetically_sorted[0]}-{alphabetically_sorted[1]}'
        # Hardcode some interaction types with its common names
        if alphabetically_sorted == ['ligand', 'protein']:
            interaction['type'] = 'protein-ligand'

    # Create interaction duplicates to avoid mutating the already processed interactions
    # Then fill these duplicates only with those fields to be uploaded to the database
    file_interactions = []
    for interaction in valid_interactions:
        file_interaction = {key: value for key, value in interaction.items() if key in UPLOAD_FIELDS}
        file_interactions.append(file_interaction)

    # Write them to disk
    save_json(file_interactions, output_analysis_filepath, indent=4)

    # Finally return the processed interactions
    return valid_interactions


def add_residues_indices(interaction: dict, structure: 'Structure'):
    """Add residue indices to an interactions object."""
    # Iterate interaction agents
    for agent in ['1', '2']:
        # Get interaction atom indices
        atom_indices = interaction[f'atom_indices_{agent}']
        # Now parse atom indices to residue indices for those analysis which work with residues
        residue_indices = sorted(list(set([structure.atoms[atom_index].residue_index for atom_index in atom_indices])))
        interaction[f'residue_indices_{agent}'] = residue_indices
        # Get interaction interface atom indices
        interface_atom_indices = interaction[f'interface_atom_indices_{agent}']
        # Then with interface atoms/residues
        interface_residue_indices = sorted(list(set([structure.atoms[atom_index].residue_index for atom_index in interface_atom_indices])))
        interaction[f'interface_residue_indices_{agent}'] = interface_residue_indices
