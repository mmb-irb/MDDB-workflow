# Energies
#
# The energies analysis computes the electrostatic and van der Waals interaction energy of every atom of an agent with
# the other agent of each interaction, frame by frame. It uses CMIP, a locally developed tool.
# CMIP stands for Classical Molecular Interaction Potential and it is useful to predict electrostatic and Van der Waals
# potentials. Electrostatics are solved with the Poisson-Boltzmann equation on a grid. Atom charges come from the topology.
#
# The frame loop, the backup and the output files are independent of the method: an engine computes the energies of
# both agents in one frame (see EnergiesEngine). CmipEngine is the only engine for now.
#
# CMIP:
# Gelpí, J.L., Kalko, S.G., Barril, X., Cirera, J., de la Cruz, X., Luque, F.J. and Orozco, M. (2001), Classical molecular interaction potentials: Improved setup procedure in molecular dynamics simulations of proteins. Proteins, 45: 428-437. doi:10.1002/prot.1159

# Imports libraries
from os import mkdir, remove
from os.path import exists
from shutil import copyfile
from glob import glob
from subprocess import run, PIPE, CompletedProcess
from typing import NamedTuple, Protocol
import re
import math
import numpy as np
from scipy.spatial import cKDTree

from mddb_workflow.tools.get_pdb_frames import get_pdb_frames
from mddb_workflow.utils.auxiliar import (
    load_json,
    save_json,
    warn,
    numerate_filename,
    get_analysis_name,
)
from mddb_workflow.utils.auxiliar import ForcedStop
from mddb_workflow.utils.constants import (
    OUTPUT_ENERGIES_FILENAME,
    ENERGIES_ANALYSIS_VERSION,
)
from mddb_workflow.utils.constants import (
    PROTEIN_RESIDUE_NAME_LETTERS,
    NUCLEIC_RESIDUE_NAME_LETTERS,
)
from mddb_workflow.utils.constants import (
    CMIP_INPUTS_CHECKONLY_SOURCE,
    CMIP_INPUTS_SOURCE,
)
from mddb_workflow.utils.constants import CMIP_VDW_SOURCE, ENERGIES_DEBUG_SCRIPT_SOURCE
from mddb_workflow.utils.structures import Structure
from mddb_workflow.utils.file import File
from mddb_workflow.utils.type_hints import *

CURSOR_UP_ONE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'
ERASE_2_PREVIOUS_LINES = (
    CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE
)
ERASE_3_PREVIOUS_LINES = (
    CURSOR_UP_ONE
    + ERASE_LINE
    + CURSOR_UP_ONE
    + ERASE_LINE
    + CURSOR_UP_ONE
    + ERASE_LINE
    + CURSOR_UP_ONE
)
# Set the minimum distance (Å) between agents beyond which CMIP is not run and energies are assumed to be 0
# Otherwise the CMIP grid, which must include both agents, becomes huge and the calculation becomes too slow and memory demanding
FAR_AGENTS_CUTOFF = 30
# CMIP fails with agents of this number of atoms or more, so these interactions are skipped
# This limit can be changed in CMIP, but CMIP must be recompiled then, which is not trivial in a conda environment
CMIP_ATOM_LIMIT = 80000
# Set the limit number of grid points
# Otherwise CMIP could return the following error:
# ERROR when trying to allocate           -2092351636  bytes of memory
# It has been found experimentally: 65945880 Grid points -> OK, 68903412 Grid points -> ERROR
GRID_POINTS_LIMIT = 65000000
# Set the default grid 'cells' size (Å) for all three dimensions
GRID_UNIT_SIZE = 0.5
# Atom energies (vdw + es, kcal/mol) above this value are set to 0
# This step is performed to filter 'infinity' values
MAX_ATOM_ENERGY = 100
# Set an exception to handle internal issues
OUTPUT_MINING_EXCEPTION = Exception('CMIP output mining failed')
EXPECTED_CMIP_PDB_LINE_LENGTH = 65


class Agent(NamedTuple):
    """One of the two agents of an interaction, in a frame."""

    # Agent name, with no spaces or slashes, so it can be used in filenames
    name: str
    # Atom indices of the agent in the structure
    atom_indices: list[int]
    # Selection of the agent atoms in the frame structure
    selection: 'Selection'


class EnergiesEngine(Protocol):
    """A method to calculate the interaction energies of both agents of an interaction in one frame."""

    def __call__(
        self,
        frame_structure: 'Structure',
        agent1: Agent,
        agent2: Agent,
        strong_bonds: Optional[list],
    ) -> tuple[list, list]:
        """Get the energies of every atom of each agent with the other agent.

        For each agent, return one (vdw, es, both) tuple per agent atom, in the order of the agent atom indices.
        Atoms whose energies are not calculated (atoms in a covalent bond between both agents) get None instead.
        """
        ...


def energies(
    trajectory_file: File,
    # Set a folder to be created in order to store residual output files from this analysis
    output_directory: str,
    structure: 'Structure',
    interactions: list,
    charges: list,
    cg_selection: 'Selection',
    dummy_selection: 'Selection',
    snapshots: int,
    frames_limit: int = 100,
    verbose: bool = False,
    debug: bool = False,
):
    """Perform the electrostatic and vdw energies analysis for each pair of interaction agents."""
    # Make sure we have interactions
    if not interactions or len(interactions) == 0:
        print('No interactions were specified')
        return

    # Skip interactions with coarse grain regions
    cg_atom_indices = set(cg_selection.atom_indices)
    supported_interactions = []
    for interaction in interactions:
        has_cg = False
        for agent in ['1', '2']:
            atom_indices = interaction[f'atom_indices_{agent}']
            if set(atom_indices).intersection(cg_atom_indices):
                warn(
                    f'Interaction {interaction["name"]} will be skipped since it has Coarse Grain (CG) regions'
                )
                has_cg = True
                break
        if has_cg:
            continue
        supported_interactions.append(interaction)

    # If there are not atomistic interactions then we are done
    if len(supported_interactions) == 0:
        print('No supported interactions')
        return

    # Make sure we have charges
    if not charges or type(charges) == Exception:
        print('Atom charges are not available -> This analysis will be skipped')
        return

    # Set the main output filepath
    output_analysis_filepath = f'{output_directory}/{OUTPUT_ENERGIES_FILENAME}'

    # Set a backup file to store some results on the fly
    # This is useful to restore these values in case the analysis is disrupt since it is a long analysis
    energies_backup = File(output_directory + '/backup.json')

    # Check the number of atoms on each interacting agent
    # If there is any agent with too many atoms CMIP will fail so we must skip this specific energies analysis by now
    for interaction in supported_interactions:
        exceeds = False
        for agent in ['1', '2']:
            atoms = interaction[f'atom_indices_{agent}']
            atom_count = len(atoms)
            if atom_count >= CMIP_ATOM_LIMIT:
                exceeds = True
                break
        if exceeds:
            warn(
                f'{interaction["name"]} is exceeding the CMIP atom count limit of {CMIP_ATOM_LIMIT} and it will be skipped for this analysis'
            )
            interaction['exceeds'] = True

    # This analysis produces many residual output files
    # Create a new folder to store all output files so they do not overcrowd the main directory
    if not exists(output_directory):
        mkdir(output_directory)

    # Adapt the structure for cmip
    energies_structure = structure.copy()

    # Rename residues according to if they are terminals
    name_terminal_residues(energies_structure)

    # Set each atom element in CMIP format
    set_cmip_elements(energies_structure)

    # Set the method to calculate the energies of each frame
    engine = CmipEngine(
        output_directory, charges, energies_structure, dummy_selection, debug
    )

    # Extract the energies for each frame in a reduced trajectory
    frames, frame_step, frame_count = get_pdb_frames(
        energies_structure, trajectory_file.path, snapshots, frames_limit
    )
    non_exceeding_interactions = [
        interaction
        for interaction in supported_interactions
        if not interaction.get('exceeds', False)
    ]
    interactions_data = get_interactions_data(
        frames, non_exceeding_interactions, engine, energies_backup, verbose
    )

    # Cleanup here some residual files since the engine is not to be run again
    engine.cleanup()

    # Now calculated atom average values through all frames for each pair of interaction agents
    save_outputs(
        interactions_data,
        non_exceeding_interactions,
        frame_step,
        output_analysis_filepath,
    )

    # Remove the backup to avoid accidentally reusing it when the output file is deleted
    energies_backup.remove()


def get_interactions_data(
    frames: Generator[str, None, None],
    interactions: list[dict],
    engine: EnergiesEngine,
    energies_backup: File,
    verbose: bool = False,
) -> list[list[dict]]:
    """Get the atom energies of every interaction in every frame.

    Return one list per interaction with one dict per frame, with the atom energies of each agent.
    A backup is saved after every frame, and frames already in the backup are not calculated again.
    """
    # Load backup data in case there is a backup file
    if energies_backup.exists and energies_backup.size > 0:
        print(' Recovering energies backup')
        interactions_data = load_json(energies_backup.path)
    else:
        interactions_data = [[] for interaction in interactions]
    for frame_number, current_frame_pdb in enumerate(frames):
        # If we already have this frame in the backup then skip it
        if frame_number < len(interactions_data[0]):
            continue
        # Run the main analysis over the current frame
        # Append the result data for each interaction
        # WARNING: At this point structure should be corrected
        # WARNING: Repeated atoms will make the analysis fail
        current_frame_structure = Structure.from_pdb_file(current_frame_pdb)
        frame_energies_data = [
            get_interaction_frame_energies(
                engine, current_frame_structure, interaction, verbose
            )
            for interaction in interactions
        ]
        for i, data in enumerate(frame_energies_data):
            interactions_data[i].append(data)
        # Save a backup just in case the process is interrupted further
        save_json(interactions_data, energies_backup.path)
    return interactions_data


def get_interaction_frame_energies(
    engine: EnergiesEngine,
    frame_structure: 'Structure',
    interaction: dict,
    verbose: bool = False,
) -> dict:
    """Get the atom energies of both agents of an interaction in one frame."""
    # Get covalent bonds between both agents, if any
    # They will be not taken in count during the calculation
    # Otherwise we would have a huge energy peak in this atom since they are very close
    strong_bonds = interaction.get('strong_bonds', None)

    # Get interaction and agent names, just for the logs
    interaction_name = interaction['name']
    print(f' Processing {interaction_name}')
    agents = []
    for agent_number in ['1', '2']:
        agent_name = (
            interaction[f'agent_{agent_number}'].replace(' ', '_').replace('/', '_')
        )
        # Get agent atom selections
        agent_atom_indices = interaction[f'atom_indices_{agent_number}']
        agent_selection = frame_structure.select_atom_indices(agent_atom_indices)
        if not agent_selection:
            raise ValueError(
                f'Empty agent {agent_number} "{agent_name}" from interaction "{interaction_name}"'
            )
        agents.append(Agent(agent_name, agent_atom_indices, agent_selection))
    agent1, agent2 = agents

    # If both agents are far from each other then skip the engine and assume all energies are 0
    if are_far_apart(frame_structure, agent1.atom_indices, agent2.atom_indices):
        print(
            f'  Agents are further than {FAR_AGENTS_CUTOFF} Å -> Energies are set to 0'
        )
        agent1_atom_energies = [(0, 0, 0) for index in agent1.atom_indices]
        agent2_atom_energies = [(0, 0, 0) for index in agent2.atom_indices]
        print(ERASE_2_PREVIOUS_LINES)
        return {'agent1': agent1_atom_energies, 'agent2': agent2_atom_energies}

    agent1_atom_energies, agent2_atom_energies = engine(
        frame_structure, agent1, agent2, strong_bonds
    )

    # Make sure the output size matches the number of atoms in the interaction
    for agent_number, agent, agent_energies in [
        ('1', agent1, agent1_atom_energies),
        ('2', agent2, agent2_atom_energies),
    ]:
        if len(agent_energies) != len(agent.atom_indices):
            print(f'Number of values: {len(agent_energies)}')
            print(f'Number of atoms: {len(agent.atom_indices)}')
            raise ValueError(
                f'Mismatch in the number values and atoms in agent {agent_number}'
            )

    # Print total energies at the end for every agent if the verbose flag has been passed
    if verbose:
        for agent_name, agent_energies in zip(
            [agent1.name, agent2.name], [agent1_atom_energies, agent2_atom_energies]
        ):
            total_vdw = total_es = total_both = 0
            for atom_energies in agent_energies:
                total_vdw += atom_energies[0]
                total_es += atom_energies[1]
                total_both += atom_energies[2]
            print(
                f' Total energies for {agent_name}: vdw {total_vdw}, es {total_es}, both {total_both}'
            )

    # Erase the 3 previous log lines
    print(ERASE_3_PREVIOUS_LINES)

    return {'agent1': agent1_atom_energies, 'agent2': agent2_atom_energies}


def are_far_apart(
    frame_structure: 'Structure',
    atom_indices_1: list[int],
    atom_indices_2: list[int],
    cutoff: float = FAR_AGENTS_CUTOFF,
) -> bool:
    """Check if every atom of one agent is further than the cutoff (Å) from every atom of the other agent."""
    agent1_coords = np.array(
        [frame_structure.atoms[index].coords for index in atom_indices_1]
    )
    agent2_coords = np.array(
        [frame_structure.atoms[index].coords for index in atom_indices_2]
    )
    distances, _ = cKDTree(agent1_coords).query(
        agent2_coords, distance_upper_bound=cutoff
    )
    return bool(np.isinf(distances).all())


def save_outputs(
    interactions_data: list[list[dict]],
    interactions: list[dict],
    frame_step: int,
    output_analysis_filepath: str,
):
    """Save the atom average energies of each interaction in its own file, and a summary listing them."""
    output_summary = []
    for i, interaction in enumerate(interactions):
        # Get the numerated output filepath for this specific interaction
        numbered_output_analysis_filepath = numerate_filename(
            output_analysis_filepath, i
        )
        # Add the root of the output analysis filename to the run data
        analysis_name = get_analysis_name(numbered_output_analysis_filepath)
        # Get the interaction name
        name = interaction['name']
        # Add this interaction to the final summary
        output_summary.append({'name': name, 'analysis': analysis_name})
        # If the output file already exists then skip this iteration
        if exists(numbered_output_analysis_filepath):
            continue
        # Get the main data
        data = interactions_data[i]
        # Format data
        agent1_output = format_data([frame['agent1'] for frame in data], frame_step)
        agent2_output = format_data([frame['agent2'] for frame in data], frame_step)
        # Format the results data and append it to the output data
        output = {
            'name': name,
            'agent1': agent1_output,
            'agent2': agent2_output,
            'version': ENERGIES_ANALYSIS_VERSION,
        }
        # Export the current interaction analysis in json format
        save_json(output, numbered_output_analysis_filepath)

    # Finally, export the summary in json format
    save_json(output_summary, output_analysis_filepath)


class CmipEngine:
    """Interaction energies with CMIP: electrostatics with the Poisson-Boltzmann equation on a grid and van der Waals.

    Each agent is calculated as a guest in the potential of the other agent (the host).
    CMIP input and output files are written in the output directory.
    """

    def __init__(
        self,
        output_directory: str,
        charges: list,
        cmip_structure: 'Structure',
        dummy_selection: 'Selection',
        debug: bool = False,
    ):
        self.output_directory = output_directory
        self.charges = charges
        # Structure with the CMIP residue names and elements
        self.cmip_structure = cmip_structure
        self.dummy_selection = dummy_selection
        # If the debug flag is passed then, instead of calculating energies, leave it all ready and stop
        self.debug = debug
        # Set the auxiliar files further required as CMIP inputs
        self.cmip_inputs_source = File(CMIP_INPUTS_SOURCE)
        # Define an additional file generated by CMIP which must be handled (written in the energies folder and deleted)
        self.restart_file = File(output_directory + '/restart')

    def __call__(
        self,
        frame_structure: 'Structure',
        agent1: Agent,
        agent2: Agent,
        strong_bonds: Optional[list],
    ) -> tuple[list, list]:
        """Get the energies of every atom of each agent with the other agent."""
        # Prepare the CMIP friendly input pdb structures for every calculation
        # First prepare the host files and the prepare the guest files
        # There is only a difference between host and guest files:
        # Host files include both agent atoms but the guest atoms are all marked as dummy
        # Guest files include only guest atoms while host atoms are removed
        cmip_pdb_inputs = dict(
            frame_structure=frame_structure,
            strong_bonds=strong_bonds,
            dummy_selection=self.dummy_selection,
            charges=self.charges,
            cmip_structure=self.cmip_structure,
            output_directory=self.output_directory,
        )
        agent1_cmip_host = pdb2cmip(
            agent_name=agent1.name,
            host_file=True,
            host_selection=agent1.selection,
            guest_selection=agent2.selection,
            **cmip_pdb_inputs,
        )
        agent2_cmip_host = pdb2cmip(
            agent_name=agent2.name,
            host_file=True,
            host_selection=agent2.selection,
            guest_selection=agent1.selection,
            **cmip_pdb_inputs,
        )
        agent1_cmip_guest = pdb2cmip(
            agent_name=agent1.name,
            host_file=False,
            host_selection=agent2.selection,
            guest_selection=agent1.selection,
            **cmip_pdb_inputs,
        )
        agent2_cmip_guest = pdb2cmip(
            agent_name=agent2.name,
            host_file=False,
            host_selection=agent1.selection,
            guest_selection=agent2.selection,
            **cmip_pdb_inputs,
        )

        # Copy the source cmip inputs file in the local directory
        # Inputs will be modified to adapt the cmip grid to both agents together
        # Note than modified inputs file is not conserved along frames
        # Structures may change along trajectory thus requiring a different grid size
        cmip_inputs = File(
            f'{self.output_directory}/{self.cmip_inputs_source.filename}'
        )
        copyfile(self.cmip_inputs_source.path, cmip_inputs.path)

        # Set the CMIP box dimensions and densities to fit both the host and the guest
        # Box origin and size are modified in the cmip inputs
        # Values returned are only used for display / debug purposes
        box_origin, box_size = adapt_cmip_grid(
            agent1_cmip_guest,
            agent2_cmip_guest,
            cmip_inputs,
            self.output_directory,
            self.restart_file,
        )

        # If the debug flag is passed then, instead of calculating energies, leave it all ready and stop here
        if self.debug:
            write_cmip_debug_files(
                self.output_directory,
                cmip_inputs,
                agent1_cmip_guest,
                agent2_cmip_host,
                agent2_cmip_guest,
                agent1_cmip_host,
            )

        # Run the CMIP software to get the desired energies
        print(
            f'  Calculating energies for {agent1.name} as guest and {agent2.name} as host'
        )
        agent1_atom_energies = get_cmip_energies(
            cmip_inputs,
            agent1_cmip_guest,
            agent2_cmip_host,
            self.output_directory,
            self.restart_file,
        )
        print(
            f'  Calculating energies for {agent2.name} as guest and {agent1.name} as host'
        )
        agent2_atom_energies = get_cmip_energies(
            cmip_inputs,
            agent2_cmip_guest,
            agent1_cmip_host,
            self.output_directory,
            self.restart_file,
        )

        # Add None in those places where there is a strong bond
        strong_bond_atom_indices = set(sum(strong_bonds, []))
        for i, atom_index in enumerate(agent1.atom_indices):
            if atom_index in strong_bond_atom_indices:
                agent1_atom_energies.insert(i, None)
        for i, atom_index in enumerate(agent2.atom_indices):
            if atom_index in strong_bond_atom_indices:
                agent2_atom_energies.insert(i, None)

        return agent1_atom_energies, agent2_atom_energies

    def cleanup(self):
        """Remove the residual CMIP files which are not needed once all frames are done."""
        # Remove the restart file since we do not use it and it may be heavy sometimes
        if self.restart_file.exists:
            self.restart_file.remove()
        # Remove fortran unit files generated when running CMIP
        # Note that we can not define where these files are written but they appear where CMIP is run
        fortran_unit_files = glob('fort.*')
        for filepath in fortran_unit_files:
            remove(filepath)


def pdb2cmip(
    agent_name: str,
    host_file: bool,
    frame_structure: 'Structure',
    host_selection: 'Selection',
    guest_selection: 'Selection',
    strong_bonds: Optional[list],
    dummy_selection: 'Selection',
    charges: list,
    cmip_structure: 'Structure',
    output_directory: str,
) -> File:
    """Transform an agent structure to a cmip input pdb, which includes charges and cmip-friendly elements.

    If this is to be a host file then leave only atoms involved in the interaction: both host and guest agents.
    If this is to be a host file, all guest atoms are set dummy* as well.
    If this is to be a guest file then remove all host atoms as well.
    Note that the frame structure is not the standard structure, which would already include charges.
    For this reason, charges are included in the frame structure before.
    Also flag some atoms as 'dummy' by adding a 'X' before the element.
    *Dummy atoms are not considered in the calculation but they stand for a region with low dielectric.
    If removed, the void left by these atoms would be considered to be filled with solvent, which has a high dielectric.
    These dielectric differences have a strong impact on the calculation.
    By default we set as dummy host atoms involved in a covalent bond with guest atoms.
    """
    # Get atom indices for host atoms involved in a covalent bond with guest atoms
    # If this is the host, they will be further marked as dummy for CMIP to ignore them in the calculations
    # If this is the guest, they will be removed
    strong_bond_atom_indices = []
    if strong_bonds:
        for bond in strong_bonds:
            strong_bond_atom_indices += bond
    # Set atoms to be flagged as dummy
    # Start including atoms which are actually "dummy", not in the sense of CMIP, but in general
    dummy_atom_indices = set(dummy_selection.atom_indices)
    # If this is to be the host file then set guest atoms and strong bond atoms dummy
    if host_file:
        dummy_atom_indices |= set(strong_bond_atom_indices)
        dummy_atom_indices |= set(guest_selection.atom_indices)
    # Set which atoms must be kept in the output file
    # If this is the host file then keep both host and guest atoms
    # If this is the guest file then keep only guest atoms
    # Always remove strong bond atoms in the guest
    selection = host_selection + guest_selection if host_file else guest_selection
    strong_bond_selection = frame_structure.select_atom_indices(
        strong_bond_atom_indices
    )
    guest_strong_bonds = guest_selection & strong_bond_selection
    selection -= guest_strong_bonds
    # Filter the selected atoms in the structure
    selected_structure = frame_structure.filter(selection)
    # Set atom charges (non standard attribute), which are used further to write the adapted cmip pdb
    # Set also the elements to match the original structure, since the frame generator messes the elements
    for a, atom in enumerate(selected_structure.atoms):
        atom_index = selection.atom_indices[a]
        charge = charges[atom_index]
        setattr(atom, 'charge', charge)
        cmip_element = cmip_structure.atoms[atom_index].element
        atom.element = cmip_element
    # Write a special pdb which contains charges as CMIP expects to find them and dummy atoms flagged
    output_filepath = (
        f'{output_directory}/{agent_name}_{"host" if host_file else "guest"}.cmip.pdb'
    )
    # Chain names must fit in 1 character
    chain_names_map = selected_structure.get_pdb_chain_names_map(
        show_warnings=False, accept_limitations=True
    )
    with open(output_filepath, 'w') as file:
        # Write a line for each atom
        for a, atom in enumerate(selected_structure.atoms):
            index = str(a + 1).rjust(5)
            atom_name = atom.name
            name = ' ' + atom_name.ljust(3) if len(atom_name) < 4 else atom_name
            residue = atom.residue
            residue_name = residue.name.ljust(3)
            if len(residue_name) > 3:
                residue_name = residue_name[0:3]
            chain = atom.chain
            chain_name = chain_names_map[chain.name].rjust(1)
            residue_number = str(residue.number).rjust(4)
            icode = residue.icode.rjust(1)
            coords = atom.coords
            x_coord, y_coord, z_coord = [
                '{:.3f}'.format(coord).rjust(8) for coord in coords
            ]
            charge = '{:.4f}'.format(atom.charge).rjust(7)
            # In case this atom is making an strong bond between both interacting agents we add an 'X' before the element
            # This way CMIP will ignore the atom. Otherwise it would return high non-sense Van der Waals values
            real_index = selection.atom_indices[a]
            # Set if atom is to be ignored
            is_dummy = real_index in dummy_atom_indices
            cmip_dummy_flag = 'X' if is_dummy else ''
            # Set dummy atoms as hydrogen, so their VDW is minimal
            element = 'H' if is_dummy else atom.element
            # Set the tail
            tail = cmip_dummy_flag + element
            atom_line = (
                f'ATOM  {index} {name} {residue_name} {chain_name}{residue_number}{icode}'
                f'   {x_coord}{y_coord}{z_coord} {charge}  {tail}\n'
            )
            # IMPORTANT: Make sure every line has the expected length, apart from the tail
            # Otherwise it means one of its parts is longer or shorter than expected
            # This will offset the position of the numbers which CMIP expects to read
            # Thus leading to missleading fortran errors
            # e.g. Fortran runtime error: Bad value during floating point read
            line_length = len(atom_line) - len(tail)
            if line_length != EXPECTED_CMIP_PDB_LINE_LENGTH:
                clean_line = atom_line.replace('\n', '')
                print(
                    f'Line "{clean_line}" has an unexpected length ({line_length}/{EXPECTED_CMIP_PDB_LINE_LENGTH}).'
                )
                raise RuntimeError(
                    'Wrong formatted CMIP input PDB file. Please report this problem to the MDDB developers.'
                )
            file.write(atom_line)

    return File(output_filepath)


def run_cmip(
    cmip_inputs: File,
    probe: File,
    host: File,
    output_file: File,
    restart_file: File,
) -> CompletedProcess:
    """Run CMIP with the given inputs, probe (guest) and host, and the CMIP van der Waals parameters."""
    vdw_source = File(CMIP_VDW_SOURCE)
    return run(
        [
            'cmip',
            '-i',
            cmip_inputs.path,
            '-pr',
            probe.path,
            '-vdw',
            vdw_source.path,
            '-hs',
            host.path,
            '-byat',
            output_file.path,
            '-rst',
            restart_file.path,
        ],
        stdout=PIPE,
        stderr=PIPE,
    )


def get_cmip_grid(
    probe: File, host: File, output_file: File, restart_file: File
) -> tuple[tuple, tuple, tuple]:
    """Run CMIP in 'checkonly' mode (i.e. start and stop) and mine the center, density and units of its grid."""
    cmip_process = run_cmip(
        File(CMIP_INPUTS_CHECKONLY_SOURCE), probe, host, output_file, restart_file
    )
    cmip_logs = cmip_process.stdout.decode()
    # Mine the grid dimensions from CMIP logs
    mining_results = mine_cmip_output(cmip_logs)
    if mining_results is OUTPUT_MINING_EXCEPTION:
        print(cmip_logs)
        cmip_error_logs = cmip_process.stderr.decode()
        print(cmip_error_logs)
        raise OUTPUT_MINING_EXCEPTION
    return mining_results


def adapt_cmip_grid(
    agent1_cmip_pdb: File,
    agent2_cmip_pdb: File,
    cmip_inputs: File,
    output_directory: str,
    restart_file: File,
) -> tuple:
    """Set a CMIP grid which includes both agents in the CMIP inputs.

    Run CMIP in 'checkonly' mode (i.e. start and stop) for both agents.
    Mine the grid generated by CMIP for both agents and calculate a new grid which would include both.
    Finally, modify the provided 'cmip_inputs' with the new grid parameters.
    Return the origin and the size of the new grid, which are used for display / debug purposes only.
    """
    # Set a name for the checkonly CMIP outputs
    # This name is not important, since the data we want is in the CMIP logs
    cmip_checkonly_output = File(output_directory + '/checkonly.energy.pdb')

    # Run CMIP in 'checkonly' mode and save the grid dimensions output, for each agent
    agent1_center, agent1_density, agent1_units = get_cmip_grid(
        agent1_cmip_pdb, agent2_cmip_pdb, cmip_checkonly_output, restart_file
    )
    agent2_center, agent2_density, agent2_units = get_cmip_grid(
        agent2_cmip_pdb, agent1_cmip_pdb, cmip_checkonly_output, restart_file
    )

    # Calculate grid dimensions for a new grid which contains both previous grids
    new_center, new_density = compute_new_grid(
        agent1_center,
        agent1_density,
        agent1_units,
        agent2_center,
        agent2_density,
        agent2_units,
    )

    # In case the number of grid points exceeds a safe limit we must reduce it
    new_density, grid_unit_size = cap_grid_points(new_density, GRID_UNIT_SIZE)

    # Add the new grid parameters to the local inputs file
    write_cmip_grid(cmip_inputs, new_center, new_density, grid_unit_size)

    # Delete the 'ckeckonly' file
    cmip_checkonly_output.remove()

    # Calculate the resulting box origin and size and return both values
    new_size = (
        new_density[0] * grid_unit_size,
        new_density[1] * grid_unit_size,
        new_density[2] * grid_unit_size,
    )
    new_origin = (
        new_center[0] - new_size[0] / 2,
        new_center[1] - new_size[1] / 2,
        new_center[2] - new_size[2] / 2,
    )
    return new_origin, new_size


def cap_grid_points(density: list, grid_unit_size: float) -> tuple[list, float]:
    """Reduce the grid density in case the number of grid points exceeds the GRID_POINTS_LIMIT.

    Reduce all dimension densities in proportion and expand the grid units size to compensate.
    Return the new density and grid unit size.
    """
    # Calculate the amount of total points with current parameters
    current_grid_points = (density[0] + 1) * (density[1] + 1) * (density[2] + 1)
    if current_grid_points <= GRID_POINTS_LIMIT:
        return density, grid_unit_size
    print(f'WARNING: Grid points limit is exceeded ({current_grid_points})')
    proportion = GRID_POINTS_LIMIT / current_grid_points
    density = [math.ceil(dimension * proportion) for dimension in density]
    grid_unit_size = math.ceil(grid_unit_size / proportion * 1000) / 1000
    print(f'WARNING: Grid resolution has been reduced -> unit size = {grid_unit_size}')
    return density, grid_unit_size


def write_cmip_grid(
    cmip_inputs: File, center: list, density: list, grid_unit_size: float
):
    """Add the grid center, density and unit size to a CMIP inputs file."""
    # Set the new lines to be written in the local CMIP inputs file
    grid_inputs = [
        f' cenx={center[0]:.1f} \n',
        f' ceny={center[1]:.1f} \n',
        f' cenz={center[2]:.1f} \n',
        f' dimx={int(density[0])} \n',
        f' dimy={int(density[1])} \n',
        f' dimz={int(density[2])} \n',
        f' intx={grid_unit_size} \n',
        f' inty={grid_unit_size} \n',
        f' intz={grid_unit_size} \n',
    ]
    # Add previous lines to the local inputs file
    with open(cmip_inputs.path, 'r+') as file:
        lines = file.readlines()
        file.seek(0)
        for line in lines:
            if line == '&end \n':
                for grid_input in grid_inputs:
                    file.write(grid_input)
            file.write(line)


def write_cmip_debug_files(
    output_directory: str,
    cmip_inputs: File,
    agent1_cmip_guest: File,
    agent2_cmip_host: File,
    agent2_cmip_guest: File,
    agent1_cmip_host: File,
):
    """Leave everything ready to run CMIP by hand in the output directory, with a README, and stop the workflow."""
    vdw_source = File(CMIP_VDW_SOURCE)
    debug_script_source = File(ENERGIES_DEBUG_SCRIPT_SOURCE)
    # Copy in the energies folder a small python script used to sum output energies
    debug_script_filename = debug_script_source.filename
    debug_script = File(f'{output_directory}/{debug_script_filename}')
    copyfile(debug_script_source.path, debug_script.path)
    # Copy in the energies folder the VDM parameters input file
    debug_vdw_filename = vdw_source.filename
    debug_vdw_parameters = File(f'{output_directory}/{debug_vdw_filename}')
    copyfile(vdw_source.path, debug_vdw_parameters.path)
    # Set auxiliar output filenames
    debug_output_1 = 'first_output.pdb'
    debug_output_2 = 'second_output.pdb'
    # Write a README file explaining what to do to manually debug CMIP
    debug_readme = File(f'{output_directory}/README')
    with open(debug_readme.path, 'w') as file:
        file.write(
            '# Energies debug\n\n'
            '# Run CMIP\n'
            f'cmip -i {cmip_inputs.filename} -pr {agent1_cmip_guest.filename} -vdw {debug_vdw_filename} -hs {agent2_cmip_host.filename} -byat {debug_output_1}\n\n'
            '# Run CMIP inverting host and guest\n'
            f'cmip -i {cmip_inputs.filename} -pr {agent2_cmip_guest.filename} -vdw {debug_vdw_filename} -hs {agent1_cmip_host.filename} -byat {debug_output_2}\n\n'
            '# Sum both energies and compare\n'
            f'python {debug_script_filename} {debug_output_1} {debug_output_2}\n'
        )
    raise ForcedStop(
        ' READY TO DEBUG -> Please go to the corresponding replica "energies" directory and follow the README instructions'
    )


def get_cmip_energies(
    cmip_inputs: File,
    guest: File,
    host: File,
    output_directory: str,
    restart_file: File,
) -> list[tuple]:
    """Run the CMIP software to get the energies of every guest atom in the potential of the host."""
    # Set the cmip output filename, which is to be read rigth after it is generated
    cmip_output_file = File(output_directory + '/cmip_output.pdb')
    # Run cmip
    cmip_process = run_cmip(cmip_inputs, guest, host, cmip_output_file, restart_file)
    cmip_logs = cmip_process.stdout.decode()
    # If the output does note exist at this poin then it means something went wrong with CMIP
    if not exists(cmip_output_file.path):
        print(cmip_logs)
        cmip_error_logs = cmip_process.stderr.decode()
        print(cmip_error_logs)
        raise RuntimeError('Something went wrong with CMIP!')
    with open(cmip_output_file.path, 'r') as file:
        lines = list(file)
    # If this file is empty it means something went wrong with CMIP
    # We print its logs and exit
    if len(lines) == 0:
        print(cmip_logs)
        raise RuntimeError('Something went wrong with CMIP!')
    # Mine energies line by line (i.e. atom by atom)
    return [mine_cmip_atom_energies(line) for line in lines]


def mine_cmip_atom_energies(line: str) -> tuple:
    """Mine the Van der Waals (vdw), electrostatic (es) and total (both) energies of an atom from a CMIP output line.

    Atoms with a total energy above MAX_ATOM_ENERGY get 0 for all three energies.
    """
    # Numbers do not always fit in fixed columns, since numbers close to 0 are expressed in E notation
    # e.g. 0.1075E-04
    # For this reason we must mine these last values relying on whitespaces between them
    line_splits = line.split()
    vdw = float(line_splits[-3])
    es = float(line_splits[-2])
    both = float(line_splits[-1])
    # Values greater than the limit are represented as 0
    if both > MAX_ATOM_ENERGY:
        return (0, 0, 0)
    return (vdw, es, both)


def set_cmip_elements(structure: 'Structure'):
    """Set the atom elements in the CMIP format.

    Hydrogens bonded to carbons remain as 'H', hydrogens bonded to oxygen are renamed as 'HO' and hydrogens bonded to
    nitrogen or sulfur are renamed as 'HN'.
    """
    # Iterate over each atom to fix their element according to CMIP standards
    for a, atom in enumerate(structure.atoms):
        element = atom.element
        # Adapt hydrogens element to CMIP requirements
        if element == 'H':
            # We must find the element of the heavy atom this hydrogen is bonded to
            atom_bonds = atom.get_bonds()
            # There should be always only 1 bond
            if len(atom_bonds) != 1:
                print(f'Atom {atom.name} ({atom.index}) has {len(atom_bonds)} bonds')
                raise ValueError('An hydrogen should always have one and only one bond')
            bonded_atom_index = atom_bonds[0]
            bonded_atom_element = structure.atoms[bonded_atom_index].element
            # Hydrogens bonded to carbons remain as 'H'
            if bonded_atom_element == 'C':
                pass
            # Hydrogens bonded to oxygen are renamed as 'HO'
            elif bonded_atom_element == 'O':
                element = 'HO'
            # Hydrogens bonded to nitrogen or sulfur are renamed as 'HN'
            elif bonded_atom_element == 'N' or bonded_atom_element == 'S':
                element = 'HN'
            else:
                raise RuntimeError(
                    f'Hydrogen bonded to not supported heavy atom: {bonded_atom_element}'
                )
        atom.element = element


# Set residue names with no endings
protein_mid_residue_names = {
    name for name in PROTEIN_RESIDUE_NAME_LETTERS if name[-1] not in {'N', 'C'}
}
nucleic_mid_residue_names = {
    name for name in NUCLEIC_RESIDUE_NAME_LETTERS if name[-1] not in {'5', '3'}
}


def name_terminal_residues(structure: 'Structure'):
    """Rename the terminal residues of protein and nucleic acid chains as CMIP expects.

    Protein terminal residues get an 'N' or 'C' at the end of their names, and nucleic acid ones a '5' or '3'.
    """
    # Iterate chains
    for chain in structure.chains:
        # Check if the first residue is tagged as a terminal residue
        # If not, rename it
        first_residue = chain.residues[0]
        # In case it is a protein
        if first_residue.name in protein_mid_residue_names:
            first_residue.name += 'N'
        # In case it is RNA
        elif first_residue.name in nucleic_mid_residue_names:
            first_residue.name += '5'
        # Check if the last residue is tagged as 'C' terminal
        # If not, rename it
        last_residue = chain.residues[-1]
        # In case it is a protein
        if last_residue.name in protein_mid_residue_names:
            last_residue.name += 'C'
        # In case it is RNA
        elif last_residue.name in nucleic_mid_residue_names:
            last_residue.name += '3'


def mine_cmip_output(logs: str) -> tuple[tuple, tuple, tuple] | Exception:
    """Mine the center, density and units of the grid from CMIP output logs.

    If any of them is missing then return the OUTPUT_MINING_EXCEPTION.
    """
    center, density, units = (), (), ()
    # Set the regular expressions to mine the relevant data
    grid_density_exp = r'^\s*Grid density:\s+(\d+)\s+(\d+)\s+(\d+)'
    grid_center_exp = r'^\s*Grid center:([- ]+\d+.\d+)([- ]+\d+.\d+)([- ]+\d+.\d+)'
    grid_units_exp = r'^\s*Grid units:\s+(\d+.\d+)\s+(\d+.\d+)\s+(\d+.\d+)'
    # Iterate log lines
    for line in logs.split('\n'):
        grid_center_groups = re.match(grid_center_exp, line)
        grid_density_groups = re.match(grid_density_exp, line)
        grid_units_groups = re.match(grid_units_exp, line)
        if grid_density_groups:
            density = tuple(float(grid_density_groups.group(i)) for i in (1, 2, 3))
        if grid_center_groups:
            center = tuple(float(grid_center_groups.group(i)) for i in (1, 2, 3))
        if grid_units_groups:
            units = tuple(float(grid_units_groups.group(i)) for i in (1, 2, 3))
    # If data mining fails there must be something wrong with the CMIP output
    if center == () or density == () or units == ():
        return OUTPUT_MINING_EXCEPTION
    return center, density, units


def compute_new_grid(
    agent1_center,
    agent1_density,
    agent1_units,
    agent2_center,
    agent2_density,
    agent2_units,
    extra_density=1,
):
    """Compute the center and density of a new grid, which is expected to cover the grids of both agents."""
    new_center = []
    new_density = []

    for k in range(3):
        min_prot = agent1_center[k] - agent1_density[k] * agent1_units[k]
        min_lig = agent2_center[k] - agent2_density[k] * agent2_units[k]
        min_new = min(min_prot, min_lig)
        max_prot = agent1_center[k] + agent1_density[k] * agent1_units[k]
        max_lig = agent2_center[k] + agent2_density[k] * agent2_units[k]
        max_new = max(max_prot, max_lig)

        dnew = int(abs(max_new - min_new) + extra_density)
        cnew = min_new + (dnew / 2)
        new_density.append(dnew)
        new_center.append(cnew)
    return new_center, new_density


def format_data(data: list, frame_step: int) -> dict:
    """Average the atom energies over all frames, and over the first and the last 20% of frames.

    Data has one list per frame with one (vdw, es, both) tuple or None per atom.
    Averages skip None values, and atoms with None in all the averaged frames get None.
    """
    # Set an array of energies by frame, atom and energy (vdw, es, both), with NaN where values are None
    frame_count = len(data)
    atom_count = len(data[0])
    values = np.full((frame_count, atom_count, 3), np.nan)
    for f, frame in enumerate(data):
        for a, atom_values in enumerate(frame):
            if atom_values is not None:
                values[f, a] = atom_values

    # Calculate the atom averages from each energy at the beginning and end of the trajectory
    # We take the initial 20% and the final 20% of frames to calculate each respectively
    p20 = round(frame_count * 0.2)
    if p20 == 0:
        p20 = 1

    # Format the results data and append it to the output data
    output = {}
    for prefix, frame_values in [
        ('', values),
        ('i', values[:p20]),
        ('f', values[-p20:]),
    ]:
        averages = average_frames(frame_values)
        for e, energy in enumerate(['vdw', 'es', 'both']):
            output[prefix + energy] = averages[:, e].tolist()
    output['step'] = frame_step
    return output


def average_frames(frame_values: np.ndarray) -> np.ndarray:
    """Average values along frames (first axis) skipping NaN values; columns with no values get None.

    Values are added frame after frame, so the averages are the same as summing them in a loop.
    """
    counts = np.sum(~np.isnan(frame_values), axis=0)
    sums = np.nansum(frame_values, axis=0)
    averages = np.divide(
        sums, counts, out=np.full(sums.shape, np.nan), where=counts > 0
    )
    return np.where(counts > 0, averages, None)
