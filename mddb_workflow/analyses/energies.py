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
from os.path import exists, abspath
from shutil import copyfile
from glob import glob
from subprocess import run, PIPE, CompletedProcess
from typing import BinaryIO, NamedTuple, Protocol
import re
import math
import json
import numpy as np
from scipy.spatial import cKDTree

from mddb_workflow.tools.get_pdb_frames import get_pdb_frames, get_coordinate_frames
from mddb_workflow.utils.auxiliar import (
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
# Set the margin (Å) added on every side of the automatic CMIP grids of both agents to make the grid of the interaction
# This value was measured on 2026-11-12 (benchmark T2, CMIP variant 'fixed') and it was not tuned
GRID_MARGIN = 5.0
# Atom energies (vdw + es, kcal/mol) above this value are set to 0, counted and reported
# This step is performed to filter 'infinity' values, which come from overlapping atoms (clashes)
# Large negative values are kept since they happen in real interactions
# e.g. the Zn2+ of carbonic anhydrase coordinated by sulfonamide ligands, -100 to -110 kcal/mol (PL-REX)
MAX_ATOM_ENERGY = 100
# Grid-free engine: Coulomb constant (kcal Å / mol e²), bulk dielectric of the Mehler-Solmajer sigmoidal dielectric,
# the van der Waals radii scale (fvdw in the CMIP inputs) and cutoff (Å, CMIP's default cutvdw)
COULOMB_CONSTANT = 332.0637
SIGMOIDAL_BULK_DIELECTRIC = 78.4
VDW_RADII_SCALE = 0.8
VDW_CUTOFF = 8.0
ENERGIES_ENGINES = ['cmip', 'gridfree']
# The backup lets an interrupted analysis resume. Its first line describes the run (version, engine, interactions and
# frames) and each further line holds the energies of one frame for every interaction, appended when the frame is done.
# Up to version 1.1.0 the backup was one JSON rewritten after every frame, so its cost grew with the square of the
# number of frames. A backup in that format is removed, not reused.
ENERGIES_BACKUP_FILENAME = 'backup.jsonl'
OLD_ENERGIES_BACKUP_FILENAME = 'backup.json'
# Set an exception to handle internal issues
OUTPUT_MINING_EXCEPTION = Exception('CMIP output mining failed')
EXPECTED_CMIP_PDB_LINE_LENGTH = 65
# CMIP runs in the output directory, with its inputs copied there under these names
# Paths are short this way: CMIP cuts file paths longer than 150 characters
CMIP_VDW_FILENAME = 'vdwprm'
CMIP_INPUTS_CHECKONLY_FILENAME = 'cmip_check.in'
# Van der Waals parameters for elements missing from the CMIP table (resources/vdwprm): R*/2 (Å), epsilon (kcal/mol),
# mass and atomic number
# CMIP matches atom types to the table as prefixes, in the order of the table, so a missing element whose first letter
# is in the table silently gets that type (e.g. Ca would be carbon) and other missing elements make CMIP fail
# The table is the Amber parm99 nonbonded set, so missing elements get parameters from the same parm99 file when it
# has them (Aqvist ions), and otherwise from the Li/Merz 12-6 ion set for TIP3P water in AmberTools
# (frcmod.ions234lm_126_tip3p), with the 2+ ion when there is one and else the ion of lowest charge
CMIP_EXTRA_VDW_PARAMETERS = {
    'Li': (1.137, 0.0183, 6.94, 3),  # parm99 Li+
    'Be': (0.956, 3.95e-06, 9.01, 4),  # Li/Merz 12-6 Be2+
    'Al': (1.297, 0.00471279, 26.98, 13),  # Li/Merz 12-6 Al3+
    'Ca': (1.7131, 0.459789, 40.08, 20),  # parm99 Ca2+
    'V': (1.364, 0.01067299, 50.94, 23),  # Li/Merz 12-6 V2+
    'Cr': (1.346, 0.00868178, 52.0, 24),  # Li/Merz 12-6 Cr2+
    'Mn': (1.407, 0.0168671, 54.94, 25),  # Li/Merz 12-6 Mn2+
    'Fe': (1.353, 0.00941798, 55.84, 26),  # Li/Merz 12-6 Fe2+
    'Co': (1.299, 0.00483892, 58.93, 27),  # Li/Merz 12-6 Co2+
    'Ni': (1.255, 0.0026232, 58.69, 28),  # Li/Merz 12-6 Ni2+
    'Cu': (1.218, 0.00148497, 63.55, 29),  # Li/Merz 12-6 Cu2+
    'Rb': (2.956, 0.00017, 85.47, 37),  # parm99 Rb+
    'Sr': (1.81, 0.22132374, 87.62, 38),  # Li/Merz 12-6 Sr2+
    'Y': (1.602, 0.08034231, 88.91, 39),  # Li/Merz 12-6 Y3+
    'Zr': (1.519, 0.04525501, 91.22, 40),  # Li/Merz 12-6 Zr4+
    'Pd': (1.303, 0.00509941, 106.42, 46),  # Li/Merz 12-6 Pd2+
    'Ag': (1.336, 0.00770969, 107.87, 47),  # Li/Merz 12-6 Ag2+
    'Cd': (1.412, 0.01773416, 112.41, 48),  # Li/Merz 12-6 Cd2+
    'In': (1.461, 0.02808726, 114.82, 49),  # Li/Merz 12-6 In3+
    'Sn': (1.666, 0.11617738, 118.71, 50),  # Li/Merz 12-6 Sn2+
    'Cs': (3.395, 8.06e-05, 132.91, 55),  # parm99 Cs+
    'Ba': (2.019, 0.40664608, 137.33, 56),  # Li/Merz 12-6 Ba2+
    'La': (1.718, 0.15060822, 138.91, 57),  # Li/Merz 12-6 La3+
    'Ce': (1.741, 0.16721338, 140.12, 58),  # Li/Merz 12-6 Ce3+
    'Pr': (1.733, 0.16134811, 140.91, 59),  # Li/Merz 12-6 Pr3+
    'Nd': (1.681, 0.12564307, 144.24, 60),  # Li/Merz 12-6 Nd3+
    'Sm': (1.819, 0.22878796, 150.36, 62),  # Li/Merz 12-6 Sm2+
    'Eu': (1.802, 0.21475916, 151.96, 63),  # Li/Merz 12-6 Eu2+
    'Gd': (1.623, 0.09126804, 157.25, 64),  # Li/Merz 12-6 Gd3+
    'Tb': (1.63, 0.09509276, 158.93, 65),  # Li/Merz 12-6 Tb3+
    'Dy': (1.609, 0.0838924, 162.5, 66),  # Li/Merz 12-6 Dy3+
    'Er': (1.602, 0.08034231, 167.26, 68),  # Li/Merz 12-6 Er3+
    'Tm': (1.602, 0.08034231, 168.93, 69),  # Li/Merz 12-6 Tm3+
    'Yb': (1.642, 0.10185975, 173.04, 70),  # Li/Merz 12-6 Yb2+
    'Lu': (1.588, 0.07351892, 174.97, 71),  # Li/Merz 12-6 Lu3+
    'Hf': (1.499, 0.03868661, 178.49, 72),  # Li/Merz 12-6 Hf4+
    'Pt': (1.266, 0.00307642, 195.08, 78),  # Li/Merz 12-6 Pt2+
    'Hg': (1.407, 0.0168671, 200.59, 80),  # Li/Merz 12-6 Hg2+
    'Tl': (1.513, 0.04321029, 204.38, 81),  # Li/Merz 12-6 Tl3+
    'Pb': (1.745, 0.17018074, 207.2, 82),  # Li/Merz 12-6 Pb2+
    'Ra': (2.019, 0.40664608, 226.0, 88),  # Li/Merz 12-6 Ra2+
    'Th': (1.708, 0.1436416, 232.04, 90),  # Li/Merz 12-6 Th4+
    'U': (1.684, 0.12758274, 238.03, 92),  # Li/Merz 12-6 U4+
    'Pu': (1.662, 0.11371963, 244.0, 94),  # Li/Merz 12-6 Pu4+
}
# Elements missing from the CMIP table which take the CMIP type of a chemically similar element
# Selenium (selenomethionine, selenocysteine) takes the type of sulfur, as CMIP gave it by prefix up to version 1.1.0
CMIP_ELEMENT_ALIASES = {'Se': 'S'}


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
    energies_engine: str = 'cmip',
):
    """Perform the electrostatic and vdw energies analysis for each pair of interaction agents.

    The energies engine is 'cmip' (Poisson-Boltzmann electrostatics, the default) or 'gridfree' (sigmoidal dielectric).
    """
    if energies_engine not in ENERGIES_ENGINES:
        raise ValueError(
            f'Unknown energies engine "{energies_engine}". Use one of: {", ".join(ENERGIES_ENGINES)}'
        )
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

    # Set a backup file to store the results of every frame on the fly, so an interrupted analysis can resume
    energies_backup = File(f'{output_directory}/{ENERGIES_BACKUP_FILENAME}')

    # Check the number of atoms on each interacting agent
    # If there is any agent with too many atoms CMIP will fail so we must skip this specific energies analysis by now
    # Note that the input interactions are not modified
    non_exceeding_interactions = []
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
            continue
        non_exceeding_interactions.append(interaction)

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

    # Check every element has CMIP van der Waals parameters
    # Elements missing from the CMIP table get parameters from AmberTools, and interactions with other elements are skipped
    cmip_vdw_types = read_cmip_vdw_types(CMIP_VDW_SOURCE)
    dummy_atom_indices = set(dummy_selection.atom_indices)
    extra_vdw_elements = set()
    computed_interactions = []
    for interaction in non_exceeding_interactions:
        elements = {
            energies_structure.atoms[index].element
            for index in interaction['atom_indices_1'] + interaction['atom_indices_2']
            if index not in dummy_atom_indices
        }
        missing_elements = elements - cmip_vdw_types
        unsupported_elements = missing_elements - set(CMIP_EXTRA_VDW_PARAMETERS)
        if unsupported_elements:
            warn(
                f'{interaction["name"]} has atoms of elements with no van der Waals parameters for CMIP '
                f'({", ".join(sorted(unsupported_elements))}) and it will be skipped for this analysis'
            )
            continue
        extra_vdw_elements |= missing_elements
        computed_interactions.append(interaction)
    if len(computed_interactions) == 0:
        print('No interactions can be calculated')
        return
    if extra_vdw_elements:
        warn(
            f'Elements {", ".join(sorted(extra_vdw_elements))} are not in the CMIP van der Waals table. '
            'They get the parameters of their ions in AmberTools (Aqvist ions in parm99, or Li/Merz 12-6 ions)'
        )

    # Set the method to calculate the energies of each frame
    if energies_engine == 'gridfree':
        engine = GridFreeEngine(
            charges, energies_structure, dummy_selection, extra_vdw_elements
        )
    else:
        engine = CmipEngine(
            output_directory,
            charges,
            energies_structure,
            dummy_selection,
            debug,
            extra_vdw_elements,
        )

    # Extract the energies for each frame in a reduced trajectory
    # CMIP's input files take the atom and residue names, residue numbers and chains of the frame structure, which
    # writing and reading each frame as a PDB file fits to the PDB format. The grid-free engine only needs coordinates.
    get_frame_structures = (
        get_coordinate_frame_structures if energies_engine == 'gridfree' else get_pdb_frame_structures
    )
    frame_structures, frame_step, frame_count = get_frame_structures(
        energies_structure, trajectory_file.path, snapshots, frames_limit
    )
    # A backup of an older version is not reused, so its frames are never mixed with the new ones
    old_backup_filepath = f'{output_directory}/{OLD_ENERGIES_BACKUP_FILENAME}'
    if exists(old_backup_filepath):
        print(' Removing an energies backup of an older version')
        remove(old_backup_filepath)
    # The backup is reused only by a run with the same setup
    backup_header = {
        'version': ENERGIES_ANALYSIS_VERSION,
        'engine': energies_engine,
        'interactions': [interaction['name'] for interaction in computed_interactions],
        'frame_step': int(frame_step),
        'frame_count': int(frame_count),
    }
    interactions_data = get_interactions_data(
        frame_structures, computed_interactions, engine, energies_backup, backup_header, verbose
    )

    # Report the atoms set to 0 since their energies are too high
    warn_zeroed_atoms(computed_interactions, interactions_data)

    # Cleanup here some residual files since the engine is not to be run again
    engine.cleanup()

    # Now calculated atom average values through all frames for each pair of interaction agents
    save_outputs(
        interactions_data,
        computed_interactions,
        frame_step,
        output_analysis_filepath,
    )

    # Remove the backup to avoid accidentally reusing it when the output file is deleted
    energies_backup.remove()


def get_pdb_frame_structures(
    structure: 'Structure', trajectory_filepath: str, snapshots: int, frames_limit: int
) -> tuple[Generator['Structure', None, None], int, int]:
    """Get the structure of every frame by writing and reading it as a PDB file.

    Return a generator of frame structures, the step between frames and the number of frames.
    """
    pdb_frames, frame_step, frame_count = get_pdb_frames(
        structure, trajectory_filepath, snapshots, frames_limit
    )
    # Each PDB file is read before the next frame is requested, which deletes it
    frame_structures = (Structure.from_pdb_file(pdb_frame) for pdb_frame in pdb_frames)
    return frame_structures, frame_step, frame_count


def get_coordinate_frame_structures(
    structure: 'Structure', trajectory_filepath: str, snapshots: int, frames_limit: int
) -> tuple[Generator['Structure', None, None], int, int]:
    """Get the structure of every frame by setting the frame coordinates in one copy of the structure.

    Return a generator of frame structures (the same copy, with the coordinates of each frame), the step between frames
    and the number of frames. Coordinates are rounded to 3 decimals (Å) as in a PDB file, since the frames were PDB files
    before and so the energies do not change.
    """
    coordinate_frames, frame_step, frame_count = get_coordinate_frames(
        structure, trajectory_filepath, snapshots, frames_limit
    )
    frame_structure = structure.copy()

    def frame_structures():
        for coordinates in coordinate_frames:
            frame_structure.set_new_coordinates(np.round(coordinates, 3).tolist())
            yield frame_structure
    return frame_structures(), frame_step, frame_count


def get_interactions_data(
    frame_structures: Generator['Structure', None, None],
    interactions: list[dict],
    engine: EnergiesEngine,
    energies_backup: File,
    backup_header: dict,
    verbose: bool = False,
) -> list[list[dict]]:
    """Get the atom energies of every interaction in every frame.

    Return one list per interaction with one dict per frame, with the atom energies of each agent.
    Each frame is appended to the backup when it is done, and frames already in the backup are not calculated again.
    """
    interactions_data = [[] for interaction in interactions]
    backup_size = 0
    if energies_backup.exists:
        interactions_data, backup_size = read_backup(energies_backup.path, backup_header, len(interactions))
        recovered_frames = len(interactions_data[0])
        if recovered_frames > 0:
            print(f' Recovering {recovered_frames} frames from the energies backup')
    # Drop anything after the last complete frame (e.g. a line cut by the interruption) before appending
    with open(energies_backup.path, 'r+b' if energies_backup.exists else 'wb') as backup:
        backup.truncate(backup_size)
        backup.seek(backup_size)
        if backup_size == 0:
            append_backup_line(backup, backup_header)
        for frame_number, current_frame_structure in enumerate(frame_structures):
            # If we already have this frame in the backup then skip it
            if frame_number < len(interactions_data[0]):
                continue
            # Run the main analysis over the current frame
            # Append the result data for each interaction
            # WARNING: At this point structure should be corrected
            # WARNING: Repeated atoms will make the analysis fail
            frame_energies_data = [
                get_interaction_frame_energies(
                    engine, current_frame_structure, interaction, verbose
                )
                for interaction in interactions
            ]
            for i, data in enumerate(frame_energies_data):
                interactions_data[i].append(data)
            # Save a backup just in case the process is interrupted further
            append_backup_frame(backup, frame_number, frame_energies_data)
    return interactions_data


def read_backup(
    backup_filepath: str, backup_header: dict, interaction_count: int
) -> tuple[list[list[dict]], int]:
    """Read the frames of an energies backup, up to the first incomplete or unexpected line.

    Return one list per interaction with one dict per frame, and the size in bytes of the lines read.
    A backup of a run with another setup (its first line differs from the header) gives no frames and size 0.
    """
    interactions_data = [[] for _ in range(interaction_count)]
    backup_size = 0
    with open(backup_filepath, 'rb') as backup:
        for line_number, line in enumerate(backup):
            # A line with no end was cut when the analysis was interrupted
            if not line.endswith(b'\n'):
                break
            try:
                record = json.loads(line)
            except ValueError:
                break
            if line_number == 0:
                if record != backup_header:
                    print(' The energies backup is from a run with another setup and it is not reused')
                    break
            elif (
                not isinstance(record, dict)
                or record.get('frame') != line_number - 1
                or len(record.get('interactions', [])) != interaction_count
            ):
                break
            else:
                for i, data in enumerate(record['interactions']):
                    interactions_data[i].append(data)
            backup_size += len(line)
    return interactions_data, backup_size


def append_backup_frame(backup: BinaryIO, frame_number: int, frame_energies_data: list[dict]):
    """Append the energies of one frame for every interaction to the backup."""
    append_backup_line(backup, {'frame': frame_number, 'interactions': frame_energies_data})


def append_backup_line(backup: BinaryIO, record: dict):
    """Append one record to the backup as a line of JSON, and flush it so it survives an interruption."""
    backup.write(json.dumps(record).encode() + b'\n')
    backup.flush()


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
    strong_bonds = interaction.get('strong_bonds', None) or []

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
        return {
            'agent1': agent1_atom_energies,
            'agent2': agent2_atom_energies,
            'agent1_zeroed': 0,
            'agent2_zeroed': 0,
        }

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

    # Set the energies of clashing atoms to 0, and count them
    agent1_atom_energies, agent1_zeroed = zero_clashes(agent1_atom_energies)
    agent2_atom_energies, agent2_zeroed = zero_clashes(agent2_atom_energies)

    # Print total energies at the end for every agent if the verbose flag has been passed
    if verbose:
        for agent_name, agent_energies in zip(
            [agent1.name, agent2.name], [agent1_atom_energies, agent2_atom_energies]
        ):
            total_vdw = total_es = total_both = 0
            for atom_energies in agent_energies:
                if atom_energies is None:
                    continue
                total_vdw += atom_energies[0]
                total_es += atom_energies[1]
                total_both += atom_energies[2]
            print(
                f' Total energies for {agent_name}: vdw {total_vdw}, es {total_es}, both {total_both}'
            )

    # Erase the 3 previous log lines
    print(ERASE_3_PREVIOUS_LINES)

    return {
        'agent1': agent1_atom_energies,
        'agent2': agent2_atom_energies,
        'agent1_zeroed': agent1_zeroed,
        'agent2_zeroed': agent2_zeroed,
    }


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
        # Add the number of atoms set to 0 in each frame (clashes)
        agent1_output['zeroed'] = [frame.get('agent1_zeroed', 0) for frame in data]
        agent2_output['zeroed'] = [frame.get('agent2_zeroed', 0) for frame in data]
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
    CMIP runs in the output directory, where its input and output files are written.
    """

    def __init__(
        self,
        output_directory: str,
        charges: list,
        cmip_structure: 'Structure',
        dummy_selection: 'Selection',
        debug: bool = False,
        extra_vdw_elements: set = set(),
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
        # Copy the CMIP inputs which do not change along the analysis in the output directory
        # Add the parameters of the elements missing from the van der Waals table
        write_cmip_vdw_table(
            f'{output_directory}/{CMIP_VDW_FILENAME}', extra_vdw_elements
        )
        copyfile(
            CMIP_INPUTS_CHECKONLY_SOURCE,
            f'{output_directory}/{CMIP_INPUTS_CHECKONLY_FILENAME}',
        )
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
        fortran_unit_files = glob(f'{self.output_directory}/fort.*')
        for filepath in fortran_unit_files:
            remove(filepath)


class GridFreeEngine:
    """Interaction energies with no grid: Coulomb with a sigmoidal dielectric and the van der Waals of CMIP.

    The electrostatics use the Mehler-Solmajer distance-dependent dielectric (bulk 78.4) with no cutoff: they have no
    low-dielectric cavity, so the shape of the agents only enters through distances.
    The van der Waals is CMIP's model: 4 eps [(s/r)^12 - (s/r)^6] with the parameters of the CMIP table, s scaled by 0.8
    and an 8 Å cutoff. As in CMIP, every atom gets the whole energy of its pairs with the other agent, atoms which are
    dummy in general contribute nothing and get 0, and atoms in strong bonds contribute nothing and get None.
    """

    def __init__(
        self,
        charges: list,
        cmip_structure: 'Structure',
        dummy_selection: 'Selection',
        extra_vdw_elements: set = set(),
    ):
        self.charges = np.array(charges, dtype=float)
        parameters = read_cmip_vdw_parameters(extra_vdw_elements)
        # Atoms with no parameters are never in an interaction, since interactions are checked first
        atom_parameters = [
            parameters.get(atom.element, (0.0, 0.0)) for atom in cmip_structure.atoms
        ]
        self.epsilons = np.array([epsilon for epsilon, _ in atom_parameters])
        self.sigmas = np.array([sigma for _, sigma in atom_parameters])
        self.dummy_atom_indices = set(dummy_selection.atom_indices)

    def __call__(
        self,
        frame_structure: 'Structure',
        agent1: Agent,
        agent2: Agent,
        strong_bonds: Optional[list],
    ) -> tuple[list, list]:
        """Get the energies of every atom of each agent with the other agent."""
        strong_bond_atom_indices = set(sum(strong_bonds or [], []))
        excluded_atom_indices = strong_bond_atom_indices | self.dummy_atom_indices
        indices_1 = np.array(agent1.atom_indices)
        indices_2 = np.array(agent2.atom_indices)
        coords_1 = np.array([frame_structure.atoms[i].coords for i in indices_1])
        coords_2 = np.array([frame_structure.atoms[i].coords for i in indices_2])
        active_1 = np.array([i not in excluded_atom_indices for i in indices_1])
        active_2 = np.array([i not in excluded_atom_indices for i in indices_2])
        vdw_1, es_1, vdw_2, es_2 = get_pair_energy_sums(
            coords_1[active_1],
            coords_2[active_2],
            self.charges[indices_1[active_1]],
            self.charges[indices_2[active_2]],
            self.epsilons[indices_1[active_1]],
            self.epsilons[indices_2[active_2]],
            self.sigmas[indices_1[active_1]],
            self.sigmas[indices_2[active_2]],
        )
        agent_energies = []
        for indices, active, vdw, es in [
            (indices_1, active_1, vdw_1, es_1),
            (indices_2, active_2, vdw_2, es_2),
        ]:
            energies = []
            values = iter(zip(vdw.tolist(), es.tolist()))
            for index, is_active in zip(indices, active):
                if is_active:
                    atom_vdw, atom_es = next(values)
                    energies.append((atom_vdw, atom_es, atom_vdw + atom_es))
                elif index in strong_bond_atom_indices:
                    energies.append(None)
                else:
                    energies.append((0.0, 0.0, 0.0))
            agent_energies.append(energies)
        return agent_energies[0], agent_energies[1]

    def cleanup(self):
        """Do nothing, since this engine leaves no files."""
        pass


def get_pair_energy_sums(
    coords_1: np.ndarray,
    coords_2: np.ndarray,
    charges_1: np.ndarray,
    charges_2: np.ndarray,
    epsilons_1: np.ndarray,
    epsilons_2: np.ndarray,
    sigmas_1: np.ndarray,
    sigmas_2: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Sum the van der Waals and electrostatic energies of the pairs of each atom with the atoms of the other group.

    Return the vdw and es sums of the atoms of group 1 and then those of group 2.
    Electrostatics use the Mehler-Solmajer sigmoidal dielectric, eps(r) = a + (bulk - a) / (1 + 7.7839 exp(-0.003627
    (bulk - a) r)) with a = -8.5525, for every pair. The sums run in a compiled numba kernel (energies_kernel), with
    threads as set by OMP_NUM_THREADS; numba is imported here, so only grid-free runs load it.
    """
    if len(coords_1) == 0 or len(coords_2) == 0:
        return np.zeros(len(coords_1)), np.zeros(len(coords_1)), np.zeros(len(coords_2)), np.zeros(len(coords_2))
    from mddb_workflow.analyses.energies_kernel import pair_energy_sums
    return pair_energy_sums(
        coords_1, coords_2, charges_1, charges_2, epsilons_1, epsilons_2, sigmas_1, sigmas_2,
        COULOMB_CONSTANT, SIGMOIDAL_BULK_DIELECTRIC, VDW_RADII_SCALE, VDW_CUTOFF,
        # Overlapping atoms get huge energies instead of infinite ones, so the clash rule sets them to 0
        minimum_distance=0.01,
    )


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
    For this reason dummy atoms keep their element, which sets their volume, so the cavity has the shape of both agents.
    By default we set as dummy host atoms involved in a covalent bond with guest atoms.
    Atoms which are dummy in general (e.g. virtual sites) have no volume, so they are written as hydrogens.
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
    general_dummy_atom_indices = set(dummy_selection.atom_indices)
    dummy_atom_indices = set(general_dummy_atom_indices)
    # If this is to be the host file then set guest atoms and strong bond atoms dummy
    if host_file:
        dummy_atom_indices |= set(strong_bond_atom_indices)
        dummy_atom_indices |= set(guest_selection.atom_indices)
    # Set which atoms must be kept in the output file
    # If this is the host file then keep both host and guest atoms
    # Host files keep the strong bond atoms of both agents as dummies, so the cavity is the same in both directions
    # If this is the guest file then keep only guest atoms, but remove its strong bond atoms
    selection = host_selection + guest_selection if host_file else guest_selection
    if not host_file:
        strong_bond_selection = frame_structure.select_atom_indices(
            strong_bond_atom_indices
        )
        selection -= guest_selection & strong_bond_selection
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
            # Dummy atoms (the other agent in host files, and atoms making a strong bond between both interacting agents)
            # get an 'X' before the element, so CMIP ignores their charge and their Van der Waals energy
            # Otherwise atoms in a strong bond would return high non-sense Van der Waals values
            real_index = selection.atom_indices[a]
            # Set if atom is to be ignored
            is_dummy = real_index in dummy_atom_indices
            cmip_dummy_flag = 'X' if is_dummy else ''
            # Keep the element of dummy atoms, which sets their volume in the low dielectric cavity
            # Atoms which are dummy in general have no volume, so they are set as hydrogens
            element = 'H' if real_index in general_dummy_atom_indices else atom.element
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
    directory: str,
) -> CompletedProcess:
    """Run CMIP with the given inputs, probe (guest) and host, and the CMIP van der Waals parameters.

    CMIP runs in the given directory, where all files must be, so the Fortran unit files it leaves (fort.*) stay there.
    """
    files = [cmip_inputs, probe, host, output_file, restart_file]
    for file in files:
        if abspath(file.basepath) != abspath(directory):
            raise ValueError(
                f'CMIP file {file.path} is not in the CMIP directory {directory}'
            )
    cmip_inputs, probe, host, output_file, restart_file = [
        file.filename for file in files
    ]
    return run(
        [
            'cmip',
            '-i',
            cmip_inputs,
            '-pr',
            probe,
            '-vdw',
            CMIP_VDW_FILENAME,
            '-hs',
            host,
            '-byat',
            output_file,
            '-rst',
            restart_file,
        ],
        stdout=PIPE,
        stderr=PIPE,
        cwd=directory,
    )


def get_cmip_grid(
    probe: File, host: File, output_file: File, restart_file: File, directory: str
) -> tuple[tuple, tuple, tuple]:
    """Run CMIP in 'checkonly' mode (i.e. start and stop) and mine the center, density and units of its grid."""
    cmip_inputs = File(f'{directory}/{CMIP_INPUTS_CHECKONLY_FILENAME}')
    cmip_process = run_cmip(
        cmip_inputs, probe, host, output_file, restart_file, directory
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
        agent1_cmip_pdb,
        agent2_cmip_pdb,
        cmip_checkonly_output,
        restart_file,
        output_directory,
    )
    agent2_center, agent2_density, agent2_units = get_cmip_grid(
        agent2_cmip_pdb,
        agent1_cmip_pdb,
        cmip_checkonly_output,
        restart_file,
        output_directory,
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
    """Coarsen the grid in case the number of grid points exceeds the GRID_POINTS_LIMIT, keeping its box.

    The grid unit size grows by the cube root of the excess, rounded up to 0.001 Å, and then by 0.001 Å until the grid
    fits, so the grid is only as coarse as needed. Return the new density and grid unit size.
    """

    def count_points(density: list) -> int:
        return (density[0] + 1) * (density[1] + 1) * (density[2] + 1)

    # Calculate the amount of total points with current parameters
    current_grid_points = count_points(density)
    if current_grid_points <= GRID_POINTS_LIMIT:
        return density, grid_unit_size
    print(f'WARNING: Grid points limit is exceeded ({current_grid_points})')
    box_size = [dimension * grid_unit_size for dimension in density]
    # Work in thousandths of Å to avoid rounding errors
    excess = current_grid_points / GRID_POINTS_LIMIT
    unit_size_milli = math.ceil(grid_unit_size * 1000 * excess ** (1 / 3))
    while True:
        new_density = [math.ceil(size * 1000 / unit_size_milli) for size in box_size]
        if count_points(new_density) <= GRID_POINTS_LIMIT:
            break
        unit_size_milli += 1
    new_grid_unit_size = unit_size_milli / 1000
    print(
        f'WARNING: Grid resolution has been reduced -> unit size = {new_grid_unit_size}'
    )
    return new_density, new_grid_unit_size


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
    """Leave everything ready to run CMIP by hand in the output directory, with a README, and stop the workflow.

    The van der Waals parameters are already in the output directory, where CMIP runs.
    """
    debug_script_source = File(ENERGIES_DEBUG_SCRIPT_SOURCE)
    # Copy in the energies folder a small python script used to sum output energies
    debug_script_filename = debug_script_source.filename
    debug_script = File(f'{output_directory}/{debug_script_filename}')
    copyfile(debug_script_source.path, debug_script.path)
    debug_vdw_filename = CMIP_VDW_FILENAME
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
    # Remove the output of a previous run, so it is never read if CMIP fails
    if cmip_output_file.exists:
        cmip_output_file.remove()
    # Run cmip
    cmip_process = run_cmip(
        cmip_inputs, guest, host, cmip_output_file, restart_file, output_directory
    )
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
    """Mine the Van der Waals (vdw), electrostatic (es) and total (both) energies of an atom from a CMIP output line."""
    # Numbers do not always fit in fixed columns, since numbers close to 0 are expressed in E notation
    # e.g. 0.1075E-04
    # For this reason we must mine these last values relying on whitespaces between them
    line_splits = line.split()
    vdw = float(line_splits[-3])
    es = float(line_splits[-2])
    both = float(line_splits[-1])
    return (vdw, es, both)


def zero_clashes(atom_energies: list) -> tuple[list, int]:
    """Set to 0 the energies of atoms whose total energy is above MAX_ATOM_ENERGY, which come from clashes.

    Return the new atom energies and the number of atoms set to 0. Atoms with no energies (None) are kept.
    """
    zeroed_count = 0
    new_atom_energies = []
    for atom_energies in atom_energies:
        if atom_energies is not None and atom_energies[2] > MAX_ATOM_ENERGY:
            atom_energies = (0, 0, 0)
            zeroed_count += 1
        new_atom_energies.append(atom_energies)
    return new_atom_energies, zeroed_count


def warn_zeroed_atoms(interactions: list[dict], interactions_data: list[list[dict]]):
    """Warn about the atoms set to 0 in each interaction, which mean overlapping atoms in a frame."""
    for interaction, data in zip(interactions, interactions_data):
        counts = [
            frame.get('agent1_zeroed', 0) + frame.get('agent2_zeroed', 0)
            for frame in data
        ]
        if sum(counts) == 0:
            continue
        frame_count = sum(count > 0 for count in counts)
        warn(
            f'{interaction["name"]}: {sum(counts)} atom energies above {MAX_ATOM_ENERGY} kcal/mol were set to 0 in '
            f'{frame_count} of {len(counts)} frames. This means overlapping atoms, which a relaxed MD frame should not have'
        )


def read_cmip_vdw_types(vdw_filepath: str) -> set[str]:
    """Read the atom types of a CMIP van der Waals parameters table."""
    with open(vdw_filepath, 'r') as file:
        lines = file.read().splitlines()[1:]
    return {line[:3].strip() for line in lines if line.strip()}


def write_cmip_vdw_table(output_filepath: str, extra_elements: set):
    """Write the CMIP van der Waals parameters table with the parameters of extra elements first.

    CMIP matches atom types as prefixes in the order of the table, so extra elements must go before the one-letter types
    (e.g. Ca before C). Parameters are written as in the table: epsilon and sigma = 2 R*/2 / 2^(1/6), with 2 decimals.
    """
    with open(output_filepath, 'w') as file:
        file.write('\n'.join(get_cmip_vdw_table_lines(extra_elements)))


def get_cmip_vdw_table_lines(extra_elements: set) -> list[str]:
    """Get the lines of the CMIP van der Waals parameters table, with the extra elements after the header."""
    with open(CMIP_VDW_SOURCE, 'r') as file:
        header, *lines = file.read().splitlines()
    extra_lines = []
    for element in sorted(extra_elements):
        half_rmin, epsilon, mass, atomic_number = CMIP_EXTRA_VDW_PARAMETERS[element]
        sigma = 2 * half_rmin / 2 ** (1 / 6)
        extra_lines.append(
            f'{element:<3}{epsilon:4.2f} {sigma:4.2f}{mass:6.2f}{atomic_number:4d}{0:7.3f}'
        )
    return [header] + extra_lines + lines


def read_cmip_vdw_parameters(extra_elements: set) -> dict[str, tuple[float, float]]:
    """Read the epsilon (kcal/mol) and sigma (Å) of each atom type in the CMIP table, as CMIP reads them."""
    parameters = {}
    for line in get_cmip_vdw_table_lines(extra_elements)[1:]:
        if line.strip():
            parameters.setdefault(
                line[:3].strip(), (float(line[3:7]), float(line[8:12]))
            )
    return parameters


def set_cmip_elements(structure: 'Structure'):
    """Set the atom elements in the CMIP format.

    Hydrogens bonded to carbons remain as 'H', hydrogens bonded to oxygen are renamed as 'HO' and hydrogens bonded to
    nitrogen or sulfur are renamed as 'HN'.
    Hydrogens bonded to other elements (e.g. P, B, Si or metals) are not hydrogen bond donors, so they remain as 'H', as
    hydrogens bonded to carbons do, with a warning.
    Elements in CMIP_ELEMENT_ALIASES (Se) take the element of their alias (S), with a warning.
    """
    other_bonded_elements = {}
    aliased_elements = {}
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
            bonded_atom_element = CMIP_ELEMENT_ALIASES.get(bonded_atom_element, bonded_atom_element)
            # Hydrogens bonded to oxygen are renamed as 'HO'
            if bonded_atom_element == 'O':
                element = 'HO'
            # Hydrogens bonded to nitrogen or sulfur are renamed as 'HN'
            elif bonded_atom_element == 'N' or bonded_atom_element == 'S':
                element = 'HN'
            # Hydrogens bonded to carbons or other elements remain as 'H'
            elif bonded_atom_element != 'C':
                other_bonded_elements[bonded_atom_element] = (
                    other_bonded_elements.get(bonded_atom_element, 0) + 1
                )
        elif element in CMIP_ELEMENT_ALIASES:
            aliased_elements[element] = aliased_elements.get(element, 0) + 1
            element = CMIP_ELEMENT_ALIASES[element]
        atom.element = element
    for element, count in sorted(aliased_elements.items()):
        warn(f'{count} {element} atoms take the CMIP type of {CMIP_ELEMENT_ALIASES[element]}')
    if other_bonded_elements:
        counts = ', '.join(
            f'{count} bonded to {element}'
            for element, count in sorted(other_bonded_elements.items())
        )
        warn(
            f'Hydrogens bonded to atoms other than C, N, O or S ({counts}) get the CMIP type of hydrogens bonded to carbon'
        )


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
    margin: float = GRID_MARGIN,
    grid_unit_size: float = GRID_UNIT_SIZE,
) -> tuple[list, list]:
    """Compute the center and density of a new grid covering the automatic grids of both agents.

    CMIP automatic grids span their center ± density × units / 2 (Å), where the density is the number of grid intervals.
    The new grid spans the union of both grids plus a margin on every side, with the given grid unit size.
    """
    new_center = []
    new_density = []
    for k in range(3):
        low = min(
            agent1_center[k] - agent1_density[k] * agent1_units[k] / 2,
            agent2_center[k] - agent2_density[k] * agent2_units[k] / 2,
        )
        high = max(
            agent1_center[k] + agent1_density[k] * agent1_units[k] / 2,
            agent2_center[k] + agent2_density[k] * agent2_units[k] / 2,
        )
        low, high = low - margin, high + margin
        new_center.append((low + high) / 2)
        new_density.append(math.ceil((high - low) / grid_unit_size))
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
