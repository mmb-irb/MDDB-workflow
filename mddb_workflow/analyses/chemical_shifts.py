# Chemical shifts analysis
#
# The chemical shifts analysis is carried by the 'LEGOLAS' software. LEGOLAS predicts the NMR chemical
# shifts of protein backbone atoms (H, HA, N, CA, CB and C) using an ensemble of 5 neural networks
# trained over atomic environment vectors (AEVs). The chemical shift is predicted for every frame and
# the standard deviation between the 5 models is returned as well as a measure of the uncertainty.
#
# https://github.com/roitberg-group/legolas
#
# This analysis was integrated by Claude

from os import mkdir
from os.path import exists
from shutil import rmtree
from subprocess import run, PIPE, STDOUT
from pathlib import Path
import csv
import json
import re
import sys
from statistics import mean, pstdev

import mdtraj as mdt

from mddb_workflow.tools.get_reduced_trajectory import get_reduced_trajectory
from mddb_workflow.utils.auxiliar import ToolError, save_json, warn
from mddb_workflow.utils.constants import GREY_HEADER, COLOR_END
from mddb_workflow.utils.constants import OUTPUT_CHEMICAL_SHIFTS_FILENAME
from mddb_workflow.utils.type_hints import *

# DANI: Esto es provisional hasta que LEGOLAS esté en conda
LEGOLAS_COMMAND = 'source ~/miniforge3/bin/activate ~/miniforge3/envs/legolas; python ~/libraries/legolas/test/legolas.py {} -t {} -o csv'

# Atom types predicted by LEGOLAS
# LEGOLAS sets the atom type by removing digits from the atom name (e.g. HA2 -> HA, H1 -> H)
LEGOLAS_ATOM_TYPES = {'H', 'HA', 'N', 'CA', 'CB', 'C'}

# Residue names supported by LEGOLAS
# Water (HOH) is also supported but we do not include it since it would make LEGOLAS much slower
LEGOLAS_RESIDUE_NAMES = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'}
# Non-standard residue names which may be translated to a supported residue name
RESIDUE_NAME_ALIASES = {
    'HID': 'HIS', 'HIE': 'HIS', 'HIP': 'HIS', 'HSD': 'HIS', 'HSE': 'HIS', 'HSP': 'HIS',
    'CYX': 'CYS', 'CYM': 'CYS', 'ASH': 'ASP', 'GLH': 'GLU', 'LYN': 'LYS',
}

# Number of decimals to keep in the output
OUTPUT_DECIMALS = 3


def chemical_shifts(
    structure_file: 'File',
    trajectory_file: 'File',
    structure: 'Structure',
    pbc_selection: 'Selection',
    cg_selection: 'Selection',
    snapshots: int,
    output_directory: str,
    frames_limit: int = 1000,
):
    """Perform the chemical shifts analysis using LEGOLAS."""
    warn('The "chemical shifts" analysis is not yet fully integrated. LEGOLAS is not installed in the conda enviornment')

    # Set the residues to be analyzed: protein residues with a residue name supported by LEGOLAS
    # Residues in PBC and coarse grain residues are excluded
    excluded_residue_indices = set(structure.get_selection_residue_indices(pbc_selection + cg_selection))
    protein_residue_indices = sorted(structure.get_selection_residue_indices(structure.select_protein()))
    legolas_residue_indices = []
    for residue_index in protein_residue_indices:
        if residue_index in excluded_residue_indices: continue
        residue = structure.residues[residue_index]
        if get_legolas_residue_name(residue.name) is None:
            warn(f'Residue {residue} is not supported by LEGOLAS and it will be skipped')
            continue
        legolas_residue_indices.append(residue_index)
    if len(legolas_residue_indices) == 0:
        print(' No protein residues to analyze')
        return

    # Set the LEGOLAS working directory, where all LEGOLAS inputs and outputs will be
    legolas_directory = f'{output_directory}/legolas'
    if not exists(legolas_directory): mkdir(legolas_directory)

    # Write a structure with only the analyzed residues
    # Residues are renamed to the names supported by LEGOLAS
    # Residues are renumbered so every residue has a unique number in LEGOLAS outputs
    # Note that LEGOLAS sorts the output by residue number, so this also guarantees the order is kept
    legolas_selection = structure.select_residue_indices(legolas_residue_indices)
    legolas_structure = structure.filter(legolas_selection)
    for r, residue in enumerate(legolas_structure.residues):
        residue.name = get_legolas_residue_name(residue.name)
        residue.number = r + 1
        residue.icode = ''
    legolas_structure_filepath = f'{legolas_directory}/{structure_file.filename}'
    legolas_structure.generate_pdb_file(legolas_structure_filepath)

    # Write a reduced trajectory with only the analyzed residues atoms
    reduced_trajectory_filepath, _, _ = get_reduced_trajectory(
        trajectory_file,
        snapshots,
        frames_limit,
    )
    print(' Filtering trajectory for LEGOLAS')
    legolas_trajectory = mdt.load(reduced_trajectory_filepath, top=structure_file.path,
        atom_indices=legolas_selection.atom_indices)
    legolas_trajectory_filepath = f'{legolas_directory}/{trajectory_file.filename}'
    legolas_trajectory.save(legolas_trajectory_filepath)

    # Run LEGOLAS
    print(' Running LEGOLAS')
    print(GREY_HEADER, end='')
    command = LEGOLAS_COMMAND.format(trajectory_file.filename, structure_file.filename)
    process = run(['bash', '-c', command], cwd=legolas_directory, stdout=PIPE, stderr=STDOUT)
    logs = process.stdout.decode()
    print(COLOR_END, end='')
    # LEGOLAS writes the output in the working directory and names it after the trajectory
    legolas_output_filepath = f'{legolas_directory}/{Path(trajectory_file.filename).stem}_cs.csv'
    if process.returncode != 0 or not exists(legolas_output_filepath):
        print(logs)
        raise ToolError('Something went wrong with LEGOLAS')

    # Set the atoms we expect LEGOLAS to predict, grouped by atom type
    # Atoms are sorted by residue number and then by atom index, just like in the LEGOLAS output
    expected_atoms = { atom_type: [] for atom_type in LEGOLAS_ATOM_TYPES }
    for atom in sorted(legolas_structure.atoms, key=lambda atom: (atom.residue.number, atom.index)):
        atom_type = get_legolas_atom_type(atom.name)
        if atom_type not in LEGOLAS_ATOM_TYPES: continue
        expected_atoms[atom_type].append(atom)

    # Read the LEGOLAS output
    # There is a row per atom and chemical shifts are lists with a value per frame
    # Rows are sorted by residue number and atom type (alphabetically)
    # Atoms with the same residue and atom type (e.g. H1, H2 and H3) are sorted by atom index
    # Note that cells are huge so we must increase the csv field size limit
    csv.field_size_limit(sys.maxsize)
    with open(legolas_output_filepath, 'r') as file:
        legolas_output = list(csv.DictReader(file))
    # Note that atom indices here are still the indices in the LEGOLAS structure
    legolas_atom_indices = []
    atom_type_counters = { atom_type: 0 for atom_type in LEGOLAS_ATOM_TYPES }
    for row in legolas_output:
        atom_type, residue_name = row['ATOM_TYPE'], row['RES_TYPE']
        residue_number = int(row['SEQ_ID'])
        atom_type_atoms = expected_atoms[atom_type]
        atom_type_counter = atom_type_counters[atom_type]
        if atom_type_counter >= len(atom_type_atoms):
            raise ToolError(f'LEGOLAS returned more {atom_type} atoms than expected')
        expected_atom = atom_type_atoms[atom_type_counter]
        expected_residue = expected_atom.residue
        if residue_number != expected_residue.number or residue_name != expected_residue.name:
            raise ToolError(f'LEGOLAS atom {atom_type} in residue {residue_name} {residue_number} does not'
                f' match atom {expected_atom.name} in residue {expected_residue.name} {expected_residue.number}')
        legolas_atom_indices.append(expected_atom.index)
        atom_type_counters[atom_type] += 1
    for atom_type, atom_type_atoms in expected_atoms.items():
        if atom_type_counters[atom_type] != len(atom_type_atoms):
            raise ToolError(f'LEGOLAS returned {atom_type_counters[atom_type]} {atom_type} atoms'
                f' but we expected {len(atom_type_atoms)}')
    # Map LEGOLAS structure atom indices to the original structure atom indices
    # Sort atoms by atom index
    atom_indices = [legolas_selection.atom_indices[atom_index] for atom_index in legolas_atom_indices]
    atom_order = sorted(range(len(atom_indices)), key=lambda row: atom_indices[row])

    # Parse chemical shifts and summarize them along frames
    # Note that LEGOLAS returns a single value instead of a list when there is only one frame
    def parse_values(text: str) -> list[float]:
        values = json.loads(text)
        if type(values) != list: values = [values]
        return values
    def round_values(values: list[float]) -> list[float]:
        return [round(value, OUTPUT_DECIMALS) for value in values]
    chemical_shift_averages = []
    chemical_shift_deviations = []
    model_deviation_averages = []
    for row in legolas_output:
        chemical_shifts = parse_values(row['CHEMICAL_SHIFT'])
        chemical_shift_averages.append(mean(chemical_shifts))
        chemical_shift_deviations.append(pstdev(chemical_shifts))
        # Note that the deviation returned by LEGOLAS is the standard deviation between its 5 models
        # i.e. it is a measure of the prediction uncertainty, not of the fluctuation along the trajectory
        model_deviations = parse_values(row['CHEMICAL_SHIFT_STD'])
        model_deviation_averages.append(mean(model_deviations))

    # Set the output analysis
    # Note that csdv is the deviation along frames while csmdv is the average deviation between LEGOLAS models
    output_analysis = {
        'atom_indices': [atom_indices[row] for row in atom_order],
        'csav': round_values([chemical_shift_averages[row] for row in atom_order]),
        'csdv': round_values([chemical_shift_deviations[row] for row in atom_order]),
        'csmdv': round_values([model_deviation_averages[row] for row in atom_order]),
    }
    save_json(output_analysis, f'{output_directory}/{OUTPUT_CHEMICAL_SHIFTS_FILENAME}')

    # Remove LEGOLAS files since they are heavy and we already have what we need
    rmtree(legolas_directory)


def get_legolas_residue_name(residue_name: str) -> Optional[str]:
    """Get the residue name as LEGOLAS expects it, or None if it is not supported."""
    residue_name = RESIDUE_NAME_ALIASES.get(residue_name, residue_name)
    if residue_name not in LEGOLAS_RESIDUE_NAMES: return None
    return residue_name


def get_legolas_atom_type(atom_name: str) -> str:
    """Get the atom type as LEGOLAS sets it: the atom name without digits."""
    return re.sub(r'\d', '', atom_name)
