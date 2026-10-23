import re
import json
from os import listdir
from os.path import isabs, isfile, dirname, basename
from glob import glob
from subprocess import call
from tempfile import TemporaryDirectory
from mddb_workflow.mwf import Project, MD
from mddb_workflow.utils.auxiliar import InputError, parse_glob, unique, warn
from mddb_workflow.utils.constants import INCOMPLETE_PREFIX, STRUCTURE_FILENAME, TRAJECTORY_FILENAME, TRAJECTORY_SAMPLE_FILENAME
from mddb_workflow.utils.constants import CORRECTED_STRUCTURE, IMAGED_STRUCTURE, FILTERED_STRUCTURE, CONVERTED_STRUCTURE
from mddb_workflow.utils.constants import CORRECTED_TRAJECTORY, IMAGED_TRAJECTORY, FILTERED_TRAJECTORY, CONVERTED_TRAJECTORY
from mddb_workflow.utils.constants import MD_INPUT_STRUCTURE_FILEPATH, MD_INPUT_TRAJECTORY_FILEPATHS, MD_INPUT_TOPOLOGY_FILEPATH, MD_REMOVED_FLAG
from mddb_workflow.utils.constants import ACCEPTED_INPUT_FILENAMES, CACHE_FILENAME
from mddb_workflow.utils.file import File
from mddb_workflow.utils.gmx_spells import run_gromacs
from mddb_workflow.utils.vmd_spells import VMD_FORMAT, vmd_supported_trajectory_formats
from mddb_workflow.tools.atom_counters import get_topology_and_trajectory_atoms, get_structure_and_trajectory_atoms
from mddb_workflow.tools.guess_and_filter import guess_and_filter_topology, guess_and_filter_structure

# Structure formats which VMD can load, with the VMD name of each format
# Topologies have no coordinates but they are fine since a trajectory is always loaded
VMD_STRUCTURE_FORMATS = {
    'pdb': 'pdb',
    'gro': 'gro',
    'cif': 'pdbx',
    'prmtop': 'parm7',
    'psf': 'psf',
}

# Topology formats which VMD can load but have no coordinates
NO_COORDINATES_FORMATS = {'prmtop', 'psf'}

# Topology formats which VMD can not load but we can convert to gro
CONVERTIBLE_TOPOLOGY_FORMATS = {'tpr'}

# Formats which can be filtered when they have more atoms than the trajectory
FILTERABLE_TOPOLOGY_FORMATS = {'tpr', 'prmtop'}
FILTERABLE_STRUCTURE_FORMATS = {'pdb', 'gro'}

# Minimum size in Angstroms of the focused region, to avoid zooming too much on a few atoms
MIN_FOCUS_SIZE = 10

# Number of frames to be loaded when no frames are specified and the trajectory is big
DEFAULT_FRAMES_COUNT = 10

# Trajectories smaller than this size (in bytes) are fully loaded when no frames are specified
FULL_LOAD_SIZE_LIMIT = 500 * 1024 * 1024

# Directory where intermediate files are written while input files are processed
INPRO_DIRECTORY = INCOMPLETE_PREFIX + 'inpro'

# Files to be searched in the MD directory, from the most processed to the least processed
STRUCTURE_CANDIDATES = [
    STRUCTURE_FILENAME,
    f'{INPRO_DIRECTORY}/{CORRECTED_STRUCTURE}',
    f'{INPRO_DIRECTORY}/{IMAGED_STRUCTURE}',
    f'{INPRO_DIRECTORY}/{FILTERED_STRUCTURE}',
    f'{INPRO_DIRECTORY}/{CONVERTED_STRUCTURE}',
]
TRAJECTORY_CANDIDATES = [
    TRAJECTORY_FILENAME,
    TRAJECTORY_SAMPLE_FILENAME,
    f'{INPRO_DIRECTORY}/{CORRECTED_TRAJECTORY}',
    f'{INPRO_DIRECTORY}/{IMAGED_TRAJECTORY}',
    f'{INPRO_DIRECTORY}/{FILTERED_TRAJECTORY}',
    f'{INPRO_DIRECTORY}/{CONVERTED_TRAJECTORY}',
]

# Cksums of files are the file size and a hash of its content (e.g. 45652-096635543c4ba329)
FILE_CKSUM_PATTERN = re.compile(r'^(\d+)-')


def find_input_file(md: MD, input_filepath: str) -> File | None:
    """Find an input file path from the inputs file, either relative to the MD or the project directory."""
    if isabs(input_filepath):
        candidate_filepaths = [input_filepath]
    else:
        candidate_filepaths = [md.pathify(input_filepath), md.project.pathify(input_filepath)]
    for candidate_filepath in candidate_filepaths:
        for parsed_filepath in parse_glob(candidate_filepath):
            parsed_file = File(parsed_filepath)
            if parsed_file.exists: return parsed_file
    return None


def find_cached_file(md: MD, cksum: str | None) -> File | None:
    """Find a file in the MD or the project directory with the same cksum."""
    if not cksum: return None
    cksum_match = FILE_CKSUM_PATTERN.match(cksum)
    if not cksum_match: return None
    size = int(cksum_match.group(1))
    for directory in unique([md.directory, md.project.directory]):
        for filename in sorted(listdir(directory)):
            candidate_file = File(f'{directory}/{filename}')
            if not isfile(candidate_file.path) or candidate_file.size != size: continue
            if candidate_file.get_cksum(unsafe=True) == cksum: return candidate_file
    return None


def get_input_files(md: MD, md_inputs: dict) -> tuple[File | None, list[File | None], File | None]:
    """Get the input structure, trajectories and topology of a MD.

    Paths are read from the inputs file.
    If there is no inputs file then files are found by comparing their cksums with the ones in the cache.
    """
    if md.project.inputs.is_file_available:
        structure_filepath = md_inputs.get(MD_INPUT_STRUCTURE_FILEPATH)
        trajectory_filepaths = md_inputs.get(MD_INPUT_TRAJECTORY_FILEPATHS) or []
        if type(trajectory_filepaths) is str: trajectory_filepaths = [trajectory_filepaths]
        # The topology may be set for the MD or for the whole project
        # Note that a 'no' value is read as False
        topology_filepath = md.get_file_input(MD_INPUT_TOPOLOGY_FILEPATH)
        return (
            find_input_file(md, structure_filepath) if structure_filepath else None,
            [find_input_file(md, filepath) for filepath in trajectory_filepaths],
            find_input_file(md, topology_filepath) if topology_filepath else None,
        )
    cksums = md.cache.retrieve(MD.input_files_processing.cache_arg_cksums_key, {})
    # Trajectory cksums are saved as a JSON list, unless the list is too long
    trajectory_cksums = cksums.get('input_trajectory_files')
    trajectory_cksums = json.loads(trajectory_cksums) if trajectory_cksums and trajectory_cksums.startswith('[') else []
    return (
        find_cached_file(md, cksums.get('input_structure_file')),
        [find_cached_file(md, cksum) for cksum in trajectory_cksums],
        find_cached_file(md, cksums.get('input_topology_file')),
    )


def find_structure(md: MD, md_inputs: dict) -> File | None:
    """Find the most processed structure in a MD directory."""
    for candidate in STRUCTURE_CANDIDATES:
        candidate_file = File(md.pathify(candidate))
        if candidate_file.exists: return candidate_file
    input_structure_file, _, input_topology_file = get_input_files(md, md_inputs)
    if input_structure_file and input_structure_file.format in VMD_STRUCTURE_FORMATS:
        return input_structure_file
    supported_topology_formats = set(VMD_STRUCTURE_FORMATS) | CONVERTIBLE_TOPOLOGY_FORMATS
    if input_topology_file and input_topology_file.format in supported_topology_formats:
        return input_topology_file
    return None


def find_trajectory(md: MD, md_inputs: dict) -> File | None:
    """Find the most processed trajectory in a MD directory."""
    for candidate in TRAJECTORY_CANDIDATES:
        candidate_file = File(md.pathify(candidate))
        if candidate_file.exists: return candidate_file
    _, input_trajectory_files, _ = get_input_files(md, md_inputs)
    if not input_trajectory_files or None in input_trajectory_files: return None
    # When there are several trajectories only the first one is loaded
    if len(input_trajectory_files) > 1:
        print(f'There are {len(input_trajectory_files)} input trajectories. Only the first one will be loaded')
    input_trajectory_file = input_trajectory_files[0]
    if input_trajectory_file.format not in vmd_supported_trajectory_formats: return None
    return input_trajectory_file


def find_in_mds(project: Project, finder: callable) -> File | None:
    """Search a file in the reference MD first and then in the rest of MDs."""
    md_indices = list(range(len(project.mds)))
    md_indices.remove(project.reference_md_index)
    md_indices = [project.reference_md_index] + md_indices
    for md_index in md_indices:
        md_inputs = project.md_config[md_index]
        if md_inputs.get(MD_REMOVED_FLAG): continue
        found_file = finder(project.mds[md_index], md_inputs)
        if found_file: return found_file
    return None


def get_project(directory: str) -> Project:
    """Set the project in the given directory.
    If there is no inputs file then MD directories are those with a cache file.
    """
    inputs_filepaths = [filepath for filename in ACCEPTED_INPUT_FILENAMES for filepath in glob(f'{directory}/*{filename}')]
    if inputs_filepaths: return Project(directory=directory)
    md_directories = sorted(basename(dirname(filepath)) for filepath in glob(f'{directory}/*/{CACHE_FILENAME}'))
    if not md_directories:
        raise InputError(f'There is no inputs file and no MD directory with a cache file in {directory}')
    print(f'There is no inputs file. Using MD directories with a cache: {", ".join(md_directories)}')
    return Project(directory=directory, input_md_directories=md_directories)


def sort_input_files(input_filepaths: list[str]) -> tuple[File | None, File | None]:
    """Find out which input file is the structure and which is the trajectory from their formats.

    A single file is always the structure.
    If both files could be either the structure or the trajectory (e.g. two pdb files) then the input order is kept.
    """
    if len(input_filepaths) > 2:
        raise InputError(f'Up to 2 input files are expected (structure and trajectory) but {len(input_filepaths)} were passed')
    input_files = [File(filepath) for filepath in input_filepaths]
    for input_file in input_files:
        if not input_file.exists:
            raise InputError(f'File {input_file.path} does not exist')
    if len(input_files) == 0: return None, None
    supported_structure_formats = set(VMD_STRUCTURE_FORMATS) | CONVERTIBLE_TOPOLOGY_FORMATS
    if len(input_files) == 1:
        structure_file = input_files[0]
        if structure_file.format not in supported_structure_formats and structure_file.format in vmd_supported_trajectory_formats:
            raise InputError(f'{structure_file.path} is a trajectory. Pass also a structure, e.g. mwf view structure.pdb {structure_file.path}')
        return structure_file, None
    first_file, second_file = input_files
    def is_valid(structure_file: File, trajectory_file: File) -> bool:
        return structure_file.format in supported_structure_formats and trajectory_file.format in vmd_supported_trajectory_formats
    if is_valid(first_file, second_file): return first_file, second_file
    if is_valid(second_file, first_file): return second_file, first_file
    raise InputError(f'Could not find out which file is the structure and which is the trajectory: {first_file.path}, {second_file.path}\n'
        f'  Supported structure formats: {", ".join(sorted(supported_structure_formats))}\n'
        f'  Supported trajectory formats: {", ".join(sorted(vmd_supported_trajectory_formats))}')


def view_project(
    directory: str = '.',
    input_filepaths: list[str] = [],
    start: int | None = None,
    end: int | None = None,
    step: int | None = None,
    frames: list[int] = [],
    selection: str | None = None,
    all_frames: bool = False,
):
    """Open VMD with a structure and some frames of a trajectory.

    Input files are the structure and the trajectory, in any order.
    If no file is passed then the most processed structure and trajectory in the project directory are used.
    If only the structure is passed then it is loaded alone.
    If all frames are requested then the whole trajectory is loaded and other frame options are ignored.
    If no frames are specified then the whole trajectory is loaded if it is small, or only the first 10 frames otherwise.
    If only some of start/end/step are specified then the end is the last frame.
    The end frame is not included, as in the subset command.
    If a selection is passed then the view is focused on it and it is also shown in CPK.
    """
    structure_file, trajectory_file = sort_input_files(input_filepaths)
    # If only the structure is passed then it is loaded alone
    if not structure_file:
        project = get_project(directory)
        structure_file = find_in_mds(project, find_structure)
        if not structure_file:
            raise InputError('No structure was found in any MD directory')
        trajectory_file = find_in_mds(project, find_trajectory)
        if not trajectory_file:
            raise InputError('No trajectory was found in any MD directory')
    print(f'Structure: {structure_file.path}')
    if trajectory_file: print(f'Trajectory: {trajectory_file.path}')
    # Temporary files are removed once VMD is closed
    temporary_directory = TemporaryDirectory()
    # Remove atoms missing in the trajectory (e.g. water) from the structure, as it is done when checking the inputs
    if trajectory_file:
        structure_file = match_trajectory_atoms(structure_file, trajectory_file, temporary_directory.name)
    # Convert topologies which VMD can not load to gro
    if structure_file.format in CONVERTIBLE_TOPOLOGY_FORMATS:
        converted_structure_file = structure_file.get_reformated_file('gro').get_relocated_file(temporary_directory.name)
        print(f'Converting {structure_file.path} to gro')
        run_gromacs(f'editconf -f {structure_file.path} -o {converted_structure_file.path}')
        structure_file = converted_structure_file
    if structure_file.format not in VMD_STRUCTURE_FORMATS:
        raise InputError(f'Structure format "{structure_file.format}" is not supported')
    if not trajectory_file and structure_file.format in NO_COORDINATES_FORMATS:
        raise InputError(f'{structure_file.path} is a topology and it has no atom coordinates, so you would see nothing in VMD.\n'
            '  Pass also a trajectory to see the system, e.g. mwf view topology.prmtop trajectory.xtc')
    if trajectory_file and trajectory_file.format not in vmd_supported_trajectory_formats:
        raise InputError(f'Trajectory format "{trajectory_file.format}" is not supported')
    # Set the VMD commands to load the structure and the trajectory frames
    structure_vmd_format = VMD_STRUCTURE_FORMATS[structure_file.format]
    commands = [f'mol new {{{structure_file.path}}} type {structure_vmd_format} waitfor all']
    if trajectory_file:
        # Frames from the structure file are removed so frame numbers match the trajectory
        commands.append('animate delete all')
        trajectory_vmd_format = VMD_FORMAT[trajectory_file.format]
        load_command = f'mol addfile {{{trajectory_file.path}}} type {trajectory_vmd_format}'
        if all_frames:
            commands.append(f'{load_command} waitfor all')
        elif frames:
            for frame in sorted(frames):
                commands.append(f'{load_command} first {frame} last {frame} waitfor all')
        else:
            no_frame_options = start is None and end is None and step is None
            if no_frame_options and trajectory_file.size >= FULL_LOAD_SIZE_LIMIT:
                print(f'The trajectory is bigger than {FULL_LOAD_SIZE_LIMIT // (1024 * 1024)} MB. '
                    f'Only the first {DEFAULT_FRAMES_COUNT} frames are loaded. Use --all to load all frames')
                end = DEFAULT_FRAMES_COUNT
            if start is None: start = 0
            if step is None: step = 1
            if end is not None and end <= start:
                raise InputError('End frame must be posterior to start frame')
            # VMD includes the last frame while we do not
            last = -1 if end is None else end - 1
            commands.append(f'{load_command} first {start} last {last} step {step} waitfor all')
    # Show the PBC box, if there is any
    commands += [
        'package require pbctools',
        'if {[molinfo top get a] > 0} { pbc box }',
    ]
    # VMD does not return an error code, so selection errors are written in a file and VMD is closed
    selection_error_filepath = f'{temporary_directory.name}/selection_error.txt'
    if selection:
        commands += get_focus_commands(selection, selection_error_filepath)
    # Write the commands in a file and open VMD
    commands_filepath = f'{temporary_directory.name}/commands.vmd'
    with open(commands_filepath, 'w') as commands_file:
        commands_file.write('\n'.join(commands) + '\n')
    call(['vmd', '-e', commands_filepath])
    selection_error = None
    if File(selection_error_filepath).exists:
        with open(selection_error_filepath) as file:
            selection_error = file.read().strip()
    temporary_directory.cleanup()
    if selection_error: raise InputError(selection_error)


def match_trajectory_atoms(structure_file: File, trajectory_file: File, output_directory: str) -> File:
    """Filter the structure when it has more atoms than the trajectory.

    Typical atom selections (e.g. no water) are tried until one matches the trajectory atom count.
    The filtered structure is written in the output directory.
    """
    if structure_file.format in FILTERABLE_TOPOLOGY_FORMATS:
        structure_atom_count, trajectory_atom_count = get_topology_and_trajectory_atoms(structure_file, trajectory_file)
    elif structure_file.format in FILTERABLE_STRUCTURE_FORMATS:
        structure_atom_count, trajectory_atom_count = get_structure_and_trajectory_atoms(structure_file, trajectory_file)
    # Other formats are loaded as they are and VMD complains if atoms do not match
    else: return structure_file
    if structure_atom_count == trajectory_atom_count: return structure_file
    mismatch_message = ('Mismatch in the number of atoms between input files:\n' +
        f' Structure "{structure_file.path}" -> {structure_atom_count} atoms\n' +
        f' Trajectory "{trajectory_file.path}" -> {trajectory_atom_count} atoms')
    if structure_atom_count < trajectory_atom_count: raise InputError(mismatch_message)
    warn(mismatch_message)
    if structure_file.format in FILTERABLE_TOPOLOGY_FORMATS:
        filtered_structure_file = structure_file.get_prefixed_file('prefiltered.').get_relocated_file(output_directory)
        guessed = guess_and_filter_topology(structure_file, filtered_structure_file, trajectory_atom_count)
    else:
        # Gro files are converted to pdb since the structure filter can not read them
        if structure_file.format == 'gro':
            pdb_structure_file = structure_file.get_reformated_file('pdb').get_relocated_file(output_directory)
            run_gromacs(f'editconf -f {structure_file.path} -o {pdb_structure_file.path}')
            structure_file = pdb_structure_file
        filtered_structure_file = structure_file.get_prefixed_file('prefiltered.').get_relocated_file(output_directory)
        guessed = guess_and_filter_structure(structure_file, filtered_structure_file, trajectory_atom_count)
    if not guessed:
        raise InputError(mismatch_message + '\nCould not guess structure atom selection to match trajectory atoms count')
    print(f'Structure filtered to match the trajectory: {filtered_structure_file.path}')
    return filtered_structure_file


def get_focus_commands(selection: str, selection_error_filepath: str) -> list[str]:
    """Set the VMD commands to focus the view on a selection and show it in CPK.

    The view is centered in the selection and zoomed in as much as the selection is smaller than the whole system.
    """
    return [
        f'set focus_text {{{selection}}}',
        'if {[catch {atomselect top $focus_text} focus]} {',
        f'    set error_file [open {{{selection_error_filepath}}} w]',
        '    puts $error_file "Wrong selection \\"$focus_text\\": $focus"',
        '    close $error_file',
        '    quit',
        '}',
        'if {[$focus num] == 0} {',
        f'    set error_file [open {{{selection_error_filepath}}} w]',
        '    puts $error_file "Selection \\"$focus_text\\" matches no atoms"',
        '    close $error_file',
        '    quit',
        '}',
        'proc get_size {atoms} {',
        '    lassign [measure minmax $atoms] lowest highest',
        '    return [tcl::mathfunc::max {*}[vecsub $highest $lowest]]',
        '}',
        'display resetview',
        'molinfo top set center_matrix [list [transoffset [vecscale -1 [measure center $focus]]]]',
        f'set zoom [expr {{[get_size [atomselect top all]] / max([get_size $focus], {MIN_FOCUS_SIZE})}}]',
        'if {$zoom > 1} { scale by $zoom }',
        'mol representation CPK',
        'mol color Name',
        'mol selection $focus_text',
        'mol material Opaque',
        'mol addrep top',
    ]
