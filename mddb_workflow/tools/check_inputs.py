from mddb_workflow.utils.auxiliar import InputError
from mddb_workflow.utils.auxiliar import warn, MISSING_TOPOLOGY
from mddb_workflow.utils.auxiliar import is_standard_topology
from mddb_workflow.utils.pyt_spells import find_all_corrupted_frames
from mddb_workflow.utils.file import File
from mddb_workflow.utils.cache import get_cached_function
from mddb_workflow.utils.constants import TOPOLOGY_SUPPORTED_FORMATS, TRAJECTORY_SUPPORTED_FORMATS
from mddb_workflow.utils.constants import STRUCTURE_SUPPORTED_FORMATS
from mddb_workflow.utils.type_hints import *

from mddb_workflow.tools.atom_counters import get_topology_and_trajectory_atoms
from mddb_workflow.tools.atom_counters import get_structure_atoms
from mddb_workflow.tools.atom_counters import get_structure_and_trajectory_atoms
from mddb_workflow.tools.guess_and_filter import guess_and_filter_topology

import re
from scipy.io import netcdf_file


# Set some known message errors
NETCDF_DTYPE_ERROR = 'When changing to a larger dtype, its size must be a divisor of the total size in bytes of the last axis of the array.'
PSF_HEADER_PATTERN = r'^\s*[0-9]* ![A-Z]*'
EMPTY_LINE = '\n'

# Set exceptions for fixes applied from here
FIXED_TOPOLOGY_EXCEPTION = Exception('Fixed topology')


def check_input_files(
    input_structure_file: 'File',
    input_trajectory_files: list['File'],
    input_topology_file: Union['File', Exception],
    cache : 'Cache') -> dict:
    """Check input files coherence and integrity.
    If there is any problem then raises an input error.
    Some exceptional problems may be fixed from here.
    In these cases, both the exception and the modified file are returned in a final dict.
    """
    # Set the exceptions dict to be returned at the end
    exceptions = {}

    # Set the topology file to be checked
    # Note that this variable may be reassigned as fixes are applied
    topology_file = input_topology_file

    # Set the structure file to be checked
    # Note that the structure file may be the topology file so it may be reassigned as well
    structure_file = input_structure_file

    # Get a sample trajectory file and then check its format
    # All input trajectory files must have the same format
    trajectory_sample = input_trajectory_files[0]

    # Check input files are supported by the workflow
    if topology_file != MISSING_TOPOLOGY and not is_standard_topology(topology_file) and topology_file.format not in TOPOLOGY_SUPPORTED_FORMATS:
        if topology_file.format in {'pdb', 'gro'}:
            raise InputError('A structure file is not supported as topology anymore. If there is no topology then use the argument "-top no"')
        raise InputError(f'Topology {topology_file.path} has a not supported format. Try one of these: {", ".join(TOPOLOGY_SUPPORTED_FORMATS)}')
    if trajectory_sample.format not in TRAJECTORY_SUPPORTED_FORMATS:
        raise InputError(f'Trajectory {trajectory_sample.path} has a not supported format. Try one of these: {", ".join(TRAJECTORY_SUPPORTED_FORMATS)}')
    if structure_file.format not in STRUCTURE_SUPPORTED_FORMATS:
        raise InputError(f'Structure {structure_file.path} has a not supported format. Try one of these: {", ".join(STRUCTURE_SUPPORTED_FORMATS)}')

    # Make sure the trajectory file is not corrupted
    # NetCDF files may have a variety of problems which are specific to this format
    if trajectory_sample.format == 'nc':
        # Check if reading the trajectory raises the following error
        # ValueError: When changing to a larger dtype, its size must be a divisor of the total size in bytes of the last axis of the array.
        # This error may happen with NetCDF files and it is a bit shady
        # Some tools may be able to read the first frames of the corrupted file: VMD and pytraj
        # Some other tools will instantly fail to read it: MDtraj and MDAnalysis
        try:
            # Iterate trajectory files
            for trajectory_file in input_trajectory_files:
                # This does not read the whole trajectory
                netcdf_file(trajectory_file.path, 'r')
        except Exception as error:
            # If the error message matches with a known error then report the problem
            error_message = str(error)
            if error_message == NETCDF_DTYPE_ERROR:
                warn(f'Corrupted trajectory file {trajectory_file.path}')
        # Get a cached version of the function to check for frame corruption
        cached_find_all_corrupted_frames = get_cached_function(find_all_corrupted_frames, cache)
        for trajectory_file in input_trajectory_files:
            pytraj_input_topology = topology_file if topology_file != MISSING_TOPOLOGY else structure_file
            # Get corrupted frames in the trajectory
            print(f'Checking NetCDF integrity in {trajectory_file.path}')
            corrupted_frames = cached_find_all_corrupted_frames(
                pytraj_input_topology.path, trajectory_file.path)
            corrupted_frames_count = len(corrupted_frames)
            # Add some extra spaces to hide the previous log
            print(f' Found {corrupted_frames_count} corrupted frames.                          ')
            if corrupted_frames_count > 0:
                print(f' First corrupted frame: {corrupted_frames[0]}')
                print(f' However some tools may still be able to read a part or the total of it.')
                raise InputError(f'Corrupted input trajectory file {trajectory_file.path}.')

    # Make sure the topology file is well formated
    # Check a specific problem affecting some PSF topologies
    if topology_file != MISSING_TOPOLOGY and topology_file.format == 'psf':
        # Set the output fixed topology file, in case it is to be created
        fixed_topology_file = topology_file.get_prefixed_file('fixed.')
        had_problem = check_and_fix_psf(topology_file, fixed_topology_file)
        # If a problem was found then report the problem and save the exception
        if had_problem:
            print(f'The input topology had format problem but it has been fixed in {fixed_topology_file.path}')
            exceptions[FIXED_TOPOLOGY_EXCEPTION] = fixed_topology_file
            # From now on use the fixed topology as the topology
            if structure_file == topology_file:
                structure_file = fixed_topology_file
            topology_file = fixed_topology_file

    # Get topology and trajectory atom counts
    topology_atom_count, trajectory_atom_count = get_topology_and_trajectory_atoms(topology_file, trajectory_sample)

    # If we have the trajectory atom count then it means we had a valid topology
    if trajectory_atom_count is not None:

        # Make sure their atom counts match
        if topology_atom_count != trajectory_atom_count:
            warn('Mismatch in the number of atoms between input files:\n' +
                f' Topology "{topology_file.path}" -> {topology_atom_count} atoms\n' +
                f' Trajectory "{trajectory_sample.path}" -> {trajectory_atom_count} atoms')
            if topology_atom_count < trajectory_atom_count:
                raise InputError('Trajectory has more atoms than topology, there is no way to fix this.')
            # If the topology has more atoms than the trajectory however we may attempt to guess
            # If we guess which atoms are the ones in the trajectory then we can filter the topology
            else:
                prefiltered_topology_file = topology_file.get_prefixed_file('prefiltered.')
                guessed = guess_and_filter_topology(
                    topology_file,
                    prefiltered_topology_file,
                    trajectory_atom_count)
                # Save the new topology file in the exceptions
                if guessed:
                    exceptions[FIXED_TOPOLOGY_EXCEPTION] = prefiltered_topology_file
                    # From now on use the prefiltered topology as the topology
                    if structure_file == topology_file:
                        structure_file = prefiltered_topology_file
                    topology_file = prefiltered_topology_file
                    topology_atom_count = trajectory_atom_count
                else: raise InputError('Could not guess topology atom selection to match trajectory atoms count')

        # If the topology file is already the structure file then there is no need to check it
        if structure_file == topology_file:
            print(f'Topology and trajectory files match in number of atoms: {trajectory_atom_count}')
            return exceptions

        # If the counts match then also get the structure atom count and compare
        structure_atom_count = get_structure_atoms(structure_file)

        # Make sure it matches the topology and trajectory atom count
        if topology_atom_count != structure_atom_count:
            raise InputError('Mismatch in the structure input file number of atoms:\n'+
                f' Topology and trajectory -> {topology_atom_count} atoms\n' +
                f' Structure "{structure_file.path}" -> {structure_atom_count} atoms')

        # If we reached this point then it means everything is matching
        print(f'All input files match in number of atoms: {trajectory_atom_count}')
        return exceptions

    # Otherwise it means we had not a valid topology file
    # We must use the structure to find trajectory atoms
    structure_atom_count, trajectory_atom_count = get_structure_and_trajectory_atoms(structure_file, trajectory_sample)

    # Make sure their atom counts match
    if structure_atom_count != trajectory_atom_count:
        raise InputError('Mismatch in the number of atoms between input files:\n' +
            f' Structure "{structure_file.path}" -> {structure_atom_count} atoms\n' +
            f' Trajectory "{trajectory_sample.path}" -> {trajectory_atom_count} atoms')

    # If we have a number of topology atoms then make sure it matches the structure and trajectory atoms
    # This may happen if the topology is our standard topology file instead of a valid topology
    if topology_atom_count is not None and topology_atom_count != trajectory_atom_count:
        raise InputError('Mismatch in the number of atoms between input files:\n' +
            f' Structure and trajectory -> {trajectory_atom_count} atoms\n' +
            f' Topology "{topology_file.path}" -> {topology_atom_count} atoms')

    # If we made it this far it means all checkings are good
    print(f'Input files match in number of atoms: {trajectory_atom_count}')
    return exceptions

def check_and_fix_psf(input_topology_file: 'File', output_topology_file: 'File') -> bool:
    """Check if a PSF topology is properly formatted in the headers.
    Wrong-formmated topologies may still valid for NAMD tools and pytraj.
    However they failed to be read by both MDtraj and MDAnalysis.
    Pytraj may write PSF files with this problem, for instance.
    If so, an output topology file will be created with the problem fixed.
    """
    # Read the content of the PSF file
    psf_content = None
    with open(input_topology_file.path, 'r') as file:
        psf_content = file.readlines()
    # Keep track of if a problem was found
    problem = False
    # Iterate lines in the file conetent
    for l, line in enumerate(psf_content):
        # We only care about headers
        if not re.search(PSF_HEADER_PATTERN, line): continue
        # Make sure the line previous to a header is an empty line
        last_line = psf_content[l-1]
        if last_line == EMPTY_LINE: continue
        # Otherwise we must fix it
        problem = True
        psf_content.insert(l, EMPTY_LINE)
    # If a problem was found then write a new fixed topology
    if not problem: return False
    # Write the fixed content to  a new topology file
    with open(output_topology_file.path, 'w') as file:
        file.writelines(psf_content)
