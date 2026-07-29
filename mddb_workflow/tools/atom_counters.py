import re

import mdtraj as mdt
import pytraj as pyt
import MDAnalysis as mda

from mddb_workflow.utils.auxiliar import InputError, ToolError
from mddb_workflow.utils.auxiliar import warn, CaptureOutput, load_json, MISSING_TOPOLOGY
from mddb_workflow.utils.auxiliar import is_standard_topology, get_auxiliar_filepath
from mddb_workflow.utils.gmx_spells import run_gromacs, mine_system_atoms_count
from mddb_workflow.utils.gmx_spells import get_tpr_atom_count
from mddb_workflow.utils.gmx_spells import get_gmx_atom_count, get_gmx_trajectory_atom_count
from mddb_workflow.utils.gmx_spells import GROMACS_SUPPORTED_TRAJECTORY_FORMATS
from mddb_workflow.utils.vmd_spells import vmd_to_pdb, get_rst7_atom_count
from mddb_workflow.utils.structures import Structure
from mddb_workflow.utils.file import File
from mddb_workflow.utils.constants import GROMACS_TRAJECTORY_SUPPORTED_FORMATS
from mddb_workflow.utils.type_hints import *

# Set some known message errors usefult to mine the actual atom counts
MDTRAJ_ATOM_MISMATCH_ERROR = r'xyz must be shape \(Any, ([0-9]*), 3\). You supplied  \(1, ([0-9]*), 3\)'
MDTRAJ_INSERTION_CODES_ERROR = r'^Could not convert residue number \[[0-9]*[a-zA-Z]\]$'
PYTRAJ_XTC_ATOM_MISMATCH_ERROR = r'Error: # atoms in XTC file \(([0-9]*)\) does not match # atoms in (topology|parm) [\w.-]* \(([0-9]*)\)'
GROMACS_ATOM_MISMATCH_ERROR = r'is larger than the number of atoms in the\ntrajectory file \(([0-9]*)\). There is a mismatch in the contents'
GROMACS_ATOM_COUNT_CHECK = r'# Atoms  ([0-9]*)'

def get_topology_atoms(topology_file: 'File') -> Optional[int]:
    """Get the atoms count from a topology alone."""
    # To do so rely on different tools depending on the topology format
    # If there is no topology file then just return None
    if topology_file == MISSING_TOPOLOGY: return None
    # If it is our standard topology then simply count the atom names
    if is_standard_topology(topology_file):
        # Parse the json and count atoms
        parsed_topology = load_json(topology_file.path)
        topology_atom_count = len(parsed_topology['atom_names'])
        return topology_atom_count
    # For a TPR use Gromacs, which is its native tool
    if topology_file.format == 'tpr':
        return get_tpr_atom_count(topology_file.path)
    # For .top files we use PyTraj since MDtraj can not handle it
    if topology_file.format == 'top':
        return get_topology_atoms_pytraj(topology_file)
    # If the topology is a restart file then MDtraj will not be able to read it
    # Use VMD to count its atoms
    if topology_file.format == 'rst7':
        return get_rst7_atom_count(topology_file.path)
    # At this point the topology should be supported by MDtraj
    try:
        topology = mdt.load_topology(topology_file.path)
        return topology.n_atoms
    except Exception as error:
        # If the error message matches with a known error then report the problem
        error_message = str(error)
        error_match = re.match(MDTRAJ_INSERTION_CODES_ERROR, error_message)
        if error_match:
            warn('The input topology has insertion codes.\n'
            ' Some tools may crash when reading the topology (MDtraj).\n'
            ' Some tools may ignore insertion codes when reading the topology (MDAnlysis, PyTraj, VMD).')
            # Use other tool to read the topology
            # Other tools could ignore the inserion codes
            # However this is not a problem here, where we only care about the number of atoms
            return get_topology_atoms_pytraj(topology_file)
        # If we do not know the error then raise it as is
        raise error

def get_topology_atoms_pytraj(topology_file: 'File') -> int:
    """Get the atoms count from a topology alone using PyTraj"""
    topology = pyt.load_topology(topology_file.path)
    return topology.n_atoms

def get_trajectory_atoms(trajectory_file: 'File') -> int:
    """Get the atoms count from a trajectory alone"""
    # We will used different tools depending on the trajectory format
    # Use Gromacs for Gromacs trajectories
    if trajectory_file.format in GROMACS_SUPPORTED_TRAJECTORY_FORMATS:
        return get_gmx_trajectory_atom_count(trajectory_file)
    # Use MDAnalysis for anything else
    # Note that MDAnalysis could also handle Gromacs trajectories actually
    universe = mda.Universe(trajectory_file.path)
    return universe.trajectory.n_atoms

def get_topology_and_trajectory_atoms_pytraj(topology_file: 'File', trajectory_file: 'File') -> tuple[int, int]:
    """Get the atoms count from topology and trajectory together using pytraj.
    This is an altermative method used when MDtraj can not handle it.
    """
    # Note that calling iterload will print a error log when atoms do not match but will not raise a proper error
    # To capture the error log we must throw this command wrapped in a stdout redirect
    trajectory = None
    with CaptureOutput('stderr') as output:
        trajectory = pyt.iterload(trajectory_file.path, top=topology_file.path)
    logs = output.captured_text
    error_match = re.match(PYTRAJ_XTC_ATOM_MISMATCH_ERROR, logs)
    if error_match:
        topology_atom_count = int(error_match[3])
        trajectory_atom_count = int(error_match[1])
    # Now obtain the number of atoms from the frame we just read
    else:
        topology_atom_count = trajectory_atom_count = trajectory.n_atoms
    return topology_atom_count, trajectory_atom_count


def get_topology_and_trajectory_atoms(topology_file: 'File', trajectory_file: 'File') -> tuple[int, int]:
    """Get the atoms count from topology and trajectory together."""
    # To do so rely on different tools depending on the topology format
    # If there is no topology file then just compare strucutre and trajectory an exit
    if topology_file == MISSING_TOPOLOGY:
        # We do not have a topology atom count to return
        # Without a valid topology we can not count trajectory atoms either
        return None, None
    # If it is our standard topology then simply count the atom names
    # Get trajectory atoms using the structure instead
    if is_standard_topology(topology_file):
        # Parse the json and count atoms
        parsed_topology = load_json(topology_file.path)
        topology_atom_count = len(parsed_topology['atom_names'])
        # Without a valid topology we can not count trajectory atoms
        return topology_atom_count, None
    # For a TPR use Gromacs, which is its native tool
    if topology_file.format == 'tpr':
        # Make sure the trajectory is compatible with gromacs
        if trajectory_file.format not in GROMACS_TRAJECTORY_SUPPORTED_FORMATS:
            raise InputError('Why loading a TPR topology with a non-gromacs trajectory?')
        # Run Gromacs just to generate a structure using all atoms in the topology and coordinates in the first frame
        # If trajectory atoms are fewer than topology atoms then we will see a specific error
        output_sample_gro_file = File('.sample.gro')
        output_logs, error_logs = run_gromacs(f'trjconv -s {topology_file.path} \
            -f {trajectory_file.path} -o {output_sample_gro_file.path} -dump 0',
            user_input='System', expected_output_filepath=None)
        # Always get error logs and mine topology atoms
        # Note that these logs include the output selection request from Gromacs
        # This log should be always there, even if there was a mismatch and then Gromacs failed
        topology_atom_count = mine_system_atoms_count(error_logs)
        # If the output does not exist at this point it means something went wrong with gromacs
        if not output_sample_gro_file.exists:
            # Check if we know the error
            error_match = re.search(GROMACS_ATOM_MISMATCH_ERROR, error_logs)
            if error_match:
                # Get the trajectory atom count
                trajectory_atom_count = int(error_match[1])
                return topology_atom_count, trajectory_atom_count
            # Otherwise just print the whole error logs and stop here anyway
            print(output_logs)
            print(error_logs)
            raise ToolError('Something went wrong with GROMACS during the checking')
        # If we had an output then it means both topology and trajectory match in the number of atoms
        # Cleanup the file we just created and proceed
        output_sample_gro_file.remove()
        # Now make sure trajectory atoms are not more than topology atoms
        # Easiest way to print trajectory atoms is using gmx check
        # However if we feed this command with the whole trajectory it will read it all
        # To prevent this we must create a single frame before
        output_sample_xtc_file = File('.sample.xtc')
        # Note that we do NOT pass the -s argument here
        # Otherwise the structure/topology would eclipse the actual number of atoms in the trajectory
        run_gromacs(f'trjconv -f {trajectory_file.path} -o {output_sample_xtc_file.path} -dump 0',
            user_input='System', expected_output_filepath=output_sample_xtc_file.path)
        # Now read the number of atoms
        output_logs, error_logs = run_gromacs(f'check -f {output_sample_xtc_file.path}')
        search_results = re.search(GROMACS_ATOM_COUNT_CHECK, error_logs)
        if not search_results:
            print(error_logs)
            raise RuntimeError('Something went wrong when reading trajectory atoms')
        # Get the trajectory atom count
        trajectory_atom_count = int(search_results[1])
        # Cleanup the file we just created and proceed
        output_sample_xtc_file.remove()
        return topology_atom_count, trajectory_atom_count
    # For .top files we use PyTraj since MDtraj can not handle it
    if topology_file.format == 'top':
        return get_topology_and_trajectory_atoms_pytraj(topology_file, trajectory_file)
    # If the trajectory is a restart file MDtraj will not be able to read it
    # Make the conversion here, since restart files are single-frame trajectories this should be fast
    use_auxiliar_pdb = False
    auxiliar_pdb_file = get_auxiliar_filepath('.auxiliar.pdb')
    if trajectory_file.format == 'rst7':
        # Generate the auxiliar PDB file
        vmd_to_pdb(topology_file.path, trajectory_file.path, auxiliar_pdb_file)
        use_auxiliar_pdb = True
    # At this point the topology should be supported by MDtraj
    try:
        trajectory_path = auxiliar_pdb_file if use_auxiliar_pdb else trajectory_file.path
        if trajectory_file.format != 'pdb':
            # Note that declaring the iterator will not fail even when there is a mismatch
            trajectory = mdt.iterload(trajectory_path, top=topology_file.path, chunk=1)
            # We must consume the generator first value to make the error raise
            frame = next(trajectory)
            # Now obtain the number of atoms from the frame we just read
            trajectory_atom_count = frame.n_atoms
        else:
            # For PDB trajectories takes a lot of memory in MDtraj, 14GB for a 950MB (100000 frames) trajectory
            # Known issue marked as wont-fix (https://github.com/mdtraj/mdtraj/issues/1172), so we use PyTraj instead
            trajectory = pyt.iterload(trajectory_path)
            trajectory_atom_count = trajectory.n_atoms
        # And still, it may happen that the topology has more atoms than the trajectory but it loads
        # MDtraj may silently load as many coordinates as possible and discard the rest of atoms in topology
        # This behaviour has been observed with a gromacs .top topology and a PDB used as trajectory
        # Two double check the match, load the topology alone with PyTraj
        topology = pyt.load_topology(topology_file.path)
        topology_atom_count = topology.n_atoms
        return topology_atom_count, trajectory_atom_count
    except Exception as error:
        # If the error message matches with a known error then report the problem
        error_message = str(error)
        error_match = re.match(MDTRAJ_ATOM_MISMATCH_ERROR, error_message)
        if error_match:
            topology_atom_count = int(error_match[1])
            trajectory_atom_count = int(error_match[2])
            return topology_atom_count, trajectory_atom_count
        error_match = re.match(MDTRAJ_INSERTION_CODES_ERROR, error_message)
        if error_match:
            warn('The input topology has insertion codes.\n'
            ' Some tools may crash when reading the topology (MDtraj).\n'
            ' Some tools may ignore insertion codes when reading the topology (MDAnlysis, PyTraj, VMD).')
            # Use other tool to read the topology
            # Other tools could ignore the inserion codes
            # However this is not a problem here, where we only care about the number of atoms
            return get_topology_and_trajectory_atoms_pytraj(topology_file, trajectory_file)
        # If we do not know the error then raise it as is
        raise error


def get_structure_atoms(structure_file: 'File') -> int:
    """Get the atoms count from a structure alone."""
    # If this is not a Structure supported file then use an alternative function
    if structure_file.format == 'gro':
        return get_gmx_atom_count(structure_file)
    # Get the number of atoms in the input structure
    structure = Structure.from_file(structure_file.path)
    return structure.atom_count


def get_structure_and_trajectory_atoms(structure_file: 'File', trajectory_file: 'File') -> tuple[int, int]:
    """Get the atoms count from structure and trajectory together."""
    try:
        # Note that declaring the iterator will not fail even when there is a mismatch
        trajectory = mdt.iterload(trajectory_file.path, top=structure_file.path, chunk=1)
        # We must consume the generator first value to make the error raise
        frame = next(trajectory)
        # Now obtain the number of atoms from the frame we just read
        trajectory_atom_count = frame.n_atoms
        # And still, it may happen that the structure has more atoms than the trajectory but it loads
        # MDtraj may silently load as many coordinates as possible and discard the rest of atoms in topology
        # This behaviour has been observed with a gromacs .top topology and a PDB used as trajectory
        # Two double check the match, load the topology alone with PyTraj if it supports the format
        if structure_file.is_pytraj_supported():
            topology = pyt.load_topology(structure_file.path)
            structure_atom_count = topology.n_atoms
        # If pytraj does not support the format (e.g. .gro files) then try other tools
        elif structure_file.format == 'gro':
            structure_atom_count = get_gmx_atom_count(structure_file)
        # If no tool supports the format then we surrender
        else:
            raise RuntimeError(f'Unsupported structure format: {structure_file.format}')
        return structure_atom_count, trajectory_atom_count
    except Exception as error:
        # If the error message matches with a known error then report the problem
        error_message = str(error)
        error_match = re.match(MDTRAJ_ATOM_MISMATCH_ERROR, error_message)
        if error_match:
            structure_atom_count = int(error_match[1])
            trajectory_atom_count = int(error_match[2])
            return structure_atom_count, trajectory_atom_count
        # DANI: Esto de aquí abajo no he probado nunca que pase con PDBs
        error_match = re.match(MDTRAJ_INSERTION_CODES_ERROR, error_message)
        if error_match:
            warn('The input structure has insertion codes.\n'
            ' Some tools may crash when reading the structure (MDtraj).\n'
            ' Some tools may ignore insertion codes when reading the topology (MDAnlysis, PyTraj, VMD).')
            # Use other tool to read the topology
            # Other tools could ignore the inserion codes
            # However this is not a problem here, where we only care about the number of atoms
            return get_topology_and_trajectory_atoms_pytraj(structure_file, trajectory_file)
        # If we do not know the error then raise it as is
        raise error