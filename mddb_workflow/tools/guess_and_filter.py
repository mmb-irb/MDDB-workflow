import pytraj as pt

from mddb_workflow.utils.auxiliar import all_cases, all_charges
from mddb_workflow.utils.constants import STANDARD_SOLVENT_RESIDUE_NAMES, STANDARD_COUNTER_ION_ATOM_NAMES
from mddb_workflow.utils.file import File
from mddb_workflow.utils.gmx_spells import tpr_filter, get_tpr_atom_count
from mddb_workflow.utils.structures import Structure
from mddb_workflow.utils.selections import Selection
from mddb_workflow.utils.type_hints import *

# Set the atom selections to be tried
ALL_COUNTER_IONS = all_charges(all_cases(STANDARD_COUNTER_ION_ATOM_NAMES))
AMBER_MASKS = {
    'No water': f'!:{",".join(STANDARD_SOLVENT_RESIDUE_NAMES)}',
    'No water or counter ions': f'!:{",".join(STANDARD_SOLVENT_RESIDUE_NAMES)},{",".join(ALL_COUNTER_IONS)}'
}

DRY_SELECTION_NAME = 'dry' # Not water
DRIER_SELECTION_NAME = 'drier' # Not water or counter ions
GROMACS_MASKS = {
    'No water': DRY_SELECTION_NAME,
    'No water or counter ions': DRIER_SELECTION_NAME,
}

# Given a topology file and a specific atom number, find an atom selection which matches the number
# Try different typical atom selections and if none matches then surrender
# Note that we have not a structure yet at this point to play with so we must use topology-specific tools
def guess_and_filter_topology (
    input_topology_file : 'File',
    output_topology_file : 'File',
    target_atom_count : int,
    verbose : bool = True) -> bool:
    if verbose:
        print(f'Trying to guess topology selection to match {target_atom_count} atoms')
    # Amber topologies and others
    if input_topology_file.format == 'prmtop':
        # Load the topology
        pt_topology = pt.load_topology(filename=input_topology_file.path)
        # Iterate possible atom selections
        for mask_name, mask in AMBER_MASKS.items():
            filtered_pt_topology = pt_topology[mask]
            filtered_atoms_count = filtered_pt_topology.n_atoms
            if verbose: print(f' {mask_name} -> {filtered_atoms_count} atoms')
            # If the filtered atom count does not match the target atom count then continue
            if filtered_atoms_count != target_atom_count: continue
            if verbose: print('  Got it!')
            # If we have a match in the atom count then save the filtered topology
            # WARNING: If the output topology is a symlink it will try to overwrite the origin
            # Remove it to avoid overwriting input data
            if output_topology_file.is_symlink(): output_topology_file.remove()
            # Now write the filtered topology
            pt.write_parm(
                filename=output_topology_file.path,
                top=filtered_pt_topology,
                format=input_topology_file.get_pytraj_parm_format(),
                overwrite=True
            )
            return True
        return False
    # Gromacs topologies
    elif input_topology_file.format == 'tpr':
        # We will not rely in Gromacs automatic atom selections
        # Instead we will make our own index file with our custom selections
        # The index file is written next to the output so input directories are not touched
        index_filepath = f'{output_topology_file.basepath}/.auxiliar.ndx'
        index_file = File(index_filepath)
        # Read and parse the TPR file
        early_structure = Structure.from_tpr_file(input_topology_file.path)
        # Set the dry selections
        # Water is found by residue name only, since atom elements are not known for coarse grained atoms yet
        water_selection = early_structure.select_water_by_residue_name()
        dry_selection = early_structure.invert_selection(water_selection)
        water_and_ion_selection = water_selection + early_structure.select_counter_ions()
        drier_selection = early_structure.invert_selection(water_and_ion_selection)
        # Write the index file
        with open(index_file.path, 'w') as file:
            file.write(dry_selection.to_ndx(DRY_SELECTION_NAME))
            file.write(drier_selection.to_ndx(DRIER_SELECTION_NAME))
        # Iterate possible atom selections
        for mask_name, mask in GROMACS_MASKS.items():
            # Filter the tpr
            tpr_filter(
                input_topology_file.path,
                output_topology_file.path,
                index_file.path, mask)
            # Count atoms in the new filtered tpr
            filtered_atoms_count = get_tpr_atom_count(output_topology_file.path)
            if verbose: print(f' {mask_name} -> {filtered_atoms_count}')
            # If the filtered atom count does not match the target atom count then continue
            if filtered_atoms_count != target_atom_count: continue
            if verbose: print('  Got it!')
            # If we have a match in the atom count then return True
            return True
        # Small cleanup
        if output_topology_file.exists: output_topology_file.remove()
        return False
    # If the format is not recognized then we surrender
    else:
        raise ValueError(f'The "guess and filter" tool does not support {input_topology_file.format} format yet')


def _select_water_and_counter_ions(pt_topology) -> tuple[set[int], list[int]]:
    """Find water atoms and counter ion atoms in a pytraj topology.
    Counter ions are single atom residues with a standard counter ion atom name.
    Note that elements are not used, so this also works for coarse grained structures.
    """
    water_atoms = set()
    counter_ion_atoms = []
    for residue in pt_topology.residues:
        if residue.name in STANDARD_SOLVENT_RESIDUE_NAMES:
            water_atoms.update(range(residue.first, residue.last))
        elif residue.n_atoms == 1:
            atom_name = pt_topology.atom(residue.first).name
            # Remove non-letter characters (e.g. '+' and '-' symbols)
            simple_atom_name = ''.join(filter(str.isalpha, atom_name.upper()))
            if simple_atom_name in STANDARD_COUNTER_ION_ATOM_NAMES:
                counter_ion_atoms.append(residue.first)
    return water_atoms, counter_ion_atoms


# Given a structure file with more atoms than the trajectory, find an atom selection which matches the number
# Try different typical atom selections and if none matches then surrender
def guess_and_filter_structure (
    input_structure_file : 'File',
    output_structure_file : 'File',
    target_atom_count : int,
    verbose : bool = True) -> bool:
    if not input_structure_file.is_pytraj_supported():
        raise ValueError(f'The "guess and filter" tool does not support {input_structure_file.format} structures yet')
    if verbose:
        print(f'Trying to guess structure selection to match {target_atom_count} atoms')
    structure = pt.load(input_structure_file.path)
    pt_topology = structure.top
    water_atoms, counter_ion_atoms = _select_water_and_counter_ions(pt_topology)
    # Ions added after solvating (e.g. with Gromacs genion) replace water molecules so they come after the water
    # Counter ions placed before the water (e.g. those of a lipid nanoparticle) are usually kept in dry trajectories
    first_water_atom = min(water_atoms) if water_atoms else pt_topology.n_atoms
    added_ion_atoms = { atom for atom in counter_ion_atoms if atom > first_water_atom }
    excluded_selections = {
        'No water': water_atoms,
        'No water or counter ions': water_atoms.union(counter_ion_atoms),
        'No water or ions added after the water': water_atoms.union(added_ion_atoms),
    }
    for selection_name, excluded_atoms in excluded_selections.items():
        filtered_atoms_count = pt_topology.n_atoms - len(excluded_atoms)
        if verbose: print(f' {selection_name} -> {filtered_atoms_count} atoms')
        if filtered_atoms_count != target_atom_count: continue
        if verbose: print('  Got it!')
        kept_atoms = [ atom for atom in range(pt_topology.n_atoms) if atom not in excluded_atoms ]
        # WARNING: If the output structure is a symlink it will try to overwrite the origin
        # Remove it to avoid overwriting input data
        if output_structure_file.is_symlink(): output_structure_file.remove()
        # Pytraj can not write mmCIF files
        if output_structure_file.format == 'cif':
            filtered_structure = Structure.from_file(input_structure_file.path).filter(Selection(kept_atoms))
            filtered_structure.generate_mmcif_file(output_structure_file.path)
        else:
            pt.write_traj(output_structure_file.path, structure[pt.AtomMask(kept_atoms)], overwrite=True)
        return True
    return False
