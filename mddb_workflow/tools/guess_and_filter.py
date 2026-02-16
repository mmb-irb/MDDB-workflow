import pytraj as pt

from mddb_workflow.utils.auxiliar import all_cases, all_charges
from mddb_workflow.utils.constants import STANDARD_SOLVENT_RESIDUE_NAMES, STANDARD_COUNTER_ION_ATOM_NAMES
from mddb_workflow.utils.file import File
from mddb_workflow.utils.gmx_spells import make_index, tpr_filter, get_tpr_atom_count
from mddb_workflow.utils.type_hints import *

# Set the atom selections to be tried
ALL_COUNTER_IONS = all_charges(all_cases(STANDARD_COUNTER_ION_ATOM_NAMES))
AMBER_MASKS = {
    'No water': f'!:{",".join(STANDARD_SOLVENT_RESIDUE_NAMES)}',
    'No water or counter ions': f'!:{",".join(STANDARD_SOLVENT_RESIDUE_NAMES)},{",".join(ALL_COUNTER_IONS)}'
}

GROMACS_MASKS = {
    'No water': '!"Water"',
    'No water or counter ions': '!"Water"&!"Ion"'
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
        # Iterate possible atom selections
        for mask_name, mask in GROMACS_MASKS.items():
            # Create the index file with the current atom selection
            index_filepath = f'{input_topology_file.basepath}/.auxiliar.ndx'
            index_file = File(index_filepath)
            filter_group_name, group_exists = make_index(input_topology_file, index_file, mask)
            # If the filter was not created then it means there was no water or ions to begin with
            if not group_exists: continue
            # Filter the tpr
            tpr_filter(
                input_topology_file.path,
                output_topology_file.path,
                index_file.path,
                filter_group_name)
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
