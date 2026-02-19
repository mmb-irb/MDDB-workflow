import os
import numpy as np
from biobb_mem.fatslim.fatslim_membranes import fatslim_membranes, parse_index
from mddb_workflow.utils.auxiliar import save_json, load_json
from mddb_workflow.utils.constants import LIPIDS_RESIDUE_NAMES
from mddb_workflow.utils.type_hints import *
from contextlib import redirect_stdout


def generate_membrane_mapping(
    inchikeys: dict[str, 'InChIKeyData'],
    lipid_references: dict[str, dict],
    structure_file: 'File',
    universe: 'Universe',
    output_file: 'File',
) -> dict:
    """Generate a list of residue numbers of membrane components from a given structure and topology file.

    Args:
        inchikeys: Dictionary mapping InChI keys to InChIKeyData objects.
        lipid_references: Dictionary mapping InChI keys to lipid database references.
        structure_file: File object representing the structure file (e.g., PDB, GRO).
        universe: MDAnalysis Universe object containing the structure and topology.
        output_file: Output JSON file to save membrane mapping.

    Returns:
        dict: A dictionary containing the membrane mapping.
            - n_mems (int): Number of membranes detected.
            - mems (dict): Dictionary where each key is a membrane index (str) and value contains:
              - leaflets (dict): Contains 'top' and 'bot' keys, each with a list of atom indices.
              - polar_atoms (dict): Contains 'top' and 'bot' keys, each with a list of polar atom indices.
            - no_mem_lipid (list): List of atom indices for lipids not assigned to any membrane.

    Notes:
        - The function identifies lipid and non-lipid residues based on InChI keys.
        - It classifies residues and checks for potential misclassifications.
        - Lipid residues are selected and neighboring lipids are found.
        - Clusters of lipids are identified, and clusters with more than 30 lipids are considered as membranes.
        - If debug is enabled, the function returns additional information including lipid residues, neighbors, counts, and clusters.

    """
    if not universe: raise RuntimeError('Missing universe')
    # Select only the lipids from the InChIKeyData
    lipid_map = {key: inchikeys[key] for key in lipid_references.keys()}
    if 'Cg' in universe.atoms.types:
        membrane_map = coarse_grain_membranes(structure_file, universe)
    else:
        membrane_map = all_atom_membranes(lipid_map, structure_file, universe)
    if membrane_map['n_mems'] > 0:
        save_json(membrane_map, output_file.path)
    # Print leaflets stats
    for mem_idx, mem_data in membrane_map['mems'].items():
        n_top = len(mem_data['leaflets']['top'])
        n_bot = len(mem_data['leaflets']['bot'])
        print(f"Membrane {mem_idx}: {n_top} lipids in top leaflet, {n_bot} lipids in bottom leaflet.")
    print(f"{len(membrane_map['no_mem_lipid'])} lipid atoms not assigned to any membrane.")
    return membrane_map


def all_atom_membranes(
    lipid_map: dict[str, 'InChIKeyData'],
    structure_file: 'File',
    universe: 'Universe',
) -> dict:
    """Generate membrane mapping for all-atom systems using FATSLiM."""
    membrane_map = {'n_mems': 0, 'mems': {}, 'no_mem_lipid': []}

    if len(lipid_map) == 0:
        # no lipids found in the structure.
        print("No lipid residues found in the structure.")
        return membrane_map

    # Select only the lipids and potential membrane members
    lipid_ridx = []
    glclipid_ridx = []
    for ref in lipid_map.values():
        lipid_ridx.extend(ref.resindices)
        if ref.moltype == 'fragment':
            glclipid_ridx.extend(ref.fragments)
    # if no lipids are found, we save the empty mapping and return
    if len(lipid_ridx) == 0:
        # no lipids found in the structure.
        return membrane_map
    mem_candidates = universe.select_atoms(f'(resindex {" ".join(map(str, (lipid_ridx)))})')

    # for better leaflet assignation we only use polar atoms
    charges = abs(np.array([atom.charge for atom in universe.atoms]))
    polar_atoms = []
    for ridx in mem_candidates.residues.resindices:
        res = universe.residues[ridx]
        res_ch = charges[res.atoms.ix]
        max_ch_idx = np.argmax(res_ch)
        polar_atoms.append(res.atoms[max_ch_idx].index)
    polar_atoms = np.array(polar_atoms)
    headgroup_sel = f'(index {" ".join(map(str, (polar_atoms)))})'
    # Run FATSLiM to find the membranes
    prop = {
        'selection': headgroup_sel,
        'cutoff': 2.2,
        'ignore_no_box': True,
        'disable_logs': True,
        'return_hydrogen': True,
        # Do not copy the input file to the sandbox
        'disable_sandbox': True,
    }
    output_ndx_path = "tmp_mem_map.ndx"
    print(' Running BioBB FATSLiM Membranes')
    with redirect_stdout(None):
        fatslim_membranes(input_top_path=structure_file.absolute_path,
                          output_ndx_path=output_ndx_path,
                          properties=prop)
    # Parse the output to get the membrane mapping
    mem_map = parse_index(output_ndx_path)
    os.remove(output_ndx_path)
    # Save the membrane mapping as a JSON file
    n_mems = len(mem_map) // 2
    membrane_map['n_mems'] = n_mems
    no_mem_lipids = set(mem_candidates.atoms.indices)
    for i in range(n_mems):
        # and they are not assigned to any membrane. FATSLiM indexes start at 1
        bot = (np.array(mem_map[f'membrane_{i + 1}_leaflet_1']) - 1).tolist()  # JSON does not support numpy arrays
        top = (np.array(mem_map[f'membrane_{i + 1}_leaflet_2']) - 1).tolist()
        top_ridx = universe.atoms[top].residues.resindices
        bot_ridx = universe.atoms[bot].residues.resindices
        # Some polar atoms from the glucids are to far from the polar atoms of the lipids
        top = set(top)
        bot = set(bot)
        remove_ridx = []
        for grp in glclipid_ridx:
            if np.in1d(grp, top_ridx).any():
                top.update(universe.residues[grp].atoms.indices)
                remove_ridx.append(grp)
                continue
            if np.in1d(grp, bot_ridx).any():
                bot.update(universe.residues[grp].atoms.indices)
                remove_ridx.append(grp)
        for grp in remove_ridx:
            glclipid_ridx.remove(grp)
        # Remove lipids not assigned to any membrane
        no_mem_lipids -= top
        no_mem_lipids -= bot
        # Check in which leaflets each of the polar atoms is and save them
        membrane_map['mems'][str(i)] = {
            'leaflets': {
                'bot': list(map(int, bot)),
                'top': list(map(int, top))
            },
            'polar_atoms': {
                'bot': polar_atoms[np.in1d(polar_atoms, list(bot))].tolist(),
                'top': polar_atoms[np.in1d(polar_atoms, list(top))].tolist()
            }
        }

    membrane_map['no_mem_lipid'] = list(map(int, no_mem_lipids))
    # Print the results and save the membrane mapping
    no_mem_lipids_str = f'{len(glclipid_ridx)} lipid/s not assigned to any membrane.' if len(glclipid_ridx) > 0 else ''
    print(f'{n_mems} membrane/s found. ' + no_mem_lipids_str)
    return membrane_map


def coarse_grain_membranes(structure_file: 'File', universe: 'Universe') -> dict:
    """Generate membrane mapping for coarse-grained systems using residue names and P atoms as headgroups."""
    membrane_map = {'n_mems': 0, 'mems': {}, 'no_mem_lipid': []}

    # Find all lipid residues in the system
    lipid_residues = []
    for resname in LIPIDS_RESIDUE_NAMES:
        lipids = universe.select_atoms(f'resname {resname}')
        if len(lipids) > 0:
            lipid_residues.append(resname)

    if len(lipid_residues) == 0:
        print("No coarse-grained lipid residues found in the structure.")
        return membrane_map

    # Select P atoms (headgroups) from lipid residues
    headgroup_atoms = f'resname {" ".join(lipid_residues)} and name P*'

    # Run FATSLiM to find the membranes
    prop = {
        'selection': headgroup_atoms,
        'cutoff': 2.5,  # Slightly larger cutoff for CG systems
        'ignore_no_box': True,
        'disable_logs': True,
        'disable_sandbox': True,
    }

    output_ndx_path = "tmp_cg_mem_map.ndx"
    print(' Running BioBB FATSLiM Membranes for coarse-grained system')

    # with redirect_stdout(None):
    fatslim_membranes(input_top_path=structure_file.absolute_path,
                      output_ndx_path=output_ndx_path,
                      properties=prop)

    # Parse the output to get the membrane mapping
    mem_map = parse_index(output_ndx_path)
    os.remove(output_ndx_path)

    # Process the membrane mapping
    n_mems = len(mem_map) // 2
    membrane_map['n_mems'] = n_mems

    all_lipid_atoms = universe.select_atoms(f'resname {" ".join(lipid_residues)}')
    no_mem_lipids = set(all_lipid_atoms.atoms.indices)

    for i in range(n_mems):
        # FATSLiM indexes start at 1, convert to 0-based
        bot_atoms = (np.array(mem_map[f'membrane_{i + 1}_leaflet_1']) - 1).tolist()
        top_atoms = (np.array(mem_map[f'membrane_{i + 1}_leaflet_2']) - 1).tolist()

        # Remove assigned atoms from no_mem_lipids
        no_mem_lipids -= set(bot_atoms)
        no_mem_lipids -= set(top_atoms)

        membrane_map['mems'][str(i)] = {
            'leaflets': {
                'bot': bot_atoms,
                'top': top_atoms
            },
            'polar_atoms': {
                'bot': universe.atoms[bot_atoms].select_atoms('name P*').indices.tolist(),
                'top': universe.atoms[top_atoms].select_atoms('name P*').indices.tolist()
            }
        }

    membrane_map['no_mem_lipid'] = list(map(int, no_mem_lipids))

    # Print results
    no_mem_lipids_str = f'{len(no_mem_lipids)} lipid atoms not assigned to any membrane.' if len(no_mem_lipids) > 0 else ''
    print(f'{n_mems} membrane/s found in coarse-grained system. ' + no_mem_lipids_str)

    return membrane_map


def display_membrane_mapping(mem_map: str, pdb: str):
    """Display the membrane mapping using nglview."""
    try:
        import nglview as nv  # type: ignore
    except ImportError:
        raise ImportError("nglview is required for displaying membrane mapping. Please install it using 'pip install nglview'.")

    mem_map = load_json(mem_map)
    top = f"@{','.join(map(str, mem_map['mems']['0']['leaflets']['top']))}"
    bot = f"@{','.join(map(str, mem_map['mems']['0']['leaflets']['bot']))}"
    no_mem = f"@{','.join(map(str, mem_map['no_mem_lipid']))}"
    polar_top = f"@{','.join(map(str, mem_map['mems']['0']['polar_atoms']['top']))}"
    polar_bot = f"@{','.join(map(str, mem_map['mems']['0']['polar_atoms']['bot']))}"

    view = nv.show_file(pdb)
    view.clear(0)
    # view.clear(1)
    view.add_point(selection='not protein', color='green', scale=1.5)
    view.add_point(selection=f'{top}', color='blue')
    view.add_point(selection=f'{bot}', color='yellow')
    view.add_point(selection=no_mem, color='red', scale=0.5)
    view.add_point(selection='protein', color='black', scale=0.5)
    view.add_ball_and_stick(selection=polar_top, color='cyan', aspectRatio=4)
    view.add_ball_and_stick(selection=polar_bot, color='orange', aspectRatio=4)
    return view
