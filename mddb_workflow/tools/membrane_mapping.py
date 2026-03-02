import numpy as np
from MDAnalysis.analysis.leaflet import LeafletFinder
from mddb_workflow.utils.auxiliar import save_json, load_json, warn
from mddb_workflow.utils.constants import LIPIDS_RESIDUE_NAMES
from mddb_workflow.utils.type_hints import *


def generate_membrane_mapping(
    inchikeys: dict[str, 'InChIKeyData'],
    lipid_references: dict[str, dict],
    universe: 'Universe',
    cg_residues: list[int],
    output_file: 'File',
) -> dict:
    """Generate a list of residue numbers of membrane components from a given structure and topology file.

    Args:
        inchikeys: Dictionary mapping InChI keys to InChIKeyData objects.
        lipid_references: Dictionary mapping InChI keys to lipid database references.
        structure_file: File object representing the structure file (e.g., PDB, GRO).
        universe: MDAnalysis Universe object containing the structure and topology.
        cg_residues: List of residue indices for coarse-grained residues.
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
    if len(cg_residues) > 0:
        headgroup_ag = coarse_grain_headgroups(universe)
        membrane_map = find_leaflets(universe, headgroup_ag)
    else:
        lipid_map = {key: inchikeys[key] for key in lipid_references.keys()}
        headgroup_ag = all_atom_headgroups(lipid_map, universe)
        membrane_map = find_leaflets(universe, headgroup_ag)

    if membrane_map['n_mems'] > 0:
        save_json(membrane_map, output_file.path)
    return membrane_map


def all_atom_headgroups(
    lipid_map: dict[str, 'InChIKeyData'],
    universe: 'Universe',
) -> dict:
    """Extract the headgroup atoms for all lipids in the system using the most charged atom.
    In case of missing charges, it will try to guess the headgroup based on the residue name.
    """
    if len(lipid_map) == 0:
        # no lipids found in the structure.
        print("No lipid residues found in the structure.")
        return None

    # Select only the lipids and potential membrane members
    lipid_resindices = []
    for ref in lipid_map.values():
        lipid_resindices.extend(ref.resindices)
    # if no lipids are found, we save the empty mapping and return
    if len(lipid_resindices) == 0:
        # no lipids found in the structure.
        return None
    all_lipids = universe.select_atoms(f'(resindex {" ".join(map(str, (lipid_resindices)))})')

    # For better leaflet assignation we only use polar atoms
    if not hasattr(universe.atoms, 'charges'):
        warn("Atom charges not found, guessing headgroups by name.")
        polar_atoms = []
        for residue in all_lipids.residues:
            if residue.resname == 'CHL':
                polar_atoms.append(residue.atoms.select_atoms('name O*').indices[0])
            else:
                # '(resname CHL and element O*) or (resname PA and name P*)'
                PO_candidates = residue.atoms.select_atoms('name P*')
                # Check if only bonded to oxygen atoms, if so, it is likely a phosphate group and we select the P atom as polar atom
                if len(PO_candidates) == 1:
                    bonded_atoms = np.unique(PO_candidates[0].bonds.to_indices().flatten())
                    bond_elements = universe.atoms[bonded_atoms].elements
                    if set(bond_elements) == {'O', 'P'}:
                        polar_atoms.append(PO_candidates[0].index)
                        continue
                warn(f"Could not find polar atom for residue {residue.resname} {residue.resid}. Skipping analysis")
                return None
    else:
        print("Using atom charges to find polar atoms for leaflet assignment.")
        charges = abs(np.array([atom.charge for atom in universe.atoms]))
        polar_atoms = []
        for ridx in all_lipids.residues.resindices:  # take only one per fragment
            res = universe.residues[ridx]
            res_ch = charges[res.atoms.ix]
            max_ch_idx = np.argmax(res_ch)
            polar_atoms.append(res.atoms[max_ch_idx].index)
    polar_atoms = np.array(polar_atoms)
    headgroup_ag = universe.atoms[polar_atoms]
    return headgroup_ag


def coarse_grain_headgroups(universe: 'Universe') -> dict:
    """Generate membrane mapping for coarse-grained systems using MDAnalysis LeafletFinder with P atoms as headgroups."""
    # Find all lipid residues in the system
    lipid_atoms = headgroup_ag = universe.atoms[[]]
    for resname in LIPIDS_RESIDUE_NAMES:
        lipid = universe.select_atoms(f'resname {resname}')
        if len(lipid) > 0:
            lipid_atoms += lipid.atoms

    if len(lipid_atoms) == 0:
        print("No coarse-grained lipid atoms found in the structure.")
        return None
    print(f"Found {len(lipid_atoms.residues)} lipid residues by name.")
    # Use P atoms as headgroups for LeafletFinder
    headgroup_ag = lipid_atoms.select_atoms('name P*')
    return headgroup_ag


def find_leaflets(universe: 'Universe', headgroup_ag: list[int] | None) -> dict:
    """Find leaflets using MDAnalysis LeafletFinder and generate membrane mapping."""
    membrane_map = {'n_mems': 0, 'mems': {}, 'no_mem_lipid': [], 'version': '0.2.0'}
    if not headgroup_ag:
        print("No headgroup atoms found, cannot assign leaflets.")
        return membrane_map
    print(' Running MDAnalysis LeafletFinder')
    L = LeafletFinder(universe, headgroup_ag)
    leaflets = []
    no_mem_lipids = []
    for group in L.groups():
        # 30 lipids like in FATSLiM
        if len(group) > 30:
            leaflets.append(group)
        else:
            no_mem_lipids.extend(group.residues.resindices)
    n_mems = len(leaflets) // 2
    membrane_map['n_mems'] = n_mems

    for i in range(n_mems):
        # Check in which leaflets each of the polar atoms is and save them
        memb_data = {
            'leaflets': {
                'bot': leaflets[i * 2].residues.atoms.indices.tolist(),
                'top': leaflets[i * 2 + 1].residues.atoms.indices.tolist()
            },
            'polar_atoms': {
                'bot': leaflets[i * 2].atoms.indices.tolist(),
                'top': leaflets[i * 2 + 1].atoms.indices.tolist()
            }
        }
        membrane_map['mems'][str(i)] = memb_data
        # Print leaflets stats
        print(f"Membrane {i}:\n"
              f"    - Top leaflet: {len(memb_data['polar_atoms']['top'])} lipids\n"
              f"    - Bottom leaflet: {len(memb_data['polar_atoms']['bot'])} lipids")

    if len(no_mem_lipids) > 0:
        print(f"Unassigned lipids: {len(no_mem_lipids)}.")
        no_mem_lipids = universe.residues[no_mem_lipids].atoms.indices.tolist()
    membrane_map['no_mem_lipid'] = no_mem_lipids

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
