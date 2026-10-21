import numpy as np
from collections import Counter
from MDAnalysis.analysis.leaflet import LeafletFinder
from mddb_workflow.utils.auxiliar import save_json, load_json, warn
from mddb_workflow.utils.constants import LIPIDS_RESIDUE_NAMES
from mddb_workflow.utils.type_hints import *

MEMBRANE_MAPPING_VERSION = '0.3.0'
# Amphiphilic aggregate types
BILAYER = 'bilayer'
VESICLE = 'vesicle'
MICELLE = 'micelle'
# Minimum number of lipids in a leaflet, 24 lipids like in FATSLiM
MIN_LEAFLET_LIPIDS = 24
# A leaflet is flat when sqrt(smallest / largest eigenvalue) of its polar atom positions is below this
FLATNESS_THRESHOLD = 0.3
# Absolute cosine between the head to tail vector and the reference direction for a lipid to be clearly oriented
STRONG_ORIENTATION = 0.5
# Fraction of lipids clearly oriented against the majority from which a group is considered mixed (merged leaflets)
MIXED_ORIENTATION = 0.25
# Minimum absolute cosine between normals of the two leaflets of a bilayer
PARALLEL_NORMALS = 0.9
# Maximum distance (Å) between the centroids of the two leaflets of a bilayer along the normal
MAX_BILAYER_SEPARATION = 80
# Maximum distance between the centers of vesicle leaflets, as a fraction of the inner leaflet radius
VESICLE_CENTER_TOLERANCE = 0.5


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
            - n_mems (int): Number of amphiphilic aggregates (membranes) detected.
            - mems (dict): Dictionary where each key is a membrane index (str) and value contains:
              - type (str): Aggregate type: 'bilayer' (flat), 'vesicle' or 'micelle'.
              - leaflets (dict): Atom indices of each leaflet. Leaflet names are 'top' and 'bot' for bilayers,
                'outer' and 'inner' for vesicles and 'outer' for micelles.
              - polar_atoms (dict): Polar atom indices (one per lipid) of each leaflet, with the same names.
            - no_mem_lipid (list): List of atom indices for lipids not assigned to any membrane.
            - version (str): Version of the membrane mapping format.

    Notes:
        - The function identifies lipid and non-lipid residues based on InChI keys.
        - It classifies residues and checks for potential misclassifications.
        - Lipid residues are selected and neighboring lipids are found.
        - Leaflets are clusters of at least 24 lipid polar atoms found by LeafletFinder.
        - Leaflets are classified into aggregates from their shape and lipid orientation (see classify_leaflets).

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
    In case of missing charges, it will try to guess the headgroup from the atom elements and bonds.
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
        warn("Atom charges not found, guessing headgroups by element.")
        polar_atoms = []
        skipped = Counter()
        for residue in all_lipids.residues:
            polar_atom = guess_polar_atom(residue)
            if polar_atom is None:
                skipped[residue.resname] += 1
                continue
            polar_atoms.append(polar_atom)
        if len(skipped) > 0:
            warn(f"Could not find a unique polar atom for some lipids, they will not be used for leaflet assignment: {dict(skipped)}")
    else:
        print("Using atom charges to find polar atoms for leaflet assignment.")
        charges = abs(np.array([atom.charge for atom in universe.atoms]))
        polar_atoms = []
        for ridx in all_lipids.residues.resindices:  # take only one per fragment
            res = universe.residues[ridx]
            res_ch = charges[res.atoms.ix]
            max_ch_idx = np.argmax(res_ch)
            polar_atoms.append(res.atoms[max_ch_idx].index)
    polar_atoms = np.array(polar_atoms, dtype=int)
    headgroup_ag = universe.atoms[polar_atoms]
    return headgroup_ag


def guess_polar_atom(residue: 'Residue') -> int | None:
    """Guess the headgroup atom of a lipid residue from elements and bonds when charges are missing.

    Candidates are tried by priority: phosphate P (bonded only to O), N (e.g. amine or ammonium),
    hydroxyl O (bonded to H, e.g. sterols) and finally any O.
    The first kind with candidates must have exactly one, otherwise the headgroup is ambiguous and None is returned.
    """
    atoms = residue.atoms
    elements = atoms.elements if hasattr(atoms, 'elements') else np.array([name[0] for name in atoms.names])
    has_bonds = hasattr(atoms, 'bonds') and hasattr(atoms, 'elements')

    def bonded_elements(atom) -> set[str]:
        return {other.element for bond in atom.bonds for other in bond.atoms if other.index != atom.index}

    phosphorus = atoms[elements == 'P']
    oxygens = atoms[elements == 'O']
    candidate_groups = [
        # Check if only bonded to oxygen atoms, if so, it is likely a phosphate group
        [atom for atom in phosphorus if not has_bonds or bonded_elements(atom) == {'O'}],
        list(atoms[elements == 'N']),
        [atom for atom in oxygens if has_bonds and 'H' in bonded_elements(atom)],
        list(oxygens),
    ]
    for candidates in candidate_groups:
        if len(candidates) == 0: continue
        if len(candidates) == 1: return candidates[0].index
        return None
    return None


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


def split_leaflet_groups(universe: 'Universe', headgroup_ag: 'AtomGroup', cutoff: float) -> tuple[list, list]:
    """Run LeafletFinder and split its groups into leaflets and residue indices of unassigned lipids."""
    L = LeafletFinder(universe, headgroup_ag, cutoff=cutoff)
    leaflets = []
    no_mem_lipids = []
    for group in L.groups():
        # 24 lipids like in FATSLiM
        if len(group) >= MIN_LEAFLET_LIPIDS:
            leaflets.append(group)
        else:
            no_mem_lipids.extend(group.residues.resindices)
    return leaflets, no_mem_lipids


def pbc_center(positions: np.ndarray, box: np.ndarray | None) -> np.ndarray:
    """Get the center of a set of positions in a periodic box using the circular mean of each dimension."""
    center = positions.mean(axis=0)
    if box is None: return center
    for dim in range(3):
        if box[dim] <= 0: continue
        theta = positions[:, dim] / box[dim] * 2 * np.pi
        center[dim] = (np.arctan2(-np.sin(theta).mean(), -np.cos(theta).mean()) + np.pi) / (2 * np.pi) * box[dim]
    return center


def minimum_image(vectors: np.ndarray, box: np.ndarray | None) -> np.ndarray:
    """Apply the minimum image convention to vectors in an orthorhombic box."""
    if box is None: return vectors
    return vectors - box * np.round(vectors / box)


def describe_leaflet(leaflet: 'AtomGroup', box: np.ndarray | None) -> dict:
    """Get the shape and lipid orientation of a LeafletFinder group of polar atoms.

    The group is flat when its thinnest principal axis is much shorter than the longest one.
    Lipid orientation is the head to tail vector (polar atom to residue center) projected on
    the leaflet normal when it is flat or on the radial direction from the leaflet center otherwise.
    Lipids in a real leaflet are oriented in the same direction, while a group with many lipids
    oriented in both directions is likely two leaflets merged by LeafletFinder.
    """
    positions = leaflet.positions
    center = pbc_center(positions, box)
    relative = minimum_image(positions - center, box)
    shift = relative.mean(axis=0)
    relative -= shift
    eigenvalues, eigenvectors = np.linalg.eigh(np.cov(relative.T))
    flat = np.sqrt(max(eigenvalues[0], 0) / eigenvalues[2]) < FLATNESS_THRESHOLD
    normal = eigenvectors[:, 0]
    # Make the normal point to the positive side of its main axis so top and bottom are consistent
    if normal[np.argmax(np.abs(normal))] < 0: normal = -normal
    radii = np.linalg.norm(relative, axis=1)
    residue_centers = np.array([atom.residue.atoms.positions.mean(axis=0) for atom in leaflet])
    tails = minimum_image(residue_centers - positions, box)
    tails /= np.maximum(np.linalg.norm(tails, axis=1), 1e-6)[:, None]
    reference = normal if flat else relative / np.maximum(radii, 1e-6)[:, None]
    cosines = (tails * reference).sum(axis=1)
    positive = np.mean(cosines > STRONG_ORIENTATION)
    negative = np.mean(cosines < -STRONG_ORIENTATION)
    return {
        'group': leaflet,
        'flat': flat,
        'centroid': center + shift,
        'normal': normal,
        'radius': radii.mean(),
        'coherent': min(positive, negative) < MIXED_ORIENTATION,
        # Flat: tails along (+1) or against (-1) the normal. Curved: tails outwards (+1) or inwards (-1)
        'direction': 1 if positive > negative else -1,
    }


def classify_leaflets(leaflets: list, box: np.ndarray | None) -> tuple[list[tuple[str, dict]], list]:
    """Group leaflets into amphiphilic aggregates: flat bilayers, vesicles and micelles.

    - Vesicle: a curved leaflet with tails inwards (outer) and a smaller concentric curved leaflet with tails outwards (inner).
    - Bilayer: two flat leaflets with parallel normals whose tails face each other.
    - Micelle: a single curved leaflet with tails inwards.

    Returns the aggregates as (type, leaflets by name) and the leaflets which could not be classified.
    """
    descriptions = [describe_leaflet(leaflet, box) for leaflet in leaflets]
    coherent = [desc for desc in descriptions if desc['coherent']]
    unclassified = [desc['group'] for desc in descriptions if not desc['coherent']]
    used = set()
    aggregates = []

    curved_in = [desc for desc in coherent if not desc['flat'] and desc['direction'] < 0]
    curved_out = [desc for desc in coherent if not desc['flat'] and desc['direction'] > 0]
    for outer in curved_in:
        candidates = [inner for inner in curved_out if id(inner) not in used and inner['radius'] < outer['radius']
            and np.linalg.norm(minimum_image(inner['centroid'] - outer['centroid'], box)) < VESICLE_CENTER_TOLERANCE * inner['radius']]
        if len(candidates) == 0: continue
        inner = min(candidates, key=lambda desc: np.linalg.norm(minimum_image(desc['centroid'] - outer['centroid'], box)))
        used.update([id(outer), id(inner)])
        aggregates.append((VESICLE, {'outer': outer['group'], 'inner': inner['group']}))

    flats = [desc for desc in coherent if desc['flat']]
    for first in flats:
        if id(first) in used: continue
        best, best_distance = None, None
        for second in flats:
            if second is first or id(second) in used: continue
            if abs(np.dot(first['normal'], second['normal'])) < PARALLEL_NORMALS: continue
            separation = minimum_image(second['centroid'] - first['centroid'], box)
            distance = np.dot(separation, first['normal'])
            # Tails of each leaflet must point towards the other leaflet
            first_faces = np.sign(distance) == first['direction']
            second_faces = np.sign(np.dot(-separation, second['normal'])) == second['direction']
            if not (first_faces and second_faces) or abs(distance) > MAX_BILAYER_SEPARATION: continue
            if best is None or abs(distance) < best_distance:
                best, best_distance = second, abs(distance)
        if best is None: continue
        used.update([id(first), id(best)])
        # The top leaflet is the one further along the normal
        separation = minimum_image(best['centroid'] - first['centroid'], box)
        top, bot = (best, first) if np.dot(separation, first['normal']) > 0 else (first, best)
        aggregates.append((BILAYER, {'top': top['group'], 'bot': bot['group']}))

    for desc in curved_in:
        if id(desc) in used: continue
        used.add(id(desc))
        aggregates.append((MICELLE, {'outer': desc['group']}))

    unclassified += [desc['group'] for desc in coherent if id(desc) not in used]
    return aggregates, unclassified


def find_leaflets(universe: 'Universe', headgroup_ag: list[int] | None, cutoffs: tuple[float, ...] = (15.0, 12.0, 10.0)) -> dict:
    """Find leaflets using MDAnalysis LeafletFinder and classify them into amphiphilic aggregates.

    When some leaflets can not be classified (e.g. two leaflets too close at a fusion stalk are
    merged by LeafletFinder) it is retried with the next smaller cutoff. The first cutoff where all
    leaflets are classified is used or, if none, the one classifying the most lipids.
    """
    membrane_map = {'n_mems': 0, 'mems': {}, 'no_mem_lipid': [], 'version': MEMBRANE_MAPPING_VERSION}
    if not headgroup_ag:
        print("No headgroup atoms found, cannot assign leaflets.")
        return membrane_map
    box = universe.dimensions[:3].copy() if universe.dimensions is not None else None
    if box is not None and not np.all(box > 0): box = None
    print(' Running MDAnalysis LeafletFinder')
    best = None
    for cutoff in cutoffs:
        if best is not None:
            print(f' {len(best[2])} leaflets could not be classified, retrying LeafletFinder with a {cutoff} Å cutoff')
        leaflets, no_mem_lipids = split_leaflet_groups(universe, headgroup_ag, cutoff)
        aggregates, unclassified = classify_leaflets(leaflets, box)
        n_classified = sum(len(group) for _, groups in aggregates for group in groups.values())
        if best is None or n_classified > best[3]:
            best = (aggregates, no_mem_lipids, unclassified, n_classified)
        if len(unclassified) == 0: break
    aggregates, no_mem_lipids, unclassified, _ = best
    if len(unclassified) > 0:
        warn(f"Could not classify {len(unclassified)} leaflet(s) with {', '.join(str(len(group)) for group in unclassified)} lipids, they will be unassigned")
        for group in unclassified:
            no_mem_lipids.extend(group.residues.resindices)
    membrane_map['n_mems'] = len(aggregates)

    for i, (aggregate_type, groups) in enumerate(aggregates):
        memb_data = {
            'type': aggregate_type,
            'leaflets': {name: group.residues.atoms.indices.tolist() for name, group in groups.items()},
            'polar_atoms': {name: group.atoms.indices.tolist() for name, group in groups.items()},
        }
        membrane_map['mems'][str(i)] = memb_data
        # Print leaflets stats
        print(f"Membrane {i} ({aggregate_type}):\n" + '\n'.join(
            f"    - {name.capitalize()} leaflet: {len(group)} lipids" for name, group in groups.items()))

    if len(no_mem_lipids) > 0:
        print(f"Unassigned lipids: {len(no_mem_lipids)}.")
        no_mem_lipids = universe.residues[no_mem_lipids].atoms.indices.tolist()
    membrane_map['no_mem_lipid'] = no_mem_lipids

    return membrane_map


def get_bilayers(membrane_map: dict | None) -> list[dict]:
    """Get the flat bilayers of a membrane mapping. Mappings before version 0.3.0 only have bilayers and no type."""
    if membrane_map is None: return []
    return [mem for mem in membrane_map['mems'].values() if mem.get('type', BILAYER) == BILAYER]


def display_membrane_mapping(mem_map: str, pdb: str, mem: str = '0'):
    """Display the leaflets of a membrane from the membrane mapping using nglview."""
    try:
        import nglview as nv  # type: ignore
    except ImportError:
        raise ImportError("nglview is required for displaying membrane mapping. Please install it using 'pip install nglview'.")

    mem_map = load_json(mem_map)
    membrane = mem_map['mems'][mem]
    no_mem = f"@{','.join(map(str, mem_map['no_mem_lipid']))}"

    view = nv.show_file(pdb)
    view.clear(0)
    view.add_point(selection='not protein', color='green', scale=1.5)
    leaflet_colors = [('blue', 'cyan'), ('yellow', 'orange')]
    for (name, atoms), (color, polar_color) in zip(membrane['leaflets'].items(), leaflet_colors):
        view.add_point(selection=f"@{','.join(map(str, atoms))}", color=color)
        polar = f"@{','.join(map(str, membrane['polar_atoms'][name]))}"
        view.add_ball_and_stick(selection=polar, color=polar_color, aspectRatio=4)
    view.add_point(selection=no_mem, color='red', scale=0.5)
    view.add_point(selection='protein', color='black', scale=0.5)
    return view
