import numpy as np
from collections import Counter
from MDAnalysis.analysis.leaflet import LeafletFinder
from MDAnalysis.lib.distances import capped_distance
from mddb_workflow.utils.auxiliar import save_json, load_json, warn
from mddb_workflow.utils.constants import LIPIDS_RESIDUE_NAMES, STEROID_RESIDUE_NAMES
from mddb_workflow.utils.type_hints import *

MEMBRANE_MAPPING_VERSION = '0.3.0'
# Amphiphilic aggregate types
BILAYER = 'bilayer'
VESICLE = 'vesicle'
MICELLE = 'micelle'
# Leaflet of vesicles for lipids between the outer and the inner leaflets
MIDPLANE = 'midplane'
# Minimum number of lipids in a leaflet, 24 lipids like in FATSLiM
MIN_LEAFLET_LIPIDS = 24
# A leaflet is flat when sqrt(smallest / largest eigenvalue) of its polar atom positions is below this
FLATNESS_THRESHOLD = 0.3
# Absolute cosine between the head to tail vector and the reference direction for a lipid to be clearly oriented
STRONG_ORIENTATION = 0.5
# Fraction of lipids clearly oriented against the majority from which a group is considered mixed (merged leaflets)
MIXED_ORIENTATION = 0.25
# Minimum fraction of lipids clearly oriented in the same direction to decide the orientation of a leaflet
CLEAR_ORIENTATION = 0.5
# Minimum absolute cosine between normals of the two leaflets of a bilayer
PARALLEL_NORMALS = 0.9
# Maximum distance (Å) between the centroids of the two leaflets of a bilayer along the normal
MAX_BILAYER_SEPARATION = 80
# Maximum distance between the centers of vesicle leaflets, as a fraction of the inner leaflet radius
VESICLE_CENTER_TOLERANCE = 0.5
# Lipids with atoms closer than this distance (Å) are in contact
CONTACT_DISTANCE = 8
# An unassigned lipid joins the leaflet it touches if it has this many times the contacts with the next leaflet
CONTACT_DOMINANCE = 2


def generate_membrane_mapping(
    inchikeys: dict[str, 'InChIKeyData'],
    lipid_references: dict[str, dict],
    universe: 'Universe',
    cg_residues: list[int],
    output_file: 'File',
    input_forced_lipids: Optional[list[str]] = None,
) -> dict:
    """Generate a list of residue numbers of membrane components from a given structure and topology file.

    Args:
        inchikeys: Dictionary mapping InChI keys to InChIKeyData objects.
        lipid_references: Dictionary mapping InChI keys to lipid database references.
        structure_file: File object representing the structure file (e.g., PDB, GRO).
        universe: MDAnalysis Universe object containing the structure and topology.
        cg_residues: List of residue indices for coarse-grained residues.
        output_file: Output JSON file to save membrane mapping.
        input_forced_lipids: Residue names to be treated as lipids in coarse grained systems,
            in addition to the known ones (all atom forced lipids are already in the lipid references).

    Returns:
        dict: A dictionary containing the membrane mapping.
            - n_mems (int): Number of amphiphilic aggregates (membranes) detected.
            - mems (dict): Dictionary where each key is a membrane index (str) and value contains:
              - type (str): Aggregate type: 'bilayer' (flat), 'vesicle' or 'micelle'.
              - leaflets (dict): Atom indices of each leaflet. Leaflet names are 'top' and 'bot' for bilayers,
                'outer', 'inner' and 'midplane' (lipids between both) for vesicles and 'outer' for micelles.
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
        headgroup_ag = coarse_grain_headgroups(universe, input_forced_lipids)
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


def coarse_grain_headgroups(universe: 'Universe', forced_lipids: Optional[list[str]] = None) -> dict:
    """Generate membrane mapping for coarse-grained systems using MDAnalysis LeafletFinder with P atoms as headgroups."""
    # Find all lipid residues in the system
    lipid_atoms = headgroup_ag = universe.atoms[[]]
    for resname in set(LIPIDS_RESIDUE_NAMES) | set(forced_lipids or []):
        lipid = universe.select_atoms(f'resname {resname}')
        if len(lipid) > 0:
            lipid_atoms += lipid.atoms

    if len(lipid_atoms) == 0:
        print("No coarse-grained lipid atoms found in the structure.")
        return None
    print(f"Found {len(lipid_atoms.residues)} lipid residues by name.")
    # Use the phosphate bead as headgroup or, if there is not a single one (e.g. cholesterol), the first bead
    # Martini lipid beads are ordered from the head to the tails
    headgroup_indices = []
    for residue in lipid_atoms.residues:
        phosphates = residue.atoms.select_atoms('name P*')
        headgroup_indices.append(phosphates[0].index if len(phosphates) == 1 else residue.atoms[0].index)
    headgroup_ag = universe.atoms[headgroup_indices]
    return headgroup_ag


def split_leaflet_groups(universe: 'Universe', headgroup_ag: 'AtomGroup', cutoff: float, pbc: bool) -> tuple[list, list]:
    """Run LeafletFinder and split its groups into leaflets and residue indices of unassigned lipids.
    With periodic boundary conditions leaflets crossing the box boundary (e.g. vesicles) are not split.
    """
    L = LeafletFinder(universe, headgroup_ag, cutoff=cutoff, pbc=pbc)
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
    The direction is undecided (0) when not enough lipids are clearly oriented in the same direction
    (e.g. small or rough bilayer patches which are not flat, whose tails are not radial).
    Sterols are ignored for the orientation since they are short and they may be not oriented
    (e.g. cholesterol in the core of lipid nanoparticles), unless there are only sterols.
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
    not_sterols = ~np.isin(leaflet.resnames, list(STEROID_RESIDUE_NAMES))
    if np.any(not_sterols): cosines = cosines[not_sterols]
    positive = np.mean(cosines > STRONG_ORIENTATION)
    negative = np.mean(cosines < -STRONG_ORIENTATION)
    direction = 0
    if max(positive, negative) >= CLEAR_ORIENTATION: direction = 1 if positive > negative else -1
    return {
        'group': leaflet,
        'flat': flat,
        'centroid': center + shift,
        'normal': normal,
        'radius': radii.mean(),
        'coherent': min(positive, negative) < MIXED_ORIENTATION,
        # Flat: tails along (+1) or against (-1) the normal. Curved: tails outwards (+1) or inwards (-1). Undecided: 0
        'direction': direction,
    }


def facing_separation(first: dict, second: dict, box: np.ndarray | None) -> np.ndarray | None:
    """Get the separation from a flat leaflet to the closest periodic image of another one where their tails face each other.

    In-plane components follow the minimum image convention. Along the normal, the closest image may be back to back
    (e.g. leaflets 40 Å apart in a 70 Å high box are 30 Å apart through the box boundary), so the images of the box
    axis closest to the normal are checked. None is returned if no image faces or they are too far apart.
    """
    separation = minimum_image(second['centroid'] - first['centroid'], box)
    images = [separation]
    if box is not None:
        axis = np.argmax(np.abs(first['normal']))
        images = [separation + shift * box[axis] * np.eye(3)[axis] for shift in (-1, 0, 1)]
    best, best_distance = None, None
    for image in images:
        distance = np.dot(image, first['normal'])
        # Tails of each leaflet must point towards the other leaflet
        first_faces = np.sign(distance) == first['direction']
        second_faces = np.sign(np.dot(-image, second['normal'])) == second['direction']
        if not (first_faces and second_faces) or abs(distance) > MAX_BILAYER_SEPARATION: continue
        if best is None or abs(distance) < best_distance:
            best, best_distance = image, abs(distance)
    return best


def classify_leaflets(leaflets: list, box: np.ndarray | None) -> tuple[list[tuple[str, dict]], list, list]:
    """Group leaflets into amphiphilic aggregates: flat bilayers and vesicles.

    - Vesicle: a curved leaflet with tails inwards (outer) and a smaller concentric curved leaflet with tails outwards (inner).
    - Bilayer: two flat leaflets with parallel normals whose tails face each other.

    Leaflets with an undecided lipid orientation are not classified.
    Returns the aggregates as (type, leaflets by name), the leaflets which could not be classified and the curved
    leaflets with tails inwards and no partner. These are micelles, unless their partner is found later (see find_leaflets).
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
        empty = outer['group'].universe.atoms[[]]
        aggregates.append((VESICLE, {'outer': outer['group'], 'inner': inner['group'], MIDPLANE: empty}))

    flats = [desc for desc in coherent if desc['flat'] and desc['direction'] != 0]
    for first in flats:
        if id(first) in used: continue
        best, best_distance, best_separation = None, None, None
        for second in flats:
            if second is first or id(second) in used: continue
            if abs(np.dot(first['normal'], second['normal'])) < PARALLEL_NORMALS: continue
            separation = facing_separation(first, second, box)
            if separation is None: continue
            distance = abs(np.dot(separation, first['normal']))
            if best is None or distance < best_distance:
                best, best_distance, best_separation = second, distance, separation
        if best is None: continue
        used.update([id(first), id(best)])
        # The top leaflet is the one further along the normal
        top, bot = (best, first) if np.dot(best_separation, first['normal']) > 0 else (first, best)
        aggregates.append((BILAYER, {'top': top['group'], 'bot': bot['group']}))

    lone_curved_in = [desc['group'] for desc in curved_in if id(desc) not in used]
    used.update(id(desc) for desc in curved_in)
    unclassified += [desc['group'] for desc in coherent if id(desc) not in used]
    return aggregates, unclassified, lone_curved_in


def find_leaflets(universe: 'Universe', headgroup_ag: list[int] | None, cutoffs: tuple[float, ...] = (15.0, 12.0, 10.0)) -> dict:
    """Find leaflets using MDAnalysis LeafletFinder and classify them into amphiphilic aggregates.

    When some leaflets can not be classified (e.g. two leaflets too close at a fusion stalk are
    merged by LeafletFinder) only their lipids are retried with the next smaller cutoff,
    so aggregates already classified with a larger cutoff are kept.
    Curved leaflets with tails inwards and no partner are micelles once there is nothing left to retry,
    since their partner may be found with a smaller cutoff (e.g. the inner leaflet of a vesicle merged
    with lipids in its lumen).
    """
    membrane_map = {'n_mems': 0, 'mems': {}, 'no_mem_lipid': [], 'version': MEMBRANE_MAPPING_VERSION}
    if not headgroup_ag:
        print("No headgroup atoms found, cannot assign leaflets.")
        return membrane_map
    box = universe.dimensions[:3].copy() if universe.dimensions is not None else None
    if box is not None and not np.all(box > 0): box = None
    print(' Running MDAnalysis LeafletFinder')
    aggregates = []
    no_mem_lipids = []
    lone_curved_in = []
    pending_headgroups = headgroup_ag
    for i, cutoff in enumerate(cutoffs):
        if i > 0:
            print(f' {len(unclassified)} leaflets could not be classified, retrying their lipids with a {cutoff} Å cutoff')
        leaflets, small_group_lipids = split_leaflet_groups(universe, pending_headgroups, cutoff, pbc=box is not None)
        no_mem_lipids.extend(small_group_lipids)
        found_aggregates, unclassified, lone_curved_in = classify_leaflets(leaflets + lone_curved_in, box)
        aggregates.extend(found_aggregates)
        if len(unclassified) == 0: break
        pending_headgroups = sum(unclassified[1:], unclassified[0])
    aggregates.extend((MICELLE, {'outer': group}) for group in lone_curved_in)
    dimensions = universe.dimensions if box is not None else None
    if len(unclassified) > 0:
        warn(f"Could not classify {len(unclassified)} leaflet(s) with {', '.join(str(len(group)) for group in unclassified)} lipids, they will be unassigned")
        for group in unclassified:
            no_mem_lipids.extend(group.residues.resindices)
    no_mem_lipids = assign_lipids_by_contact(aggregates, no_mem_lipids, headgroup_ag, dimensions)
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


def assign_lipids_by_contact(aggregates: list[tuple[str, dict]], lipid_resindices: list[int],
        headgroup_ag: 'AtomGroup', dimensions: np.ndarray | None) -> list[int]:
    """Add unassigned lipids to the leaflet they are in contact with.

    Lipids out of LeafletFinder groups (e.g. sterols in the core of lipid nanoparticles) are added to a leaflet
    when they have clearly more atom contacts with it than with any other leaflet. It is repeated until no more
    lipids are added, so lipids only touching other added lipids (e.g. sterol clusters) are also added.
    Lipids in contact with both leaflets of a vesicle alike (e.g. cholesterol between leaflets) are added to its
    midplane leaflet. Lipids in contact with several leaflets of different aggregates or with none remain unassigned.
    Leaflet groups are updated in place and the residue indices of lipids which remain unassigned are returned.
    """
    if len(aggregates) == 0 or len(lipid_resindices) == 0: return lipid_resindices
    universe = headgroup_ag.universe
    # Unassigned lipids come from LeafletFinder groups, so all of them have a polar atom
    polar_atoms = np.full(len(universe.residues), -1)
    polar_atoms[headgroup_ag.resindices] = headgroup_ag.indices
    remaining = np.array(lipid_resindices, dtype=int)
    leaflets = [(groups, name) for _, groups in aggregates for name in groups]
    # Leaflet index of every atom, -1 for atoms out of leaflets
    owners = np.full(len(universe.atoms), -1)
    for leaflet_index, (groups, name) in enumerate(leaflets):
        owners[groups[name].residues.atoms.indices] = leaflet_index
    # Midplane leaflet index of the aggregate of every leaflet, -1 for aggregates without midplane
    midplane_indices = {id(groups): index for index, (groups, name) in enumerate(leaflets) if name == MIDPLANE}
    midplanes = np.array([midplane_indices.get(id(groups), -1) for groups, _ in leaflets] + [-1])
    # Atom contacts of every remaining lipid with every leaflet, with an extra empty leaflet to rank at least two
    contacts = np.zeros((len(remaining), len(leaflets) + 1), dtype=int)
    rows = np.full(len(universe.residues), -1)
    # Contacts are counted only once: in the first round with all leaflet atoms, then with the atoms of the added lipids
    new_atoms = np.where(owners >= 0)[0]
    while len(remaining) > 0 and len(new_atoms) > 0:
        rows[remaining] = np.arange(len(remaining))
        unassigned_atoms = universe.residues[remaining].atoms
        pairs = capped_distance(unassigned_atoms.positions, universe.atoms[new_atoms].positions, CONTACT_DISTANCE,
            box=dimensions, return_distances=False)
        np.add.at(contacts, (rows[unassigned_atoms.resindices[pairs[:, 0]]], owners[new_atoms[pairs[:, 1]]]), 1)
        # Rank the leaflets of every lipid by contacts
        ranking = np.argsort(-contacts, axis=1, kind='stable')
        first, second = ranking[:, 0], ranking[:, 1]
        first_contacts = np.take_along_axis(contacts, ranking[:, :1], axis=1)[:, 0]
        second_contacts = np.take_along_axis(contacts, ranking[:, 1:2], axis=1)[:, 0]
        dominant = first_contacts >= CONTACT_DOMINANCE * second_contacts
        # Lipids touching two leaflets of the same vesicle alike are in its midplane
        same_vesicle = (midplanes[first] == midplanes[second]) & (midplanes[first] >= 0)
        targets = np.where(dominant, first, np.where(same_vesicle, midplanes[first], -1))
        targets[first_contacts == 0] = -1
        added = targets >= 0
        for leaflet_index in np.unique(targets[added]):
            groups, name = leaflets[leaflet_index]
            resindices = remaining[targets == leaflet_index]
            merged = groups[name] + universe.atoms[polar_atoms[resindices]]
            groups[name] = merged[np.argsort(merged.indices)]
            owners[universe.residues[resindices].atoms.indices] = leaflet_index
        new_atoms = universe.residues[remaining[added]].atoms.indices
        remaining, contacts = remaining[~added], contacts[~added]
    n_assigned = len(lipid_resindices) - len(remaining)
    if n_assigned > 0:
        print(f' Assigned {n_assigned} unassigned lipids to the leaflets they are in contact with')
    return remaining.tolist()


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
