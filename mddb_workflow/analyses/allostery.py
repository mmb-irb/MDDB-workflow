# Allostery analysis
#
# The allostery analysis is carried by the 'ComPASS' software. ComPASS builds a network where nodes
# are residues (proteins and nucleic acids) and edges connect residues which are close and whose
# dynamics are coupled. Couplings are the result of a PCA over generalized correlations,
# interactions (non-bonded contacts, salt bridges and hydrogen bonds) and communication propensity.
#
# From this network we keep:
#   - The network itself, including the centrality of every node and edge
#   - The allosteric hotspots: nodes with the highest centralities
#   - The communities: a partition of the network found with the Leiden algorithm
#   - The cliques: non-overlapping groups of residues where all residues are connected to each other
#     Cliques are found over a second network built with a larger distance cutoff
#
# Bheemireddy S, González Alemán R, Bignon E, Karami Y. Communication pathway analysis within
# protein-nucleic acid complexes. bioRxiv. 2025:2025-02.
# https://github.com/yasamankarami/compass
#
# This analysis was integrated by Claude

from os import mkdir, cpu_count
from os.path import exists
from shutil import rmtree
from subprocess import run, PIPE, STDOUT
import json
import re

import mdtraj as mdt

from mddb_workflow.tools.get_reduced_trajectory import get_reduced_trajectory
from mddb_workflow.utils.auxiliar import ToolError, save_json, warn
from mddb_workflow.utils.constants import GREY_HEADER, COLOR_END
from mddb_workflow.utils.constants import OUTPUT_ALLOSTERY_FILENAME
from mddb_workflow.utils.type_hints import *

# DANI: Esto es provisional hasta que ComPASS esté en conda
COMPASS_COMMAND = 'source ~/miniforge3/bin/activate ~/miniforge3/envs/compass; compass {}'

# Set the residues which ComPASS considers as nodes
# Note that these are the same criteria used by ComPASS internally
# Proteins are residues including a CA atom and nucleic acids are residues including a C5' atom
PROTEIN_NODE_ATOM = 'CA'
NUCLEIC_NODE_ATOM = "C5'"
DNA_PATTERN = re.compile(r'(5|3)?D([ATGC]){1}(3|5)?$')
RNA_PATTERN = re.compile(r'(3|5)?R?([AUGC]){1}(3|5)?$')

# Set the ComPASS input and output filenames
COMPASS_CONFIG_FILENAME = 'compass.cfg'
COMPASS_OUTPUT_DIRECTORY = 'output'
COMPASS_JOB_NAME = 'allostery'


def allostery(
    structure_file: 'File',
    trajectory_file: 'File',
    structure: 'Structure',
    pbc_selection: 'Selection',
    cg_selection: 'Selection',
    snapshots: int,
    output_directory: str,
    # ComPASS time scales linearly with the number of frames, but most of the time is spent in
    # steps which depend only in the number of residues, so we can afford many frames
    frames_limit: int = 1000,
    # Maximum distance (Å) between two residues to be connected in the network
    # Communities, centralities and hotspots are calculated over this network
    graph_cutoff: int = 5,
    # Maximum distance (Å) between two residues to be connected in the cliques network
    clique_cutoff: int = 10,
):
    """Perform the allostery analysis using ComPASS."""
    warn('The "allostery" analysis is not yet fully integrated. ComPASS is not installed in the conda enviornment')

    # Set the residues to be nodes of the network
    # Residues in PBC and coarse grain residues are excluded
    excluded_residue_indices = set(structure.get_selection_residue_indices(pbc_selection + cg_selection))
    candidate_selection = structure.select_protein() + structure.select_nucleic()
    candidate_residue_indices = sorted(structure.get_selection_residue_indices(candidate_selection))
    node_residue_indices = []
    for residue_index in candidate_residue_indices:
        if residue_index in excluded_residue_indices: continue
        residue = structure.residues[residue_index]
        if not is_compass_node(residue): continue
        node_residue_indices.append(residue_index)
    # ComPASS needs at least a few residues to make a network
    if len(node_residue_indices) < 3:
        print(' No protein or nucleic acid residues to analyze')
        return

    # Set the ComPASS working directory, where all ComPASS inputs and outputs will be
    compass_directory = f'{output_directory}/compass'
    if not exists(compass_directory): mkdir(compass_directory)

    # Write a structure with only the node residues
    # Residues are renumbered so every residue has a unique number in ComPASS outputs
    # Note that ComPASS sorts residues by chain and number, so this also guarantees the order is kept
    node_selection = structure.select_residue_indices(node_residue_indices)
    node_structure = structure.filter(node_selection)
    for r, residue in enumerate(node_structure.residues):
        residue.number = r + 1
        residue.icode = ''
    compass_structure_filepath = f'{compass_directory}/{structure_file.filename}'
    node_structure.generate_pdb_file(compass_structure_filepath)

    # Write a reduced trajectory with only the node residues atoms
    reduced_trajectory_filepath, step, _ = get_reduced_trajectory(
        structure_file,
        trajectory_file,
        snapshots,
        frames_limit,
    )
    print(' Filtering trajectory for ComPASS')
    compass_trajectory = mdt.load(reduced_trajectory_filepath, top=structure_file.path,
        atom_indices=node_selection.atom_indices)
    compass_trajectory_filepath = f'{compass_directory}/{trajectory_file.filename}'
    compass_trajectory.save(compass_trajectory_filepath)

    # Write the ComPASS configuration file
    # Paths are relative to the configuration file
    # Paths between residues are not calculated since they require custom sources and targets
    compass_config_filepath = f'{compass_directory}/{COMPASS_CONFIG_FILENAME}'
    with open(compass_config_filepath, 'w') as file:
        file.write(f"""[generals]
topology = {structure_file.filename}
trajectory = {trajectory_file.filename}
output_dir = {COMPASS_OUTPUT_DIRECTORY}
n_cores = {cpu_count()}
job_name = {COMPASS_JOB_NAME}

[non_bond]
non_bond_cut = 0.39

[salt_bridges]
NO_cut = 0.32

[hbonds]
DA_cut = 0.39
HA_cut = 0.25
DHA_cut = 90
heavy = S N O

[distance cutoffs]
Graph = {graph_cutoff}
Cliques = {clique_cutoff}

[paths]
find_path = False
sources =
targets =
""")

    # Run ComPASS
    print(' Running ComPASS')
    print(GREY_HEADER, end='')
    process = run(['bash', '-c', COMPASS_COMMAND.format(COMPASS_CONFIG_FILENAME)],
        cwd=compass_directory, stdout=PIPE, stderr=STDOUT)
    logs = process.stdout.decode()
    print(COLOR_END, end='')
    if process.returncode != 0 or 'Normal Termination' not in logs:
        print(logs)
        raise ToolError('Something went wrong with ComPASS')

    # Mine the ComPASS outputs
    network_directory = f'{compass_directory}/{COMPASS_OUTPUT_DIRECTORY}/network'
    graph_prefix = f'{network_directory}/graph_cutoff_{graph_cutoff}'
    clique_prefix = f'{network_directory}/graph_cutoff_{clique_cutoff}'

    # Read the network
    # Nodes are numbered from 0 in the same order that residues in the ComPASS structure
    with open(f'{graph_prefix}.json', 'r') as file:
        compass_graph = json.load(file)
    graph = compass_graph['graph']
    atom_mapping = compass_graph['atom_mapping']
    node_count = len(graph['nodes'])
    if node_count != len(node_residue_indices):
        raise ToolError(f'ComPASS found {node_count} nodes but we expected {len(node_residue_indices)}')
    # Make sure every node is the residue we expect
    # Residue numbers are enough since we made them unique
    # Note that residue names may not match since MDtraj standardizes them (e.g. HIE -> HIS)
    # Text outputs refer to residues by chain and number, so map these labels to nodes as well
    label_to_node = {}
    for node in range(node_count):
        residue_name, _, residue_number, chain_name = atom_mapping[str(node)]
        expected_residue = node_structure.residues[node]
        if residue_number != expected_residue.number:
            raise ToolError(f'ComPASS node {node} ({residue_name} {residue_number}) does not match'
                f' residue {expected_residue.name} {expected_residue.number}')
        label_to_node[(chain_name.strip(), residue_number)] = node

    def label_node(chain_name: str, residue_number: str) -> int:
        return label_to_node[(chain_name.strip(), int(residue_number))]

    # Read node centralities
    # Lines look like "Node (12,A)\t0.0031\t0.1562\t7"
    betweenness = [None] * node_count
    closeness = [None] * node_count
    degree = [None] * node_count
    with open(f'{graph_prefix}_centralities.txt', 'r') as file:
        next(file)
        for line in file:
            node_label, node_betweenness, node_closeness, node_degree = line.rstrip('\n').split('\t')
            residue_number, chain_name = re.match(r'Node \((\d+),(.*)\)', node_label).groups()
            node = label_node(chain_name, residue_number)
            betweenness[node] = float(node_betweenness)
            closeness[node] = float(node_closeness)
            degree[node] = int(node_degree)

    # Read edge centralities
    # Lines look like "(12,A)-(13,A)\t0.0044"
    edge_betweenness = {}
    with open(f'{graph_prefix}_edge_betweenness.txt', 'r') as file:
        next(file)
        for line in file:
            edge_label, value = line.rstrip('\n').split('\t')
            number_1, chain_1, number_2, chain_2 = re.match(r'\((\d+),(.*)\)-\((\d+),(.*)\)', edge_label).groups()
            edge = tuple(sorted([label_node(chain_1, number_1), label_node(chain_2, number_2)]))
            edge_betweenness[edge] = float(value)

    # Set the edges
    # Note that the edges key depends on the networkx version: 'links' in old versions and 'edges' in new ones
    links = graph['edges'] if 'edges' in graph else graph['links']
    edges = sorted((*sorted([link['source'], link['target']]), link['weight']) for link in links)

    # Read the hotspots
    # Lines look like "Node (12, A)"
    hotspots = []
    with open(f'{graph_prefix}_top_5_percent_nodes.txt', 'r') as file:
        next(file)
        for line in file:
            residue_number, chain_name = re.match(r'Node \((\d+),(.*)\)', line.strip()).groups()
            hotspots.append(label_node(chain_name, residue_number))

    # Read communities and cliques
    communities = read_groups(f'{graph_prefix}_communities_leiden.txt', label_node)
    cliques = read_groups(f'{clique_prefix}_cliques.txt', label_node)

    # Get the modularity from the logs
    modularity_match = re.search(r'modularity (-?\d+\.\d+)', logs)
    modularity = float(modularity_match.group(1)) if modularity_match else None

    # Set the output analysis
    # Note that edges, hotspots, communities and cliques refer to nodes by their position in 'residues'
    output_analysis = {
        'graph_cutoff': graph_cutoff,
        'clique_cutoff': clique_cutoff,
        'residues': node_residue_indices,
        'betweenness': betweenness,
        'closeness': closeness,
        'degree': degree,
        'hotspots': sorted(hotspots),
        'edges': {
            'source': [edge[0] for edge in edges],
            'target': [edge[1] for edge in edges],
            'weight': [edge[2] for edge in edges],
            'betweenness': [edge_betweenness[(edge[0], edge[1])] for edge in edges],
        },
        'communities': communities,
        'modularity': modularity,
        'cliques': cliques,
        'start': 0,
        'step': step,
    }
    save_json(output_analysis, f'{output_directory}/{OUTPUT_ALLOSTERY_FILENAME}')

    # Remove ComPASS files since they are heavy and we already have what we need
    rmtree(compass_directory)


def is_compass_node(residue: 'Residue') -> bool:
    """Check if a residue will be considered a node by ComPASS."""
    atom_names = set(atom.name for atom in residue.atoms)
    if PROTEIN_NODE_ATOM in atom_names: return True
    is_nucleic_name = DNA_PATTERN.match(residue.name) or RNA_PATTERN.match(residue.name)
    return bool(is_nucleic_name) and NUCLEIC_NODE_ATOM in atom_names


def read_groups(filepath: str, label_node: Callable) -> list[list[int]]:
    """Read a ComPASS communities or cliques file.
    Lines look like "Community 0: A_86, A_87, A_88".
    """
    groups = []
    with open(filepath, 'r') as file:
        next(file)
        for line in file:
            _, members = line.split(':', 1)
            nodes = []
            for member in members.split(','):
                chain_name, residue_number = member.strip().rsplit('_', 1)
                nodes.append(label_node(chain_name, residue_number))
            groups.append(sorted(nodes))
    return groups
