import os
from os.path import exists
import json
from model_workflow.utils.auxiliar import save_json
from model_workflow.utils.type_hints import *
from model_workflow.tools.get_pdb_frames import get_pdb_frames
from model_workflow.utils.structures import Structure
import subprocess
from model_workflow.utils.file import File
from collections import defaultdict
import networkx as nx
import matplotlib.pyplot as plt
from model_workflow.utils.constants import OUTPUT_LIGANDS_INTERACTIONS_FILENAME
import MDAnalysis as mda
from rdkit import DataStructs
import pandas as pd
import seaborn as sns
#### AGUS: hay que instalar en el env las librerias de arpeggio y prolif
import prolif as plf
from arpeggio.core import InteractionComplex

def ligand_interactions (
    structure: Structure,
    structure_file: str,
    output_directory: str,
    trajectory_file: str,
    topology_file: str,
    #actual_md: str,
    snapshots: int,
    frames_limit: int = 500,
    ligand_map: dict[str, str] = None,
    input_ligands: List[str] = None,
    mercy: bool = False,
) -> None:
    
    print('-> Running ligand interactions analysis')
    ### AGUS: este análisis utiliza el ligand map (que viene de ligmap) para minar los residuos identificados como ligandos en la estructura
    ### Si no hay ligand map, no se puede hacer este análisis
    # print('LIGAND MAP:',ligand_map)

    # Check if the ligand map is empty
    if not ligand_map:
        print('  No ligands to analyze')
        return
    
    '''
    ### AGUS: toda esta parte comentada incluye los análisis con arpeggio, el cual devuelve unos json que se podrían ver en el cliente de alguna forma
    ### AGUS: En total corre tres análisis: conteo de interacciones, persistencia de interacciones y red de interacciones
    ### AGUS: Todo el resto del código dentro de esta función (ligand_interactions) es un análisis alternativo que utiliza la librería prolif y MDAnalysis
    ### AGUS: Esta parte está a medias y el output son plots, habría que buscar la forma de devolver los datos en json para que se puedan ver en el cliente o viceversa

    ############# ARPEGGIO ANALYSIS #############
    # Check the ligand map in order to run the analysis for each ligand respectively
    data = {
        "data": []
    }
    
    # Iterate over the ligand map to obtain the ligands
    for ligand in ligand_map:
        # Create a dict with the ligand information
        ligand_data = {
            "name": ligand["name"],
        }
        # Run arpeggio tool and obtain the results of the different analyses to save them in a json file
        interactions_count, interactions_persistence, interactions_graph = get_arpeggio_trajectory_analysis (
            structure, 
            trajectory_file, 
            output_directory, 
            actual_md, 
            snapshots, 
            frames_limit, 
            ligand )
        # Add the results of the analyses to the ligand data
        ligand_data["interactions_count"] = interactions_count
        ligand_data["interactions_persistence"] = interactions_persistence
        ligand_data["interactions_graph"] = interactions_graph
        # Now, add this ligand dict to the general data dict to be exported and saved
        data["data"].append(ligand_data)
    
    # Save the data in a json file in the MD folder (not the ligands folder)
    json_path = os.path.join(actual_md, OUTPUT_LIGANDS_INTERACTIONS_FILENAME)
    print(f"Saving ligand interactions data to {json_path}")
    with open(json_path, "w") as f:
        json.dump(data, f, indent=2)

    # Check if the json file was created
    if not os.path.exists(json_path):
        raise InputError(f"Error: {json_path} was not created")
    # Remove the ligands folder to save space
    if os.path.exists(output_directory):
        shutil.rmtree(output_directory)
    ########## END ARPEGGIO ANALYSIS #############
    '''
    
    # load topology and trajectory
    u = mda.Universe(topology_file.path, trajectory_file.path)
    #### AGUS: solo he hecho las pruebas para el primero de los ligandos (debería hacerse para todos si hay más de un ligando en la estructura)
    ligand = ligand_map[0]  # Assuming we are only analyzing the first ligand in the map
    ligand_residue = ligand['residue_indices'][0]  # Get the first residue index of the ligand
    # create selections for the ligand and protein
    ligand_selection = u.select_atoms(f"resid {ligand_residue}")
    protein_selection = u.select_atoms("protein")
    print(f"Selected ligand: {ligand_selection}, Protein: {protein_selection}")
    protein_selection = u.select_atoms(
        "protein and byres around 20.0 group ligand", ligand=ligand_selection
    )
    print(f"Protein selection after applying ligand proximity: {protein_selection}")
    # create a molecule from the MDAnalysis selection
    ligand_mol = plf.Molecule.from_mda(ligand_selection, NoImplicit=True, force=True)
    protein_mol = plf.Molecule.from_mda(protein_selection, NoImplicit=True, force=True)
    # IF the bonds are not defined, we need to use this
    # ligand_selection.guess_bonds()
    # protein_selection.guess_bonds()

    # use default interactions
    fp = plf.Fingerprint()
    # AGUS: por defecto selecciona automaticamente los residuos vecinos con un radio  de 6 angstroms, 
    # AGUS: se puede cambiar con: plf.Fingerprint(vicinity_cutoff=7.0)

    # run on a slice of the trajectory frames: u.trajectory[::10] --> from begining to end with a step of 10
    fp.run(u.trajectory[::10], ligand_selection, protein_selection)
    fp.to_pickle("fingerprint.pkl")
    fp = plf.Fingerprint.from_pickle("fingerprint.pkl")
    
    df = fp.to_dataframe()
    # show only the 10 first frames
    print(df.head(10))

    # percentage of the trajectory where each interaction is present
    print(df.mean().sort_values(ascending=False).to_frame(name="%").T * 100)

    # same but we regroup all interaction types
    print(
        df.T.groupby(level=["ligand", "protein"])
        .sum()
        .T.astype(bool)
        .mean()
        .sort_values(ascending=False)
        .to_frame(name="%")
        .T
        * 100
    )

    # percentage of the trajectory where PiStacking interactions are present, by residue
    # (
    #     df.xs("PiStacking", level="interaction", axis=1)
    #     .mean()
    #     .sort_values(ascending=False)
    #     .to_frame(name="%")
    #     .T
    #     * 100
    # )

    # percentage of the trajectory where interactions with SER212 occur, by interaction type
    # print(
    #     df.xs("ARG588", level="protein", axis=1)
    #     .mean()
    #     .sort_values(ascending=False)
    #     .to_frame(name="%")
    #     .T
    #     * 100
    # )

    # percentage of the trajectory where each interaction type is present
    print(
        df.T.groupby(level="interaction")
        .sum()
        .T.astype(bool)
        .mean()
        .sort_values(ascending=False)
        .to_frame(name="%")
        .T
        * 100
    )

    # 10 residues most frequently interacting with the ligand
    print(
        df.T.groupby(level=["ligand", "protein"])
        .sum()
        .T.astype(bool)
        .mean()
        .sort_values(ascending=False)
        .head(10)
        .to_frame("%")
        .T
        * 100
    )

    bitvectors = fp.to_bitvectors()
    tanimoto_sims = DataStructs.BulkTanimotoSimilarity(bitvectors[0], bitvectors)
    #print(tanimoto_sims)

    # Tanimoto similarity matrix
    bitvectors = fp.to_bitvectors()
    similarity_matrix = []
    for bv in bitvectors:
        similarity_matrix.append(DataStructs.BulkTanimotoSimilarity(bv, bitvectors))
    similarity_matrix = pd.DataFrame(similarity_matrix, index=df.index, columns=df.index)
    # display heatmap
    fig, ax = plt.subplots(figsize=(3, 3), dpi=200)
    colormap = sns.diverging_palette(
        300, 145, s=90, l=80, sep=30, center="dark", as_cmap=True
    )
    sns.heatmap(
        similarity_matrix,
        ax=ax,
        square=True,
        cmap=colormap,
        vmin=0,
        vmax=1,
        center=0.5,
        xticklabels=5,
        yticklabels=5,
    )
    ax.invert_yaxis()
    plt.yticks(rotation="horizontal")
    #plt.patch.set_facecolor("white")
    plt.savefig(
        "ligand_interactions_heatmap.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()
    # %matplotlib ipympl
    # Plot and save barcode
    fp.plot_barcode()
    plt.savefig("barcode_plot.png", dpi=300, bbox_inches="tight")
    plt.close()
    # Plot and save ligand network
    fp.plot_lignetwork(ligand_mol)
    plt.savefig("ligand_network.png", dpi=300, bbox_inches="tight")
    plt.close()

    fp_count = plf.Fingerprint(count=True)
    fp_count.run(u.trajectory[::10], ligand_selection, protein_selection)
    fp_count.plot_lignetwork(ligand_mol, frame=0, display_all=True)
    plt.savefig("ligand_count.png", dpi=300, bbox_inches="tight")
    plt.close()

    print('-> Finished ligand interactions analysis')

    return


# Obtain a .cif file from the structure pdb file
def obtain_cif_from_structure_pdb (structure_file : str) -> None:
    # Set the name of the cif file to be the same as the structure file 
    output_cif = structure_file.replace('.pdb', '.cif')
    # Check if the cif file already exists to avoid unnecessary conversion
    if os.path.exists(output_cif):
        return output_cif
    # Use the gemmi tool to convert the pdb file to a cif file (structure.pdb -> structure.cif)
    try:
        result = subprocess.run(["gemmi", "convert", structure_file, output_cif], 
                   check=True,
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   text=True)
    except subprocess.CalledProcessError as e:
        print("ERROR: Falló la conversión con gemmi.")
        print("Código de salida:", e.returncode)
        print("Comando:", e.cmd)
        print("stdout:", e.stdout)
        print("stderr:", e.stderr)

    print(f"Converted: {structure_file} -> {output_cif}")
    #raise InputError(f"Error during the conversion: {e}")
    # The output cif path is returned to be used in the next steps
    return output_cif

# Run arpeggio to get the interactions of the ligand with the structure
def run_pdbe_arpeggio (structure_cif : str, selection : List[str], ligand_name : str, frame_number: str, output_directory: str) -> None:
    # Create the interaction complex and import the arpeggio tool to analyze the structure
    complex = InteractionComplex(structure_cif)
    complex.structure_checks()
    complex.address_ambiguities()
    complex.initialize()
    complex.run_arpeggio([selection],
                         interacting_cutoff=5, # cutoff for 'proximal' interactions
                         vdw_comp=0.1, # 'compensation' factor to address structural inconsistencies
                         include_sequence_adjacent=False
                         )
    contacts = complex.get_contacts()

    # Save the contacts in a json file
    save_json(contacts, f'{output_directory}/frame{frame_number}_{ligand_name}_contacts.json')

# Get pdb frames along the trajectory
def get_arpeggio_trajectory_analysis (
        structure : 'Structure',
        trajectory_file : File,
        output_directory : str,
        md_folder : str,
        snapshots : int,
        frames_limit : int,
        ligand: dict,
        ):
    
    # First, create the ligands folder if it does not exist to save the intermediary files
    if not exists(output_directory):
        os.mkdir(output_directory)
    # Copy the structure file in case we need to modify it (AGUS: puede que en un futuro necesite modicar esto)
    ligands_structure = structure.copy()
    ligands_structure_file = File(output_directory + '/ligands.pdb')
    ligands_structure.generate_pdb_file(ligands_structure_file.path)
    frames, step, count = get_pdb_frames(ligands_structure_file.path, trajectory_file.path, snapshots, frames_limit)
    for frame_number, current_frame_pdb in enumerate(frames):
        current_frame_cif = obtain_cif_from_structure_pdb(current_frame_pdb)
        # Obtain ligand name
        ligand_name = ligand['name']
        # Select chain and residue from the ligand map
        chain, residue = ligand["chain_residue_selection"][0]
        arpeggio_selection = f'/{chain}/{residue + 1}/'
        print("ligand name: ", ligand_name)
        print("arpeggio selection: ", arpeggio_selection)
        # Run arpeggio over the current frame
        run_pdbe_arpeggio(current_frame_cif, arpeggio_selection, ligand_name, frame_number, output_directory)
        
        # Remove the current cif file to save space
        os.remove(current_frame_cif)

    # Read all json files generated by arpeggio to obtain the interactions along the trajectory
    frames_data = []
    files = [
        fname for fname in os.listdir(output_directory)
            if fname.endswith("_contacts.json")
        ]

    # Sort the files by frame number
    files_sorted = sorted(files, key=lambda x: int(x.split('_')[0].replace("frame", "")))
    # Obtain the information of the json files and delete them to save space
    for fname in files_sorted:
        parts = fname.split('_')
        frame_number = int(parts[0].replace("frame", ""))
        ligand_name = parts[1]
        json_path = os.path.join(output_directory, fname)
        # Read the json file generated by arpeggio
        with open(json_path, "r") as f:
            data = json.load(f)
        # Add the data to the frame data list
        frames_data.append({
            "frame": frame_number,
            "ligand": ligand_name,
            "interactions": data
        })
        # Remove the json file
        os.remove(json_path)

    # First analaysis: Interaction types and subtypes count
    print(' |---> Analyzing ligand interactions count')
    interactions_count = analyze_interactions_count(frames_data, output_directory)

    # Second analysis: Ligand persistence interactions
    print(' |---> Analyzing ligand persistence interactions')
    interactions_persistence = analyze_ligand_persistence(frames_data, output_directory)

    # Third analysis: 
    #print(' |---> Analyzing ligand hydrogen bond interactions')
    #analyze_hbond_interactions(frames_data, output_directory)

    # Forth analysis:
    print(' |---> Analyzing ligand interactions web')
    interactions_graph = analyze_ligands_interactions_web(frames_data, output_directory)

    return interactions_count, interactions_persistence, interactions_graph


# Analysis of the ligand interactions web
def analyze_ligands_interactions_web (frame_data: List, output_directory: str) -> None:
    # Interacciones acumuladas
    interaction_counts = defaultdict(int)
    interaction_types = defaultdict(set)
    interaction_frames = defaultdict(set)
    for frame_idx, frame in enumerate(frame_data):
        interactions = frame.get("interactions", [])
        seen_in_this_frame = set()  # evita contar la misma interacción más de una vez por frame
        for interaction in interactions:
            # Ligand - Protein, not any other type of interaction (intramolecular, etc)
            if interaction.get("interacting_entities") != "INTER":
                continue
            # Begging and end of the interaction
            bgn = interaction["bgn"]
            end = interaction["end"]
            # Canonic name of the nodes: VAL_313_A (residue name, residue number, chain id)
            node1 = f"{bgn['label_comp_id']}_{bgn['auth_seq_id']}_{bgn['auth_asym_id']}"
            node2 = f"{end['label_comp_id']}_{end['auth_seq_id']}_{end['auth_asym_id']}"
            # Sort the nodes to avoid duplicates
            # We use a tuple to store the nodes in a hashable way
            pair = tuple(sorted([node1, node2]))
            interaction_frames[pair].add(frame_idx)
            # Multiple interactions can be present in the same frame, so we need to check if we have already counted this interaction
            if pair not in seen_in_this_frame:
                interaction_counts[pair] += 1
                for ctype in interaction.get("contact", []):
                    interaction_types[pair].add(ctype)
                seen_in_this_frame.add(pair)

    # Create a graph from the interaction counts
    G = nx.Graph()
    # Add nodes and edges with weights and labels
    for (n1, n2), weight in interaction_counts.items():
        G.add_node(n1)
        G.add_node(n2)
        G.add_edge(n1, n2, weight=weight, contact_types=list(interaction_types[(n1, n2)]))

    
    # Call the function to export the graph to a json file
    graph_data = export_graph_with_frames(G, interaction_counts, interaction_types, interaction_frames, f"{output_directory}/red_interacciones.json")
    return graph_data
    # Layout
    # pos = nx.spring_layout(G, seed=42)
    # edge_weights = [G[u][v]['weight'] for u, v in G.edges()]
    # max_weight = max(edge_weights)
    # normalized_weights = [1 + 4*(w / max_weight) for w in edge_weights] 
    # plt.figure(figsize=(12, 10))
    # nx.draw(
    #     G, pos,
    #     with_labels=True,
    #     node_size=500,
    #     width=normalized_weights,
    #     edge_color='gray',
    #     font_size=8
    # )

    # plt.title("Red de Interacciones (Proteína - Ligando)")
    # plt.savefig(f"{output_directory}/red_interacciones.png", dpi=300, bbox_inches='tight')


# Create a function to export the graph to a json file 
def export_graph_with_frames(G, interaction_counts, interaction_types, interaction_frames, output_path):
    data = {
        "nodes": [],
        "edges": []
    }

    for node in G.nodes():
        node_type = "ligand" if node.endswith('_B') else "residue"
        data["nodes"].append({
            "id": node,
            "type": node_type
        })

    for u, v in G.edges():
        pair = tuple(sorted([u, v]))
        frames = interaction_frames.get(pair, [])
        frame_ranges = group_consecutive_frames(frames) if frames else []
        data["edges"].append({
            "source": u,
            "target": v,
            "weight": interaction_counts.get(pair, 1),
            "contacts": list(interaction_types.get(pair, [])),
            #"frames": sorted(list(frames)),
            "frame_ranges": frame_ranges
        })

    with open(output_path, 'w') as f:
        json.dump(data, f, indent=2)

    return data

def group_consecutive_frames(frames):
    frames = sorted(frames)
    ranges = []
    start = prev = frames[0]

    for frame in frames[1:]:
        if frame == prev + 1:
            prev = frame
        else:
            ranges.append((start, prev))
            start = prev = frame
    ranges.append((start, prev))
    return [f"{s}" if s == e else f"{s}-{e}" for s, e in ranges]


'''
# Analysis of the ligand interactions web : SOLAMENTE PARA CONTAR LAS INTERACCIONES 1 VEZ (NO GENERA FRECUENCIAS)
def analyze_ligands_interactions_web_not_frecuency (frame_data: List, output_directory: str) -> None:
    # Interacciones acumuladas
    interaction_counts = defaultdict(int)
    interaction_types = defaultdict(set)
    seen = set()
    for frame in frame_data:
        interactions = frame.get("interactions", [])
        for interaction in interactions:
            # Ligand - Protein, not any other type of interaction (intramolecular, etc)
            if interaction.get("interacting_entities") != "INTER":
                continue
            # Begging and end of the interaction
            bgn = interaction["bgn"]
            end = interaction["end"]
            # Canonic name of the nodes: VAL_313_A (residue name, residue number, chain id)
            node1 = f"{bgn['label_comp_id']}_{bgn['auth_seq_id']}_{bgn['auth_asym_id']}"
            node2 = f"{end['label_comp_id']}_{end['auth_seq_id']}_{end['auth_asym_id']}"
            # Sort the nodes to avoid duplicates
            # We use a tuple to store the nodes in a hashable way

            pair = tuple(sorted([node1, node2]))
            print("pair: ",pair)
            if pair not in seen:
                seen.add(pair)
                interaction_counts[pair] += 1
                # Add the interaction type to the set of interaction types
                for ctype in interaction.get("contact", []):
                    interaction_types[pair].add(ctype)

    print("interaction counts: ",interaction_counts)
    G = nx.Graph()

    # Agregar nodos y aristas con pesos y etiquetas
    for (n1, n2), weight in interaction_counts.items():
        G.add_node(n1)
        G.add_node(n2)
        G.add_edge(n1, n2, weight=weight, contact_types=list(interaction_types[(n1, n2)]))

    # Layout
    pos = nx.spring_layout(G, seed=42)
    edge_weights = [G[u][v]['weight'] for u, v in G.edges()]

    plt.figure(figsize=(12, 10))
    nx.draw(
        G, pos,
        with_labels=True,
        node_size=500,
        width=edge_weights,
        edge_color='gray',
        font_size=8
    )
    plt.title("Red de Interacciones (Proteína - Ligando)")
    plt.savefig(f"{output_directory}/red_interacciones.png", dpi=300, bbox_inches='tight')
'''

# Analysis of the ligand persistence interactions
def analyze_ligand_persistence (frame_data: List, output_directory: str) -> None:
    # Define the interaction signature with the fields that arpeggio returns: begging and end, the interaction type and the contact
    def interaction_signature(interaction):
        bgn = interaction["bgn"]
        end = interaction["end"]
        contact = ",".join(sorted(interaction.get("contact", [])))
        return (
            interaction["type"],
            bgn["auth_asym_id"], bgn["auth_seq_id"], bgn["auth_atom_id"],
            end["auth_asym_id"], end["auth_seq_id"], end["auth_atom_id"],
            contact
        )

    # Find consecutive sequences of frames
    # We set 3 frames as the minimum length of a sequence
    # AGUS: 3 es representativo? pregnutar a los expertos
    def find_consecutive_sequences(frames, min_length=3):
        frames = sorted(frames)
        sequences = []
        temp = [frames[0]]
        for i in range(1, len(frames)):
            if frames[i] == frames[i-1] + 1:
                temp.append(frames[i])
            else:
                if len(temp) >= min_length:
                    sequences.append(temp.copy())
                temp = [frames[i]]
        if len(temp) >= min_length:
            sequences.append(temp)
        return sequences
    # Create a dictionary to store the presence of each interaction signature across frames
    interaction_presence = defaultdict(list)  # key: signature, value: list of frames
    for frame in frame_data:
        frame_num = frame["frame"]
        for interaction in frame["interactions"]:
            # Save the interaction with the important fields
            sig = interaction_signature(interaction)
            # Add the frame number to the list of frames for this interaction signature
            interaction_presence[sig].append(frame_num)

    # Create the dictionary to store the persistent interactions
    # The same interactions may have multiple blocks of frames because they could be interrupted
    # by other interactions, so we need to group them e.g. 1-3, 5-10
    persistent_interactions_grouped = defaultdict(list)
    for sig, frames in interaction_presence.items():
        # Find the consecutive sequences of frames for this interaction signature
        sequences = find_consecutive_sequences(frames)
        # If there are consecutive sequences, save the interaction and the frames
        if sequences:
            interaction_example = next(
                (intx for f in frame_data for intx in f["interactions"] if interaction_signature(intx) == sig),
                None
            )
            # If we have an interaction example, save it with the sequences of frames
            if interaction_example:
                persistent_interactions_grouped[sig].append({
                    "interaction": interaction_example,
                    "frames_blocks": sequences
                })

    # Save the list with the persistent interactions that are repeated in the trajectory
    # A list of repeated frames is saved for each interaction if it is repeated in the trajectory
    final_grouped = []
    for group in persistent_interactions_grouped.values():
        if group:
            interaction_data = group[0]["interaction"]
            all_blocks = []
            for g in group:
                all_blocks.extend(g["frames_blocks"])
            # Convert the list of lists of frames to a list of ranges
            all_blocks_ranges = format_frame_ranges_from_lists(all_blocks)            
            # Save the blocks of frames for this interaction
            final_grouped.append({
                "interaction": interaction_data,
                "frames": all_blocks_ranges,
            })
    # Save the persistent interactions to a JSON file
    with open(f"{output_directory}/persistent_interactions.json", "w") as f:
        json.dump(final_grouped, f, indent=2)

    return final_grouped

# Format the frame ranges from a list of frames
def format_frame_ranges_from_lists(frame_lists):
    ranges = []
    for group in frame_lists:
        if not group:
            continue
        start, end = group[0], group[-1]
        if start == end:
            ranges.append(str(start))
        else:
            ranges.append(f"{start}-{end}")
    return ranges

# Analysis of the hydrogen bond interactions
def analyze_hbond_interactions (frame_data: List, output_directory: str) -> None:
    # Initialize a dictionary to store the hydrogen bond interactions
    hbond_interactions = []
    # Iterate over each frame and count interactions
    for frame in frame_data:
        # Obtain frame number
        frame_number = frame["frame"]
        # Iterate over the interactions in each frame
        for interaction in frame["interactions"]:
            i_contact = interaction.get("contact", None)
            if i_contact:
                for contact in i_contact:
                    if contact == "hbond":
                        # Count the hydrogen bond interaction type
                        hbond_interactions.append(interaction)
    
    print(hbond_interactions)
    raise
    # Save the hydrogen bond interactions to JSON files
    with open(f"{output_directory}/hbond_interactions.json", "w") as f:
        json.dump(hbond_interactions, f, indent=2)

# Analysis of the interactions count along the trajectory
def analyze_interactions_count(frame_data: List, output_directory: str) -> None:
    known_types = [
        "atom-atom",
        "atom-plane",
        "plane-plane",
        "group-group",
        "group-plane"
    ]

    combined_data = []

    for frame in frame_data:
        frame_number = frame["frame"]
        interaction_counts = defaultdict(int)
        subtype_counts = {itype: defaultdict(int) for itype in known_types}

        for interaction in frame["interactions"]:
            i_type = interaction.get("type", "")
            if i_type not in known_types:
                continue

            interaction_counts[i_type] += 1

            contact = interaction.get("contact", [])
            if contact:
                for subtype in contact:
                    subtype_counts[i_type][subtype] += 1
            else:
                # Contador especial para interacciones sin contacto definido
                subtype_counts[i_type]["__empty__"] += 1

        # Formatear los subtipos como dicts simples
        formatted_subtypes = {
            i_type: dict(subtype_counts[i_type])
            for i_type in known_types if subtype_counts[i_type]
        }

        combined_data.append({
            "frame": frame_number,
            "types": dict(interaction_counts),
            "subtypes": formatted_subtypes
        })

    with open(f"{output_directory}/interactions_over_time.json", "w") as f:
        json.dump(combined_data, f, indent=2)
    
    return combined_data