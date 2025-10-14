import os
from os import chdir, getcwd 
import glob
import json
from pathlib import Path 
from typing import List
from model_workflow.utils.constants import NASSA_ANALYSES_CANALS

def generate_nassa_config ( 
        folder_path: list[str],
        seq_path: str,
        output_path: str,
        unit_len: int,
        n_sequences: int,
 ):
    nassa_config = {
                    'unit_name': 'hexamer',
                    'unit_len': 6,
                    'n_lines': 5000,
                    'tail': True,
                    'bimod': True,
                    'save_tables': True,
                    'save_plots': True,
                    'save_path': None,
                    'sequence_files': [],
                    'coordinate_info': {}
                }
                # If the sequence path is given, we will use it as the base path to search the sequences
    if seq_path:
        # As the base path of sequences is given, we will search for the sequences in the given path and create a list of them
        seq_path = os.path.abspath(seq_path)
        for sequence_file in os.listdir(seq_path):
            # We assume that the sequences are in fasta format 
            # AGUS: I think we should allow for other formats as: .seq, .txt
            if sequence_file.endswith('.fasta') or sequence_file.endswith('.fa'):
                nassa_config['sequence_files'].append(str(os.path.join(seq_path,sequence_file)))
                if n_sequences and len(nassa_config['sequence_files']) == n_sequences:
                    break
    # If the output path is given, we will use it as the save path
    if output_path:
        output_nassa_path = os.path.abspath(output_path)
        nassa_config['save_path'] = output_nassa_path
    # If the output path is not given, we will use the current path as the save path
    else:
        nassa_config['save_path'] = os.path.abspath('.')
    # If the unit length is given, we will use it as the unit length
    if unit_len:
        nassa_config['unit_len'] = unit_len
    # We will create the configuration file asuming that the coordinates are in the same path as the configuration file (helical_parameters folder)
    # The path given as argument -m must be the base path of the helical_parameters folder or the folder where the coordinates are
    # AGUS: habría que explorar más casos, no sé si será siempre así
    #actual_path = os.path.abspath(folder_path)
    #print('actual_path: ', actual_path)
    actual_path = getcwd()
    
    for path in folder_path:
        md_path = os.path.join(actual_path, path)
        if os.path.exists(os.path.join(md_path, 'helical')):
            # If canals + curves have previously been calculated, we will use these outputs
            # We will create a list of the different .ser archives that we have interest in according to the values of NASSA_ANALYSES_CANALS
            coordinates = []
            [coordinates.extend(value) for value in NASSA_ANALYSES_CANALS.values()]
            coordinates = list(set(coordinates))
            for seq_file in os.listdir(os.path.join(md_path, 'helical')):
                if seq_file.endswith('.ser'):
                    seq_file_coordinate = seq_file.split('_')[2].replace('.ser', '')
                    # Filter the archives with the correct coordinate
                    if seq_file_coordinate in coordinates:
                        if seq_file_coordinate not in nassa_config["coordinate_info"]:
                            nassa_config["coordinate_info"][seq_file_coordinate] = []
                        nassa_config["coordinate_info"][seq_file_coordinate].append(os.path.join(md_path, 'helical', seq_file))
                        if n_sequences:
                            if len(nassa_config['coordinate_info'][seq_file_coordinate]) == n_sequences:
                                continue
        else:
            # If the helical folder does not exist, we will search for the sequences in the given path
            # In this case, the sequences are in different folders, each folder is a coordinate
            # AGUS: pueden existir estos archivos en diferentes carpetas, por lo que habría que buscar en todas las carpetas ¿?
            folders = [f for f in os.listdir(actual_path) if os.path.isdir(os.path.join(actual_path, f))]
            all_coordinates = []
            for coordinate in NASSA_ANALYSES_CANALS.values():
                all_coordinates.extend(coordinate)
            all_coordinates = list(set(all_coordinates))
            for folder in folders:
                for coordinate in all_coordinates:
                    if coordinate == folder:
                        nassa_config["coordinate_info"][coordinate] = []
                        for seq_file in os.listdir(os.path.join(actual_path, folder)):
                            if seq_file.endswith('.ser'):
                                nassa_config["coordinate_info"][coordinate].append(os.path.join(actual_path, folder, seq_file))
                                if n_sequences:
                                    if len(nassa_config['coordinate_info'][coordinate]) == n_sequences:
                                        break
    
    # Sometimes, the number of .ser archives could be less than the whole sequence files, so we will filter the sequence files that have a .ser archive and sort them
    num_archives = None
    for coordinate, archives in nassa_config['coordinate_info'].items():
        count = len(archives)
        if num_archives is None:
            num_archives = count
        # If the number of archives is not the same, we will raise an error and print the number of archives for each coordinate
        elif count != num_archives:
            for coordinate, archives in nassa_config['coordinate_info'].items():
                count = len(archives)
                print(f"Number of coordinate archives in {coordinate}: {count}")
            raise ValueError("Not all coordinate archives have the same number")
    # If the number of archives for each coordinate is the same, we will sort and select the sequence files
    if len(nassa_config['sequence_files']) != num_archives:
        nassa_config['sequence_files'].sort()
        nassa_config['sequence_files'] = nassa_config['sequence_files'][:num_archives]

    # At this point is strange but If something wrong could happen, we will raise an error    
    if len(nassa_config['sequence_files']) != num_archives:
        print(f"Number of sequence files: {len(nassa_config['sequence_files'])}")
        print(f"Number of coordinate archives: {num_archives}")
        raise ValueError("The number of sequence files is not the same as the number of coordinate archives")
    # Save nassa_config as a JSON file
    if output_path:
        save_path = os.path.join(output_path, 'nassa.yml')
    else:
        save_path = os.path.join(os.path.abspath('.'), 'nassa.yml')
    
    with open(save_path, 'w') as f:
        json.dump(nassa_config, f)
    
    return save_path