import numpy as np
import os
from os.path import exists
import pandas as pd
import math
import subprocess
from shutil import move
import glob

from model_workflow.utils.auxiliar import save_json, store_binary_data
from model_workflow.utils.type_hints import *

from model_workflow.utils.nucleicacid import NucleicAcid
from model_workflow.tools.nassa_loaders import load_sequence2
from model_workflow.tools.nassa_loaders import load_serfile
conda_prefix = os.environ['CONDA_PREFIX']
# If this path does not exist then it means curves is not installed
curves_path = conda_prefix + '/.curvesplus'
# Note that this is not a file, but the prefix to 3 different files
standard_prefix = curves_path + '/standard'

# TAKE INTO ACCOUNT THAT EL CODE WAS DEVELOP AND IMPLEMENTED IN THIS WORKFLOW USING AS A TEMPLATE CODE USED IN IRBBARCELONA BIOBB 
# HERE THERE IS THE LINK RELATED TO THEIR WEBPAGE WITH MORE WORKFLOWS COMPUTING OTHER STUFF
# https://mmb.irbbarcelona.org/biobb/workflows

# HERE THERE IS THE LINK RELATED TO THE GITHUB WEBPAGE WITH THE PYTHON CODE
# https://github.com/bioexcel/biobb_wf_dna_helparms

# HERE THERE IS THE LINK RELATED TO THE INFORMATION OF HOW TO EXECUTE CURVES+ AND CANALS AND INTERPRET THEIR OUTPUTS AS WELL AS THE INPUT NEEDED FOR BOTH
# http://gbio-pbil.ibcp.fr/tools/curves_plus/canal-user-guide.html

# List with the different files names that correspond to single bases
hp_singlebases = [
    "shear",
    "stagger",
    "stretch",
    "buckle",
    "opening",
    "propel",
    "inclin",
    "tip",
    "xdisp",
    "ydisp"
]

# List with the different files names that correspond to base pairs
hp_basepairs = [
    "majw",
    "majd",
    "minw",
    "mind",
    "rise",
    "shift",
    "slide",
    "roll",
    "tilt",
    "twist"
]

# List with the different files names which we will use degrees units
hp_angular = [
    "roll",
    "tilt",
    "twist",
    "buckle",
    "opening",
    "propel",
    "inclin",
    "tip"
]

# List with the different files names which we will use amstrongs units
hp_translational = [
    "shear",
    "stagger",
    "stretch",
    "xdisp",
    "ydisp",
    "majw",
    "majd",
    "minw",
    "mind",
    "rise",
    "shift",
    "slide"
]

# List with the different files names corresponding to Backbone block of Helical Parameters
hp_backbone = [
    "alphac",
    "alphaw",
    "betac",
    "betaw",
    "gammac",
    "gammaw",
    "deltac",
    "deltaw",
    "chic",
    "chiw",
    "phasec",
    "phasew",
    "epsilc",
    "epsilw",
    "zetac",
    "zetaw"
]

helical_parameters = hp_basepairs + hp_singlebases + hp_backbone
baselen = 0
hp_unit = ""

def helical_parameters (
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_analysis_filename : str,
    structure : 'Structure',
    structure_filename : str,
    frames_limit : int,
    dna_selection : str = None
):
    actual_path = os.getcwd()
    print('-> Running helical analysis')

    # Get a selection from both chains in the B-DNA
    selection = None
    # In case a selection was provided then use it
    if dna_selection:
        selection = structure.select(dna_selection, syntax='vmd')
    # Otherwise, get all nucleic acids
    else:
        selection = structure.select_nucleic()

    # If there is nothing to analyze then we return here
    if not selection:
        print(' There are no nucleic acids')
        return

    # Check curves is installed
    if not exists(curves_path):
        raise SystemExit(' Cannot find curves path. Is Curves+ installed?')

    # Get the sequence from the selected chains
    chain_indices = structure.get_selection_chain_indices(selection)
    if len(chain_indices) != 2:
        raise SystemExit('We have a different number of chains than 2. Is this canonical B-DNA? Make sure every strand has an independent chain')
    sequences = []
    residue_index_ranges = []
    for c, chain_index in enumerate(chain_indices):
        chain = structure.chains[chain_index] # Obtain Chain of nucleic acids
        sequences.append(chain.get_sequence()) # Retrieve its sequence
        first_residue_index = chain.residue_indices[0] +1
        last_residue_index = chain.residue_indices[-1] +1
        residue_index_range = (first_residue_index, last_residue_index) if c == 0 else (last_residue_index, first_residue_index)
        residue_index_ranges.append(residue_index_range)
    
    # Call function to execute Curves+ and Canals software to generate the desired output in order to do different calculations
    folder_path = False
    if '/' in input_trajectory_filename:
        folder_path = input_trajectory_filename.split('/')[-2]

    # Run the hydrogen bond analysis to generate .dat file
    # hydrogen_bonds(
    #     input_topology_filename,
    #     input_trajectory_filename,
    #     output_analysis_filename,
    #     folder_path,
    #     structure_filename,
    # )

    # terminal_execution(input_trajectory_filename, structure_filename, residue_index_ranges, sequences[0],folder_path)
    # Save in a dictionary all the computations done by the different functions called by send_files function
    dictionary_information = send_files(sequences[0], frames_limit, folder_path)
    # Set the path into the original directory outside the folder helicalparameters
    os.chdir(actual_path)
    # Convert the dictionary into a json file
    # DANI: Estamos teniendo NaNs que no son soportados más adelante, hay que eliminarlos
    save_json(dictionary_information, output_analysis_filename)
    
# Function to execute cpptraj and calculate .dat file (Hydrogen Bonds)
def hydrogen_bonds(
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_analysis_filename : str,
    folder_path : str,
    input_structure_filename : str,
):
    # This analysis is done using cpptraj and its done by two steps
    # 1. Crate the cpptraj.in file to generate the dry trajectory (cpptraj_dry.inpcrd)
    # parm_path = f"/orozco/projects/ABCix/ProductionFiles" 
    # setup_path = f"/orozco/projects/ABCix/ProductionFiles/SETUP" 

    # # Config cpptraj
    # cpptraj_in_file = f"{folder_path}/cpptraj.in"
    # cpptraj_out_file = f"{folder_path}/cpptraj_dry.inpcrd"

    # if not os.path.exists(cpptraj_out_file):


    #     # Exec cpptraj
    #     cpptraj_command = f"cpptraj {cpptraj_in_file}" 
    #     res = subprocess.run(cpptraj_command, shell=True)


    ################
    # 2. From the dry.inpcrd file we will generate the .dat file
    # Create the helical parameters folder
    helical_parameters_folder = f"{folder_path}/helical_parameters"
    # If the folder already exists dont create it again
    if not os.path.exists(helical_parameters_folder):
        os.mkdir(helical_parameters_folder)
    output_dat_file = os.path.join(helical_parameters_folder, "mdf.nahbonds.dat")

    print(f"Saving .dat file to: {output_dat_file}")

    # Crear el contenido del input de cpptraj usando las variables
    cpptraj_script = f"""
    # Load the topology file
    parm {input_topology_filename}

    # Load the initial reference structure
    trajin {input_structure_filename}

    # Load the trajectory
    trajin {input_trajectory_filename}

    # Run the NA structure analysis
    nastruct NA

    # Run the command to process the trajectory
    run

    # Write out hydrogen bond information to a file
    writedata {output_dat_file} NA[hb]

    # Quit
    quit
    """
    print("Executing cpptraj...")
    # Execute cpptraj with the generated script
    subprocess.run(["cpptraj"], input=cpptraj_script, text=True, check=True)

    # Check if the output file was created
    if not os.path.exists(output_dat_file):
        raise SystemExit(f"Error: {output_dat_file} was not created. Check the cpptraj command and input files.")

    # Renumber the dat file because the first line is repeated
    # AGUS: con esto tuvimos muchos problemas en el proyecto de ABCix y al final optamos por eliminar la primera línea y renumerar todo 
    # AGUS: además, también sirvió para comprimir el archivo y leerlo mejor
    def renumber_dat_file(file_path):
        with open(file_path, 'r') as file:
            lines = file.readlines()  

        counter = 1
        first = True

        renumbered_lines = []
        for line in lines:
            if line.strip().startswith('#'):
                renumbered_lines.append(line)
            elif first:
                first = False
                continue
            else:
                renumbered_lines.append(f"{counter} {' '.join(map(str, line.split()[1:]))}\n")
                counter += 1

        with open(file_path, 'w') as file:
            file.writelines(renumbered_lines)
        print("Renumbering the dat file...")
    # Renumber the dat file
    renumber_dat_file(output_dat_file)
    
    # Set the function to compress the dat file to .bin file
    def compress_dat_to_bin(output_dat_file):
        data = []
        with open(output_dat_file, 'r') as file:
            # Skip the first line since it is only the headers
            next(file)
            for line in file:
                # Skip the first line since it is only the row number
                numbers = [ int(s) for s in line.strip().split()[1:] ]
                data += numbers

        # Data is a list of numeric values
        # Bit size is the number of bits for each value in data to be occupied
        def store_bits_dat (data : list, bit_size : int, filepath : str):
            from struct import pack
            # Check bit size to make sense
            if bit_size <= 0:
                raise ValueError('Bit size must a number greater than 0')
            if bit_size % 8 == 0:
                raise ValueError('Bit size is multiple of 8 so bytes must be used instead of bits')
            # Start writting the output file
            with open(filepath, 'wb') as file:
                bit_count = 0
                current_byte = ''
                # Iterate over data list values
                for value in data:
                    # Parse the value to binary and make sure the binary is as long as the bit size
                    bits = format(value, 'b').zfill(bit_size)
                    if len(bits) != bit_size:
                        raise ValueError(f'Value {value} cannot be stored in {bit_size} bits')
                    # Add bits one by one to the current byte to be written
                    for bit in bits:
                        current_byte += bit
                        bit_count += 1
                        # If the current byte is full then write it to the output file
                        if bit_count == 8:
                            #print(current_byte + ' -> ' + str(int(current_byte, 2)))
                            file.write(pack('<B', int(current_byte, 2)))
                            current_byte = ''
                            bit_count = 0
                # If last byte is truncated then fill it with 0s and write it
                if bit_count != 0:
                    last_byte = current_byte.ljust(8, '0')
                    file.write(pack('<B', int(last_byte, 2)))
        
        # Set the name of the binary dat file
        bin_file_path = output_dat_file.replace('.dat', '.bin')
        # Call the function to store the bits in the binary file
        store_bits_dat(data, 2, bin_file_path)
        # Set the name of the meta file and create the meta data corresponding to the binary dat file
        meta_file_path = bin_file_path + '.meta.json'
        meta_data = {
            'x': {
                'name': 'bases',
                'length': 20
            },
            'y': {
                'name': 'frames',
                'length': 500000
            },
            'bitsize': 2,
        }
        save_json(meta_data, meta_file_path)
    
    # Compress the dat file to .bin file
    compress_dat_to_bin(output_dat_file)
    print("Dat file compressed to bin file, removing original dat file...")
    # Remove the original dat file
    os.remove(output_dat_file)


# Function to execute Curves+ and Canals software to generate the output needed
def terminal_execution(trajectory_input,topology_input,strand_indexes,sequence,folder_path):
    helical_parameters_folder = f"{folder_path}/helical_parameters"
    # If the folder already exists dont create it again
    if not os.path.exists(helical_parameters_folder):
        os.mkdir(helical_parameters_folder)
    # Purge residual files from previous runs
    possible_residual_filenames = glob.glob(helical_parameters_folder+'/*.cda') + glob.glob(helical_parameters_folder+'/*.lis') + glob.glob(helical_parameters_folder+'/*.ser')
    for filename in possible_residual_filenames:
        os.remove(filename)
        
    # Indicate all the instructions with the desired commands, keep in mind that there could be more commands to include and obtain other results and files
    instructions = [
        "Cur+ <<!",
        " &inp",
        f"  file={trajectory_input},ftop={topology_input},lis={helical_parameters_folder}/test,lib={standard_prefix}",
        " &end",
        "2 1 -1 0 0",
        f"{strand_indexes[0][0]}:{strand_indexes[0][1]}",
        f"{strand_indexes[1][0]}:{strand_indexes[1][1]}",
        "!"
    ]
    instructions = ["\n".join(instructions)]
    cmd = " ".join(instructions)
    print(' Running curves')
    # Execute the software using the instructions written previously
    process = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,executable=os.getenv('SHELL', '/bin/sh'))

    out, err = process.communicate()
    process.wait()
    out.decode("utf-8")

    # Enter inside helicalparameters folder to execute Canals software as it needs the Curves+ output as input
    os.chdir(helical_parameters_folder)
    
    # If output has not been generated then warn the user
    if not os.path.exists('test.cda'):
        raise SystemExit('Something went wrong with Curves+ software')
    # Change the file name as Canals just take the name of the .cda file to differentiate between .lis and .cda
    move("test.cda","cinput.cda")
    # We have set up level1 and level2 to 0 as if lev1=lev2=0, lev1 is set to 1 and lev2 is set to the length of the oligmer
    level1 = 0
    level2 = 0
    # We just want to compute series files, but if you want more output just use the link written before 
    # And include them in the instructions similarly to series input
    
    series = ".t."
    instructions2 = [
        "Canal <<! ",
        "&inp",
        "  lis=canal_output,",
        f" lev1={level1},lev2={level2},",
        f" series={series},",
        "&end",
        f"cinput {sequence}",
        "!"
    ]
    instructions2 = ["\n".join(instructions2)]
    cmd2 = " ".join(instructions2)
    #logs = subprocess.run(instructions,stderr=subprocess.PIPE).stderr.decode()
    print(' Running canals')
    process = subprocess.Popen(cmd2,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,executable=os.getenv('SHELL', '/bin/sh'))
    
    out, err = process.communicate()
    process.wait()
    out.decode("utf-8")

    # If output has not been generated then warn the user
    if not os.path.exists('canal_output_phaseC.ser'):
        raise SystemExit('Something went wrong with Canals software')


def send_files(sequence,frames_limit, folder_path):
    files_averages = []
    files_backbones = []
    files_allbackbones = []
    # Iterate over all files in the directory
    helical_parameters_folder = f"{folder_path}/helical_parameters"
    os.chdir(helical_parameters_folder)
    for file in os.listdir():
        if not file.endswith(".ser"): # Check that we are going to analyze the correct files, that have .ser extension
            continue
        # data = []
        # with open(file, 'r') as ser_file:
        #     # Skip the first line since it is only the headers
        #     n_lines = 0
        #     for line in ser_file:
        #         # Skip the first line since it is only the row number
        #         numbers = [ s for s in line.strip().split()[1:] ]
        #         data += numbers
        #         n_lines += 1
        #         # if len(data) > 1000:
        #         #     break
        # # Save the file with the same name but with .bin extension and mdf prefix
        # file_name = 'mdf.' + file
        # file_path = os.path.join(os.getcwd(), file_name.replace('.ser', '.bin'))
        # # Store data in a binary format
        # store_binary_data(
        #     data=data,
        #     byte_size=4, # 32 bits
        #     filepath=file_path
        # )
        # # We need the meta file to store the information about the binary file (with the same name and .meta.json extension)
        # name_meta_file = file_path + '.meta.json'
        # meta_data = {
        #     'x': {
        #         'name': 'bases',
        #         'length': len(numbers) # Number of columns.  It is the last value of the bucle for line in ser_file (maybe it causes problems)
        #     },
        #     'y': {
        #         'name': 'frames',
        #         'length': n_lines # Total number of lines in the file os the number of frames
        #     },
        #     'bitsize': 32, # 32 bits
        # }
        # save_json(meta_data, name_meta_file)


        #  THIS PART OF THE CODE IS TO ANALYZE THE FILES AND COMPUTE THE AVERAGES, STANDARD DEVIATIONS AND TIME SERIES
        # BUT IT SHOULD BE REMOVED
        words = file.split('.')
        helpar = words[0].split('_')[-1].lower()
        # Check to which block the file came from
        if helpar in ["shear","stagger","stretch","buckle","opening","propel","inclin","tip","xdisp","ydisp", 
                        "majw","majd","minw","mind","rise","shift","slide","roll","tilt","twist"]:
            
            files_averages.append(file)

        elif helpar in ["alphac","alphaw","betac","betaw","gammac","gammaw","phasec","phasew","epsilc","epsilw","zetac","zetaw"]:
            
            files_backbones.append(file)

        if helpar in ['alphac', 'alphaw', 'betac', 'betaw', 'gammac', 'gammaw', 'deltac', 'deltaw', 'epsilc', 'epsilw', 
                        'zetac', 'zetaw', 'chic', 'chiw', 'phasec', 'phasew']:
            
            files_allbackbones.append(file)

    

    information_dictionary = {'avg_res':{'backbone':{},'grooves':{},'axisbp':{},'intrabp':{},'interbp':{}, 'stiffness':{}}}

    information_dictionary0 = get_stiffness(sequence,files_averages,information_dictionary,frames_limit)
    # 'ts':{'backbone':{},'grooves':{},'axisbp':{},'intrabp':{},'interbp':{}}
    # Before there was another field called 'ts' or time series but it was removed because now it is included in the .ser and .bin files

    # Work with files that correspond to the different blocks of Helical Parameters and perform the different computations to each block 

    # Call function to distribute the different files to compute their averages
    information_dictionary1 = flow_files_average(files_averages,information_dictionary0,sequence)

    # Call function to distribute the files related to backbone torsions and perform different calculcations to each
    information_dictionary2 = flow_files_backbone(files_backbones,information_dictionary1) 

    # Call function to distribut all files and compute Time Series in all of them
    #info_dictionary = flow_files_timeseries(files_averages,files_allbackbones,information_dictionary2,frames_limit) 
    return information_dictionary2 # Return the dictionary to convert it to a json

# Function to compute the stiffness 
def get_stiffness(sequence,files,info_dict,frames_limit):
    # Create and object (as NASSA does) to store the sequence and the corresponding .ser info so then we can calculate the stiffness
    extracted = {}
    # First parse the sequence to obtain the NucleicAcid object as NASSA wf does
    seq = []
    seq.append(load_sequence2(sequence, unit_len=6))
    extracted['sequences'] = seq
    # Set the hexamer coordinates
    # AGUS: si son pentámeros habría que añadir las coordenadas:  shear stagger stretch chic chiw buckle opening propel
    cordinate_list = ['shift','slide','rise','roll','tilt','twist']
    for coord in cordinate_list:
        crd_data = []
        for ser_file in files:
            if ser_file.split('_')[-1].split('.')[0].lower() == coord:
                ser = load_serfile(
                    ser_file)
                crd_data.append(ser)
        extracted[coord.lower()] = crd_data
        print(f"loaded {len(crd_data)} files for coordinate <{coord}>")
    results = transform(extracted)
    # Add the results stiffness to the dictionary
    info_dict['avg_res']['stiffness'] = results  
    return info_dict

# AGUS: estas tres funciones siguientes son una copia de NASSA pero modificada para este caso ya que NASSA 
# AGUS: funciona para todas las secuencias o MD al mismo tiempo y el WF para cada proyecto/MD
# AGUS: además, ahora se define por defecto el nombre de la unidad como hexámero
def transform(data, unit_name='hexamer'):
    sequences = data.pop('sequences')
    results = {"stiffness": [], "covariances": {}, "constants": {}}
    for traj, seq in enumerate(sequences):
        traj_series = {coord.lower(): data[coord][traj]
                        for coord in data.keys()}
        traj_results = get_stiffness_seq(
            seq,
            traj_series)
        results["stiffness"].append(traj_results["stiffness"])
        covariances_serializable = {
            key: df.to_dict(orient="split")  # 'index', 'columns', 'data'
            for key, df in traj_results["covariances"].items()
        }
        results["covariances"].update(covariances_serializable)
        constants_serializable = {
            key: df.to_dict(orient="split")  # 'index', 'columns', 'data'
            for key, df in traj_results["constants"].items()
        }
        results["constants"].update(constants_serializable)

    # AGUS: aquí se hace un cambio en la estructura del diccionario de stiffness para que sea más fácil de manejar en front end
    stiffness_by_coord = {}
    for item in results["stiffness"]:
        for key, value in item.items():
            if key == "hexamer":
                continue  
            if key not in stiffness_by_coord:
                stiffness_by_coord[key] = []
            stiffness_by_coord[key].append(value)
    results["stiffness"] = stiffness_by_coord
    for key, value in results["stiffness"].items():
        # Si es una lista con Series adentro, concatenamos y convertimos a list
        if isinstance(value, list) and all(hasattr(v, "tolist") for v in value):
            combined = pd.concat(value)
            results["stiffness"][key] = combined.tolist()
        # Si es una sola Series
        elif hasattr(value, "tolist"):
            results["stiffness"][key] = value.tolist()
    return results
def get_stiffness_seq(
        sequence,
        series_dict):
    # get stiffness table for a given trajectory
    coordinates = list(series_dict.keys())
    results = {"stiffness": {}, "covariances": {}, "constants": {}}
    diagonals = {}
    start = sequence.flanksize # + 2
    end = sequence.size - (sequence.baselen + sequence.flanksize - 1) # +2
    print(f"start: {start}, end: {end}")
    for i in range(start, end):
        hexamer = sequence.get_subunit(i)
        #ic_hexamer = sequence.inverse_complement(hexamer)
        cols_dict = {coord: series_dict[coord][i+1]
                        for coord in series_dict.keys()}
        #print(f"cols_dict: {cols_dict}")
        stiffness_diag, cte, cov_df = get_subunit_stiffness(
            cols_dict,
            coordinates)
        diagonals[hexamer] = np.append(
            stiffness_diag,
            [np.product(stiffness_diag), np.sum(stiffness_diag)])
        # diagonals[ic_hexamer] = np.append(
        #     stiffness_diag,
        #     [np.product(stiffness_diag), np.sum(stiffness_diag)])
        results["covariances"][hexamer] = cov_df
        #results["covariances"][ic_hexamer] = cov_df
        results["constants"][hexamer] = cte
        #results["constants"][ic_hexamer] = cte
    # build stiffness table
    columns = [sequence.unit_name] + coordinates + ["product", "sum"]
    results["stiffness"] = pd.DataFrame.from_dict(
        diagonals,
        orient="index").reset_index()
    results["stiffness"].columns = columns
    return results
def get_subunit_stiffness(
        cols_dict,
        coordinates,
        scaling=[1, 1, 1, 10.6, 10.6, 10.6],
        KT=0.592186827):
    import numpy.ma as ma
    def circ_avg(xarr, degrees=True):
        n = len(xarr)
        if degrees:
            # convert to radians
            xarr = xarr * np.pi / 180
        av = np.arctan2(
            (np.sum(np.sin(xarr)))/n,
            (np.sum(np.cos(xarr)))/n) * 180 / np.pi
        return av
    # AGUS: insisto, está adaptado solamente para hexámeros
    # if (self.unit_len % 2) == 0:
    SH_av = cols_dict["shift"].mean()
    SL_av = cols_dict["slide"].mean()
    RS_av = cols_dict["rise"].mean()
    TL_av = circ_avg(cols_dict["tilt"])
    RL_av = circ_avg(cols_dict["roll"])
    TW_av = circ_avg(cols_dict["twist"])
    # elif (self.unit_len % 2) == 1:
    #     SH_av = cols_dict["shear"].mean()
    #     SL_av = cols_dict["stretch"].mean()
    #     RS_av = cols_dict["stagger"].mean()
    #     CW_av = cols_dict["chiw"].mean()
    #     CC_av = cols_dict["chic"].mean()
    #     TL_av = self.circ_avg(cols_dict["buckle"])
    #     RL_av = self.circ_avg(cols_dict["propel"])
    #     TW_av = self.circ_avg(cols_dict["opening"])
    cols_arr = [cols_dict[coord] for coord in coordinates]
    cols_arr = np.array(cols_arr).T

    cv = ma.cov(ma.masked_invalid(cols_arr), rowvar=False)
    cv.filled(np.nan)

    cov_df = pd.DataFrame(cv, columns=coordinates, index=coordinates)
    stiff = np.linalg.inv(cv) * KT
    #print(stiff)
    # Added two new variables: ChiC and ChiW -> 8 (for PENTAMERS) => in NASSA (here only for hexamers)
    if (6 % 2) == 0:
        last_row = [SH_av, SL_av, RS_av, TL_av, RL_av, TW_av] #, CW_av, CC_av]
        stiff = np.append(stiff, last_row).reshape(7, 6)
    # elif (self.unit_len % 2) == 1:
    #     last_row = [SH_av, SL_av, RS_av, TL_av, RL_av, TW_av, CW_av, CC_av]
    #     print(last_row)
    #     stiff = np.append(stiff, last_row).reshape(9, 8)
    #     scaling=[1, 1, 1, 10.6, 10.6, 10.6, 1, 1]
    
    stiff = stiff.round(6)
    stiff_diag = np.diagonal(stiff) * np.array(scaling)

    cte = pd.DataFrame(stiff)
    cte.columns = coordinates
    cte.index = coordinates + ["avg"]
    return stiff_diag, cte, cov_df

# Function to check whether the file is from one block or another regarding Helical parameters, the value of baselen could change
def checking(helpar_name):
    # get base length and unit from helical parameter name
    if helpar_name.lower() in hp_basepairs:
        baselen = 1
    elif helpar_name.lower() in hp_singlebases:
        baselen = 0
    return baselen

# Function to check whether the file is from one block or another regarding Helical parameters, the value of baselen could change
def checking2(helpar_name):
    # get base length and unit from helical parameter name
    if helpar_name.lower() in hp_singlebases:
        baselen = 0
    else:
        baselen = 1
    if helpar_name.lower() in hp_angular:
        hp_unit = "Degrees"
    else:
        hp_unit = "Angstroms"
    return baselen,hp_unit

# Function to compute the inverse complement DNA or RNA sequence it is comented becuase it is not used but in case it must be used it is here
'''
# Compute inverse complement DNA or RNA sequence 
def reverse_sequence(sequence,DNA=True):
    if DNA: # If DNA flag is tru we want to compute the inverse using T instead of U
        A_base = "T"
    else: # Now it is RNA so we want to convert T to U
        A_base = "U"
    inverse = {"A":A_base,"G":"C","C":"G",A_base:"A"} # Dictionary to convert easily the sequence 
    inv_seq = ""
    for i in sequence[::-1]: # Traverse the sequence from the end to the beginning
        inv_seq += inverse[i] # Obtain the complementary base 
    return inv_seq # Return the inverse complement
'''
    
# Function to read a given .ser file and transform it into pandas dataframe
def read_series(input_serfile, usecols=None):
    """Read .ser file"""
    extra_kwargs = dict(
        header=None,
        sep='\\s+',
        index_col=0)
    ser_data = pd.read_csv(input_serfile, **extra_kwargs)
    if usecols is not None:
        if 0 in usecols:
            usecols.pop(usecols.index(0))
        ser_data = ser_data[[i+1 for i in usecols]]
    return ser_data


# AVERAGES FUNCTION TO COMPUTE THE AVERAGE OF THE FOLLOWING BLOCKS REGARDING THE HELICAL PARAMETERS
# GROOVES, INTRA BASE PAIR, INTER BASE PAIR, AXIS BASE PAIR

# FUNCTION TO CONTROL FLOW OF FILES THAT ENTERS THE FUNCTION TO BE ANALYZED IN ORDER TO COMPUTE AVERAGE (MEAN AND STANDARD DEVIATION)
def flow_files_average(files,info_dict,sequence):
    for file in files: # Iterate over all files that must be analyzed
        word = file.split('.')
        helpword = word[0].split('_')[-1].lower()
        # Obtain baselen value depending on the file
        baselen = checking(helpword) 
        # Convert the files into Pandas dataframe to make the computations and data manipulation easily
        dataframe = read_series(file) 
        df1,df2 = average_std(dataframe,baselen,sequence)
        df1i,df2i = average_std_intra(dataframe,baselen,sequence)
        if helpword in ["roll","tilt","twist","rise","shift","slide"]: 
            # INTER BASEPAIR BLOCK
            info_dict['avg_res']['interbp'][helpword] = {'avg':{},'std':{}}
            #  we want to store all the information regarding the averages 
            info_dict['avg_res']['interbp'][helpword]['avg'] = df1.T.values.tolist() 

            #  we want to store all the information regarding the standard deviations
            info_dict['avg_res']['interbp'][helpword]['std'] = df2.T.values.tolist() 

        elif helpword in ["shear","stagger","stretch","buckle","opening","propel"]: 
            # INTRA BASEPAIR BLOCK
            info_dict['avg_res']['intrabp'][helpword] = {'avg':{},'std':{}}
            #  we want to store all the information regarding the averages 
            info_dict['avg_res']['intrabp'][helpword]['avg'] = df1i.T.values.tolist() 

            #  we want to store all the information regarding the standard deviations
            info_dict['avg_res']['intrabp'][helpword]['std'] = df2i.T.values.tolist() 

        elif helpword in ["xdisp","ydisp","inclin","tip"]:
            # AXIS BASEPAIR BLOCK 
            info_dict['avg_res']['axisbp'][helpword] = {'avg':{},'std':{}} 

            #  we want to store all the information regarding the averages
            info_dict['avg_res']['axisbp'][helpword]['avg'] = df1i.T.values.tolist()  

            #  we want to store all the information regarding the standard deviations
            info_dict['avg_res']['axisbp'][helpword]['std'] = df2i.T.values.tolist() 

        else:
            # GROOVES BLOCK
            info_dict['avg_res']['grooves'][helpword] = {'avg':{},'std':{}} 

            #  we want to store all the information regarding the averages 
            info_dict['avg_res']['grooves'][helpword]['avg'] = df1.T.values.tolist() 

            #  we want to store all the information regarding the standard deviations
            info_dict['avg_res']['grooves'][helpword]['std'] = df2.T.values.tolist() 
    
    return info_dict # Return dictionary with all information computed (averages and standard deviation)

# Compute average and standard deviation
def average_std(dataf,baselen,sequence):
    # For hexamers, baselen has to be 2
    baselen = 2 
    # discard first and last base(pairs) from sequence
    dataf = dataf[dataf.columns[2:17]]
    # sequence = sequence[1:]
    # print("sequence: ",sequence)
    xlabels = [
        f"{sequence[i:i+1+baselen]}"
        for i in range(len(dataf.columns))] # - baselen
    means = dataf.mean(axis=0).iloc[:len(xlabels)]
    stds = dataf.std(axis=0).iloc[:len(xlabels)]
    return means,stds 

def average_std_intra(dataf,baselen,sequence):
    # For hexamers, baselen has to be 2
    baselen = 2 
    # discard first and last base(pairs) from sequence
    dataf = dataf[dataf.columns[2:18]]
    # sequence = sequence[1:]
    # print("sequence: ",sequence)
    xlabels = [
        f"{sequence[i:i+1+baselen]}"
        for i in range(len(dataf.columns))] # - baselen
    means = dataf.mean(axis=0).iloc[:len(xlabels)]
    stds = dataf.std(axis=0).iloc[:len(xlabels)]
    return means,stds 

# OTHER FUNCTIONS THAT MUST PERFORM OTHER CALCULATIONS REGARDING TO THE LAST BLOCK OF HELICAL PARAMETERS (BACKBONE TORSIONS)
# THERE ARE 3 DIFFERENT FUNCTIONS EACH ONE COMPUTING DIFFERENT ASPECTS OF THE FLEXIBILITY IN THE BACKBONE 
# SUGAR PUCKERING, CANONICAL ALPHA/GAMMA AND BI/BII POPULATIONS

# FUNCTION TO CONTROL 
def flow_files_backbone(files,info_dict):
    # Store files related to Puckering computations
    puckering_files = [] 
    # Store files related to Canonical Alpha Gamma computations
    canonicalg_files = [] 
    # Store files related to BI and BII populations
    bipopulations_files = [] 
    # Iterate over all files that must be analyzed
    for file in files: 
        word = file.split('.')
        helpword = word[0].split('_')[-1].lower()
        # Sugar Puckering
        if helpword in ['phasec','phasew']: 
            puckering_files.append(file)
        # Canonical Alpha/Gamma
        elif helpword in ['alphac','alphaw','gammac','gammaw']: 
            canonicalg_files.append(file)
        # BI/BII Population
        elif helpword in ['epsilc','epsilw','zetac','zetaw']: 
            bipopulations_files.append(file)
    info_dict['avg_res']['backbone'] = {'puckering':{'north':{},'east':{},'west':{},'south':{}},'bi':{},'bii':{},'canonical_alphagamma':{}}
    # Call function to perform computations regarding Sugar Puckering
    npop, epop, wpop, spop = sugar_puckering(puckering_files) 
    info_dict['avg_res']['backbone']['puckering']['north'] = npop.T.values.tolist()
    info_dict['avg_res']['backbone']['puckering']['east'] = epop.T.values.tolist()
    info_dict['avg_res']['backbone']['puckering']['west'] =  wpop.T.values.tolist()
    info_dict['avg_res']['backbone']['puckering']['south'] = spop.T.values.tolist()
    # Call function to perform computations regarding Canonical Alpha Gamma
    canon = canonical_alphagamma(canonicalg_files) 
    info_dict['avg_res']['backbone']['canonical_alphagamma'] = canon.T.values.tolist()
    # Call function to perform computations regarding BI and BII populations
    b1,b2 = bi_populations(bipopulations_files) 
    info_dict['avg_res']['backbone']['bi'] = b1.T.values.tolist()
    info_dict['avg_res']['backbone']['bii'] = b2.T.values.tolist()
    return info_dict

######## SUGAR PUCKERING 

# READ FILES AND CORRECT ANGLES 
def sugar_puckering(pcukfiles):
    for file in pcukfiles: # Iterate over all files 
        if file[-5:] == "C.ser":
            phasec = file
        else:
            phasew = file
    
    # Convert the files into Pandas dataframe to make the computations and data manipulation easily
    phasec = read_series(phasec) 
    # Convert the files into Pandas dataframe to make the computations and data manipulation easily
    phasew = read_series(phasew) 
    # Fix angles that are negative
    phasec = fix_angles(phasec)  
    # Fix angles that are negative
    phasew = fix_angles(phasew) 

    return cpuck(phasec,phasew)

# Correct angles values as we dont want negative angles measurements
def fix_angles(dataset):
    values = np.where(dataset < 0, dataset + 360, dataset)
    dataset = pd.DataFrame(values)
    return dataset

# Compute sugar puckering, using 2 files phaseC and phaseW
def cpuck(phaseC,phaseW):
    separator_df = pd.DataFrame({"-": np.nan}, index=range(1, len(phaseC)))
    phase = pd.concat([
        phaseW,
        separator_df,
        phaseC[phaseC.columns[::-1]]],
        axis=1)

    Npop = np.logical_or(phase > 315, phase < 45).mean() * 100
    Epop = np.logical_and(phase > 45, phase < 135).mean() * 100
    Wpop = np.logical_and(phase > 225, phase < 315).mean() * 100
    Spop = np.logical_and(phase > 135, phase < 225).mean() * 100

    return Npop, Epop, Wpop, Spop

######## CANONICAL ALPHA/GAMMA

# READ FILES AND CORRECT ANGLES REGARDING CANICAL ALPHA GAMMA 
def canonical_alphagamma(canonfiles):
    for file in canonfiles:
        if file[-7:] == "haC.ser":
            alphac = file
        elif file[-7:] == "haW.ser":
            alphaw = file
        elif file[-7:] == "maC.ser":
            gammac = file
        elif file[-7:] == "maW.ser":
            gammaw = file

    # Convert the files into Pandas dataframe to make the computations and data manipulation easily
    alphac = read_series(alphac) 
    # Convert the files into Pandas dataframe to make the computations and data manipulation easily
    alphaw = read_series(alphaw) 
    # Convert the files into Pandas dataframe to make the computations and data manipulation easily
    gammac = read_series(gammac) 
    # Convert the files into Pandas dataframe to make the computations and data manipulation easily
    gammaw = read_series(gammaw) 
    # Fix angles that are negative and over 360 degrees
    alphac = fix_angles2(alphac) 
    # Fix angles that are negative and over 360 degrees
    alphaw = fix_angles2(alphaw)
    # Fix angles that are negative and over 360 degrees
    gammac = fix_angles2(gammac)
    # Fix angles that are negative and over 360 degrees
    gammaw = fix_angles2(gammaw) 

    return check_alpgamm(alphac,gammac,alphaw,gammaw)

# Correct angles values as we dont want negative angles measurements and more than 360 degrees
def fix_angles2(dataset):
    values = np.where(dataset < 0, dataset + 360, dataset)
    values = np.where(values > 360, values - 360, values)
    dataset = pd.DataFrame(values)
    return dataset

# Compute canonical populations using 4 different files alphaC, gammaC, alphaW and gammaW
def check_alpgamm(alphaC,gammaC,alphaW,gammaW):
    separator_df = pd.DataFrame({"-": np.nan}, index=range(len(gammaW)))
    gamma = pd.concat([
        gammaW,
        separator_df,
        gammaC[gammaC.columns[::-1]]],
        axis=1)
    alpha = pd.concat([
        alphaW,
        separator_df,
        alphaC[alphaC.columns[::-1]]],
        axis=1)
    alpha_filter = np.logical_and(alpha > 240, alpha < 360)
    gamma_filter = np.logical_and(gamma > 0, gamma < 120)
    canonical_alpha_gamma = np.logical_and(
        alpha_filter, gamma_filter).mean() * 100

    return canonical_alpha_gamma

# BI / BII POPULATIONS

def bi_populations(bfiles):
    for file in bfiles:
        if file[-7:] == "ilC.ser":
            epsilc = file
        elif file[-7:] == "ilW.ser":
            epsilw = file
        elif file[-7:] == "taC.ser":
            zetac = file
        elif file[-7:] == "taW.ser":
            zetaw = file

    # Convert the files into Pandas dataframe to make the computations and data manipulation easily
    epsilc = read_series(epsilc) 
    # Convert the files into Pandas dataframe to make the computations and data manipulation easily
    epsilw = read_series(epsilw) 
    # Convert the files into Pandas dataframe to make the computations and data manipulation easily
    zetac = read_series(zetac) 
    # Convert the files into Pandas dataframe to make the computations and data manipulation easily
    zetaw = read_series(zetaw) 
    # Fix angles that are negative and over 360 degrees
    epsilc = fix_angles2(epsilc) 
    # Fix angles that are negative and over 360 degrees
    epsilw = fix_angles2(epsilw)
    # Fix angles that are negative and over 360 degrees
    zetac = fix_angles2(zetac)
    # Fix angles that are negative and over 360 degrees
    zetaw = fix_angles2(zetaw) 
    # Compute differences between Epsilon and Zeta values
    diff_epsil_zeta = angles_diff_ze(epsilc,zetac,epsilw,zetaw) 
    BI,BII = bi_pop(diff_epsil_zeta)
    return BI,BII

# Compute difference between epsil and zeta
def angles_diff_ze( epsilC, zetaC, epsilW, zetaW):
    # concatenate zeta and epsil arrays
    separator_df = pd.DataFrame({"-": np.nan}, index=range(len(zetaW)))
    zeta = pd.concat([
        zetaW,
        separator_df,
        zetaC[zetaC.columns[::-1]]],
        axis=1)
    epsil = pd.concat([
        epsilW,
        separator_df,
        epsilC[epsilC.columns[::-1]]],
        axis=1)

    # difference between epsilon and zeta coordinates
    diff_epsil_zeta = epsil - zeta
    return diff_epsil_zeta

# Compute BI and BII populations
def bi_pop(diff_epsil_zeta):
    BI = (diff_epsil_zeta < 0).sum(axis=0) * 100 / len(diff_epsil_zeta)
    BII = 100 - BI
    return BI, BII

######## TIME SERIES

# Function to distribute all .ser files in order to compute later Time Series
def flow_files_timeseries(files_average,files_backbone,info_dict,frames_limit):
    for file in files_average: # Iterate over all files that must be analyzed
        word = file.split('.')
        helpword = word[0].split('_')[-1].lower()
        df = time_series(file,helpword,frames_limit)

        # Store all computations from files that are related to the block Inter Basepair 
        if helpword in ["roll","tilt","twist","rise","shift","slide"]: 
            info_dict['ts']['interbp'][helpword] = {}
            # Store all information 
            info_dict['ts']['interbp'][helpword] = df.T.values.tolist()  

        # Store all computations from files that are related to the block Intra Basepair
        elif helpword in ["shear","stagger","stretch","buckle","opening","propel"]: 
            info_dict['ts']['intrabp'][helpword] = {}
            # Store all information 
            info_dict['ts']['intrabp'][helpword] = df.T.values.tolist() 

        # Store all computations from files that are related to the block Axis Basepair
        elif helpword in ["xdisp","ydisp","inclin","tip"]: 
            info_dict['ts']['axisbp'][helpword] = {}
            # Store all information
            info_dict['ts']['axisbp'][helpword] = df.T.values.tolist()  


        else:
            # Store all computations from files that are related to the block Grooves
            info_dict['ts']['grooves'][helpword] = {} 
            # Store all information 
            info_dict['ts']['grooves'][helpword] = df.T.values.tolist() 
            
    for file2 in files_backbone:
        word = file2.split('.')
        helpword = word[0].split('_')[-1].lower()
        df1 = time_series(file2,helpword,frames_limit)
        # Store all computations from files that are related to the block Backbone Torsions
        info_dict['ts']['backbone'][helpword] = {} 
        # Store all information
        info_dict['ts']['backbone'][helpword] = df1.T.values.tolist() 


    return info_dict

# Function to compute Time Series for each files that flow_files_timeseries() passes
def time_series(file,word,reduced_trajectory_frames_limit):
    sequence = "GCTTTTAGCGGTTGAACGGC"
    baselen,hp_unit = checking2(word)
    ser_data = read_series(file)
    ser_data = ser_data[ser_data.columns[1:-1]]
    sequence = sequence[1:]

    '''
    subunits = [
        f"{i+1}_{sequence[i:i+1+baselen]}"
        for i in range(len(ser_data.columns))]
    ser_data.columns = subunits
    '''
    snapshots = len(ser_data.index)

    step = math.ceil(snapshots)# / reduced_trajectory_frames_limit)

    ser_data = ser_data[0:len(ser_data.index):step]
    return ser_data
