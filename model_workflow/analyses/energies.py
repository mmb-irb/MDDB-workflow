# Energies
#
# The energies analysis is carried by ACPYPE and a locally developed tool: CMIP.
# ACPYPE is a tool based in Python to use Antechamber to generate topologies for chemical compounds and to interface with others python applications.
# CMIP stands for Classical Molecular Interaction Potential and it is usefull to predict electrostatic and Van der Waals potentials.
#
# ACPYPE:
# SOUSA DA SILVA, A. W. & VRANKEN, W. F. ACPYPE - AnteChamber PYthon Parser interfacE. BMC Research Notes 5 (2012), 367 doi: 10.1186/1756-0500-5-367 http://www.biomedcentral.com/1756-0500/5/367
# BATISTA, P. R.; WILTER, A.; DURHAM, E. H. A. B. & PASCUTTI, P. G. Molecular Dynamics Simulations Applied to the Study of Subtypes of HIV-1 Protease. Cell Biochemistry and Biophysics 44 (2006), 395-404. doi: 10.1385/CBB:44:3:395
#
# Antechamber:
# WANG, J., WANG, W., KOLLMAN, P. A., and CASE, D. A. Automatic atom type and bond type perception in molecular mechanical calculations. Journal of Molecular Graphics and Modelling 25, 2 (2006), 247–260. doi: 10.1016/j.jmgm.2005.12.005
# WANG, J., WOLF, R. M., CALDWELL, J. W., KOLLMAN, P. A., and CASE, D. A. Development and testing of a General Amber Force Field. Journal of Computational Chemistry 25, 9 (2004), 1157–1174. doi: 10.1002/jcc.20035
#
# CMIP:
# Gelpí, J.L., Kalko, S.G., Barril, X., Cirera, J., de la Cruz, X., Luque, F.J. and Orozco, M. (2001), Classical molecular interaction potentials: Improved setup procedure in molecular dynamics simulations of proteins. Proteins, 45: 428-437. doi:10.1002/prot.1159

# Imports libraries
import os
from pathlib import Path
import re
import numpy
import math
from subprocess import run, PIPE, Popen
import json
from typing import Optional, List, Tuple

from model_workflow.tools.get_pdb_frames import get_pdb_frames

from mdtoolbelt.structures import Structure
from mdtoolbelt.selections import Selection

CURSOR_UP_ONE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'
ERASE_2_PREVIOUS_LINES = CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE + ERASE_LINE + CURSOR_UP_ONE

# Set the path to auxiliar files required for this analysis
resources = str(Path(__file__).parent.parent / "utils" / "resources")
#preppdb_source = resources + '/preppdb.pl'
preppdb_source = resources + '/old_preppdb.pl' # DANI: Este no renombra los terminales
cmip_inputs_checkonly_source = resources + '/check.in'
cmip_inputs_source = resources + '/input.in'
vdw_source = resources + '/vdwprm'

# Set a folder to be created in order to store residual output files from this analysis
current_directory = str(Path().absolute())
energies_folder = current_directory + '/energies'

# Perform the electrostatic and vdw energies analysis for each pair of interaction agents
def energies (
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_analysis_filename : str,
    structure : 'Structure',
    interactions : list,
    charges : list,
    frames_limit : int):

    if not interactions or len(interactions) == 0:
        print('No interactions were specified')
        return

    if not charges or len(charges) == 0:
        print('No charges were passed')
        return

    # Check the number of atoms on each interacting agent
    # If there is any agent with more than 80000 atoms CMIP will fail so we must skip this specific energies analysis by now
    # DANI: Este valor límite se puede cambiar en CMIP, pero hay que recompilar y eso no es banal en un environment de conda
    cmip_atom_limit = 80000
    for interaction in interactions:
        exceeds = False
        for agent in ['1', '2']:
            residues = interaction['residues_' + agent]
            atom_count = sum([ len(residue.atom_indices) for residue in residues ])
            if atom_count >= cmip_atom_limit:
                exceeds = True
                break
        if exceeds:
            print('WARNING: ' + interaction['name'] + ' is exceeding the CMIP atom count limit of ' + str(cmip_atom_limit) + ' and it will be skipped for this analysis')
            interaction['exceeds'] = True

    # This anlaysis produces many residual output files
    # Create a new folder to store all ouput files so they do not overcrowd the main directory
    if not os.path.exists(energies_folder):
        os.mkdir(energies_folder)

    # Adapt the structure for cmip
    energies_structure = structure.copy()

    # Rename residues according to if they are terminals
    name_terminal_residues(energies_structure)

    # Set each atom element in CMIP format
    set_cmip_elements(energies_structure)

    # Save the structure back to a pdb
    energies_structure_filename = energies_folder + '/energies.pdb'
    energies_structure.generate_pdb_file(energies_structure_filename)

    # Transform an agent structure to a cmip input pdb, which includes charges
    # Charges have been previously taken from the charges topology and in the structure as atom additional atributes
    def selection2cmip (name : str, residues : List['Residue'], structure : 'Structure', strong_bonds : Optional[list]) -> str:
        # Parse the residue indices selection
        # Convert residue indices to atom indices
        atom_indices = sum([ residue.atom_indices for residue in residues ],[])
        parsed_selection = Selection(atom_indices)
        # Raise an error if the selection is empty
        if not parsed_selection:
            raise ValueError('ERROR: Agent "' + name + '" with selection "' + selection + '" has no atoms')
        # Filter the selected atom in the structure
        selected_structure = structure.filter(parsed_selection)
        # Set atom charges (non standard attribute), which are used further to write the adapted cmip pdb
        # Set also the elements to macth the original structure, since the frame generator messes the elements
        for a, atom in enumerate(selected_structure.atoms):
            atom_index = parsed_selection.atom_indices[a]
            charge = charges[atom_index]
            setattr(atom, 'charge', charge)
            cmip_element = energies_structure.atoms[atom_index].element
            atom.element = cmip_element
        # Get indices of atoms in current agent with a strong bond with an atom in the other agent
        # They will be further marked for CMIP to ignore them
        strong_bond_indexes = []
        if strong_bonds:
            for bond in strong_bonds:
                strong_bond_indexes += bond
        # Write a special pdb which contains charges as CMIP expects to find them
        output_filename = energies_folder + '/' + name + '.cmip.pdb'
        with open(output_filename, "w") as file:
            # Write a line for each atom
            for a, atom in enumerate(selected_structure.atoms):
                
                index = str(a+1).rjust(5)
                atom_name = atom.name
                name =  ' ' + atom_name.ljust(3) if len(atom_name) < 4 else atom_name
                residue = atom.residue
                residue_name = residue.name.ljust(3)
                chain = atom.chain
                chain_name = chain.name.rjust(1)
                residue_number = str(residue.number).rjust(4)
                icode = residue.icode.rjust(1)
                coords = atom.coords
                x_coord, y_coord, z_coord = [ "{:.3f}".format(coord).rjust(8) for coord in coords ]
                charge = "{:.4f}".format(atom.charge) # Charge was manually added before, it is not a standard attribute
                # In case this atom is making an strong bond between both interacting agents we add an 'X' before the element
                # This way CMIP will ignore the atom. Otherwise it would return high non-sense Van der Waals values
                real_index = atom.index
                cmip_ignore_flag = 'X' if real_index in strong_bond_indexes else ''
                element = cmip_ignore_flag + atom.element
                
                atom_line = ('ATOM  ' + index + ' ' + name + ' ' + residue_name + ' '
                    + chain_name + residue_number + icode + '   ' + x_coord + y_coord + z_coord
                    + ' ' + str(charge).rjust(7) + '  ' + element + '\n')
                file.write(atom_line)

        return output_filename

    # Given a pdb structure, use CMIP to extract energies
    # Output energies are already added by residues
    def get_frame_energy (frame_structure : 'Structure') -> List[dict]:

        # WARNING: At this point structure should be corrected
        # WARNING: Repeated atoms will make the analysis fail

        # Repeat the whole process for each interaction
        data = []
        for interaction in interactions:

            # Check if the interaction has been marked as 'exceeds', in which case we skip it
            if interaction.get('exceeds', False):
                continue

            name = interaction['name']
            strong_bonds = interaction.get('strong_bonds', None)

            # Select the first agent and extract it in a CMIP friendly pdb format
            agent1_name = interaction['agent_1'].replace(' ', '_').replace('/', '_')
            agent1_residue_indices = interaction['residues_1']
            agent1_cmip = selection2cmip(agent1_name, agent1_residue_indices, frame_structure, strong_bonds)

            # Repeat the process with agent 2
            agent2_name = interaction['agent_2'].replace(' ', '_').replace('/', '_')
            agent2_residue_indices = interaction['residues_2']
            agent2_cmip = selection2cmip(agent2_name, agent2_residue_indices, frame_structure, strong_bonds)

            # Copy the source cmip inputs file in the local directory
            # Inputs will be modified to adapt the cmip grid to both agents together
            # Note than modified inputs file is not conserved along frames
            # Structures may change along trajectory thus requiring a different grid size
            cmip_inputs = energies_folder + '/cmip.in'
            run([ "cp", cmip_inputs_source, cmip_inputs ], stdout=PIPE)

            # Select the first agent interface and extract it in a CMIP friendly pdb format
            agent1_interface_name = agent1_name + '_interface'
            agent1_interface_residue_indices = interaction['interface_1']
            agent1_interface_cmip = selection2cmip(agent1_interface_name, agent1_interface_residue_indices, frame_structure, strong_bonds)

            # Repeat the process with agent 2
            agent2_interface_name = agent2_name + '_interface'
            agent2_interface_residue_indices = interaction['interface_2']
            agent2_interface_cmip = selection2cmip(agent2_interface_name, agent2_interface_residue_indices, frame_structure, strong_bonds)

            # Set the CMIP box dimensions and densities to fit both the host and the guest
            # Box origin and size are modified in the cmip inputs
            # Values returned are only used for display / debug purposes
            box_origin, box_size = adapt_cmip_grid(agent1_interface_cmip, agent2_interface_cmip, cmip_inputs)

            # Run the CMIP software to get the desired energies
            print(' Calculating energies for ' + agent1_name + ' as host and ' + agent2_name + ' as guest')
            agent1_output_energy = energies_folder + '/' + agent1_name + '.energy.pdb'
            agent1_residue_energies, agent1_atom_energies = get_cmip_energies(cmip_inputs, agent1_cmip, agent2_cmip, agent1_output_energy)
            print(' Calculating energies for ' + agent2_name + ' as host and ' + agent1_name + ' as guest')
            agent2_output_energy = energies_folder + '/' + agent2_name + '.energy.pdb'
            agent2_residue_energies, agent2_atom_energies = get_cmip_energies(cmip_inputs, agent2_cmip, agent1_cmip, agent2_output_energy)

            # DANI: Usa esto para escribir los resultados de las energías por átomo
            sample = {
                'agent1': { 'energies': agent1_atom_energies, 'pdb': agent1_cmip },
                'agent2': { 'energies': agent2_atom_energies, 'pdb': agent2_cmip },
                'box': { 'origin': box_origin, 'size': box_size } 
            }
            with open('energies-sample.json', 'w') as file:
                json.dump(sample, file)

            data.append({ 'agent1': agent1_residue_energies, 'agent2': agent2_residue_energies })

            # Delete current agent pdb files before going for the next interaction
            # run([
            #     "rm",
            #     agent1_cmip,
            #     agent2_cmip,
            #     agent1_output_energy,
            #     agent2_output_energy,
            # ], stdout=PIPE).stdout.decode()

            # Erase the 2 previous log lines
            print(ERASE_2_PREVIOUS_LINES)

        return data

    # Extract the energies for each frame in a reduced trajectory
    frames, step, count = get_pdb_frames(energies_structure_filename, input_trajectory_filename, frames_limit)
    interactions_data = [[] for interaction in interactions if not interaction.get('exceeds', False)]
    for current_frame_pdb in frames:
        
        # Run the main analysis over the current frame
        # Append the result data for each interaction
        current_frame_structure = Structure.from_pdb_file(current_frame_pdb)
        frame_energies_data = get_frame_energy(current_frame_structure)
        for i, data in enumerate(frame_energies_data):
            interactions_data[i].append(data)

    # Now calculated residue average values through all frames for each pair of interaction agents
    output_analysis = []
    for i, interaction in enumerate(interactions):

        # Check if the interaction as been marked as 'exceeds', in which case we skip it
        if interaction.get('exceeds', False):
            continue

        # Get the main data
        data = interactions_data[i]

        # Format data
        agent1_output = format_data([ frame['agent1'] for frame in data ])
        agent2_output = format_data([ frame['agent2'] for frame in data ])

        # Format the results data and append it to the output data
        output = {
            'name': interaction['name'],
            'agent1': agent1_output,
            'agent2': agent2_output,
        }
        output_analysis.append(output)

    # Finally, export the analysis in json format
    with open(output_analysis_filename, 'w') as file:
        json.dump({'data': output_analysis}, file)

    # Finally remove the reduced topology
    logs = run([
        "rm",
        energies_structure_filename,
    ], stdout=PIPE).stdout.decode()

# Given a topology (e.g. pdb, prmtop), extract the atom elements in a CMIP friendly format
# Hydrogens bonded to carbons remain as 'H'
# Hydrogens bonded to oxygen are renamed as 'HO'
# Hydrogens bonded to nitrogen or sulfur are renamed as 'HN'
# Some heavy atom elements may be also modified (e.g. 'CL' -> 'Cl')
def set_cmip_elements (structure : 'Structure'):
    # Iterate over each atom to fix their element according to CMIP standards
    for a, atom in enumerate(structure.atoms):
        element = atom.element
        # Adapt hydrogens element to CMIP requirements
        if element == 'H':
            # We must find the element of the heavy atom this hydrogen is bonded to
            atom_bonds = atom.get_bonds()
            # There should be always only 1 bond
            if len(atom_bonds) != 1:
                print('Atom ' + atom.name + ' (' + str(atom.index) + ') has ' + str(len(atom_bonds)) + ' bonds')
                raise ValueError('An hydrogen should always have one and only one bond')
            bonded_atom_index = atom_bonds[0]
            bonded_atom_element = structure.atoms[bonded_atom_index].element
            # Hydrogens bonded to carbons remain as 'H'
            if bonded_atom_element == 'C':
                pass
            # Hydrogens bonded to oxygen are renamed as 'HO'
            elif bonded_atom_element == 'O':
                element = 'HO'
            # Hydrogens bonded to nitrogen or sulfur are renamed as 'HN'
            elif bonded_atom_element == 'N' or bonded_atom_element == 'S':
                element = 'HN'
            else:
                raise SystemExit(
                    'ERROR: Hydrogen bonded to not supported heavy atom: ' + bonded_atom_element)
        atom.element = element

protein_residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
                    'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

rna_residues = ['RA', 'RU', 'RC', 'RG']

# Change residue names in a structure to meet the CMIP requirements
# Change terminal residue names by adding an 'N' or 'C'
def name_terminal_residues (structure : 'Structure'):

    for chain in structure.chains:

        residues = chain.residues

        # Check if the first residue is tagged as a terminal residue
        # If not, rename it
        first_residue = residues[0]
        first_residue_name = first_residue.name
        # In case it is a protein
        if first_residue_name in protein_residues:
            first_residue_name += 'N'
            first_residue.name = first_residue_name
        # In case it is RNA
        elif first_residue_name in rna_residues:
            first_residue_name += '5'
            first_residue.name = first_residue_name

        # Check if the last residue is tagged as 'C' terminal
        # If not, rename it
        last_residue = residues[-1]
        last_residue_name = last_residue.name
        # In case it is a protein
        if last_residue_name in protein_residues:
            last_residue_name += 'C'
            last_residue.name = last_residue_name
        # In case it is RNA
        elif last_residue_name in rna_residues:
            last_residue_name += '3'
            last_residue.name = last_residue_name

# Run CMIP in 'checkonly' mode (i.e. start and stop) for both agents
# Mine the grid generated by CMIP for both agents and calculate a new grid which would include both
# Finally, modify the provided 'cmip_inputs' with the new grid parameters
# Keep residues in the interface only to determine the size of the box since residues which are far will never contribute
def adapt_cmip_grid (agent1_cmip_pdb : str, agent2_cmip_pdb : str, cmip_inputs : str):

    # Set a name for the checkonly CMIP outputs
    # This name is not important, since the data we want is in the CMIP logs
    cmip_checkonly_output = energies_folder + '/checkonly.energy.pdb'

    # Move to the energies folder for a moment for the residual cmip files to get generated inside of it
    os.chdir(energies_folder)

    # Run CMIP in 'checkonly' mode and save the grid dimensions output
    # First do it for the agent 1
    cmip_logs_agent1 = run([
        "cmip",
        "-i",
        cmip_inputs_checkonly_source,
        "-pr",
        agent1_cmip_pdb,
        "-vdw",
        vdw_source,
        "-hs",
        agent2_cmip_pdb,
        "-byat",
        cmip_checkonly_output,
    ], stdout=PIPE, stderr=PIPE).stdout.decode()

    # Mine the grid dimensions from CMIP logs
    agent1_center, agent1_density, agent1_units = mine_cmip_output(
        cmip_logs_agent1.split("\n"))

    # Run CMIP in 'checkonly' mode and save the grid dimensions output
    # Now do it for the agent 2
    cmip_logs_agent2 = run([
        "cmip",
        "-i",
        cmip_inputs_checkonly_source,
        "-pr",
        agent2_cmip_pdb,
        "-vdw",
        vdw_source,
        "-hs",
        agent1_cmip_pdb,
        "-byat",
        cmip_checkonly_output,
    ], stdout=PIPE, stderr=PIPE).stdout.decode()

    # Mine the grid dimensions from CMIP logs
    agent2_center, agent2_density, agent2_units = mine_cmip_output(
        cmip_logs_agent2.split("\n"))

    # Get back to the main folder
    os.chdir("..")

    # Calculate grid dimensions for a new grid which contains both previous grids
    new_center, new_density = compute_new_grid(
        agent1_center,
        agent1_density,
        agent1_units,
        agent2_center,
        agent2_density,
        agent2_units)

    # Here we must calculate how many grid points the new box would have
    # In case the number of points exceeds a safe limit we must reduce it
    # Otherwise CMIP could return the following error:
    # ERROR when trying to allocate           -2092351636  bytes of memory

    # Set the limit number of grid points
    # It has been found experimentally: 65945880 Grid points -> OK, 68903412 Grid points -> ERROR
    grid_points_limit = 65000000

    # Set the default grid 'cells' size for all three dimensions
    grid_unit_size = 0.5

    # Calculate the amount of total points with current parameters
    current_grid_points = (new_density[0] + 1) * (new_density[1] + 1) * (new_density[2] + 1)

    # In case the current number of grid points exceeds the limit...
    # Reduce all dimension densities in proportion and expand the grid units size to compensate
    if current_grid_points > grid_points_limit:
        print('WARNING: Grid points limit is exceeded (' + str(current_grid_points) + ')')
        proportion = grid_points_limit / current_grid_points
        new_density[0] = math.ceil(new_density[0] * proportion)
        new_density[1] = math.ceil(new_density[1] * proportion)
        new_density[2] = math.ceil(new_density[2] * proportion)
        grid_unit_size = math.ceil(grid_unit_size / proportion * 1000) / 1000
        print('WARNING: Grid resolution has been reduced -> unit size = ' + str(grid_unit_size))

    # Set the new lines to be written in the local CMIP inputs file
    grid_inputs = [
        f" cenx={new_center[0]:.1f} \n",
        f" ceny={new_center[1]:.1f} \n",
        f" cenz={new_center[2]:.1f} \n",
        f" dimx={int(new_density[0])} \n",
        f" dimy={int(new_density[1])} \n",
        f" dimz={int(new_density[2])} \n",
        f" intx={grid_unit_size} \n",
        f" inty={grid_unit_size} \n",
        f" intz={grid_unit_size} \n",
    ]
    # Add previous lines to the local inputs file
    with open(cmip_inputs, "r+") as file:
        lines = file.readlines()
        file.seek(0)
        for line in lines:
            if line == '&end \n':
                for grid_input in grid_inputs:
                    file.write(grid_input)
            file.write(line)

    # Delete the 'ckeckonly' file
    run([
        "rm",
        cmip_checkonly_output,
    ], stdout=PIPE).stdout.decode()

    # Calculate the resulting box origin and size and return both values
    # These values are used for display / debug purposes only
    new_size = (new_density[0] * grid_unit_size, new_density[1] * grid_unit_size, new_density[2] * grid_unit_size)
    new_origin = (new_center[0] - new_size[0] / 2, new_center[1] - new_size[1] / 2, new_center[2] - new_size[2] / 2)
    return new_origin, new_size

def mine_cmip_output (logs):
    center, density, units = (), (), ()
    grid_density_exp = r"^\s*Grid density:\s+([-]*\d+)\s+([-]*\d+)\s+([-]*\d+)"
    grid_center_exp = r"^\s*Grid center:\s+([-]*\d+.\d+)\s+([-]*\d+.\d+)\s+([-]*\d+.\d+)"
    grid_units_exp = r"^\s*Grid units:\s+(\d+.\d+)\s+(\d+.\d+)\s+(\d+.\d+)"
    for line in logs:
        grid_center_groups = re.match(grid_center_exp, line)
        grid_density_groups = re.match(grid_density_exp, line)
        grid_units_groups = re.match(grid_units_exp, line)
        if grid_density_groups:
            density = tuple(float(grid_density_groups.group(i))
                            for i in (1, 2, 3))
        if grid_center_groups:
            center = tuple(float(grid_center_groups.group(i))
                           for i in (1, 2, 3))
        if grid_units_groups:
            units = tuple(float(grid_units_groups.group(i))
                          for i in (1, 2, 3))
    # If data mining fails there must be something wrong with the CMIP output
    if center == () or density == () or units == ():
        for line in logs:
            print(line)
        print('WARNING: CMIP output mining failed')
        raise SystemExit('ERROR: Something was wrong with CMIP')
    return center, density, units

# This function is used to create new grid parameters
# The new grid is expected to cover both input grids: agent 1 grid and agent 2 grid
def compute_new_grid (
        agent1_center,
        agent1_density,
        agent1_units,
        agent2_center,
        agent2_density,
        agent2_units,
        extra_density=1):
    new_center = []  # [(i + j) / 2 for i, j in zip(agent1_center, agent2_center)]
    new_density = []

    for k in range(3):
        min_prot = agent1_center[k] - agent1_density[k] * agent1_units[k]
        min_lig = agent2_center[k] - agent2_density[k] * agent2_units[k]
        min_new = min(min_prot, min_lig)
        max_prot = agent1_center[k] + agent1_density[k] * agent1_units[k]
        max_lig = agent2_center[k] + agent2_density[k] * agent2_units[k]
        max_new = max(max_prot, max_lig)

        dnew = int(abs(max_new - min_new) + extra_density)
        cnew = min_new + (dnew / 2)
        new_density.append(dnew)
        new_center.append(cnew)
    return new_center, new_density

# Run the CMIP software to get the desired energies
def get_cmip_energies (cmip_inputs, pr, hs, cmip_output):
    # Move to the energies folder for a moment for the residual cmip files to get generated inside of it
    os.chdir(energies_folder)
    # Run cmip
    cmip_logs = run([
        "cmip",
        "-i",
        cmip_inputs,
        "-pr",
        pr,
        "-vdw",
        vdw_source,
        "-hs",
        hs,
        "-byat",
        cmip_output,
    ], stdout=PIPE, stderr=PIPE).stdout.decode()
    # Get back to the main folder
    os.chdir("..")

    # Mine the electrostatic (es) and Van der Walls (vdw) energies for each atom
    # Group the results by residues adding their values
    atom_energies = []
    residue_energies = {}
    with open(cmip_output, 'r') as file:
        lines = list(file)
        # If this file is empty it means something went wrong with CMIP
        # We print its logs and exit
        if len(lines) == 0:
            if type(cmip_logs) == str:
                print(cmip_logs)
            else:
                for line in cmip_logs:
                    print(line)
            raise SystemExit('ERROR: Something went wrong with CMIP!')
        # Mine energies line by line (i.e. atom by atom)
        for line in lines:
            chain = line[21:22]
            residue_id = line[22:28]
            residue = chain + ':' + residue_id.strip() # strip removes white spaces
            vdw = float(line[42:53])
            es = float(line[57:68])
            both = float(line[72:83])
            # Values greater than 100 are represented as 0
            # This step is performed to filter 'infinity' values
            energies = (vdw, es, both)
            if both > 100:
                #print('WARNING: We have extremly high values in energies which are beeing discarded')
                energies = (0, 0, 0)
            # Add current atom energy values to the atom energies list
            # DANI: Esto ahora mismo solo sirve para poder samplear los resultados antes de que se junten por residuos
            atom_energies.append(energies)
            # Add current atom energy values to the accumulated residue energy value
            if residue in residue_energies:
                residue_energies[residue] = tuple([a+b for a, b in zip(energies, residue_energies[residue])])
            else:
                residue_energies[residue] = energies

    return residue_energies, atom_energies
    

# Format data grouping atom energies by residue
# Calculate means of the whole
def format_data (data : list) -> dict:

    # First, reorder data by residues and energies
    residues_number = len(data[0])
    residues_labels = [residue for residue in data[0]]

    residues_vdw_values = [[] for n in range(residues_number)]
    residues_es_values = [[] for n in range(residues_number)]
    residues_both_values = [[] for n in range(residues_number)]
    for frame in data:
        for r, residue in enumerate(frame):
            values = frame[residue]
            residues_vdw_values[r].append(values[0])
            residues_es_values[r].append(values[1])
            residues_both_values[r].append(values[2])

    # Calculate the residue averages from each energy
    residues_vdw_avg = [sum(v) / len(v) for v in residues_vdw_values]
    residues_es_avg = [sum(v) / len(v) for v in residues_es_values]
    residues_both_avg = [sum(v) / len(v) for v in residues_both_values]

    # Calculate the residue averages from each energy at the beginig and end of the trajectory
    # We take the initial 20% and the final 20% of frames to calculate each respectively
    p20_frames = len(data)
    p20 = round(p20_frames*0.2)
    if p20 == 0:
        p20 = 1

    # Initials
    residues_vdw_values_initial = [[] for n in range(residues_number)]
    residues_es_values_initial = [[] for n in range(residues_number)]
    residues_both_values_initial = [[] for n in range(residues_number)]
    for frame in data[:p20]:
        for r, residue in enumerate(frame):
            values = frame[residue]
            residues_vdw_values_initial[r].append(values[0])
            residues_es_values_initial[r].append(values[1])
            residues_both_values_initial[r].append(values[2])

    residues_vdw_avg_initial = [sum(v) / len(v)
                                for v in residues_vdw_values_initial]
    residues_es_avg_initial = [sum(v) / len(v)
                                for v in residues_es_values_initial]
    residues_both_avg_initial = [sum(v) / len(v)
                                    for v in residues_both_values_initial]

    # Finals
    residues_vdw_values_final = [[] for n in range(residues_number)]
    residues_es_values_final = [[] for n in range(residues_number)]
    residues_both_values_final = [[] for n in range(residues_number)]
    for frame in data[-p20:]:
        for r, residue in enumerate(frame):
            values = frame[residue]
            residues_vdw_values_final[r].append(values[0])
            residues_es_values_final[r].append(values[1])
            residues_both_values_final[r].append(values[2])

    residues_vdw_avg_final = [sum(v) / len(v)
                                for v in residues_vdw_values_final]
    residues_es_avg_final = [sum(v) / len(v)
                                for v in residues_es_values_final]
    residues_both_avg_final = [sum(v) / len(v)
                                for v in residues_both_values_final]

    # Format the results data and append it to the output data
    output = {
        'labels': residues_labels,
        'vdw': residues_vdw_avg,
        'es': residues_es_avg,
        'both': residues_both_avg,
        'ivdw': residues_vdw_avg_initial,
        'ies': residues_es_avg_initial,
        'iboth': residues_both_avg_initial,
        'fvdw': residues_vdw_avg_final,
        'fes': residues_es_avg_final,
        'fboth': residues_both_avg_final,
    }

    return output