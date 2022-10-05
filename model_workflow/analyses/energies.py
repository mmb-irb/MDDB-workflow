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
from typing import Optional

import prody
import pytraj as pt

from model_workflow.tools.get_pdb_frames import get_pdb_frames

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
def energies(
        input_topology_filename: str,
        input_trajectory_filename: str,
        output_analysis_filename: str,
        interactions: list,
        charges: list,
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

    # Correct the topology by renaming residues according to if they are terminals
    energies_topology = energies_folder + '/energies.pdb'
    energies_topology_prody = prody.parsePDB(input_topology_filename)
    energies_topology_prody = name_terminal_residues(energies_topology_prody)
    prody.writePDB(energies_topology, energies_topology_prody)

    # Get each atom element in CMIP format
    elements = get_topology_cmip_elements(input_topology_filename)

    # Given a pdb structure, use CMIP to extract energies
    # Output energies are already added by residues
    def get_frame_energy(frame_pdb):
        # Parse the pdb file to prody format and then correct it
        frame_structure = prody.parsePDB(frame_pdb)
        # Add chains according to the reference topology, since gromacs has deleted chains
        # DANI: Esto ya no debería hacer falta ya que ahora las frames vienen de pytraj, no de Gromacs
        #frame_structure.setChids(reference.topology.getChids())

        # Set charges and elements according to the topology mined data
        frame_structure.setElements(elements)
        frame_structure.setCharges(charges)

        # WARNING: At this point topology should be corrected
        # WARNING: Repeated atoms will make the analysis fail

        # Repeat the whole process for each interaction
        data = []
        for interaction in interactions:

            # Check if the interaction as been marked as 'exceeds', in which case we skip it
            if interaction.get('exceeds', False):
                continue

            name = interaction['name']
            strong_bonds = interaction.get('strong_bonds', None)

            # Select the first agent and extract it in a CMIP friendly pdb format
            agent1_name = interaction['agent_1'].replace(' ', '_').replace('/', '_')
            agent1_selection = frame_structure.select(interaction['selection_1'])
            if not agent1_selection:
                raise SystemExit('ERROR: Agent "' + interaction['agent_1'] + '" with selection "' + 
                    interaction['selection_1'] + '" has no atoms')
            agent1_cmip = selection2cmip(agent1_name, agent1_selection, strong_bonds)

            # Repeat the process with agent 2
            agent2_name = interaction['agent_2'].replace(' ', '_').replace('/', '_')
            agent2_selection = frame_structure.select(interaction['selection_2'])
            if not agent2_selection:
                raise SystemExit('ERROR: Agent "' + interaction['agent_2'] + '" with selection "' + 
                    interaction['selection_2'] + '" has no atoms')
            agent2_cmip = selection2cmip(agent2_name, agent2_selection, strong_bonds)

            # Copy the source cmip inputs file in the local directory
            # Inputs will be modified to adapt the cmip grid to both agents together
            # Note than modified inputs file is not conserved along frames
            # Structures may change along trajectory thus requiring a different grid size
            cmip_inputs = energies_folder + '/cmip.in'
            run([ "cp", cmip_inputs_source, cmip_inputs ], stdout=PIPE)

            # Set the CMIP box dimensions and densities to fit both the host and the guest
            adapt_cmip_grid(agent1_cmip, agent2_cmip, cmip_inputs)

            # Run the CMIP software to get the desired energies
            agent1_output_energy = energies_folder + '/' + agent1_name + '.energy.pdb'
            agent1_data = get_cmip_energies(cmip_inputs, agent1_cmip, agent2_cmip, agent1_output_energy)
            agent2_output_energy = energies_folder + '/' + agent2_name + '.energy.pdb'
            agent2_data = get_cmip_energies(cmip_inputs, agent2_cmip, agent1_cmip, agent2_output_energy)

            data.append({ 'agent1': agent1_data, 'agent2': agent2_data })

            # Delete current agent pdb files before going for the next interaction
            run([
                "rm",
                agent1_cmip,
                agent2_cmip,
                agent1_output_energy,
                agent2_output_energy,
            ], stdout=PIPE).stdout.decode()

        return data

    # Extract the energies for each frame in a reduced trajectory
    frames, step, count = get_pdb_frames(input_topology_filename, input_trajectory_filename, frames_limit)
    interactions_data = [[] for interaction in interactions if not interaction.get('exceeds', False)]
    for current_frame in frames:
        
        # Run the main analysis over the current frame
        # Append the result data for each interaction
        frame_energies_data = get_frame_energy(current_frame)
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
        energies_topology,
    ], stdout=PIPE).stdout.decode()

# Transform prody selection to a cmip input pdb, which includes charges
# Charges have been previously taken from the charges topology and injected in prody
def selection2cmip(agent_name : str, agent_selection : str, strong_bonds : Optional[list]) -> str:
    
    print('Setting ' + agent_name + ' charges')
    atoms = agent_selection.iterAtoms()

    strong_bond_indexes = []
    if strong_bonds:
        for bond in strong_bonds:
            strong_bond_indexes += bond

    # Write a special pdb which contains charges as CMIP expects to find them
    output_filename = energies_folder + '/' + agent_name + '.cmip.pdb'
    with open(output_filename, "w") as file:

        for a, atom in enumerate(atoms):
            
            index = str(a+1).rjust(5)
            atomname = atom.getName()
            name =  ' ' + atomname.ljust(3) if len(atomname) < 4 else atomname
            residue_name = atom.getResname().ljust(3)
            chain = atom.getChid().rjust(1)
            residue_number = str(atom.getResnum()).rjust(4)
            icode = atom.getIcode().rjust(1)
            coords = atom.getCoords()
            x_coord, y_coord, z_coord = [ "{:.3f}".format(coord).rjust(8) for coord in coords ]
            charge = "{:.4f}".format(atom.getCharge())
            # In case this atom is making an strong bond between both interacting agents we add an 'X' before the element
            # This way CMIP will ignore the atom. Otherwise it would return high non-sense Van der Waals values
            real_index = atom.getIndex()
            cmip_ignore_flag = 'X' if real_index in strong_bond_indexes else ''
            element = cmip_ignore_flag + atom.getElement()
            
            atom_line = ('ATOM  ' + index + ' ' + name + ' ' + residue_name + ' '
                + chain + residue_number + icode + '   ' + x_coord + y_coord + z_coord
                + ' ' + str(charge).rjust(7) + '  ' + element + '\n')
            file.write(atom_line)

    return output_filename

# Given a topology (e.g. pdb, prmtop), extract the atom elements in a CMIP friendly format
# Hydrogens bonded to carbons remain as 'H'
# Hydrogens bonded to oxygen are renamed as 'HO'
# Hydrogens bonded to nitrogen or sulfur are renamed as 'HN'
# Some heavy atom elements may be also modified (e.g. 'CL' -> 'Cl')
# There are two ways to do this: the canonical (faster) and the alternative (error proof)
def get_topology_cmip_elements (input_topology_filename : str):
    try:
        elements = get_topology_cmip_elements_canonical(input_topology_filename)
    except Exception as err:
        print(err)
        print('The canonical elements mining failed. Retrying with alternative mining')
        elements = get_topology_cmip_elements_alternative(input_topology_filename)
    return elements

# Use pytraj for this task
def get_topology_cmip_elements_canonical (input_topology_filename : str):
    
    topology = pt.load_topology(filename=input_topology_filename)

    # Set all supported elements
    # This is required to transform element names (returned by pytraj) to element letters
    standard_elements = {
        'hydrogen': 'H',
        'carbon': 'C',
        'oxygen': 'O',
        'nitrogen': 'N',
        'sulfur': 'S',
        'sodium': 'Na',
        'chlorine': 'Cl',
        'zinc': 'Zn',
        'fluorine': 'F',
        'magnesium': 'Mg',
        'phosphorus': 'P',
    }
    # Iterate over each atom to save their CMIP element
    elements = []
    atoms = list(topology.atoms)
    for a, atom in enumerate(atoms):
        residue = atom.resname
        # Skip this atom if we already found it
        name = atom.name
        element = standard_elements[atom.element]
        # Adapt hydrogens element to CMIP requirements
        if element == 'H':
            # There should we always only 1 bond
            # If you have the error below you may need to updated the pytraj version or reintsall pytraj
            # ValueError: Buffer dtype mismatch, expected 'int' but got 'long'
            bonded_heavy_atom_index = atom.bonded_indices()[0]
            bonded_heavy_atom = atoms[bonded_heavy_atom_index]
            bonded_heavy_atom_element = standard_elements[bonded_heavy_atom.element]
            # Hydrogens bonded to carbons remain as 'H'
            if bonded_heavy_atom_element == 'C':
                pass
            # Hydrogens bonded to oxygen are renamed as 'HO'
            elif bonded_heavy_atom_element == 'O':
                element = 'HO'
            # Hydrogens bonded to nitrogen or sulfur are renamed as 'HN'
            elif bonded_heavy_atom_element == 'N' or bonded_heavy_atom_element == 'S':
                element = 'HN'
            else:
                raise SystemExit(
                    'ERROR: Hydrogen bonded to not supported heavy atom: ' + bonded_heavy_atom_element)
        elements.append(element)
    return elements

# Use prody for this task
def get_topology_cmip_elements_alternative (input_topology_filename : str):

    topology = prody.parsePDB(input_topology_filename)

    # Try to guess the atom element from the name of the atom
    # This is used only when element is missing
    cmip_supported_elements = ['Cl', 'Na', 'Zn', 'Mg', 'C', 'N', 'P', 'S',
                                'HW', 'HO', 'HN', 'H', 'F', 'OW', 'O', 'IP', 'IM', 'I']
    def guess_name_element(name: str) -> str:
        length = len(name)
        next_character = None
        for i, character in enumerate(name):
            # Get the next character, since element may be formed by 2 letters
            if i < length - 1:
                next_character = name[i+1]
                # If next character is not a string then ignore it
                if not next_character.isalpha():
                    next_character = None
            # Try to get all possible matches between the characters and the supported atoms
            # First letter is always caps
            character = character.upper()
            # First try to match both letters together
            if next_character:
                # Start with the second letter in caps
                next_character = next_character.upper()
                both = character + next_character
                if both in cmip_supported_elements:
                    return both
                # Continue with the second letter in lowers
                next_character = next_character.lower()
                both = character + next_character
                if both in cmip_supported_elements:
                    return both
            # Finally, try with the first character alone
            if character in cmip_supported_elements:
                return character
        raise SystemExit(
            "ERROR: Not recognized element in '" + name + "'")

    # Set functions to find the closest atom of a specified atom
    def get_distance(coords1, coords2):
        squared_distance = numpy.sum((coords1 - coords2)**2, axis=0)
        distance = numpy.sqrt(squared_distance)
        return distance

    # Find the closest atom to the provided atom
    # Since this function is only used with hydrogen atoms we only search in the current atom residue
    def find_closest_atom(atom):
        atom_resnum = atom.getResnum()
        residue = topology.select('resnum ' + str(atom_resnum))
        residue_atoms = list(residue.iterAtoms())
        residue_coords = residue.getCoords()
        current_coords = atom.getCoords()
        distances = [get_distance(current_coords, c)
                        for c in residue_coords]
        sorted_distances = [d for d in distances]
        sorted_distances.sort()
        # We take the second minimum, since the first minimum will be always 0
        smallest_distance = sorted_distances[1]
        index = distances.index(smallest_distance)
        closest_atom = residue_atoms[index]
        return closest_atom

    # Harvest the element of each atom
    # Update the element of each hydrogen according to CMIP needs
    elements = []
    for atom in topology.iterAtoms():
        # First of all, correct tha name by moving numbers from the start to the end
        # e.g. 1HD1 -> HD11
        name = atom.getName()
        characters = name
        for character in characters:
            if not character.isalpha():
                name = name[1:] + name[0]
            else:
                break
        atom.setName(name)
        # Then, correct the element
        # If element is missing try to guess it from the name
        element = atom.getElement()
        if not element:
            element = guess_name_element(name)
            atom.setElement(element)
        # Find hydrogens by element
        # WARNING: Avoid finding hydrogens by name. It is very risky
        if element == 'H':
            bonded_heavy_atom = find_closest_atom(atom)
            bonded_heavy_atom_element = bonded_heavy_atom.getElement()
            if not bonded_heavy_atom_element:
                bonded_heavy_atom_element = guess_name_element(bonded_heavy_atom.getName())
                bonded_heavy_atom.setElement(bonded_heavy_atom_element)
            # Hydrogens bonded to carbons remain as 'H'
            if bonded_heavy_atom_element == 'C':
                pass
            # Hydrogens bonded to oxygen are renamed as 'HO'
            elif bonded_heavy_atom_element == 'O':
                element = 'HO'
            # Hydrogens bonded to nitrogen or sulfur are renamed as 'HN'
            elif bonded_heavy_atom_element == 'N' or bonded_heavy_atom_element == 'S':
                element = 'HN'
            else:
                raise SystemExit(
                    'ERROR: Hydrogen bonded to not supported heavy atom: ' + bonded_heavy_atom_element)
        # Update other elements naming
        elif element == 'CL':
            element = 'Cl'
        elif element == 'BR':
            element = 'Br'
        elif element == 'ZN':
            element = 'Zn'
        elif element == 'NA':
            element = 'Na'
        elif element == 'MG':
            element = 'Mg'
        # Get the correct element
        atom.setElement(element)
        elements.append(element)
    # Return the corrected prody topology
    return elements

protein_residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
                    'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

rna_residues = ['RA', 'RU', 'RC', 'RG']

# Change residue names in a prody selection to meet the CMIP requirements
# Change terminal residue names by adding an 'N' or 'C'
def name_terminal_residues(selection):
    
    # Set a new selection
    new_selection = selection.toAtomGroup()

    # Get all atoms and atom coordinates
    chains = new_selection.iterChains()

    for chain in chains:

        residues = list(chain.iterResidues())

        # Check if the first residue is tagged as a terminal residue
        # If not, rename it
        first_residue = residues[0]
        first_residue_name = first_residue.getResname()
        # In case it is a protein
        if first_residue_name in protein_residues:
            first_residue_name += 'N'
            first_residue.setResname(first_residue_name)
        # In case it is RNA
        elif first_residue_name in rna_residues:
            first_residue_name += '5'
            first_residue.setResname(first_residue_name)

        # Check if the last residue is tagged as 'C' terminal
        # If not, rename it
        last_residue = residues[-1]
        last_residue_name = last_residue.getResname()
        # In case it is a protein
        if last_residue_name in protein_residues:
            last_residue_name += 'C'
            last_residue.setResname(last_residue_name)
        # In case it is RNA
        elif last_residue_name in rna_residues:
            last_residue_name += '3'
            last_residue.setResname(last_residue_name)

    # Return the corrected prody topology
    return new_selection

# Run CMIP in 'checkonly' mode (i.e. start and stop) for both agents
# Mine the grid generated by CMIP for both agents and calculate a new grid which would include both
# Finally, modify the provided 'cmip_inputs' with the new grid parameters
def adapt_cmip_grid(agent1, agent2, cmip_inputs):
    print('Adapting grid')

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
        agent1,
        "-vdw",
        vdw_source,
        "-hs",
        agent2,
        "-byat",
        cmip_checkonly_output,
    ], stdout=PIPE).stdout.decode()

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
        agent2,
        "-vdw",
        vdw_source,
        "-hs",
        agent1,
        "-byat",
        cmip_checkonly_output,
    ], stdout=PIPE).stdout.decode()

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

def mine_cmip_output(logs):
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
def compute_new_grid(
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
def get_cmip_energies(cmip_inputs, pr, hs, cmip_output):
    print('Calculating energies')
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
    ], stdout=PIPE).stdout.decode()
    # Get back to the main folder
    os.chdir("..")

    # Mine the electrostatic (es) and Van der Walls (vdw) energies for each atom
    # Group the results by residues adding their values
    residues = {}
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
                print('WARNING: We have extremly high values in energies which are beeing discarded')
                energies = (0, 0, 0)
            if residue in residues:
                residues[residue] = tuple(
                    [a+b for a, b in zip(energies, residues[residue])])
            else:
                residues[residue] = energies

    return residues

# Format data grouping atom energies by residue
# Calculate means of the whole
def format_data(data):

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