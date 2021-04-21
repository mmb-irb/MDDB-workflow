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
import prody
import re
import numpy
import math
from subprocess import run, PIPE, Popen
import json

from model_workflow.tools.get_reduced_trajectory import get_reduced_trajectory

# Set the path to auxiliar files required for this analysis
resources = str(Path(__file__).parent.parent / "utils" / "resources")
reslib_source = resources + '/res.lib'
#preppdb_source = resources + '/preppdb.pl'
preppdb_source = resources + '/old_preppdb.pl' # DANI: Este no renombra los terminales
cmip_inputs_checkonly_source = resources + '/check.in'
cmip_inputs_source = resources + '/input.in'
vdw_source = resources + '/vdwprm'

# Perform the electrostatic and vdw energies analysis for each pair of interaction agents
def energies(
        input_topology_filename: str,
        input_trajectory_filename: str,
        output_analysis_filename: str,
        reference,
        snapshots: int,
        interactions: list):

    if not interactions or len(interactions) == 0:
        print('No interactions were specified')
        return

    # Correct the topology by renaming residues according to if they are terminals
    energies_topology = 'energies.pdb'
    energies_topology_prody = prody.parsePDB(input_topology_filename)
    energies_topology_prody = name_terminal_residues(energies_topology_prody)
    prody.writePDB(energies_topology, energies_topology_prody)

    # Given a pdb structure, use CMIP to extract energies
    # Output energies are already added by residues
    def get_frame_energy(frame_pdb):
        # Parse the pdb file to prody format and then correct it
        frame_structure = prody.parsePDB(frame_pdb)
        # Add chains according to the reference topology, since gromacs has deleted chains
        frame_structure.setChids(reference.topology.getChids())
        # Set element names according to the reference topology
        # Gromacs may mess some element names with two letters (e.g. 'CL')
        frame_structure.setElements(reference.topology.getElements())

        # WARNING: At this point topology should be corrected
        # WARNING: Repeated atoms will make the analysis fail

        # Repeat the whole process for each interaction
        data = []
        for interaction in interactions:

            name = interaction['name']
            
            # Copy the source 'res.lib' file in the local directory if it does not exists yet
            # If ACPYPE is used then the local res.lib will be modified
            # Theorically this should only happen once:
            # - The modified res.lib is conserved, so further frames should not miss any residue
            # - The modified res.lib is shared across interactions
            reslib_filename = 'res.lib'
            if not os.path.exists(reslib_filename):
                run([
                    "cp",
                    reslib_source,
                    reslib_filename,
                ], stdout=PIPE).stdout.decode()

            # Select the first agent, extract a pdb only with its atoms and parse it to CMIP
            # Run ACPYPE if required and then modify the res.lib file
            agent1 = interaction['agent_1'].replace(' ', '_').replace('/', '_')
            agent1_pdb = agent1 + '.pdb'
            agent1_selection = frame_structure.select(interaction['selection_1'])
            agent1_cmip = selection2cmip(agent1_selection, agent1, reslib_filename)

            # Repeat the process with agent 2
            agent2 = interaction['agent_2'].replace(' ', '_').replace('/', '_')
            agent2_pdb = agent2 + '.pdb'
            agent2_selection = frame_structure.select(interaction['selection_2'])
            agent2_cmip = selection2cmip(agent2_selection, agent2, reslib_filename)

            # Copy the source cmip inputs file in the local directory
            # Inputs will be modified to adapt the cmip grid to both agents together
            # Note than modified inputs file is not conserved along frames
            # Structures may change along trajectory thus requiring a different grid size
            cmip_inputs = 'cmip.in'
            logs = run([
                "cp",
                cmip_inputs_source,
                cmip_inputs,
            ], stdout=PIPE).stdout.decode()

            adapt_cmip_grid(agent1_cmip, agent2_cmip, cmip_inputs)

            # Run the CMIP software to get the desired energies
            agent1_energy = agent1 + '.energy.pdb'
            agent1_data = get_cmip_energies(cmip_inputs, agent1_cmip, agent2_cmip, agent1_energy)
            agent2_energy = agent2 + '.energy.pdb'
            agent2_data = get_cmip_energies(cmip_inputs, agent2_cmip, agent1_cmip, agent2_energy)

            data.append({ 'agent1': agent1_data, 'agent2': agent2_data })

            # Delete current agent pdb files before going for the next interaction
            run([
                "rm",
                agent1_pdb,
                agent2_pdb,
                agent1_cmip,
                agent2_cmip,
                agent1_energy,
                agent2_energy,
            ], stdout=PIPE).stdout.decode()

        return data

    # Set the frames where we extract energies to calculate the average
    # WARNING: The gromacs '-fr' option counts frames starting at 1, not at 0
    frames = range(1, snapshots +1)

    # Set a maximum of frames
    # If trajectory has more frames than the limit create a reduced trajectory
    energies_trajectory_filename = input_trajectory_filename
    frames_number = 100
    if snapshots > frames_number:
        energies_trajectory_filename = 'energies.trajectory.xtc'
        # WARNING: The gromacs '-fr' option counts frames starting at 1, not at 0
        frames = range(1, frames_number +1)  # if frames_number > 1 else [1]
        #if not os.path.exists(energies_trajectory_filename): # DANI: Esto es temporal
        get_reduced_trajectory(
            energies_topology,
            input_trajectory_filename,
            energies_trajectory_filename,
            snapshots,
            frames_number,
        )
    else:
        frames_number = snapshots

    # Extract the energies for each frame
    frames_ndx = 'frames.ndx'
    interactions_data = [[] for i in interactions]
    for f in frames:
        print('Frame ' + str(f) + ' / ' + str(frames_number))
        # Extract the current frame
        current_frame = 'frame' + str(f) + '.pdb'
        # The frame selection input in gromacs works with a 'ndx' file
        with open(frames_ndx, 'w') as file:
            file.write('[frames]\n' + str(f))
        p = Popen([
            "echo",
            "System",
        ], stdout=PIPE)
        logs = run([
            "gmx",
            "trjconv",
            "-s",
            energies_topology,
            "-f",
            energies_trajectory_filename,
            '-o',
            current_frame,
            "-fr",
            frames_ndx,
            '-quiet'
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()
        # Run the main analysis over the current frame
        # Append the result data for each interaction
        frame_energies_data = get_frame_energy(current_frame)
        for i, data in enumerate(frame_energies_data):
            interactions_data[i].append(data)
        # Delete current frame files before going for the next frame
        run([
            "rm",
            current_frame,
            frames_ndx,
        ], stdout=PIPE).stdout.decode()

    # Now calculated residue average values through all frames for each pair of interaction agents
    output_analysis = []
    for i, interaction in enumerate(interactions):

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

    # Finally remove the reduced trajectory since it is not required anymore
    if energies_trajectory_filename == 'energies.trajectory.xtc':
        logs = run([
            "rm",
            energies_topology,
            energies_trajectory_filename,
        ], stdout=PIPE).stdout.decode()

# Save all unique residues in res.lib
def get_reslib_residues(reslib_filename : str) -> list:
    residues = []
    with open(reslib_filename, 'r') as file:
        lines = file.readlines()
        for line in lines:
            resname = line[0:4].strip()
            residues.append(resname)
    return list(set(residues))

# Transform prody selection to a cmip input pdb, which includes charges
# Charges are taken from the 'reslib' provided file
# Use ACPYPE to calculate charges which are not in the 'reslib' by default
# Modify the 'reslib' if ACPYPE is used by adding the new charges
def selection2cmip(input_selection, selection_name, reslib_filename) -> str:
    print('Setting ' + selection_name + ' charges')
    selection = input_selection
    name = selection_name

    # Set a pdb filename to write the selection when required
    selection_filename = name + '.pdb'

    # Correct the protein by renaming residues according to if they are terminals
    #selection = name_terminal_residues(selection)
    selection = selection.toAtomGroup()

    # Get a list with all unique reslib residues
    reslib_residues = get_reslib_residues(reslib_filename)

    # Check if any residue in the selection is not found in the reslib residues
    # If so, it means it is not a regular protein or nucleic acid
    # Charges must be calculated thorugh ACPYPE
    acpype_is_required = False
    for residue in selection.iterResidues():
        residue_name = residue.getResname()
        if residue_name not in reslib_residues:
            print('There is at least 1 residue not found in res.lib: ' + residue_name)
            print('ACPYPE will be run')
            acpype_is_required = True
            break
    
    # When required, calculate the energies with acpype and then mine them from the output file
    energies = []
    if acpype_is_required:

        # Remove all residue names to prevent duplicated residue names, which makes ACPYPE fail
        selection = remove_residue_names(selection)

        # Write the selection to a pdb file since it is the accepted ACPYPE format
        prody.writePDB(selection_filename, selection)

        # Run acpype
        # test = !obabel -ipdb $selection_filename -omol2 whatever
        acpype_logs = run([
            "acpype",
            "-i",
            selection_filename,
            # "-n",
            # "0",
        ], stdout=PIPE).stdout.decode()

        energies_file = name + '.acpype/' + name + '.mol2'
        # Check if the output file exists. If not, send error and print acpype logs
        if not os.path.exists(energies_file):
            print(acpype_logs)
            raise SystemExit(
                "ERROR: Something was wrong with ACPYPE")
        # Write all lines but the last line: 'END'
        energies_lines = None
        with open(energies_file, "r") as file:
            energies_lines = file.readlines()

        # Mine the energies in the acpype output file
        # These energies are used by CMIP as a reference
        # Count how many lines are energy 0
        zero_count = 0
        for line in energies_lines:
            line = line.split()
            if len(line) == 9:
                energy_value = line[-1]
                if float(energy_value) == 0:
                    zero_count += 1
                energies.append(energy_value)

        # If all lines are energy 0 return here
        if zero_count == len(energies):
            raise SystemExit(
                "ERROR: All ACPYPE output energies are zero :/")

        # Check the number of atoms matches the number of energies
        if len(list(selection.iterAtoms())) != len(energies):
            raise SystemExit(
                "Stop!! The number of atoms and energies does not match :/")

    # Adapt atom elements to match what CMIP expects to find
    # Hydrogens bonded to oxygen or nitrogen or sulfur must be renamed as HO or HN
    # CMIP uses atom names to set the atom radius so this step is critical
    # Overwrite the previous 'selection' value and the selection pdb file
    # WARNING: It is very important to make the CMIP adaptation in this specific place!!
    # WARNING: If it is done before running ACPYPE then ACPYPE will fail
    # WARNING: It is is done after editing the reslib file then elements will not match
    selection = adapt_cmip_atoms(selection)
    prody.writePDB(selection_filename, selection)

    if acpype_is_required:

        # Modify the reslib file with the new ACPYPE charges
        with open(reslib_filename, "a") as reslib:

            # Prepare the standard formatted line for each atom
            # Then write it into the 'res.lib' file
            atoms = selection.iterAtoms()
            for i, atom in enumerate(atoms):
                residue = atom.getResname().ljust(5)
                atomname = atom.getName().ljust(5)
                element = atom.getElement().ljust(6)
                energy = energies[i]
                if energy[0] != '-':
                    energy = ' ' + energy
                line = residue + atomname + element + energy
                reslib.write(line + '\n')

    # Create the pdb of the current corrected selection
    prody.writePDB(selection_filename, selection)
    # Run a perl script which assigns charges to each atom
    # This modification is what makes the pdb a valid CMIP input
    output_filename = name + '.cmip.pdb'
    run([
        "perl",
        preppdb_source,
        reslib_filename,
        selection_filename,
        output_filename,
    ], stdout=PIPE).stdout.decode()

    # It is very usual that some atoms are not found in the res.lib file
    # When this happens, the output cmip file has '????' instead of a valid charge
    # Check the output file to be clean of question marks before returning
    with open(output_filename, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if '????' in line:
                # There is an irregularity in the topology and it must be solved by hand
                raise SystemExit('ERROR: Atom not found in res.lib: \n' + line)

    # Finally, remove the pdb header since sometimes it makes cmip fail
    with open(output_filename, 'r+') as file:
        lines = file.readlines()
        file.seek(0)
        for line in lines:
            if line[0:4] == 'ATOM':
                file.write(line)
        file.truncate()

    return output_filename

# Change elements in a prody selection to meet the CMIP requirements
# Hydrogens bonded to carbons remain as 'H'
# Hydrogens bonded to oxygen are renamed as 'HO'
# Hydrogens bonded to nitrogen or sulfur are renamed as 'HN'
# Some heavy atom elements may be also modified (e.g. 'CL' -> 'Cl')
def adapt_cmip_atoms(selection):

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
        residue = selection.select('resnum ' + str(atom_resnum))
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

    # Update the element of each hydrogen according to CMIP needs
    for atom in selection.iterAtoms():
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
                continue
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
        # Set the correct element
        atom.setElement(element)
    # Return the corrected prody topology
    return selection

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
    cmip_checkonly_output = 'checkonly.energy.pdb'

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

    # Calculate grid dimensions for a new grid which contains both previous grids
    new_center, new_density = compute_new_grid(
        agent1_center,
        agent1_density,
        agent1_units,
        agent2_center,
        agent2_density,
        agent2_units)

    # Set the new lines to be written in the local inputs file
    grid_inputs = [
        f" cenx={new_center[0]:.1f} \n",
        f" ceny={new_center[1]:.1f} \n",
        f" cenz={new_center[2]:.1f} \n",
        f" dimx={int(new_density[0])} \n",
        f" dimy={int(new_density[1])} \n",
        f" dimz={int(new_density[2])} \n",
        #" intx=0.5 \n",
        #" inty=0.5 \n",
        #" intz=0.5 \n",
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


# Set all residue names as 'UNK' in a prody selection to meet the ACPYPE requirements
# ACPYPE will return error if there are more than one different residue name
def remove_residue_names(selection):

    # Set a new selection
    new_selection = selection.toAtomGroup()

    # Get all residues
    residues = new_selection.iterResidues()

    # Change tehir names one by one
    for residue in residues:
        residue.setResname('UNK')

    # Return the corrected prody topology
    return new_selection

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

    # Mine the electrostatic (es) and Van der Walls (vdw) energies for each atom
    # Group the results by reidues adding their values
    residues = {}
    with open(cmip_output, 'r') as file:
        lines = list(file)
        # If this file is empty it means something went wrong with CMIP
        # We print its logs and exit
        if len(lines) == 0:
            for line in cmip_logs:
                print(line)
            raise SystemExit('ERROR: Something went wrong with CMIP!')
        for line in lines:
            chain = line[21:22]
            residue_id = line[22:28]
            residue = chain + ':' + str(int(residue_id))
            vdw = float(line[42:53])
            es = float(line[57:68])
            both = float(line[72:83])
            # Values greater than 100 are represented as 0
            # This step is performed to filter 'infinity' values
            energies = (vdw, es, both) if both < 100 else (0, 0, 0)
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
    frames_number = len(data)
    p20 = round(frames_number*0.2)

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