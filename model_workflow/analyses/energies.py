# Energies
#
# The energies analysis is carried by ACPYPE and a locally developed tool: CMIP.
# ACPYPE is a tool based in Python to use Antechamber to generate topologies for chemical compounds and to interface with others python applications.
# CMIP stands for Classical Molecular Interaction Potential and it is usefull to predict electrostatic and Van der Waals potentials.
# Energies analysis is performed only in protein-ligand trajectories. First, ACPYPE is used to calculate the ligand (small molecule) atom charges. Then, CMIP is used to calculate the electrostatic potential and Van der Waals forces of each atom in both the protein and the ligand.
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
preppdb_source = resources + '/preppdb.pl'  # For proteins
old_preppdb_source = resources + '/old_preppdb.pl'  # For ligands
cmip_inputs_checkonly_source = resources + '/check.in'
cmip_inputs_source = resources + '/input.in'
vdw_source = resources + '/vdwprm'


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
# The new grid is expected to cover both input grids: protein grid and ligand grid
def compute_new_grid(
        prot_center,
        prot_density,
        prot_units,
        lig_center,
        lig_density,
        lig_units,
        extra_density=1):
    new_center = []  # [(i + j) / 2 for i, j in zip(prot_center, lig_center)]
    new_density = []

    for k in range(3):
        min_prot = prot_center[k] - prot_density[k] * prot_units[k]
        min_lig = lig_center[k] - lig_density[k] * lig_units[k]
        min_new = min(min_prot, min_lig)
        max_prot = prot_center[k] + prot_density[k] * prot_units[k]
        max_lig = lig_center[k] + lig_density[k] * lig_units[k]
        max_new = max(max_prot, max_lig)

        dnew = int(abs(max_new - min_new) + extra_density)
        cnew = min_new + (dnew / 2)
        new_density.append(dnew)
        new_center.append(cnew)
    return new_center, new_density


# Perform the electrostatic and vdw energies analysis for each ligand
# DANI: En principio soporta casos en que hay multiples ligandos, pero no se ha provado
def energies(
        input_topology_filename: str,
        input_trajectory_filename: str,
        output_analysis_filename: str,
        reference,
        snapshots: int,
        ligands: list):

    if not ligands or len(ligands) == 0:
        print('No ligands were specified')
        return

    # Change elements in a prody selection to meet the CMIP requirements
    # Hydrogens bonded to carbons remain as 'H'
    # Hydrogens bonded to oxygen are renamed as 'HO'
    # Hydrogens bonded to nitrogen or sulfur are renamed as 'HN'
    # Some heavy atom elements may be also modified (e.g. 'CL' -> 'Cl')
    def adapt_cmip_atoms(selection):
        # Get all atoms and atom coordinates
        selection_coords = selection.getCoords()
        selection_atoms = list(selection.iterAtoms())

        # Set functions to find the closest atom of a specified atom
        def get_distance(coords1, coords2):
            squared_distance = numpy.sum((coords1 - coords2)**2, axis=0)
            distance = numpy.sqrt(squared_distance)
            return distance

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

        def find_closest_atom(atom):
            current_coords = atom.getCoords()
            distances = [get_distance(current_coords, c)
                         for c in selection_coords]
            sorted_distances = [d for d in distances]
            sorted_distances.sort()
            # We take the second minimum, since the first minimum will be always 0
            smallest_distance = sorted_distances[1]
            index = distances.index(smallest_distance)
            closest_atom = selection_atoms[index]
            return closest_atom

        # Update the element of each hydrogen according to CMIP needs
        for atom in selection_atoms:
            # First of all, correct tha name by removing numbers at the start
            name = atom.getName()
            for character in name:
                if not character.isalpha():
                    name = name[1:]
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
                bonded_heavy_atom = find_closest_atom(atom).getElement()
                # Hydrogens bonded to carbons remain as 'H'
                if bonded_heavy_atom == 'C':
                    continue
                # Hydrogens bonded to oxygen are renamed as 'HO'
                elif bonded_heavy_atom == 'O':
                    element = 'HO'
                # Hydrogens bonded to nitrogen or sulfur are renamed as 'HN'
                elif bonded_heavy_atom == 'N' or bonded_heavy_atom == 'S':
                    element = 'HN'
                else:
                    raise SystemExit(
                        'ERROR: Hydrogen bonded to not supported heavy atom')
            # Update other elements naming
            if element == 'CL':
                element = 'Cl'
            if element == 'BR':
                element = 'Br'
            if element == 'ZN':
                element = 'Zn'
            if element == 'NA':
                element = 'Na'
            if element == 'MG':
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

    # Given a pdb structure, use CMIP to extract energies
    # Output energies are already added by residues
    def get_frame_energy(frame_pdb):
        # Parse the pdb file to prody format and then correct it
        original_topology = prody.parsePDB(frame_pdb)
        # Add chains according to the reference topology, since gromacs has deleted chains
        original_topology.setChids(reference.topology.getChids())
        # Set element names according to the reference topology
        # Gromacs may mess some element names with two letters (e.g. 'CL')
        original_topology.setElements(reference.topology.getElements())

        # WARNING: At this point topology should be corrected
        # WARNING: Repeated atoms will make the analysis fail

        # Prepare the protein files
        # Create a new pdb only with the protein
        protein_pdb = 'protein.pdb'
        selection = original_topology.select('protein')
        # Correct the protein by renaming residues according to if they are terminals
        selection = name_terminal_residues(selection)
        prody.writePDB(protein_pdb, selection)

        # Get the cmip protein input file
        protein_cmip = 'protein.cmip.pdb'
        run([
            "perl",
            preppdb_source,
            reslib_source,
            protein_pdb,
            protein_cmip,
        ], stdout=PIPE).stdout.decode()

        # Prepare the ligand files for each ligand
        ligand_data = []
        for ligand in ligands:

            # Create a new pdb only with the current ligand
            # Adapt this topology to ACPYPE by changing some element names
            # WARNING: Name spaces must be replaced by underscore or processes will fail
            name = ligand['name'].replace(' ', '_')
            ligand_pdb = name + '.pdb'
            selection = original_topology.select(ligand['prody'])
            # If acypype is required we remove all residue names to prevent duplicated residue names
            # If no acpype is run it means the ligand is known by the default cmip reslib
            # For this case, we better check the terminal residues to be well set
            acpype_is_required = ligand.get('acpype', True)
            if acpype_is_required:
                selection = remove_residue_names(selection)
            else:
                selection = name_terminal_residues(selection)
            prody.writePDB(ligand_pdb, selection)

            # Set the source reslib file as the reslib to be used further by CMIP
            # Then, if acpype is runned, we copy and modify the source reslib
            reslib_filename = reslib_source

            # Calculate the energies with acpype and then mine them
            # This process is done only if the ligand value 'acpype' is True
            energies = []

            if acpype_is_required:

                # Run acpype
                # test = !obabel -ipdb $ligand_pdb -omol2 whatever
                acpype_logs = run([
                    "acpype",
                    "-i",
                    ligand_pdb,
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
                        "All energies are zero! Check files and try again")

                # Check the number of atoms matches the number of energies
                if len(list(selection.iterAtoms())) != len(energies):
                    raise SystemExit(
                        "Stop!! The number of atoms and energies does not match :/")

            # Adapt atom elements to match what CMIP expects to find
            # Hydrogens bonded to oxygen or nitrogen or sulfur must be renamed as HO or HN
            # CMIP uses atom names to set the atom radius so this step is critical
            # Overwrite the previous 'selection' value and the ligand pdb file
            # WARNING: It is very important to make the CMIP adaptation in this specific place!!
            # WARNING: If it is done before running ACPYPE then ACPYPE will fail
            # WARNING: It is is done after editing the reslib file then elements will not match
            selection = adapt_cmip_atoms(selection)
            prody.writePDB(ligand_pdb, selection)

            if acpype_is_required:

                # Copy the 'res.lib' file in the local path and open it to 'a'ppend new text
                reslib_filename = name + '_res.lib'
                run([
                    "cp",
                    reslib_source,
                    reslib_filename,
                ], stdout=PIPE).stdout.decode()

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

            # Get the cmip ligand input file
            ligand_cmip = name + '.cmip.pdb'
            logs = run([
                "perl",
                old_preppdb_source,
                reslib_filename,
                ligand_pdb,
                ligand_cmip,
            ], stdout=PIPE).stdout.decode()
            #print (logs)

            # Remove the pdb header since sometimes it makes cmip failing
            pdb_content = None
            with open(ligand_cmip, 'r+') as file:
                pdb_content = file.readlines()
            with open(ligand_cmip, 'w') as file:
                for line in pdb_content:
                    if line[0:4] == 'ATOM':
                        file.write(line)

            cmip_output = 'protein-' + name + '.energy.pdb'

            # Run CMIP in 'checkonly' mode and save the grid dimensions output
            # First do it for the protein
            cmip_logs_protein = run([
                "cmip",
                "-i",
                cmip_inputs_checkonly_source,
                "-pr",
                ligand_cmip,
                "-vdw",
                vdw_source,
                "-hs",
                protein_cmip,
                "-byat",
                cmip_output,
            ], stdout=PIPE).stdout.decode()

            # Mine the grid dimensions from CMIP logs
            prot_center, prot_density, prot_units = mine_cmip_output(
                cmip_logs_protein.split("\n"))

            # Run CMIP in 'checkonly' mode and save the grid dimensions output
            # Now do it for the ligand
            cmip_logs_ligand = run([
                "cmip",
                "-i",
                cmip_inputs_checkonly_source,
                "-pr",
                protein_cmip,
                "-vdw",
                vdw_source,
                "-hs",
                ligand_cmip,
                "-byat",
                cmip_output,
            ], stdout=PIPE).stdout.decode()

            # Mine the grid dimensions from CMIP logs
            lig_center, lig_density, lig_units = mine_cmip_output(
                cmip_logs_ligand.split("\n"))

            # Calculate grid dimensions for a new grid which contains both previous grids
            new_center, new_density = compute_new_grid(
                prot_center,
                prot_density,
                prot_units,
                lig_center,
                lig_density,
                lig_units)

            # Copy the cmip inputs source file locally and write the new grid parameters on it
            # Set the name for the local inputs file
            cmip_inputs = 'cmip.in'
            # Copy the source
            logs = run([
                "cp",
                cmip_inputs_source,
                cmip_inputs,
            ], stdout=PIPE).stdout.decode()
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

            # Run the cmip software to get the desired energies
            cmip_output = 'protein-' + name + '.energy.pdb'
            cmip_logs = run([
                "cmip",
                "-i",
                cmip_inputs,
                "-pr",
                protein_cmip,
                "-vdw",
                vdw_source,
                "-hs",
                ligand_cmip,
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

            ligand_data.append(residues)

        return ligand_data

    # Set the frames where we extract energies to calculate the average
    frames = range(1, snapshots)

    # Set a maximum of frames
    # If trajectory has more frames than the limit create a reduced trajectory
    energies_trajectory_filename = input_trajectory_filename
    frames_number = 100
    if snapshots > frames_number:
        energies_trajectory_filename = 'energies.trajectory.xtc'
        get_reduced_trajectory(
            input_topology_filename,
            input_trajectory_filename,
            energies_trajectory_filename,
            snapshots,
            frames_number,
        )
        frames = range(1, frames_number)  # if frames_number > 1 else [1]

    # Extract the energies for each frame
    #data = []
    frames_ndx = 'frames.ndx'
    ligands_data = [[] for l in ligands]
    for f in frames:
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
            input_topology_filename,
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
        # Append the result data for each ligand
        energies_data = get_frame_energy(current_frame)
        for i, data in enumerate(energies_data):
            ligands_data[i].append(data)
        # Delete current frame files before going for the next frame
        run([
            "rm",
            current_frame,
            frames_ndx,
        ], stdout=PIPE).stdout.decode()

    # Now calculated residue average values through all frames for each ligand
    output_analysis = []
    for i, ligand in enumerate(ligands):

        # Get the main data
        data = ligands_data[i]

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
            'name': ligand['name'],
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
        output_analysis.append(output)

    # Finally, export the analysis in json format
    with open(output_analysis_filename, 'w') as file:
        json.dump({'data': output_analysis}, file)

    # Finally remove the reduced trajectory since it is not required anymore
    if energies_trajectory_filename == 'energies.trajectory.xtc':
        logs = run([
            "rm",
            energies_trajectory_filename,
        ], stdout=PIPE).stdout.decode()
