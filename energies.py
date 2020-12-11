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
from pathlib import Path
import prody
import re
import numpy
import math
from subprocess import run, PIPE, Popen

import json

# Set the path to auxiliar files required for this analysis
repo_path = str(Path(__file__).parent)
source_reslib = repo_path + '/aux/res.lib'
preppdb_source = repo_path + '/aux/preppdb.pl'
check_source = repo_path + '/aux/check.in'
test_source = repo_path + '/aux/test.in'
vdw_source = repo_path + '/aux/vdwprm'

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
    return center, density, units

def compute_new_grid(
        prot_center,
        prot_density,
        prot_units,
        lig_center,
        lig_density,
        lig_units):
    new_center = [(i + j) / 2 for i, j in zip(prot_center, lig_center)]
    new_density = []
    for k in range(3):
        min_prot = prot_center[k] - prot_density[k] * prot_units[k]
        min_lig = lig_center[k] - lig_density[k] * lig_units[k]
        min_new = min(min_prot, min_lig)
        max_prot = prot_center[k] + prot_density[k] * prot_units[k]
        max_lig = lig_center[k] + lig_density[k] * lig_units[k]
        max_new = max(max_prot, max_lig)
        dnew = int(abs(max_new - min_new))
        new_density.append(dnew)
    return new_center, new_density


def write_CMIP_input(test_input_file, weighted_center, weighted_density):
    test_in_file = [
        "Electrostatic + VdW Interaction energy calculation. \n",
        "&cntrl \n",
        " tipcalc=3,irest=0 \n",
        " ebyatom = 1 \n",
        " readgrid=0 \n",
        f" cenx={weighted_center[0]:.1f} \n",
        f" ceny={weighted_center[1]:.1f} \n",
        f" cenz={weighted_center[2]:.1f} \n",
        f" dimx={int(weighted_density[0])} \n",
        f" dimy={int(weighted_density[1])} \n",
        f" dimz={int(weighted_density[2])} \n",
        " intz=0.5 \n",
        " intz=0.5 \n",
        " intz=0.5 \n",
        " coorfmt=2 \n",
        " fvdw=0.8 \n",
        "&end \n"]
    with open(test_input_file, "w") as f:
        for line in test_in_file:
            f.write(line)


# Perform the electrostatic and vdw energies analysis for each ligand
# DANI: En principio soporta casos en que hay multiples ligandos, pero no se ha provado
def energies (
    input_topology_filename : str,
    input_trajectory_filename : str,
    output_analysis_filename : str,
    reference,
    snapshots : int,
    ligands : dict ):

    # Change elements in a prody topology to meet the CMIP requirements
    # Hydrogens bonded to carbons remain as 'H'
    # Hydrogens bonded to oxygen are renamed as 'HO'
    # Hydrogens bonded to nitrogen or sulfur are renamed as 'HN'
    def correct_orphan_atoms (prody_topology, prody_selection):
        selection = prody_topology.select(prody_selection)
        # Get all atoms and atom coordinates
        selection_coords = selection.getCoords()
        selection_atoms = list(selection.iterAtoms())
        # Set functions to find the closest atom of a specified atom
        def get_distance (coords1, coords2):
            squared_distance = numpy.sum((coords1 - coords2)**2, axis=0)
            distance = numpy.sqrt(squared_distance)
            return distance
        def find_closest_atom (atom):
            current_coords = atom.getCoords()
            distances = [get_distance(current_coords, c) for c in selection_coords]
            sorted_distances = [d for d in distances]
            sorted_distances.sort()
            # We take the second minimum, since the first minimum will be always 0
            smallest_distance = sorted_distances[1]
            index = distances.index(smallest_distance)
            closest_atom = selection_atoms[index]
            return closest_atom
        # Update the element of each hydrogen according to CMIP needs
        for atom in selection_atoms:
            if atom.getName() == 'H':
                bonded_heavy_atom = find_closest_atom(atom).getElement()
                # Hydrogens bonded to carbons remain as 'H'
                if bonded_heavy_atom == 'C':
                    continue
                # Hydrogens bonded to oxygen are renamed as 'HO'
                if bonded_heavy_atom == 'O':
                    atom.setElement('HO')
                # Hydrogens bonded to nitrogen or sulfur are renamed as 'HN'
                if bonded_heavy_atom == 'N' or bonded_heavy_atom == 'S':
                    atom.setElement('HN')
                if atom.getName() == 'CL':
                    atom.setElement('Cl')
                if atom.getName() == 'BR':
                    atom.setElement('Br')
        # Return the corrected prody topology
        return prody_topology

    # Given a pdb structure, use CMIP to extract energies
    # Output energies are already added by residues
    def get_frames_energy (frame_pdb):
        # Parse the pdb file to prody format and the correct it
        original_topology = prody.parsePDB(frame_pdb)
        # Add chains according to the reference topology, since gromacs has deleted chains
        original_topology.setChids(reference.topology.getChids())
        # Correct posible 'orphan' hydrogens in each ligand by renaming them
        # Hydrogens bonded to oxygen or nitrogen or sulfur must be renamed as HO or HN
        # CMIP uses atom names to set the atom radius so this step is important
        for ligand in ligands:
            original_topology = correct_orphan_atoms(original_topology, ligand['prody_selection'])
        # WARNING: At this point topology should be corrected
        # WARNING: Repeated atoms will make the analysis fail

        # Prepare the protein files
        # Create a new pdb only with the protein
        protein_pdb = 'protein.pdb'
        selection = original_topology.select('protein')
        prody.writePDB(protein_pdb, selection)

        # Get the 'cmip' input file
        protein_cmip = 'protein.cmip.pdb'
        run([
            "perl",
            preppdb_source,
            source_reslib,
            protein_pdb,
            protein_cmip,
        ], stdout=PIPE).stdout.decode()

        # Prepare the ligand files for each ligand
        ligand_data = []
        for ligand in ligands:

            # Create a new pdb only with the current ligand
            name = ligand['name']
            ligand_pdb = name + '.pdb'
            selection = original_topology.select(ligand['prody_selection'])
            prody.writePDB(ligand_pdb, selection)

            # Calculate the energies with acpype and then mine them
            #test = !obabel -ipdb $ligand_pdb -omol2 whatever
            acpype_logs = run([
                "acpype",
                "-i",
                ligand_pdb,
                #"-n",
                #"0",
            ], stdout=PIPE).stdout.decode()

            energies_file = name + '.acpype/' + name + '.mol2'
            # Write all lines but the last line: 'END'
            with open(energies_file, "r") as file:
                energies_lines = file.readlines()

            energies = []
            zero_count = 0
            for line in energies_lines:
                line = line.split()
                if len(line) == 9:
                    energy_value = line[-1]
                    if float(energy_value) == 0:
                        zero_count += 1
                    energies.append(energy_value)
            # print(energies)

            if zero_count == len(energies):
                raise SystemExit(
                    "All energies are zero! Check files and try again.")

            # Check the number of atoms matches the number of energies
            if len(list(selection.iterAtoms())) != len(energies):
                raise SystemExit("Stop!! The number of atoms and energies does not match :/")

            # Copy the 'res.lib' file in the local path and open it to 'a'ppend new text
            reslib_filename = name + '_res.lib'
            #!cp $source_reslib ./$reslib_filename
            run([
                "cp",
                source_reslib,
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

            # Get the cmip input file
            ligand_cmip = name + '.cmip.pdb'
            run([
                "perl",
                preppdb_source,
                reslib_filename,
                ligand_pdb,
                ligand_cmip,
            ], stdout=PIPE).stdout.decode()

            cmip_output = 'protein-' + name + '.energy.pdb'

            # check host and guest dimensions to build grid
            cmip_logs_host = run([
                "cmip",
                "-i",
                check_source,
                "-pr",
                ligand_cmip,
                "-vdw",
                vdw_source,
                "-hs",
                protein_cmip,
                "-byat",
                cmip_output,
            ], stdout=PIPE).stdout.decode()
            prot_center, prot_density, prot_units = mine_cmip_output(
                cmip_logs_host.split("\n"))

            cmip_logs_guest = run([
                "cmip",
                "-i",
                check_source,
                "-pr",
                protein_cmip,
                "-vdw",
                vdw_source,
                "-hs",
                ligand_cmip,
                "-byat",
                cmip_output,
            ], stdout=PIPE).stdout.decode()
            lig_center, lig_density, lig_units = mine_cmip_output(
                cmip_logs_guest.split("\n"))

            new_center, new_density = compute_new_grid(
                prot_center,
                prot_density,
                prot_units,
                lig_center,
                lig_density,
                lig_units)

            test_input_file = repo_path + '/aux/input.in'
            write_CMIP_input(
                test_input_file,
                new_center,
                new_density)

            # Run the cmip software to get the desired energies
            cmip_output = 'protein-' + name + '.energy.pdb'
            cmip_logs = run([
                "cmip",
                "-i",
                test_input_file,
                # test_source,
                "-pr",
                protein_cmip,
                "-vdw",
                vdw_source,
                "-hs",
                ligand_cmip,
                "-byat",
                cmip_output,
            ], stdout=PIPE).stdout.decode()
            #print(cmip_logs)

            # Mine the electrostatic (es) and Van der Walls (vdw) energies for each atom
            # Group the results by reidues adding their values
            residues = {}
            with open(cmip_output,'r') as file:
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
                    energies =  (vdw, es, both) if both < 100 else (0, 0, 0)
                    if residue in residues:
                        residues[residue] = tuple(
                            [a+b for a,b in zip(energies, residues[residue])])
                    else:
                        residues[residue] = energies
            
            ligand_data.append(residues)

        return ligand_data
        

    # Set the number of frames where we extract energies to calculate the average
    frames_number = 100
    frames = None
    if snapshots > frames_number:
        frames = [f * math.floor(snapshots / frames_number) for f in range(1, frames_number + 1)]
    else:
        frames = range(1, snapshots + 1)

    # Extract the energies for each frame
    #data = []
    ligands_data = [[] for l in ligands]
    for f in frames:
        # Extract the current frame
        current_frame = 'frame' + str(f) + '.pdb'
        # The frame selection input in gromacs works with a 'ndx' file
        frames_ndx = 'frames.ndx'
        with open(frames_ndx,'w') as file:
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
            input_trajectory_filename,
            '-o',
            current_frame,
            "-fr",
            frames_ndx,
            '-quiet'
        ], stdin=p.stdout, stdout=PIPE).stdout.decode()
        p.stdout.close()
        # Run the main analysis over the current frame
        # Append the result data for each ligand
        energies_data = get_frames_energy(current_frame)
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

        residues_vdw_avg_initial = [sum(v) / len(v) for v in residues_vdw_values_initial]
        residues_es_avg_initial = [sum(v) / len(v) for v in residues_es_values_initial]
        residues_both_avg_initial = [sum(v) / len(v) for v in residues_both_values_initial]

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

        residues_vdw_avg_final = [sum(v) / len(v) for v in residues_vdw_values_final]
        residues_es_avg_final = [sum(v) / len(v) for v in residues_es_values_final]
        residues_both_avg_final = [sum(v) / len(v) for v in residues_both_values_final]

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
        json.dump({ 'data': output_analysis }, file)