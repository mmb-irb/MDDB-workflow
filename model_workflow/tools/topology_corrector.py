from pathlib import Path
from prody import parsePDB, writePDB, calcDistance
from subprocess import run, PIPE
from collections import Counter

import itertools

# Import local tools
from model_workflow.tools.vmd_processor import vmd_chainer

# Set the path to the 'addChains' script
repo_path = str(Path(__file__).parent.parent)
add_chains_script = repo_path + '/utils/resources/addChainCV19.pl'

# Analyze the topology looking for irregularities and then modify the topology to standarize the format
# Both analysis and modifications are carried by Prody
#
# Supported cases:
#
# * Missing chains -> Chains are added through VMD
# * Repeated chains -> Chains are renamed (e.g. A, G, B, G, C, G -> A, G, B, H, C, I)
# * (Not supported yet) Hyrdogens are placed at the end ->
# * Repeated residues -> Residues are renumerated (e.g. 1, 2, 3, 1, 2, 1, 2 -> 1, 2, 3, 4, 5, 6, 7)
# * Repeated atoms -> Atoms are renamed with their numeration (e.g. C, C, C, O, O -> C1, C2, C3, O1, O2)


def topology_corrector(
    input_topology_filename: str,
    output_topology_filename: str):

    print('Correcting topology')

    # Store some messages for the final logs
    logs = []

    # Track if there has been any modification and then topology must be rewritten
    modified = False

    # Remove all bytes which can not be decoded by utf-8 codec
    purgeNonUtf8(input_topology_filename)

    # Import the topology to prody
    test = parsePDB(input_topology_filename)

    # Set a function to get a chain residue number not used yet in topology which is after the last (highest) used number
    def get_next_free_number(chain : str) -> int:
        numbers = []
        for residue in test.iterResidues():
            if residue.getChid() != chain:
                continue
            numbers.append(residue.getResnum())
        return max(numbers) + 1

    # Set a function to find out if an atom is close enought to any other atom in a list
    # Distance limit is 3 Ångstroms (Å)
    def is_close(atom, other_atoms : list, distance : float = 3) -> bool:
        for other_atom in other_atoms:
            if calcDistance(other_atom, atom) < 3:
                return True
        return False

    # ------------------------------------------------------------------------------------------
    # Missing chains --------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # Check if chains are missing. If so, create a new chainned topology and set it as the reference
    chains = list(test.iterChains())
    no_chains = len(chains) == 1 and ( chains[0].getChid() == ' ' or chains[0].getChid() == 'X' )
    if no_chains:
        modified = True
        print('WARNING: chains are missing')

        # Use VMD to set chains according to fragments
        vmd_chainer(input_topology_filename, input_topology_filename)

        # Reload the prody topology and chains
        test = parsePDB(input_topology_filename)
        chains = list(test.iterChains())

        # This system is deprecated. Use vmd processor
        #chainned_topology = 'chainned_topology.pdb'
        #run([
        #     "perl",
        #     add_chains_script,
        #     input_topology_filename,
        #     chainned_topology,
        # ], stdout=PIPE).stdout.decode()
        # test = parsePDB(chainned_topology)
        # run([
        #     "rm",
        #     chainned_topology,
        # ], stdout=PIPE).stdout.decode()
        # logs.append('- Chains have been asigned automatically')

    # ------------------------------------------------------------------------------------------
    # Repeated chains --------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # Save the string letters of each chain
    # Find repeated chains which would give problems further if topology is not modified:
    chain_letters = []
    last_chain_letter = None
    repeated_chain_letters = []
    last_repeated_chain_letter = None

    # DANI: No se ha provado con iterAtoms
    for a in test.iterAtoms():
        chain_letter = a.getChid()
        if chain_letter == last_chain_letter:
            continue
        if chain_letter == last_repeated_chain_letter:
            chain = chains[a.getChindex()]
            chain.setChid(last_chain_letter)
            continue
        # When a new chain letter is found
        last_repeated_chain_letter = None
        # If it is repeated modify the chain letter
        if chain_letter in chain_letters:
            repeated_chain_letters.append(chain_letter)
            last_repeated_chain_letter = chain_letter
            # Find the next alphabetic letter until we find a not repeated chain letter to define a new chain
            last_chain_letter = chain_letter
            while last_chain_letter in chain_letters:
                last_chain_letter = get_next_letter(last_chain_letter)
            # Set the new chain letter
            chain_letters.append(last_chain_letter)
            print('CHAINS')
            print(chains)
            print('INDEX: ' + str(a.getChindex()))
            chain = chains[a.getChindex()]
            chain.setChid(last_chain_letter)
            logs.append('- One of the repeated chains "' + last_repeated_chain_letter +
                        '" has been renamed as "' + last_chain_letter + '"')
        # If is not repeated add it to the list
        else:
            last_chain_letter = chain_letter
            chain_letters.append(last_chain_letter)

    # Check if were there repetaed chains
    reps = len(repeated_chain_letters)
    if reps > 0:
        modified = True
        print('WARNING! There are ' + str(reps) +
              ' repated chains (e.g. ' + repeated_chain_letters[0] + ')')

    # ------------------------------------------------------------------------------------------
    # Repeated residues ------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # Iter atom by atom recording each atom chain, residue number and icode
    # In case non-contiguous atoms have identical recorded values check distance between those atoms
    # If atoms are close it means they are indeed the same residue
    # In this case, we do not modify anything
    # If atoms are far it means they are different residues with the same notation
    # In this case, renumerate the second residue by raising the number

    # Splitted residues are found in some pdbs and they are supported by some tools
    # These tools consider all atoms with the same 'record' as the same residue
    # However, there are other tools which would consider the splitted residue as two different resdiues
    # This workflow has some fixes to give partiall support for splitted residues but it is not ideal
    # This situation has no easy fix since we can not change the order in atoms because trajectory data depends on it

    # Duplicated residues are usual and they are always supported but with some problems
    # For example, pytraj analysis outputs use to sort results by residues and each residue is tagged
    # If there are duplicated residues with the same tag it may be not possible to know which result belongs to each residue
    # Another example are NGL selections once in the web client
    # If you select residue ':A and 1' and there are multiple residues 1 in chain A all of them will be displayed
    
    # Residues are recorded as a tuple with format (chain, number, icode)
    records = []
    # Residue atoms of each record are saved
    record_atoms = []

    # Get names of splitted and repeated residues only for display
    splitted_residues = []
    repeated_residues = []
    
    last_record = None
    last_repeated_record = None
    for atom in test.iterAtoms():
        chain = atom.getChid()
        resnum = atom.getResnum()
        icode = atom.getIcode()
        record = chain, resnum, icode
        # If atom is in the same residue that the previous atom
        if record == last_record:
            record_atoms[-1].append(atom)
            continue
        # If atom is the first atom of a residue which has not been recorded yet
        if record not in records:
            last_record = record
            records.append(record)
            record_atoms.append([atom])
            continue
        # If the residue has been previously recorded
        index = records.index(record)
        # If atom is in the same previously recorded residue that the previous atom
        if record == last_repeated_record:
            previous_corrected_number = last_record[1]
            atom.setResnum(previous_corrected_number)
            record_atoms[index].append(atom)
            continue
        last_repeated_record = record
        recorded_atoms = record_atoms[index]
        tag = chain + ':' + str(resnum) + icode
        # If atom is close enought to the previously recorded residue (i.e. it is a splitted residue)
        if is_close(atom, recorded_atoms):
            splitted_residues.append(tag)
            last_record = record
            record_atoms[index].append(atom)
            continue
        # If atom belongs to a different residue which is duplicated
        repeated_residues.append(tag)
        new_corrected_number = get_next_free_number(chain)
        atom.setResnum(new_corrected_number)
        record = chain, new_corrected_number, icode
        last_record = record
        records.append(record)
        record_atoms.append([atom])

    # Warn the user in case we have repeated or splitted residues
    if len(splitted_residues) > 0:
        print('WARNING: There are ' + str(len(splitted_residues)) + ' splitted residues')
        print(splitted_residues)
    if len(repeated_residues) > 0:
        modified = True
        logs.append('- Some repeated residues have been re-numerated')
        print('WARNING: There are ' + str(len(repeated_residues)) + ' repeated residues')
        print(repeated_residues)

    # ------------------------------------------------------------------------------------------

    # Save the residue identifieres, which are strings composed of each chain letter, residue number and icode
    # This ids are used to find repeated residues which would give problems further if topology is not modified:
    # Pytraj would consider repeated residues as a single residue
    # The web client would have no way to refer a residue by its numeration (e.g. 'C:17A')
    # DEPRECTAED: This system relies on prody's iterResidues which may fail itself when counting residues

    # residue_ids = []
    # repeated_residues = []

    # residues = test.iterResidues()
    # for r in residues:
    #     chain = r.getChid()
    #     num = r.getResnum()
    #     icode = r.getIcode()
    #     rid = chain + str(num) + icode
    #     if rid in residue_ids:
    #         repeated_residues.append(rid)
    #         while rid in residue_ids:
    #             num += 1
    #             rid = chain + str(num) + icode
    #         r.setResnum(num)
    #     residue_ids.append(rid)

    # # Check if were there repetaed residues
    # reps = len(repeated_residues)
    # if reps > 0:
    #     modified = True
    #     logs.append('- Some repeated residues have been re-numerated')
    #     print('WARNING! There are ' + str(reps) +
    #           ' repated residues (e.g. ' + repeated_residues[0] + ')')

    # ------------------------------------------------------------------------------------------
    # Repeated atoms --------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # Atoms with identical chain, residue, icode and name are renamed with an ordinal number
    # e.g.
    # ATOM      1  C   UNK L   0   -->   ATOM      1  C1  UNK L   0
    # ATOM      1  C   UNK L   0   -->   ATOM      1  C2  UNK L   0
    # ATOM      1  C   UNK L   0   -->   ATOM      1  C3  UNK L   0

    reporter = 0  # Just count how many atoms have been renamed for the output logs
    example = None # Just save an example of a renamed atom for the output logs
    residues = test.iterResidues()
    for r in residues:
        # Iterate over the residue atoms and change all their names
        counter = Counter()
        atoms = r.iterAtoms()
        for atom in atoms:
            # We check if the atom name already exists. If not, go to the next atom
            name = atom.getName()
            if counter[name] == 0:
                counter[name] += 1
                continue
            reporter += 1
            # We set the name of the atom as the element letter + the count of this element
            # If element is missing take the first character of the atom name
            initial = atom.getElement()
            if not initial or initial == ' ':
                initial = name[0]
            number = counter[initial] + 1
            new_name = initial + str(number)
            while counter[new_name] > 0:
                number += 1
                new_name = initial + str(number)
            counter[new_name] += 1
            if not example:
                example = '(' + str(atom.getIndex() + 1) + ') ' + name + ' -> ' + new_name
            atom.setName(new_name)
            #print(str(atom.getIndex() + 1) + ' / ' + r.getResname() + ': ' + name + ' -> ' + new_name)

    # Check if were there repetaed residues
    if reporter > 0:
        modified = True
        logs.append('- Some repeated atoms have been renamed')
        print('WARNING! There are ' + str(reporter) + ' repated atoms')
        print('e.g. ' + example)

    # ------------------------------------------------------------------------------------------
    # Final output -----------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # Write a new topology if any modification was done
    if modified:
        # Print all stored logs
        print('The topology file has been modificated:')
        for log in logs:
            print(log)
        writePDB(output_topology_filename, test)


# Set a function to get the next letter from an input letter in alphabetic order
def get_next_letter(letter: str) -> str:
    if letter == 'z' or letter == 'Z':
        raise SystemExit("Limit of chain letters has been reached")
    next_letter = chr(ord(letter) + 1)
    return next_letter

# Remove all bytes which can not be decoded by utf-8 codec
# This prevents the prody parsePDB function to return the following error:
# UnicodeDecodeError: 'utf-8' codec can't decode byte ...
def purgeNonUtf8 (filename : str):
    with open(filename, mode="r+", encoding="utf-8", errors= 'ignore') as file:
        lines = file.readlines()
        file.seek(0)
        for line in lines:
            file.write(line)
        file.truncate()