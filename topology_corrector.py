import prody
from collections import Counter

# Analyze the topology looking for irregularities and then modify the topology to standarize the format
# Both analysis and modifications are carried by Prody
# 
# Supported cases:
# 
# * Repeated chains -> Chains are renamed (e.g. A, G, B, G, C, G -> A, G, B, H, C, I)
# * (Not supported yet) Hyrdogens are placed at the end ->
# * Repeated residues -> Residues are renumerated (e.g. 1, 2, 3, 1, 2, 1, 2 -> 1, 2, 3, 4, 5, 6, 7)
# * Repeated atoms -> Atoms are renamed with their numeration (e.g. C, C, C, O, O -> C1, C2, C3, O1, O2)
def topology_corrector (
    input_topology_filename : str,
    output_topology_filename : str ):

    # Store some messages for the final logs
    logs = []

    # Track if there has been any modification and then topology must be rewritten
    modified = False

    # Import the topology to prody
    test = prody.parsePDB(input_topology_filename)

    # ------------------------------------------------------------------------------------------
    # Repeated chains --------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # Save the string letters of each chain
    # Find repeated chains which would give problems further if topology is not modified:
    chain_letters = []
    last_chain = None
    repeated_chains = []
    last_repeated_chain = None

    residues = test.iterResidues()
    for r in residues:
        chain = r.getChid()
        if chain == last_chain:
            continue
        if chain == last_repeated_chain:
            r.getChain().setChid(last_chain)
            continue
        # When a new chain is found
        last_repeated_chain = None
        # If it is repeated modify the chain letter
        if chain in chain_letters:
            repeated_chains.append(chain)
            last_repeated_chain = chain
            # Find the next alphabetic letter until we find a not repeated chain letter to define a new chain
            last_chain = chain
            while last_chain in chain_letters:
                last_chain = get_next_letter(last_chain)
            # Set the new chain
            chain_letters.append(last_chain)
            r.getChain().setChid(last_chain)
            logs.append('- One of the repeated chains "' + last_repeated_chain + '" has been renamed as "' + last_chain + '"')
        # If is not repeated add it to the list
        else:
            last_chain = chain
            chain_letters.append(last_chain)

    # Check if were there repetaed chains
    reps = len(repeated_chains)
    if reps > 0:
        modified = True
        print('WARNING! There are ' + str(reps) + ' repated chains (e.g. ' + repeated_chains[0] + ')')

    # ------------------------------------------------------------------------------------------
    # Repeated residues ------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # Save the residue identifieres, which are strings composed of each chain letter, residue number and icode
    # This ids are used to find repeated residues which would give problems further if topology is not modified:
    # Pytraj would consider repeated residues as a single residue
    # The web client would have no way to refer a residue by its numeration (e.g. 'C:17A')
    residue_ids = []
    repeated_residues = []

    residues = test.iterResidues()
    for r in residues:
        chain = r.getChid()
        num = r.getResnum()
        icode = r.getIcode()
        rid = chain + str(num) + icode
        if rid in residue_ids:
            repeated_residues.append(rid)
            while rid in residue_ids:
                num += 1
                rid = chain + str(num) + icode
            r.setResnum(num)
        residue_ids.append(rid)

    # Check if were there repetaed residues
    reps = len(repeated_residues)
    if reps > 0:
        modified = True
        logs.append('- Some repeated residues have been re-numerated')
        print('WARNING! There are ' + str(reps) + ' repated residues (e.g. ' + repeated_residues[0] + ')')
        
    # ------------------------------------------------------------------------------------------
    # Repeated atoms --------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # Atoms with identical chain, residue, icode and name are renamed with an ordinal number
    # e.g.
    # ATOM      1  C   UNK L   0   -->   ATOM      1  C1  UNK L   0
    # ATOM      1  C   UNK L   0   -->   ATOM      1  C2  UNK L   0
    # ATOM      1  C   UNK L   0   -->   ATOM      1  C3  UNK L   0

    # Iterate over all atoms and change all their names
    atoms = test.iterAtoms()
    reporter = 0 # Just count how many atoms have been renamed for the report
    counter = Counter()
    current = None
    for atom in atoms:
        # Find the chain, residue number and icode where this atom belongs to
        chainLetter = atom.getChid()
        residueNumber = atom.getResnum()
        icode = atom.getIcode()
        label = chainLetter + str(residueNumber) + icode
        # If this is a new resiude reset the counter
        if label != current:
            counter = Counter()
            current = label
        # We set the name of the atom as the element letter + the count of this element
        element = atom.getElement()
        counter[element] += 1
        reporter += 1
        atom.setName(element + str(counter[element]))

    # Check if were there repetaed residues
    reps = reporter
    if reps > 0:
        modified = True
        logs.append('- Some repeated atoms have been renamed')
        print('WARNING! There are ' + str(reps) + ' repated atoms')

    # ------------------------------------------------------------------------------------------
    # Final output -----------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------

    # Write a new topology if any modification was done
    if modified:
        # Print all stored logs
        print('The topology file has been modificated:')
        for l in logs:
            print(l)
        prody.writePDB(output_topology_filename, test)


# Set a function to get the next letter from an input letter in alphabetic order
def get_next_letter (letter : str):
    if letter == 'z' or letter == 'Z':
        raise SystemExit("Limit of chain letters has been reached")
    next_letter = chr(ord(letter) + 1)
    return next_letter