import pytraj as pt

# Given a topology which includes charges
# Extract those charges and save them in a list to be returned
# This list will include also the atom element in a CMIP friendly format
# Supported formats (tested): prmtop, psf (standard psf, not from DESRES)
# Supported formats (not tested): mol2, cif, tpr, sdf
# Non supported formats: mae, pdb (has no charges)
def get_topology_charges (input_topology_filename : str):
    
    topology = pt.load_topology(filename=input_topology_filename)
    # WARNING: We must convert this numpy ndarray to a normal list
    # Otherwise the search by index is extremly ineficient
    topology_charges = list(topology.charge)

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
    }
    # Iterate over each atom to save their charge and CMIP element
    charges = []
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
            # DANI: Error aquí
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
        charge = topology_charges[a]
        charges.append(charge)
    return charges, elements
            

# Given a topology which includes charges
# Extract those charges and write them in a new file with a specific format
# Supported formats (tested): prmtop, psf (standard psf, not from DESRES)
# Supported formats (not tested): mol2, cif, tpr, sdf
# Non supported formats: mae, pdb (has no charges)
def get_topology_charges_reslib (
    input_topology_filename : str,
    output_charges_filename : str = 'charges.txt'):
    
    topology = pt.load_topology(filename=input_topology_filename)

    # Save the name of already found unknown residues so we can skip them when repeated
    already_parsed_atoms = []
    # Set all supported elements
    # This is required to transform element names (returned by pytraj) to element letters
    elements = {
        'hydrogen': 'H',
        'carbon': 'C',
        'oxygen': 'O',
        'nitrogen': 'N',
        'sulfur': 'S',
        'sodium': 'Na',
        'chlorine': 'Cl',
    }
    with open(output_charges_filename, 'w') as file:
        atoms = list(topology.atoms)
        for a, atom in enumerate(atoms):
            residue = atom.resname
            # Skip this atom if we already found it
            name = atom.name
            if next((at for at in already_parsed_atoms if at['residue'] == residue and at['name'] == name), None):
                continue
            element = elements[atom.element]
            # Make sure the first element letter is upper
            # If the element has 2 letters, make sure the second letter is lower
            first_letter = element[0].upper()
            second_letter = ''
            if len(element) == 2:
                second_letter = element[1].lower()
            element = first_letter + second_letter
            # Adapt hydrogens element to CMIP requirements
            if element == 'H':
                # There should we always only 1 bond
                # DANI: Error aquí
                # ValueError: Buffer dtype mismatch, expected 'int' but got 'long'
                bonded_heavy_atom_index = atom.bonded_indices()[0]
                bonded_heavy_atom = atoms[bonded_heavy_atom_index]
                bonded_heavy_atom_element = elements[bonded_heavy_atom.element]
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
            charge = topology.charge[a]
            # Stringify the charge
            # Then add a white space at the start if it is not negative
            # Then cut the string to 7 characters
            str_charge = str(charge)
            if charge >= 0:
                str_charge = ' ' + str_charge
            str_charge = str_charge[0:8]
            reslib_line = residue.ljust(5) + name.ljust(5) + element.ljust(6) + str_charge
            already_parsed_atoms.append({ 'residue': residue, 'name': name })
            file.write(reslib_line)