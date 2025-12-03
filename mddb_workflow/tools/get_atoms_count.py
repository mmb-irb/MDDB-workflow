# Count different type of atoms and residues in the structure
def get_atoms_count (structure : 'Structure') -> tuple:

    # Number of system atoms and residues
    system_atoms = len(structure.atoms)
    system_residues = len(structure.residues)
    # Number of protein atoms and residues
    protein_selection = structure.select_protein()
    protein_atoms = len(protein_selection)
    protein_residues = len(structure.get_selection_residue_indices(protein_selection))
    # Number of nucleic atoms and residues
    nucleic_selection = structure.select_nucleic()
    nucleic_atoms = len(nucleic_selection)
    nucleic_residues = len(structure.get_selection_residue_indices(nucleic_selection))
    # Number of lipid atoms and residues
    lipids_selection = structure.select_lipids()
    lipid_atoms = len(lipids_selection)
    lipid_residues = len(structure.get_selection_residue_indices(lipids_selection))
    # Number of carbohydrates atoms and residues
    carbohydrates_selection = structure.select_carbohydrates()
    carbohydrates_atoms = len(carbohydrates_selection)
    carbohydrates_residues = len(structure.get_selection_residue_indices(carbohydrates_selection))
    # Number of solvent atoms and residues
    solvent_selection = structure.select_water()
    solvent_atoms = len(solvent_selection)
    solvent_residues = len(structure.get_selection_residue_indices(solvent_selection))
    # Number of counter ions
    counter_cations = len(structure.select_counter_ions(charge='+'))
    counter_anions = len(structure.select_counter_ions(charge='-'))
    counter_ions = counter_cations + counter_anions

    # Display a summary of atom and residue counts
    print('Atom and residue counts:')
    print(f' System atoms: {system_atoms}')
    print(f' System residues: {system_residues}')
    print(f' Protein atoms: {protein_atoms}')
    print(f' Protein residues: {protein_residues}')
    print(f' Nucleic atoms: {nucleic_atoms}')
    print(f' Nucleic residues: {nucleic_residues}')
    print(f' Lipid atoms: {lipid_atoms}')
    print(f' Lipid residues: {lipid_residues}')
    print(f' Carbohydrate atoms: {carbohydrates_atoms}')
    print(f' Carbohydrate residues: {carbohydrates_residues}')
    print(f' Solvent atoms: {solvent_atoms}')
    print(f' Solvent residues: {solvent_residues}')
    print(f' Counter cations: {counter_cations}')
    print(f' Counter anions: {counter_anions}')
    print(f' Counter ions: {counter_ions}')

    # Return the counts
    return (system_atoms, system_residues, protein_atoms, protein_residues,
        nucleic_atoms, nucleic_residues, lipid_atoms, lipid_residues,
        carbohydrates_atoms, carbohydrates_residues, solvent_atoms, solvent_residues,
        counter_cations, counter_anions, counter_ions)
