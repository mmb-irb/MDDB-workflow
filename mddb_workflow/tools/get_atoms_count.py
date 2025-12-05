# Count different type of atoms and residues in the structure
def get_atoms_count (structure : 'Structure') -> tuple:

    # Number of system atoms and residues
    system_atoms = len(structure.atoms)
    system_residues = len(structure.residues)
    # Number of protein atoms and residues
    protein_selection = structure.select_protein()
    protein_atoms = len(protein_selection)
    protein_residues = len(structure.get_selection_residue_indices(protein_selection))
    # Number of DNA atoms and residues
    dna_selection = structure.select_by_classification('dna')
    dna_atoms = len(dna_selection)
    dna_residues = len(structure.get_selection_residue_indices(dna_selection))
    # Number of RNA atoms and residues
    rna_selection = structure.select_by_classification('rna')
    rna_atoms = len(rna_selection)
    rna_residues = len(structure.get_selection_residue_indices(rna_selection))
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
    # Number of ions
    counter_cations = len(structure.select_counter_ions(charge='+'))
    counter_anions = len(structure.select_counter_ions(charge='-'))
    counter_ions = counter_cations + counter_anions
    total_ions = len(structure.select_ions())
    non_counter_ions = total_ions - counter_ions
    # Other atoms which do not fail in any of the previous sections
    other_atoms = system_atoms - protein_atoms - dna_atoms - rna_atoms - lipid_atoms - carbohydrates_atoms - total_ions

    # Display a summary of atom and residue counts
    print('Atom and residue counts:')
    print(f' System atoms: {system_atoms}')
    print(f' System residues: {system_residues}')
    print(f' Protein atoms: {protein_atoms}')
    print(f' Protein residues: {protein_residues}')
    print(f' DNA atoms: {dna_atoms}')
    print(f' DNA residues: {dna_residues}')
    print(f' RNA atoms: {rna_atoms}')
    print(f' RNA residues: {rna_residues}')
    print(f' Lipid atoms: {lipid_atoms}')
    print(f' Lipid residues: {lipid_residues}')
    print(f' Carbohydrate atoms: {carbohydrates_atoms}')
    print(f' Carbohydrate residues: {carbohydrates_residues}')
    print(f' Solvent atoms: {solvent_atoms}')
    print(f' Solvent residues: {solvent_residues}')
    print(f' Counter cations: {counter_cations}')
    print(f' Counter anions: {counter_anions}')
    print(f' Counter ions: {counter_ions}')
    print(f' Non-counter ions: {non_counter_ions}')
    print(f' Other atoms: {other_atoms}')

    # Return the counts
    return (system_atoms, system_residues, protein_atoms, protein_residues,
        dna_atoms, dna_residues, rna_atoms, rna_residues, lipid_atoms, lipid_residues,
        carbohydrates_atoms, carbohydrates_residues, solvent_atoms, solvent_residues,
        counter_cations, counter_anions, counter_ions, non_counter_ions, other_atoms)
