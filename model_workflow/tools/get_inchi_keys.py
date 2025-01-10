from typing import Optional
import MDAnalysis
from rdkit import Chem
from model_workflow.utils.structures import Structure
import requests


def get_inchi_keys (
    u : 'MDAnalysis.Universe',
    structure : 'Structure'
) -> dict:
    """
    Generate a dictionary mapping InChI keys to residue information for non-standard residues.

    This function uses MDAnalysis to parse the input structure and topology files and identifies
    residues that are not classified as 'ion', 'solvent', 'nucleic', or 'protein'. For each
    identified residue, it converts the structure to RDKit format to obtain the InChI key
    and InChI string. The resulting data is stored in dictionaries to map InChI keys to residue
    details and residue names to InChI keys. PDB coordinates are necesary to distinguish stereoisomers.

    Args:
        input_structure_file (File): The input structure file.
        input_topology_file (Optional[File]): The input topology file.
        structure (Structure): The Structure object containing residues.

    Returns:
        dict: A dictionary where keys are InChI keys and values are dictionaries containing
            residue information such as associated residues, InChI strings, bond information,
            and classification.
    """
    key_2_name = {} # To see if different name for same residue
    name_2_key = {} # To see if different residues under same name
    # TO-DO: multiprocessing
    for residue in structure.residues:
        # Skip residues that are aminoacids, nucleics, or too small
        if residue.classification in ['ion', 'solvent', 'nucleic', 'protein']:
            continue
        # Select residues atoms with MDAnalysis
        res_sele = residue.get_selection().to_mdanalysis()
        res_mda = u.select_atoms(res_sele)
        # Convert to RDKIT and get InChi data
        res_RD = res_mda.atoms.convert_to("RDKIT") # slow step, 50% of time 
        inchikey = Chem.MolToInchiKey(res_RD)
        # If key don't existe we create the default entry with info that is only put once
        if not key_2_name.get(inchikey, None):
            inchi = Chem.MolToInchi(res_RD)
            has_bonds = residue.get_bonded_atoms()
            key_2_name[inchikey] = {'residues':[], 
                                    'inchi': inchi, 
                                    'has_bonds': has_bonds, 
                                    'classification': residue.classification}
        # We append all the residue with this key
        key_2_name[inchikey]['residues'].append(residue)

        # This dict is only to check if a residue has multiple keys:
        # Incorrect residue name, estereoisomers, lose of atoms...
        if not name_2_key.get(residue.name, None):
            name_2_key[residue.name] = []
        name_2_key[residue.name].append(inchikey)
    
    # DATA CHECKS
    # Check if there are multiple names for the same InChi key
    for inchikey, data in key_2_name.items():
        unique_names = list(set((res.name for res in data['residues'])))
        if data["has_bonds"]:
            print(f'WARNING: inChiKey {inchikey} is a substructure of type {data["classification"]}'\
                  ' and will result in inprecise search in PubChemb due lo lack of hidrogens and different charges on the bonded atoms.')
        if len(unique_names) > 1:
            print('WARNING: Same residue with different names: ' + inchikey + ' -> ' + str(unique_names))
        data['resname'] = list(unique_names) # Normally there is only one name

    # Check if there are multiple InChi keys for the same name
    for name, inchikeys in name_2_key.items():
        inchikeys = set(inchikeys)
        if len(inchikeys) > 1:
            key_counts = '  '.join([f'({len(key_2_name[key]["residues"])}): {key}' for key in inchikeys])
            print('WARNING: Same residue with different InChi keys: ' + name + ' -> ' + key_counts)
   
    return key_2_name

def is_in_LIPID_MAPS(inchikey) -> dict:
    """Search the InChi keys in LIPID MAPS"""
    headers = {'accept': 'json'}
    url = f"https://www.lipidmaps.org/rest/compound/inchi_key/{inchikey}/all"
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        js = response.text
        if js != '[]':
           return js        
    else:
        print(f"Error for {inchikey}: {response.status_code}")

def get_atoms_info(atom):
    """Get the RDkit atom information."""
    info = {
        'índice': atom.GetIdx(),
        'símbolo': atom.GetSymbol(),
        'carga_formal': atom.GetFormalCharge(),
        'carga_parcial': atom.GetDoubleProp('_GasteigerCharge') if atom.HasProp('_GasteigerCharge') else 'N/A',
        'valencia': atom.GetTotalValence(),
        'grado': atom.GetDegree(),
        'hidrógenos_explícitos': atom.GetNumExplicitHs(),
        'hidrógenos_implícitos': atom.GetNumImplicitHs(),
        'aromático': atom.GetIsAromatic()
    }
    return info


def inchi_2_mol(inchi, withHs=False, neutralize=False):
    """Converts an InChI string to a RDKIT. Usefull to display them."""
    mol = Chem.MolFromInchi(inchi)

    if neutralize:
        # Crear una copia editable de la molecula
        rwmol = Chem.RWMol(mol)
        for atom in mol.GetAtoms():
            charge = atom.GetFormalCharge()
            if charge:
                idx = atom.GetIdx()
                Hs = atom.GetNumExplicitHs()
                
                rw_atom = rwmol.GetAtomWithIdx(idx)
                rw_atom.SetFormalCharge(0)
                rw_atom.SetNumExplicitHs(Hs-charge)
        mol = rwmol

    if withHs:
        return Chem.AddHs(mol)
    else:
        return mol

def Residue_2_Mol(res):
    mol = Chem.RWMol()
    # map index in universe to index in mol
    atom_mapper = {}
    bonds = []
    non_res_bonds = res.get_bonded_atom_indices()
    for i, atom in enumerate(res.atoms):
        # create atom
        rdatom = Chem.Atom(atom.element)
        # enable/disable adding implicit H to the molecule
        rdatom.SetNoImplicit(True)
        for bond_to in atom.bonds:
            if bond_to in non_res_bonds:
                continue
            bond = sorted([atom.index, bond_to])
            if bond not in bonds:
                bonds.append(bond)
        # add atom
        index = mol.AddAtom(rdatom)
        atom_mapper[atom.index] = index
    for bond in bonds:
        mol.AddBond(bond[0], atom_mapper[bond[1]], Chem.BondType.SINGLE)

    # sanitize if possible
    err = Chem.SanitizeMol(mol, catchErrors=True)
    # infer bond orders and formal charges
    MDAnalysis.converters.RDKit._infer_bo_and_charges(mol)
    return mol
