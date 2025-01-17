import MDAnalysis
from rdkit import Chem
from model_workflow.utils.structures import Structure
from model_workflow.utils.type_hints import *
from model_workflow.utils.warnings import warn
from functools import lru_cache
import requests


def process_residue(res_atoms: 'MDAnalysis.AtomGroup', 
                    resindex : int) -> Tuple[str, str, int]:
    """
    Process a single residue to get its InChI key and related information.
    """
    # Convert to RDKIT and get InChI data
    res_RD = res_atoms.convert_to("RDKIT")
    # Calculate InChI key and string
    inchikey = Chem.MolToInchiKey(res_RD)  # slow step, 50% of time 
    inchi = Chem.MolToInchi(res_RD)
    return (inchikey, inchi, resindex)


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
    # 1) Prepare residue data for parallel processing
    results = []
    residues: List[Residue] = structure.residues
    for resindex, residue in enumerate(residues):
        # Skip residues that are aminoacids, nucleics, or too small
        if residue.classification in ['ion', 'solvent', 'nucleic', 'protein']:
            continue
        # Select residues atoms with MDAnalysis
        res_atoms = u.select_atoms(residue.get_selection().to_mdanalysis())
        results.append(process_residue(res_atoms,resindex))

    # 2) Process results and build dictionaries
    key_2_name = {} # To see if different name for same residue
    name_2_key = {} # To see if different residues under same name
    for result in results:
        inchikey, inchi, index = result

        # If key don't existe we create the default entry with info that is only put once
        if inchikey not in key_2_name:
            key_2_name[inchikey] = {'inchi': inchi,
                                    'resindices': [],
                                    'resname':[], 
                                    'classification': [],
                                    'has_bonds': residues[index].get_bonded_atoms()} # TO-DO: check for every residue not only the first
        # Add residue index to the list
        key_2_name[inchikey]['resindices'].append(index)
        # Add residue name to the list
        resname = residues[index].name
        if resname not in key_2_name[inchikey]['resname']:
            key_2_name[inchikey]['resname'].append(resname)
        # Add residue class to the list
        classification = residues[index].classification
        if classification not in key_2_name[inchikey]['classification']:
            key_2_name[inchikey]['classification'].append(classification)

        # Incorrect residue name, estereoisomers, lose of atoms...
        if resname not in name_2_key:
            name_2_key[resname] = []
        name_2_key[resname].append(inchikey)
    
    # 3) Check data coherence
    # Check if there are multiple names for the same InChI key
    for inchikey, data in key_2_name.items():
        if data["has_bonds"]: # TO-DO quitar, neutralizar bonds.
            warn(f"The InChIKey {inchikey} is a substructure ({data['classification']})\n"
                   "and will result in imprecise search in PubChemb due lo lack of hidrogens\n"
                   "and different charges on the bonded atoms.")
        if len(data['resname']) > 1:
            warn('Same residue with different names:\n'
                f'{inchikey} -> {str(data["resname"])}')
        if len(data['classification']) > 1:
            warn('Same residue with different classifications:\n'
                 f'{inchikey} + -> {str(data["classification"])} for names {str(data["resname"])}')   

    # Check if there are multiple InChI keys for the same name
    for name, inchikeys in name_2_key.items():
        inchikeys = set(inchikeys)
        if len(inchikeys) > 1:
            key_counts = '\n'.join([f'{key}: {len(key_2_name[key]["resindices"])}. {key_2_name[key]["inchi"]}' for key in inchikeys])
            warn(f'The residue {name} has more than one InChi key:\n'
                 f'{key_counts}')
    return key_2_name


@lru_cache(maxsize=None)
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
    "WiP"
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
