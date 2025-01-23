import MDAnalysis
from rdkit import Chem
from model_workflow.utils.structures import Structure
from model_workflow.utils.type_hints import *
from model_workflow.utils.warnings import warn
from functools import lru_cache
import requests


def get_connected_residues(
        residue: 'MDAnalysis.Residue',
        visited=None):
    """
    Recursively finds all residues that are connected through bonds to the given residue.
    Parameters
    ----------
    residue : MDAnalysis.core.groups.Residue
        The residue to start the search from
    visited : set, optional
        Set of already visited residue IDs to prevent infinite recursion. 
        Default is None, which initializes an empty set.
    Returns
    -------
    list
        A set containing the residue IDs of all residues that are connected 
        through bonds to the input residue, either directly or indirectly.
    Notes
    -----
    This function traverses the molecular structure recursively, following bonds 
    between residues. It keeps track of visited residues to avoid cycles in the
    traversal. Both direct bonds between residues and indirect connections through
    other residues are included in the result.
    """
    
    if visited is None:
        visited = set()
    
    # Add current residue to visited set
    visited.add(residue.resid)
    
    # Get direct external bonds
    external_bonds = set()
    for bond in residue.atoms.bonds:
        external_bonds.update(bond.atoms.resindices)
    
    # Recursively check external bonds
    all_external_bonds = external_bonds.copy()
    u = residue.universe
    for ext_resid in external_bonds:
        if ext_resid not in visited:
            # Get external bonds for the connected residue
            recursive_bonds = get_connected_residues(u.residues[ext_resid], visited)
            all_external_bonds.update(recursive_bonds)
    
    return list(all_external_bonds)


def process_residue(res_atoms: 'MDAnalysis.AtomGroup', 
                    resindex : int) -> Tuple[str, str, int]:
    """
    Process a single residue to get its InChI key and related information.
    """
    # Convert to RDKIT and get InChI data
    res_RD = res_atoms.convert_to("RDKIT")
    # Calculate InChI key and string
    inchikey = Chem.MolToInchiKey(res_RD)  # slow step, 50% of time 
    # rdinchi.MolToInchi so it doesnt print the warnings
    inchi, retcode, message, logs, aux  = Chem.rdinchi.MolToInchi(res_RD)
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
    # 1) Prepare residue data for parallel processing (not yet implemented)
    # First group residues that are bonded together
    results = []
    residues: List[Residue] = structure.residues
    visited = set()
    for resindex in range(len(residues)):
        # Skip residues that have already been visited
        if resindex in visited:
            continue
        # Residues grouped by bonded atoms
        res_grp_idx = get_connected_residues(u.residues[resindex], visited) 
        classes = [residues[grp_idx].classification for grp_idx in res_grp_idx]
        # Skip residues that are aminoacids, nucleics, or too small
        # We also skips residues connected to them: glicoprotein, lipid-anchored protein...
        if any(cls in ['ion', 'solvent', 'nucleic', 'protein'] for cls in classes):
            continue
        # Select residues atoms with MDAnalysis
        res_atoms = u.residues[res_grp_idx].atoms
        # Convert to RDKit and get InChI data
        results.append(process_residue(res_atoms,res_grp_idx))

    # 2) Process results and build dictionaries
    key_2_name = {} # To see if different name for same residue
    name_2_key = {} # To see if different residues under same name
    for result in results:
        inchikey, inchi, indexes = result

        # If key don't existe we create the default entry with info that is only put once
        if inchikey not in key_2_name:
            key_2_name[inchikey] = {'inchi': inchi,
                                    'resindices': [],
                                    'resname': set(), 
                                    'classification': set()
                                    }
        # Add residue index to the list
        key_2_name[inchikey]['resindices'].extend(indexes)
        # Add residue name to the list. For multi residues we join the names
        resname = '-'.join([residues[index].name for index in indexes])
        key_2_name[inchikey]['resname'].add(resname)
        # Add residue class to the list
        classes = set([residues[index].classification for index in indexes])
        key_2_name[inchikey]['classification'].update(classes)

        # Incorrect residue name, estereoisomers, lose of atoms...
        if resname not in name_2_key:
            name_2_key[resname] = []
        name_2_key[resname].append(inchikey)
    
    # 3) Check data coherence
    # Check if there are multiple names for the same InChI key
    for inchikey, data in key_2_name.items():
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
            key_counts = '\n'.join([f'{key}: {len(key_2_name[key]["resindices"]): >4}. '
                                    f'{key_2_name[key]["inchi"]}. '
                                    f'Sample index: {key_2_name[key]["resindices"][0]}'
                                      for key in inchikeys])
            warn(f'The residue {name} has more than one InChi key:\n'
                 f'{key_counts}')
    return key_2_name


@lru_cache(maxsize=None)
def is_in_LIPID_MAPS(inchikey) -> dict:
    """Search the InChi keys in LIPID MAPS"""
    headers = {'accept': 'json'}
    # https://www.lipidmaps.org/resources/rest
    # Output item = physchem, is the only one that returns data for the inchi key
    # for only the two first layers (main and atom connection)
    # To see InChiKey layers: https://www.inchi-trust.org/
    first_layer = inchikey[:14]
    url = f"https://www.lipidmaps.org/rest/compound/inchi_key/{first_layer}/physchem"
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


def Residue_2_Mol(res: 'Residue'):
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
