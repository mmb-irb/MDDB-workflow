import MDAnalysis
from rdkit import Chem
from model_workflow.utils.structures import Structure
from model_workflow.utils.auxiliar import warn
from model_workflow.utils.type_hints import *
from functools import lru_cache
import requests
import multiprocessing


def residue_to_inchi(task: Tuple['MDAnalysis.AtomGroup', int]) -> Tuple[str, str, int]:
    """
    Process a single residue to get its InChI key and related information.
    """
    resatoms , resindices = task
    # Convert to RDKIT and get InChI data
    res_RD = resatoms.convert_to("RDKIT")
    # Calculate InChI key and string
    inchikey = Chem.MolToInchiKey(res_RD)  # slow step, 50% of time 
    # rdinchi.MolToInchi so it doesnt print the warnings
    inchi, retcode, message, logs, aux  = Chem.rdinchi.MolToInchi(res_RD)
    return (inchikey, inchi, resindices)


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
        dict: A dictionary where keys are InChI keys and values are dictionaries containing:
                - 'inchi' (str): The InChI string for the residue
                - 'resindices' (list): A list of all residue indices with this InChI key
                - 'resgroups' (list): Lists of residue indices that are connected as a single group
                - 'resname' (set): Set of residue names associated with this InChI key
                - 'classification' (set): Set of residue classifications for this InChI key
    Notes:
        The function also performs consistency checks, warning if multiple residue names 
        map to the same InChI key or if multiple InChI keys map to the same residue name,
        which can indicate mismatched residue definitions or stereoisomers.
    """
    # 1) Prepare residue data for parallel processing (not yet implemented)
    # First group residues that are bonded together
    tasks = []
    residues: List[Residue] = structure.residues

    # Check if the residues are the same as in the universe TODO: BORRAR
    for i in range(len(residues)):
        assert u.residues.resnames[i] == residues[i].name

    # Fragment = residues that are bonded together
    fragments = u.atoms.fragments  
    skip_class = {'ion', 'solvent', 'dna', 'rna', 'protein'}
    for i, fragment in enumerate(fragments):
        resindices = list(set(fragment.resindices))
        
        # Continue to the next fragment if any of its 
        # residues are of a disallowed classification
        if any(residues[resindex].classification in skip_class
            for resindex in resindices):
            continue

        # Select residues atoms with MDAnalysis
        resatoms = u.residues[resindices].atoms
        # If you pass a residue selection to a parallel worker, you a passing a whole MDAnalysis
        # universe, slowing the process down because you have to pickle the object
        # To avoid this we create
        resatoms = MDAnalysis.Merge(resatoms).universe.atoms
        # Convert to RDKit and get InChI data
        tasks.append((resatoms, resindices))

    results = []
    # Execute tasks in parallel
    with multiprocessing.Pool() as pool:
        results = pool.map(residue_to_inchi, tasks)

    # 2) Process results and build dictionaries
    key_2_name = {} # To see if different name for same residue
    name_2_key = {} # To see if different residues under same name
    for (inchikey, inchi, resindices) in results:

        # If key don't existe we create the default entry with info that is only put once
        if inchikey not in key_2_name:
            key_2_name[inchikey] = {'inchi': inchi,
                                    'resindices': [],
                                    'resgroups': [], # For glucolipids
                                    'resname': set(), 
                                    'classification': set()
                                    }
        # Add residue index to the list
        key_2_name[inchikey]['resindices'].extend(resindices)
        # Add residue name to the list. For multi residues we join the names
        resname = '-'.join(sorted([residues[index].name for index in resindices]))
        key_2_name[inchikey]['resname'].add(resname)
        # Add residue class to the list
        if len(resindices) > 1:
            classes = tuple(set([residues[index].classification for index in resindices]))
            key_2_name[inchikey]['classification'].add(classes)
            # Glucolipids saved the groups of residues the form a 'fragment' to solve a 
            # problem with FATSLiM later (ex: A01IR, A01J5)
            key_2_name[inchikey]['resgroups'].append(list(map(int, resindices)))
        else:
            key_2_name[inchikey]['classification'].add(residues[resindices[0]].classification)

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
            key_counts = '\n'.join([f'\t{key}: {len(key_2_name[key]["resindices"]): >4}'
                                      for key in inchikeys])
            warn(f'The residue {name} has more than one InChi key:\n'
                 f'{key_counts}')
    return key_2_name


@lru_cache(maxsize=None)
def is_in_LIPID_MAPS(inchikey, only_first_layer=False) -> dict:
    """Search the InChi keys in LIPID MAPS"""
    headers = {'accept': 'json'}
    # https://www.lipidmaps.org/resources/rest
    # Output item = physchem, is the only one that returns data for the inchi key
    # for only the two first layers (main and atom connection)
    # To see InChiKey layers: 
    # https://www.inchi-trust.org/about-the-inchi-standard/
    # Or https://www.rhea-db.org/help/inchi-inchikey#What_is_an_InChIKey_
    key = inchikey[:14] if only_first_layer else inchikey
    url = f"https://www.lipidmaps.org/rest/compound/inchi_key/{key}/all"
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        js = response.json()
        if js != []:
           return js      
        else:
            return False  
    else:
        print(f"Error for {inchikey}: {response.status_code}")


@lru_cache(maxsize=None)
def is_in_swiss_lipids(inchikey, only_first_layer=False, 
                       protonation=True) -> dict:
    """Search the InChi keys in SwissLipids. Documentation: https://www.swisslipids.org/#/api"""
    key = inchikey[:14] if only_first_layer else inchikey
    headers = {'accept': 'json'}
    url = f"https://www.swisslipids.org/api/index.php/search?term={key}"
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        return response.json()[0]
    else:
        return False
    
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
