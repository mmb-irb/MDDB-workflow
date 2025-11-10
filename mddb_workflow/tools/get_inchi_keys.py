import MDAnalysis
import multiprocessing
from rdkit import Chem
from mddb_workflow.utils.structures import Structure
from mddb_workflow.utils.auxiliar import warn
from mddb_workflow.utils.type_hints import *


def get_inchikeys (
    universe : 'MDAnalysis.Universe',
    structure : 'Structure',
    skip_class = {'ion', 'solvent'},
) -> dict:
    """
    Generate a dictionary mapping InChI keys to residue information for non-standard residues.

    This function uses MDAnalysis to parse the input structure and topology files and identifies
    residues that are not classified as 'ion', 'solvent', 'nucleic', or 'protein'. For each
    identified residue, it converts the structure to RDKit format to obtain the InChI key
    and InChI string. The resulting data is stored in dictionaries to map InChI keys to residue
    details and residue names to InChI keys. PDB coordinates are necesary to distinguish stereoisomers.

    Args:
        universe (Universe): The MDAnalysis Universe object containing the structure and topology.
        structure (Structure): The Structure object containing residues.
        skip_class (set): A set of residue classifications to skip. Defaults to {'ion', 'solvent'}.

    Returns:
        dict: A dictionary where keys are InChI keys and values are dictionaries containing:
            - 'inchi' (str): The InChI string for the residue
            - 'resindices' (list): A list of all residue indices with this InChI key
            - 'fragments' (list[list]): Lists of residue indices that are connected as a single group
            - 'frag_len' (int): Length of the fragments. 1 if no fragments are present.
            - 'resname' (set): Set of residue names associated with this InChI key
            - 'classification' (set): Set of residue classifications for this InChI key
    
    Notes:
        The function also performs consistency checks, warning if multiple residue names 
        map to the same InChI key or if multiple InChI keys map to the same residue name,
        which can indicate mismatched residue definitions or stereoisomers.
    """
    try:
        universe.universe.atoms.charges
    except:
        warn('Topology file does not have charges, InChI keys may be unreliable.')

    # 1) Prepare residue data for parallel processing
    # First group residues that are bonded together
    tasks = []
    residues = structure.residues

    # Fragment = residues that are bonded together
    fragments = universe.atoms.fragments
    for i, fragment in enumerate(fragments):
        resindices = fragment.residues.resindices.tolist()

        # Continue to the next fragment if any of its
        # residues are of a disallowed classification
        classes = {residues[resindex].classification for resindex in resindices}
        if (classes.intersection(skip_class) or
            (classes.intersection({'dna', 'rna', 'protein'}) and len(resindices) > 1)):
            continue

        # Select residues atoms with MDAnalysis
        resatoms = universe.residues[resindices].atoms
        if 'Cg' in resatoms.types:
            # Skip coarse grain residues
            continue
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
    key_2_data = {} # To see if different name for same residue
    name_2_key = {} # To see if different residues under same name
    for (inchikey, inchi, resindices) in results:

        # If key don't existe we create the default entry with info that is only put once
        if inchikey not in key_2_data:
            key_2_data[inchikey] = {'inchi': inchi,
                                    'resindices': [],
                                    'fragments': [], # For glucolipids
                                    'resname': set(), 
                                    'classification': set()
                                    }
        # Add residue index to the list
        key_2_data[inchikey]['resindices'].extend(resindices)
        # Add residue name to the list. For multi residues we join the names
        resname = '-'.join(sorted([residues[index].name for index in resindices]))
        key_2_data[inchikey]['resname'].add(resname)
        # Add residue class to the list
        if len(resindices) > 1:
            classes = tuple(set([residues[index].classification for index in resindices]))
            key_2_data[inchikey]['classification'].add(classes)
            # Glucolipids saved the groups of residues the form a 'fragment' to solve a 
            # problem with FATSLiM later (ex: A01IR, A01J5)
            key_2_data[inchikey]['fragments'].append(list(map(int, resindices)))
        else:
            key_2_data[inchikey]['classification'].add(residues[resindices[0]].classification)

        # Incorrect residue name, estereoisomers, lose of atoms...
        if resname not in name_2_key:
            name_2_key[resname] = []
        name_2_key[resname].append(inchikey)
    
    # 3) Check data coherence
    for inchikey, data in key_2_data.items():
        # Check if there are multiple names for the same InChI key 
        if len(data['resname']) > 1:
            warn('Same residue with different names:\n'
                f'{inchikey} -> {str(data["resname"])}')
        # Check if there are multiple classifications for the same InChI key
        if len(data['classification']) > 1:
            warn('Same residue with different classifications:\n'
                 f'{inchikey} + -> {str(data["classification"])} for names {str(data["resname"])}')
        # Check if there are multiple fragments length for the same InChI key
        if len(data['fragments']) == 0:
            data['frag_len'] = 1
        else:
            frag_len = len(set([len(fragment) for fragment in data['fragments']]))
            assert frag_len == 1, \
                f'Fragments of different lengths for InChI key {inchikey}: {str(frag_len)}'
            data['frag_len'] = frag_len

    # Check if there are multiple InChI keys for the same name
    for name, inchikeys in name_2_key.items():
        inchikeys = list(set(inchikeys))
        if len(inchikeys) < 2: continue
        counts = {}
        # Count the number of fragments 
        for key in inchikeys:
            # If there are not fragments, we use the number of residues
            counts[key] = (key_2_data[key]["frag_len"] 
                           if key_2_data[key]["frag_len"] > 1 
                           else len(key_2_data[key]["resindices"]))
        # Format the counts for printing
        key_counts = '\n'.join([f'\t{key}: {counts[key]: >4}' for key in inchikeys])
        warn(f'The fragment {name} has more than one InChi key:\n'
                f'{key_counts}')
    
    # Convert sets to lists for JSON serialization
    for inchikey in key_2_data:
        key_2_data[inchikey]['resname'] = list(key_2_data[inchikey]['resname'])
        key_2_data[inchikey]['classification'] = list(key_2_data[inchikey]['classification'])

    return key_2_data

def residue_to_inchi(task: tuple['MDAnalysis.AtomGroup', int]) -> tuple[str, str, int]:
    """
    Process a single residue to get its InChI key and related information.
    """
    resatoms , resindices = task
    # Convert to RDKIT and get InChI data
    res_RD = resatoms.convert_to.rdkit()
    # Calculate InChI key and string
    inchikey = Chem.MolToInchiKey(res_RD)
    # rdinchi.MolToInchi so it doesnt print the warnings
    inchi, retcode, message, logs, aux  = Chem.rdinchi.MolToInchi(res_RD)
    return (inchikey, inchi, resindices)

def obtain_inchikey_from_pdb(pdb_file : str) -> Optional[str]:
    mol = Chem.MolFromPDBFile(pdb_file, sanitize=True)
    if mol is None:
        print("❌ InchiKey couldn't be constructed from PDB file.")
        return None
    return Chem.inchi.MolToInchiKey(mol)

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
