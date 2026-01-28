import MDAnalysis
import traceback
import multiprocessing
from dataclasses import dataclass, field
from mddb_workflow.tools.get_ligands import pubchem_standardization
from mddb_workflow.utils.structures import Structure
from mddb_workflow.utils.auxiliar import warn, save_json, timeout
from mddb_workflow.utils.type_hints import *
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.rdDetermineBonds import DetermineBondOrders
from MDAnalysis.converters.RDKitInferring import MDAnalysisInferrer


@dataclass
class InChIKeyData:
    """Data structure for InChI key information.

    Attributes:
        inchi (str): The InChI string for the residue.
        resindices (list[int]): List of all residue indices with this InChI key.
        fragments (list[list[int]]): Lists of residue indices that are connected as a single group.
        resnames (set[str]): Set of residue names associated with this InChI key.
        molname (str): Representative molecule name for this InChI key.
        moltype (Literal['residue', 'fragment']): Type of the molecule.
        classification (set): Set of residue classifications for this InChI key.
        frag_len (int): Length on number of residues of the fragments. 1 if no fragments are present.
        references (dict): Additional database references related to this InChI key.

    """
    inchi: str
    resindices: list[int] = field(default_factory=list)
    fragments: list[list[int]] = field(default_factory=list)
    resnames: set[str] = field(default_factory=set)
    molname: str = ''
    moltype: Literal['residue', 'fragment'] = 'residue'
    classification: set = field(default_factory=set)
    frag_len: int = 1
    references: dict = field(default_factory=dict)

    @classmethod
    def load_cache(cls, cache_dict: dict[dict | None]) -> dict[str, 'InChIKeyData']:
        """Load data from cache by converting dictionaries to InChIKeyData objects."""
        if len(cache_dict) > 0 and isinstance(next(iter(cache_dict.values())), dict):
            inchikeys = {}
            for inchikey, data in cache_dict.items():
                inchidata = cls(inchi=data['inchi'])
                for key, value in data.items():
                    setattr(inchidata, key, value)
                inchikeys[inchikey] = inchidata
            return inchikeys
        else:
            return cache_dict


def is_ferroheme(mda_atoms: 'MDAnalysis.AtomGroup') -> bool:
    """Check if the given MDAnalysis AtomGroup corresponds to a ferroheme molecule."""
    # Create a copy of the universe to remove the bonds safely
    resatoms = MDAnalysis.Merge(mda_atoms)
    resatoms.delete_bonds(resatoms.select_atoms('element Fe').bonds)
    # Convert to basic RDKit molecule (no bond order nor formal charges)
    mol = resatoms.select_atoms('not element Fe').convert_to.rdkit(inferrer=None)
    # Add hydrogen to first two nitrogens in exchange
    # for the removed Fe bonds. Which N to does not matter
    # as we standardize later
    mol_editable = Chem.RWMol(mol)
    n_atoms = [at for at in mol_editable.GetAtoms() if at.GetSymbol() == 'N']
    for n_atom in n_atoms[:2]:
        # Add H atom bonding to N
        h_idx = mol_editable.AddAtom(Chem.Atom('H'))
        mol_editable.AddBond(n_atom.GetIdx(), h_idx, Chem.BondType.SINGLE)
    # Convert back to regular molecule
    mol_with_h = mol_editable.GetMol()
    # rdDepictor.Compute2DCoords(mol_with_h)  # Optional: compute 2D coordinates for visualization
    # Use MDAnalysisInferrer to get formal charge of the molecule
    mol_mda = MDAnalysisInferrer()(mol)
    # DetermineBondOrders can crash so we first check if the standardization
    # gives us ferroheme directly without it
    standar_cid = pubchem_standardization(Chem.MolToInchi(mol_mda))
    if standar_cid and any(cid['pubchem'] in ['4971', '5353910', '25246109'] for cid in standar_cid):
        return True
    formal_charge = Chem.GetFormalCharge(mol_mda)
    DetermineBondOrders(mol_with_h, charge=formal_charge, maxIterations=1000)
    unch_mol = rdMolStandardize.ChargeParent(mol_with_h)
    inchi = Chem.MolToInchi(unch_mol)
    standar_cid = pubchem_standardization(inchi)
    return standar_cid[0]['pubchem'] == '4971'  # CID for ferroheme without Fe


@timeout(180)
def residue_to_inchi(task: tuple['MDAnalysis.AtomGroup', int]) -> tuple[str, str, int]:
    """Process a single residue to get its InChI key and related information."""
    resatoms, resindices = task
    # Convert to RDKIT and get InChI data
    error = None
    try:
        res_RD = resatoms.convert_to.rdkit(force=True)
    except KeyError as e:
        tb_str = traceback.format_exc()
        # Error that happend with accesoin A01J7 (united atoms)
        if 'bond.SetBondType(RDBONDORDER[order])' in tb_str:
            error = f'Invalid bond order {e}, failed to convert to RDKit.'
            return ('', '', resindices, error)
    except Exception:
        try:
            res_RD = resatoms.convert_to.rdkit(inferrer=None)
        except Exception as e:
            raise e
    if 'Fe' in set(resatoms.atoms.elements) and len(resatoms.atoms) > 1:
        # Metallo proteins are not well handled by RDKit/InChI so we hardcode
        # the InChI key for ferroheme that is the only case we have found
        if is_ferroheme(resatoms):
            inchikey = 'KABFMIBPWCXCRK-UHFFFAOYSA-L'
            inchi = 'InChI=1S/C34H34N4O4.Fe/c1-7-21-17(3)25-13-26-19(5)23(9-11-33(39)40)31(37-26)16-32-24(10-12-34(41)42)20(6)28(38-32)15-30-22(8-2)18(4)27(36-30)14-29(21)35-25;/h7-8,13-16H,1-2,9-12H2,3-6H3,(H4,35,36,37,38,39,40,41,42);/q;+2/p-2'
        else:
            raise NotImplementedError('Non-ferroheme residues with Fe are not supported.')
    else:
        formal_charge = (int(resatoms.atoms.charges.sum().round())
                         if hasattr(resatoms.atoms, 'charges')
                         else None)
        # For charged residues, DetermineBondOrders is better as we
        # can set the formal charge on the molecule. Step needed for A026E
        # If MDAnalysisInferrer infers the correct formal charge, we skip this step
        if formal_charge and Chem.GetFormalCharge(res_RD) != formal_charge:
            # Try/except because DetermineBondOrders can fail
            try:
                res_RD_copy = Chem.Mol(res_RD)
                DetermineBondOrders(res_RD_copy, charge=formal_charge)
                res_RD = res_RD_copy
            except Exception:
                pass
        # Calculate InChI key and string
        inchikey = Chem.MolToInchiKey(res_RD)
        # rdinchi.MolToInchi so it doesnt print the warnings
        inchi, retcode, message, logs, aux = Chem.rdinchi.MolToInchi(res_RD)
    if not inchi: raise RuntimeError('Failed to find inchi. Is the atom group too big?')
    return (inchikey, inchi, resindices, error)


def generate_inchikeys(
    universe: 'MDAnalysis.Universe',
    structure: 'Structure',
    cg_selection : 'Selection',
    parallel: bool = False,
) -> dict[str, InChIKeyData]:
    """Generate a dictionary mapping InChI keys to residue information for non-standard residues.

    This function uses MDAnalysis to parse the input structure and topology files and identifies
    residues that are not classified as 'ion', 'solvent', 'nucleic', or 'protein'. For each
    identified residue, it converts the structure to RDKit format to obtain the InChI key
    and InChI string. The resulting data is stored in dictionaries to map InChI keys to residue
    details and residue names to InChI keys. PDB coordinates are necessary to distinguish stereoisomers.

    Args:
        universe (Universe): The MDAnalysis Universe object containing the structure and topology.
        structure (Structure): The Structure object containing residues.
        parallel (bool): Whether to use multiprocessing (default False). Set to True for parallel processing (higher memory).

    Returns:
        dict: A dictionary mapping InChI keys to InChIKeyData objects.

    Notes:
        The function also performs consistency checks, warning if multiple residue names
        map to the same InChI key or if multiple InChI keys map to the same residue name,
        which can indicate mismatched residue definitions or stereoisomers.

    """
    try:
        universe.universe.atoms.charges
    except Exception:
        warn('Topology file does not have charges, InChI keys may be unreliable.')

    # 1) Prepare residue data for parallel processing
    # First group residues that are bonded together
    tasks = []
    residues = structure.residues

    # Get coarse grain atom indices
    coarse_gran_atom_indices = set(cg_selection.atom_indices)

    # Fragment = residues that are bonded together
    fragments = universe.atoms.fragments
    for i, fragment in enumerate(fragments):

        # If the fragment contains coarse grain regions then skip it
        fragment_atom_indices = set([ int(index) for index in fragment.indices ])
        if coarse_gran_atom_indices.intersection(fragment_atom_indices): continue

        resindices = fragment.residues.resindices.tolist()

        # Continue to the next fragment if any of its
        # residues are of a disallowed classification
        classes = {residues[resindex].classification for resindex in resindices}
        if (classes.intersection({'solvent'}) or
            (classes.intersection({'dna', 'rna', 'protein'}) and len(resindices) > 1)):
            # print(f'Skipping fragment {i} with classes {classes} and residues {resindices}')
            continue

        # Select residues atoms with MDAnalysis
        resatoms = universe.residues[resindices].atoms
        # If you pass a residue selection to a parallel worker, you a passing a whole MDAnalysis
        # universe, slowing the process down because you have to pickle the object
        # To avoid this we create
        resatoms = MDAnalysis.Merge(resatoms).universe.atoms
        # Convert to RDKit and get InChI data
        tasks.append((resatoms, resindices))

    results = []
    if parallel:
        # Execute tasks in parallel
        # Use 'spawn' context to avoid fork-safety issues with threading (RDKit?)
        # https://pythonspeed.com/articles/python-multiprocessing/
        ctx = multiprocessing.get_context('spawn')
        with ctx.Pool() as pool:
            try:
                print(f'Processing {len(tasks)} residues to get InChI keys (parallel mode)...')
                async_results = pool.map_async(residue_to_inchi, tasks)
                # Add timeout (e.g., 60 seconds = 1 minute per task)
                results = async_results.get(timeout=60 * len(tasks))
            except multiprocessing.TimeoutError:
                warn('Pool timeout - some tasks did not complete.')
                pool.terminate()  # Force kill all workers
                pool.join()
                raise
            except Exception:
                pool.terminate()
                pool.join()
                raise
    else:
        print(f'Processing {len(tasks)} residues to get InChI keys (sequential mode)...')
        for task in tasks:
            results.append(residue_to_inchi(task))

    # 2) Process results and build dictionaries
    inchikeys: dict[str, InChIKeyData] = {}  # To see if different name for same residue
    name_2_key = {}  # To see if different residues under same name
    errors = {}
    for (inchikey, inchi, resindices, error) in results:
        if error:
            errors.setdefault(error, []).append(*resindices)
            continue
        # Get or create the entry for this InChI key
        data = inchikeys.setdefault(inchikey, InChIKeyData(inchi=inchi))

        # Add residue index to the list
        data.resindices.extend(resindices)
        # Add residue name to the list. For multi residues we join the names
        resnames = '-'.join(sorted([residues[index].name for index in resindices]))
        data.resnames.add(resnames)
        # Add residue class to the list
        if len(resindices) > 1:
            classes = tuple(set([residues[index].classification for index in resindices]))
            data.classification.add(classes)
            # Glucolipids saved the groups of residues the form a 'fragment' to solve a
            # problem with FATSLiM later (ex: A01IR, A01J5)
            data.fragments.append(list(map(int, resindices)))
        else:
            data.classification.add(residues[resindices[0]].classification)

        # Incorrect residue name, stereoisomers, loss of atoms...
        name_2_key.setdefault(resnames, []).append(inchikey)
    if errors:
        for error, resindices in errors.items():
            warn(f'Error processing residues {resindices}: {error}')
    # 3) Check data coherence
    for inchikey, data in inchikeys.items():
        # Check if there are multiple names for the same InChI key
        if len(data.resnames) > 1:
            warn(f'Same residue with different names:\n {inchikey} -> {str(data.resnames)}')
        data.molname = list(data.resnames)[0]  # Just pick one name
        # Check if there are multiple classifications for the same InChI key
        if len(data.classification) > 1:
            warn('Same residue with different classifications:\n'
                 f'{inchikey} + -> {str(data.classification)} for names {str(data.resnames)}')
        # Check if there are multiple fragments length for the same InChI key
        if len(data.fragments) == 0:
            data.frag_len = 1
        else:
            data.moltype = 'fragment'
            frag_lens = set([len(fragment) for fragment in data.fragments])
            assert len(frag_lens) == 1, \
                f'Fragments of different lengths for InChI key {inchikey}: {str(frag_lens)}'
            data.frag_len = frag_lens.pop()

    # Check if there are multiple InChI keys for the same name
    for name, keys in name_2_key.items():
        keys = list(set(keys))
        if len(keys) < 2: continue
        counts = {}
        # Count the number of fragments
        for key in keys:
            # If there are not fragments, we use the number of residues
            counts[key] = (inchikeys[key].frag_len
                           if inchikeys[key].frag_len > 1
                           else len(inchikeys[key].resindices))
        # Format the counts for printing
        key_counts = '\n'.join([f'\t{k}: {c: >4}' for k, c in counts.items()])
        warn(f'The fragment {name} has more than one InChi key:\n'
                f'{key_counts}')

    return inchikeys


def generate_inchi_references(
    inchikeys: dict[str, 'InChIKeyData'],
    lipid_references: dict[str, dict],
    ligand_references: dict[str, dict],
    output_file: 'File',
) -> list[dict]:
    """Generate InChI references for the database."""
    inchikey_references = []
    inchikey_map = []
    for inchikey, res_data in inchikeys.items():
        # If there is force ligands, the inchikey may have changed
        ref_inchikey = ligand_references.get(inchikey, {}).get('inchikey', inchikey)
        ref_inchi = ligand_references.get(inchikey, {}).get('inchi', res_data.inchi)
        inchikey_references.append({
            'inchikey': ref_inchikey,
            'inchi': ref_inchi,
            'swisslipids': lipid_references.get(inchikey, {}).get('swisslipids', {}),
            'lipidmaps': lipid_references.get(inchikey, {}).get('lipidmaps', {}),
            'pubchem': ligand_references.get(inchikey, {}),
        })
        # Sort dictionary entries for consistency when uploading to database
        for k, v in inchikey_references[-1].items():
            if type(v) is dict:
                inchikey_references[-1][k] = dict(sorted(v.items()))

        # Get residue indices from ligand forced selections if available
        resindices = ligand_references.get(inchikey, {}).get('resindices', list(map(int, res_data.resindices)))
        inchikey_map.append({
            'inchikey': ref_inchikey,
            'name': list(res_data.resnames)[0],  # For rmsds
            # 'inchi': ref_inchi,
            # 'fragments': res_data.fragments,
            'residue_indices': resindices,
            'is_lipid': inchikey in lipid_references,
            'match': {
                'ref': {'inchikey': inchikey}
            }
        })
    save_json(inchikey_references, output_file.path)
    return inchikey_map
