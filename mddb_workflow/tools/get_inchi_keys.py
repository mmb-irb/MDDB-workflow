import MDAnalysis
import traceback
import multiprocessing
from mddb_workflow.utils.auxiliar import warn, save_json, SocketTimeout, DatabaseUnavailableError
from mddb_workflow.utils.constants import INCHI_REFERENCE_VERSION, INCOMPLETE_REFERENCE_VERSION
from mddb_workflow.core.dataclass import InChIKeyData
from mddb_workflow.utils.cache import get_cached_function
from mddb_workflow.utils.pubchem import pubchem_standardization
from mddb_workflow.utils.selections import Selection
from mddb_workflow.utils.structures import Structure
from mddb_workflow.utils.type_hints import *
from re import match
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.rdDetermineBonds import DetermineBondOrders
from MDAnalysis.converters.RDKitInferring import MDAnalysisInferrer

# Regular expression matching formulas of single ions
ION_FORMULA = r'^[A-Z]{1}[a-z]?[+-]?\d?$'


def is_ferroheme(mda_atoms: 'MDAnalysis.AtomGroup') -> bool:
    """Check if the given MDAnalysis AtomGroup corresponds to a ferroheme molecule."""
    try:
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
    except Exception:
        return False


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
            error = 'Non-ferroheme residues with Fe are not supported.'
            return ('', '', resindices, error)
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
                DetermineBondOrders(res_RD_copy, charge=formal_charge, maxIterations=1000)
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
    cg_selection: 'Selection',
    parallel: bool = True,
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
        cg_selection (Selection): Selection object for coarse-grained regions.
        parallel (bool): Whether to use multiprocessing (default True). Uses more memory, set to False to process residues serially.

    Returns:
        dict: A dictionary mapping InChI keys to InChIKeyData objects.

    Notes:
        The function also performs consistency checks, warning if multiple residue names
        map to the same InChI key or if multiple InChI keys map to the same residue name,
        which can indicate mismatched residue definitions or stereoisomers.

    """
    socket_timeout = SocketTimeout(120)
    socket_timeout.start()
    if not hasattr(universe.atoms, 'charges'):
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
        fragment_atom_indices = set([int(index) for index in fragment.indices])
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
        # Skip residues that only contain dummy atoms (Dm) as tey have no element and RDKit breaks
        if set(resatoms.types) == {'Dm'}:
            continue
        # If you pass a residue selection to a parallel worker, you a passing a whole MDAnalysis
        # universe, slowing the process down because you have to pickle the object
        # To avoid this we create
        resatoms = MDAnalysis.Merge(resatoms).universe.atoms
        # Convert to RDKit and get InChI data
        tasks.append((resatoms, resindices))

    results = []
    # Starting workers is slow, so a single residue is processed directly
    if parallel and len(tasks) > 1:
        # Execute tasks in parallel
        # Use 'spawn' context to avoid fork-safety issues with threading (RDKit?)
        # https://pythonspeed.com/articles/python-multiprocessing/
        ctx = multiprocessing.get_context('spawn')
        with ctx.Pool(processes=min(len(tasks), multiprocessing.cpu_count())) as pool:
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
            errors.setdefault(error, []).extend(resindices)
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
    socket_timeout.stop()
    return inchikeys


def get_database_inchikey_reference(
    inchikey: str,
    cache: 'Cache',
    database: Optional['Database'] = None,
) -> Optional[dict]:
    """Given an InChI key, get its reference from the MDDB database, if any.
    References generated with a different version are ignored so they are generated again.
    """
    if not database: return None
    cached_get_database_reference = get_cached_function(database.get_reference_data, cache)
    try:
        reference = cached_get_database_reference('inchikeys', inchikey)
    # If the database is not available then the reference will be generated again
    except DatabaseUnavailableError:
        return None
    if not reference: return None
    if reference.get('version') != INCHI_REFERENCE_VERSION: return None
    return reference


def make_inchikey_reference(
    inchikey: str,
    inchi: str,
    lipid_reference: Optional[dict],
    ligand_reference: Optional[dict],
) -> dict:
    """Make the database reference of an InChI key from its lipid and ligand references.

    Lipids have their LIPID MAPS category and whether they were forced by the user, so lipids missing
    in the lipid databases (e.g. synthetic ionizable lipids) can be told apart from other molecules.
    References of lipids made while the lipid databases were not available get an incomplete version,
    so they are not reused from the MDDB database and they are remade later.
    """
    lipid = None
    version = INCHI_REFERENCE_VERSION
    if lipid_reference:
        lipid = {'category': lipid_reference.get('category'), 'forced': lipid_reference.get('forced', False)}
        if lipid_reference.get('version') == INCOMPLETE_REFERENCE_VERSION:
            version = INCOMPLETE_REFERENCE_VERSION
    lipid_reference = lipid_reference or {}
    reference = {
        'inchikey': inchikey,
        'inchi': inchi,
        'lipid': lipid,
        'swisslipids': lipid_reference.get('swisslipids'),
        'lipidmaps': lipid_reference.get('lipidmaps'),
        'ligand': ligand_reference,
        'version': version
    }
    # Sort dictionary entries for consistency when uploading to database
    return {k: dict(sorted(v.items())) if type(v) is dict else v for k, v in reference.items()}


def generate_inchi_references(
    inchikeys: dict[str, 'InChIKeyData'],
    lipid_references: dict[str, dict],
    ligand_references: dict[str, dict],
    membrane_map: dict,
    structure: 'Structure',
    output_file: 'File',
) -> list[dict]:
    """Generate InChI references for the database."""
    inchikey_references = []  # For the database
    inchikey_map = []  # For the workflow
    # Add the user selected ligands to the references, if any
    for k, v in ligand_references.items():
        if v.get('forced', False):
            inchikeys[k] = InChIKeyData(inchi=v['inchi'], resnames=set(v['resnames']), resindices=v['resindices'])
    for inchikey, res_data in inchikeys.items():
        # If there is force ligands, the inchikey may have changed
        ref_inchikey = ligand_references.get(inchikey, {}).get('inchikey', inchikey)
        ref_inchi = ligand_references.get(inchikey, {}).get('inchi', res_data.inchi)
        inchikey_references.append(make_inchikey_reference(ref_inchikey, ref_inchi,
            lipid_references.get(inchikey), ligand_references.get(inchikey)))

        # Get residue indices from ligand forced selections if available
        resindices = ligand_references.get(inchikey, {}).get('resindices', list(map(int, res_data.resindices)))
        inchikey_map.append({
            'inchikey': ref_inchikey,  # InChI with the possible corrections to find the ligand reference
            'name': list(res_data.resnames)[0],  # For rmsds workflow analysis
            'residue_indices': resindices,
            'is_lipid': inchikey in lipid_references,
            'is_ligand': False,  # Updated later
            'generated_inchikey': inchikey,  # Original inchikey generated from the structure
        })
    save_json(inchikey_references, output_file.path)

    # Use inchikeys to 'classify' residues as ligands
    no_mem_lipid = Selection(membrane_map['no_mem_lipid'])
    for inchikey_data in inchikey_map:
        inchikey = inchikey_data['generated_inchikey']
        if inchikey not in ligand_references:
            continue
        # Skip ions
        formula = ligand_references[inchikey].get('formula', 'missing formula')
        is_ion = match(ION_FORMULA, formula)
        if is_ion:
            continue
        if inchikey_data['is_lipid']:
            ligand_selection = structure.select_residue_indices(inchikey_data['residue_indices'])
            # Take only the lipids that are not part of a membrane (fatty-like ligands)
            no_mem_ligand_selection = ligand_selection & no_mem_lipid
            # If most of the lipid is in a membrane, we do not mark it as a ligand
            if len(no_mem_ligand_selection) / len(ligand_selection) < 0.5:
                continue
        inchikey_data['is_ligand'] = True
    return inchikey_map
