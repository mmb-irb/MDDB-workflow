"""Classify lipids from their structure following the LIPID MAPS categories.

https://lipidmaps.org/resources/education/classification

Each category is recognised by SMARTS of its core structure and, where it makes sense,
the hydrophobic chains attached to it, so drugs and other small molecules are not lipids:

- GP glycerophospholipids: glycerol phosphate with acyl/alkyl chains.
- SP sphingolipids: sphingoid base (2-amino-1,3-diol with a long chain).
- GL glycerolipids: glycerol with acyl/alkyl chains.
- ST sterol lipids: steroid nucleus or secosteroid (vitamin D) core with an aliphatic side chain in the
  D ring, as in sterols, bile acids and vitamin D. Steroid hormones and steroid drugs are not lipids here.
- SL saccharolipids: sugars with acyl chains directly attached.
- PR prenol lipids: chains of 3+ head to tail isoprene units and retinoids/carotenoids with long chains.
- FA fatty acyls: a long hydrophobic chain making a large part of the molecule, two long chains
  (synthetic double tailed lipids, e.g. ionizable lipids) or the prostanoid skeleton.
- PK polyketides are not recognised: they share no core structure (flavonoids, macrolides,
  tetracyclines...) and they are drugs or natural products rather than membrane lipids.

Peptides with acyl chains (lipopeptide antibiotics) are not saccharolipids nor fatty acyls.

SMARTS use generic bonds and no hydrogen counts since the molecules come from MD structures,
whose InChIs may have unusual tautomers (e.g. amides as imidic acids) or charges.
"""
from collections import Counter, deque
from rdkit import Chem, rdBase
from mddb_workflow.utils.type_hints import *

FATTY_ACYLS = 'FA'
GLYCEROLIPIDS = 'GL'
GLYCEROPHOSPHOLIPIDS = 'GP'
SPHINGOLIPIDS = 'SP'
STEROL_LIPIDS = 'ST'
PRENOL_LIPIDS = 'PR'
SACCHAROLIPIDS = 'SL'
POLYKETIDES = 'PK'
LIPID_CATEGORIES = {
    FATTY_ACYLS: 'Fatty Acyls',
    GLYCEROLIPIDS: 'Glycerolipids',
    GLYCEROPHOSPHOLIPIDS: 'Glycerophospholipids',
    SPHINGOLIPIDS: 'Sphingolipids',
    STEROL_LIPIDS: 'Sterol Lipids',
    PRENOL_LIPIDS: 'Prenol Lipids',
    SACCHAROLIPIDS: 'Saccharolipids',
    POLYKETIDES: 'Polyketides',
}

# Atoms in rings larger than this are considered part of chains (e.g. archaeal macrocyclic tetraethers)
MAX_RING_SIZE = 30
# Minimum hydrophobic chain length (carbons only bonded to carbons) of the chains attached to lipid cores
# Hexanoyl chains have 5 of these carbons since the carbonyl carbon is not counted (e.g. DHPC)
MIN_CORE_CHAIN = 4
# Minimum hydrophobic chain length of the sphingoid base
MIN_SPHINGOID_CHAIN = 8
# Minimum hydrophobic chain length of the steroid side chain in the D ring, 0 to accept any steroid
# Bile acids have 4 and cholesterol 7, while steroid hormones have no side chain (estradiol, testosterone)
# or it has 1 carbon (progesterone, cortisol), so they are left out since they do not form membranes
MIN_STEROID_SIDE_CHAIN = 4
# Minimum hydrophobic chain length of fatty acyls (lauric acid has 11)
MIN_FATTY_CHAIN = 10
# Minimum fraction of heavy atoms in hydrophobic chains for a single chain fatty acyl
# Polyethylene glycol atoms are not counted since it is a polymeric head (e.g. PEG-lipids)
MIN_CHAIN_FRACTION = 0.3
# Molecules with this many peptide bonds are peptides
MIN_PEPTIDE_BONDS = 2

# Carbon not in small rings
_C = f'[#6;!r{{3-{MAX_RING_SIZE}}}]'
# Oxygen not in small rings
_O = f'[#8;!r{{3-{MAX_RING_SIZE}}}]'
# Carbons with 2 and 3 heavy atoms bonded (CH2 and CH when they are saturated)
_CH2 = f'[#6;D2;!r{{3-{MAX_RING_SIZE}}}]'
_CH = f'[#6;D3;!r{{3-{MAX_RING_SIZE}}}]'
SMARTS = {
    # Glycerol phosphate (glycerol backbone in the sn-3 position, sn-1 in archaea) and glycerol backbone
    'glycerophosphate': f'[#15](~[#8])(~[#8])(~[#8])-{_O}-{_CH2}-{_CH}(-{_O})-{_CH2}-{_O}',
    'glycerol': f'{_CH2}(-{_O})-{_CH}(-{_O})-{_CH2}-{_O}',
    # Sphingoid base: CH2(O)-CH(N)-CH(O)-C (1-hydroxy may be phosphorylated or glycosylated, 3-hydroxy may be oxo)
    'sphingoid': f'{_CH2}(-[#8])-{_CH}(-[#7;!r{{3-{MAX_RING_SIZE}}}])-{_C}(~{_O})-{_C}',
    # Steroid nucleus (gonane): fused 6-6-6-5 rings
    'steroid': '[#6]1~[#6]~[#6]~[#6]2~[#6](~[#6]~1)~[#6]~[#6]~[#6]1~[#6]~2~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]~[#6]~1~2',
    # Secosteroid (vitamin D): A ring linked by a diene to the CD rings
    'secosteroid': '[#6]1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#6;!R]~[#6;!R]~[#6]1~[#6]~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]~[#6]~1~2',
    # Pyranose and furanose rings
    'pyranose': '[#6;R]1~[#6;R]~[#6;R]~[#6;R]~[#6;R]~[#8;R]~1',
    'furanose': '[#6;R]1~[#6;R]~[#6;R]~[#6;R]~[#8;R]~1',
    # Three head to tail isoprene units: methyl branches every 4 carbons, with unbranched carbons between them
    'isoprenoid': f'{_CH}(-[#6;D1])~{_CH2}~{_CH2}~{_CH2}~{_CH}(-[#6;D1])~{_CH2}~{_CH2}~{_CH2}~{_CH}(-[#6;D1])',
    # Retinoids and carotenoids: 1,1,5-trimethylcyclohexene with a polyene chain with methyl branches
    'retinoid': '[#6;D1]-[#6;R]1(-[#6;D1])~[#6;R]~[#6;R]~[#6;R]~[#6;R](-[#6;D1])~[#6;R]~1-[#6;!R]~[#6;!R]~[#6;!R](-[#6;D1])~[#6;!R]~[#6;!R]~[#6;!R]~[#6;!R](-[#6;D1])',
    # Prostanoid skeleton: cyclopentane with two chains in adjacent carbons, one ending in a carboxyl
    'prostanoid': '[#6;r5]1(-[#6;!R]~[#6;!R]~[#6;!R]~[#6;!R]~[#6;!R]~[#6;!R]~[#6;!R](~[#8])~[#8])~[#6;r5](-[#6;!R]~[#6;!R]~[#6;!R]~[#6;!R]~[#6;!R])~[#6;r5]~[#6;r5]~[#6;r5]~1',
    # Polyethylene glycol unit
    'peg': f'{_O}-{_CH2}-{_CH2}-{_O}',
    # Peptide bond between two amino acids: N-CA-C(=O)-N
    'peptide': '[#7]-[#6;X4]-[#6;X3](~[#8])~[#7]-[#6;X4]',
}
PATTERNS = {name: Chem.MolFromSmarts(smarts) for name, smarts in SMARTS.items()}


def inchi_to_mol(inchi: str) -> Optional['Chem.Mol']:
    """Build the RDKit molecule of an InChI, keeping only its largest component (e.g. no counter ions)."""
    try:
        # Block the RDKit InChI warnings
        with rdBase.BlockLogs():
            mol = Chem.MolFromInchi(inchi)
    except Exception:
        return None
    if mol is None or mol.GetNumAtoms() == 0: return None
    fragments = Chem.GetMolFrags(mol, asMols=True)
    return max(fragments, key=lambda fragment: fragment.GetNumHeavyAtoms())


class HydrophobicChains:
    """Hydrophobic chains of a molecule: connected carbons out of small rings which are only bonded to carbons.

    Carbons bonded to heteroatoms (e.g. carbonyl, glycerol or polyethylene glycol carbons) break the chains,
    so a hexadecanoyl chain has 15 carbons.
    """

    def __init__(self, mol: 'Chem.Mol'):
        self.mol = mol
        ring_info = mol.GetRingInfo()
        self.small_ring_atoms = {atom.GetIdx() for atom in mol.GetAtoms()
            if 0 < ring_info.MinAtomRingSize(atom.GetIdx()) <= MAX_RING_SIZE}
        self.atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6
            and atom.GetIdx() not in self.small_ring_atoms
            and all(neighbor.GetAtomicNum() in (1, 6) for neighbor in atom.GetNeighbors())}
        self.neighbors = {index: [neighbor.GetIdx() for neighbor in mol.GetAtomWithIdx(index).GetNeighbors()
            if neighbor.GetIdx() in self.atoms] for index in self.atoms}
        # Chain (connected component) of each atom and the length of the longest path of each chain
        self.chain_of = {}
        self.lengths = []
        for index in self.atoms:
            if index in self.chain_of: continue
            farthest, _, component = self._bfs(index)
            _, length, _ = self._bfs(farthest)
            for atom in component: self.chain_of[atom] = len(self.lengths)
            self.lengths.append(length)
        counts = Counter(self.chain_of.values())
        self.sizes = [counts[chain] for chain in range(len(self.lengths))]

    def _bfs(self, start: int) -> tuple[int, int, set[int]]:
        """Get the farthest atom from the start atom, its distance in atoms and the visited atoms."""
        distances = {start: 1}
        queue = deque([start])
        farthest = start
        while queue:
            current = queue.popleft()
            for neighbor in self.neighbors[current]:
                if neighbor in distances: continue
                distances[neighbor] = distances[current] + 1
                if distances[neighbor] > distances[farthest]: farthest = neighbor
                queue.append(neighbor)
        return farthest, distances[farthest], set(distances)

    @property
    def longest(self) -> int:
        """Length of the longest chain."""
        return max(self.lengths, default=0)

    def count_longer_than(self, min_length: int) -> int:
        """Count the chains with at least the given length."""
        return sum(length >= min_length for length in self.lengths)

    def substituent_length(self, start: int, excluded: set[int]) -> int:
        """Get the longest chain in the substituent starting at the given atom.
        The substituent is made of the carbons out of small rings connected to the start atom, without crossing
        heteroatoms nor the excluded atoms (e.g. the lipid core), so it includes carbonyl and hydroxylated carbons.
        """
        if start in excluded: return 0
        visited = {start}
        queue = deque([start])
        length = 0
        while queue:
            current = queue.popleft()
            atom = self.mol.GetAtomWithIdx(current)
            if atom.GetAtomicNum() != 6 or current in self.small_ring_atoms: continue
            if current in self.chain_of:
                length = max(length, self.lengths[self.chain_of[current]])
            for neighbor in atom.GetNeighbors():
                index = neighbor.GetIdx()
                if index in visited or index in excluded: continue
                visited.add(index)
                queue.append(index)
        return length

    def oxygen_substituent_length(self, oxygen: int, excluded: set[int]) -> int:
        """Get the longest chain linked to an oxygen (or nitrogen), as in esters, ethers and amides."""
        atom = self.mol.GetAtomWithIdx(oxygen)
        return max((self.substituent_length(neighbor.GetIdx(), excluded) for neighbor in atom.GetNeighbors()
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in excluded), default=0)


def classify_lipid(inchi: str) -> Optional[str]:
    """Get the LIPID MAPS category of the molecule of an InChI, or None if it is not a lipid."""
    mol = inchi_to_mol(inchi)
    if mol is None: return None
    return classify_lipid_mol(mol)


def classify_lipid_mol(mol: 'Chem.Mol') -> Optional[str]:
    """Get the LIPID MAPS category of an RDKit molecule, or None if it is not a lipid."""
    chains = HydrophobicChains(mol)
    if has_acylated_glycerol(mol, chains, PATTERNS['glycerophosphate']):
        return GLYCEROPHOSPHOLIPIDS
    if has_sphingoid_base(mol, chains):
        return SPHINGOLIPIDS
    if has_acylated_glycerol(mol, chains, PATTERNS['glycerol']):
        return GLYCEROLIPIDS
    if has_sterol_core(mol, chains):
        return STEROL_LIPIDS
    if mol.HasSubstructMatch(PATTERNS['isoprenoid']) or (
        mol.HasSubstructMatch(PATTERNS['retinoid']) and chains.longest >= MIN_FATTY_CHAIN):
        return PRENOL_LIPIDS
    if len(mol.GetSubstructMatches(PATTERNS['peptide'])) >= MIN_PEPTIDE_BONDS:
        return None
    if has_acylated_sugar(mol, chains):
        return SACCHAROLIPIDS
    if is_fatty_acyl(mol, chains):
        return FATTY_ACYLS
    return None


def has_acylated_glycerol(mol: 'Chem.Mol', chains: HydrophobicChains, pattern: 'Chem.Mol') -> bool:
    """Check if a glycerol backbone has an oxygen linked to a chain (acyl, alkyl or alkenyl)."""
    for match in mol.GetSubstructMatches(pattern):
        oxygens = [index for index in match if mol.GetAtomWithIdx(index).GetAtomicNum() == 8]
        # The phosphate oxygens are not part of the glycerol
        oxygens = oxygens[-3:]
        if any(chains.oxygen_substituent_length(oxygen, set(match)) >= MIN_CORE_CHAIN for oxygen in oxygens):
            return True
    return False


def has_sphingoid_base(mol: 'Chem.Mol', chains: HydrophobicChains) -> bool:
    """Check if there is a sphingoid base with its long chain."""
    for match in mol.GetSubstructMatches(PATTERNS['sphingoid']):
        if chains.substituent_length(match[-1], set(match[:-1])) >= MIN_SPHINGOID_CHAIN:
            return True
    return False


def has_sterol_core(mol: 'Chem.Mol', chains: HydrophobicChains) -> bool:
    """Check if there is a steroid or secosteroid core with the minimum side chain in the D ring, if any."""
    for pattern in [PATTERNS['steroid'], PATTERNS['secosteroid']]:
        for match in mol.GetSubstructMatches(pattern):
            core = set(match)
            # The last 5 atoms of both patterns are the D ring
            for ring_atom in match[-5:]:
                for neighbor in mol.GetAtomWithIdx(ring_atom).GetNeighbors():
                    if chains.substituent_length(neighbor.GetIdx(), core) >= MIN_STEROID_SIDE_CHAIN:
                        return True
    return False


def has_acylated_sugar(mol: 'Chem.Mol', chains: HydrophobicChains) -> bool:
    """Check if a sugar ring has a fatty acid directly linked as ester or amide."""
    sugar_rings = mol.GetSubstructMatches(PATTERNS['pyranose']) + mol.GetSubstructMatches(PATTERNS['furanose'])
    for ring in sugar_rings:
        ring_atoms = set(ring)
        # Ring carbons and their exocyclic carbons (e.g. C6 of hexoses)
        carbons = {index for index in ring if mol.GetAtomWithIdx(index).GetAtomicNum() == 6}
        carbons |= {neighbor.GetIdx() for index in list(carbons) for neighbor in mol.GetAtomWithIdx(index).GetNeighbors()
            if neighbor.GetAtomicNum() == 6 and not neighbor.IsInRing()}
        for carbon in carbons:
            for heteroatom in mol.GetAtomWithIdx(carbon).GetNeighbors():
                if heteroatom.GetAtomicNum() not in (7, 8) or heteroatom.GetIdx() in ring_atoms: continue
                for acyl in heteroatom.GetNeighbors():
                    if acyl.GetIdx() == carbon or not is_acyl_carbon(acyl): continue
                    if chains.substituent_length(acyl.GetIdx(), ring_atoms | {carbon, heteroatom.GetIdx()}) >= MIN_CORE_CHAIN:
                        return True
    return False


def is_acyl_carbon(atom: 'Chem.Atom') -> bool:
    """Check if a carbon is a carbonyl carbon of an ester or amide (also as imidic acid tautomer)."""
    if atom.GetAtomicNum() != 6 or atom.IsInRing(): return False
    heteroatoms = [neighbor for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() in (7, 8)]
    return len(heteroatoms) >= 2 and any(neighbor.GetAtomicNum() == 8 for neighbor in heteroatoms)


def is_fatty_acyl(mol: 'Chem.Mol', chains: HydrophobicChains) -> bool:
    """Check if the molecule is a fatty acyl: a long hydrophobic chain making a large part of the molecule,
    two long chains (e.g. synthetic ionizable and cationic lipids) or the prostanoid skeleton.
    """
    if chains.count_longer_than(MIN_FATTY_CHAIN) >= 2:
        return True
    if mol.HasSubstructMatch(PATTERNS['prostanoid']):
        return True
    if chains.longest < MIN_FATTY_CHAIN:
        return False
    peg_atoms = {index for match in mol.GetSubstructMatches(PATTERNS['peg']) for index in match}
    chain_atoms = sum(size for size, length in zip(chains.sizes, chains.lengths) if length >= MIN_CORE_CHAIN)
    return chain_atoms / (mol.GetNumHeavyAtoms() - len(peg_atoms)) >= MIN_CHAIN_FRACTION
