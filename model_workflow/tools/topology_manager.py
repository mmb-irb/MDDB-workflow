# Topology management tools
# 
# Set tools for handling topology numeration conversions
# 
# There are 2 different numeration systems to handle in this workflow:
# 
#     Source: The oiginal topology numeration.It is defined by the chain letter, the residue number and, optionally, and icode. Here is defined a specific class called 'sourceResidue' used to define residues with this numeration system. The workflow input and output values use this system (e.g. the analyses labels).
#     Pytraj: The ordinal residue numbers from 1 to 'n'. It is defined just by a string number. Pytraj functions input and output values use this system (e.g. the 'hbonds' function)
# 
# In the following cell there are functions to convert values from one numeration system to the other through a topolgy reference. This reference is set using ProDy.
# 
# In addition other functions are set:
# 
#     Selector to get residues from a string selection

from subprocess import run, PIPE
from collections import OrderedDict
import re
import prody
import pytraj as pt

# DANI: Hay que comprobar que no haya residuos de número repetido sin icode
# DANI: Si esto sucede ProDy no verá los residuos y pytraj sí
# DANI: Análisis de pytraj fallarán por este baile de números (e.g. rmsd-perres)

# This is composed of chain, residue number and optionally an icode
class sourceResidue:
    
    def __init__ (self, chain_letter: str, residue_number: int, icode: str = ''):
        self.chain_letter = chain_letter
        self.residue_number = residue_number
        self.icode = self.valid_icode(icode)

    # Set the residue form the own prody residue
    @classmethod
    def from_prody(cls, residue):
        chain_letter = residue.getChid()
        residue_number = residue.getResnum()
        icode = residue.getIcode()
        return cls(chain_letter, residue_number, icode)

    # Set the residue form a tag (e.g. 'A:1')
    @classmethod
    def from_tag(cls, tag : str):
        match = re.match(r'(\w):(\d*)(\w?)', tag)
        if not match:
            raise SystemExit("ERROR: Tag '" + tag + "' is not valid")
        chain_letter = match.group(1)
        residue_number = int(match.group(2))
        icode = match.group(3)
        return cls(chain_letter, residue_number, icode)
    
    def __str__ (self):
        return self.tag()
    
    def __repr__ (self):
        return self.tag()
    
    def __eq__ (self, other):
        return (self.chain_letter == other.chain_letter and 
                self.residue_number == other.residue_number and
                self.icode == other.icode)
    
    def __hash__(self):
        return hash((self.chain_letter, self.residue_number, self.icode))

    def __lt__(self, other):
        if self.chain_letter == other.chain_letter:
            return self.residue_number < other.residue_number
        return self.chain_letter < other.chain_letter
    
    # Replace the default emtpy icode (' ') by an empty string ('')
    def valid_icode (self, icode):
        if icode == ' ':
            return ''
        return icode

    def tag (self):
        return self.chain_letter + ':' + str(self.residue_number) + self.icode

# Load the whole topology in a standarized and accessible format
class TopologyReference:
    
    def __init__ (self, topology_filename):
        self.topology_filename = topology_filename
        self.topology = prody.parsePDB(topology_filename)
        self.pytraj_topology = pt.load_topology(filename = topology_filename)
        #print(dir(self.topology))
        # Get chains
        self.chains = list(self.topology.iterChains())
        # Get residues
        self.residues = list(self.topology.iterResidues())
        # Fet atoms
        self.atoms = list(self.topology.iterAtoms())
        # Get pytraj residues
        self.pytraj_residues = list(self.pytraj_topology.residues)

    # Get the source formatted name of a specific prody residue
    def get_residue_name (self, residue):
        return sourceResidue.from_prody(residue).tag()
    
    # Transform a topology residue number to the pytraj residue numeration (1, 2, ... n)
    def source2pytraj (self, source : sourceResidue) -> int:
        residx = None
        # Search all reference residues till we find the one that matches our number
        for index, residue in enumerate(self.residues):
            #print(str(index) + ' -> ' + str(residue.getResnum()) + ':' + residue.getChid())
            if (residue.getResnum() == source.residue_number and
                residue.getChid() == source.chain_letter and
                residue.getIcode() == source.icode):
                residx = index
        
        # If there is no previous index we stop here
        if (residx == None):
            return None

        # Otherwise, get data from the recently found index residue
        residue = self.residues[residx]
        residue_number = residue.getResnum()
        residue_name = residue.getResname()[0:3]
        
        # And check that this residue data matches the pytraj residues data
        ptres = self.pytraj_residues[residx]
        if (residue_number == ptres.original_resid and residue_name == ptres.name):
            return residx + 1
            
        # If not, we must iterate over all pytraj residues to find a match
        for index, ptres in enumerate(self.pytraj_residues):
            if (residue_number == ptres.original_resid and residue_name == ptres.name):
                return index + 1
        
        # Return None if there is no match
        return None

    # Transform a pytraj residue numeration (1, 2, ... n) to the topology residue number
    def pytraj2source (self, pytraj : int) -> sourceResidue:
        index = pytraj - 1
        pytrajResidue = self.pytraj_residues[index]
        expectedNumber = pytrajResidue.original_resid
        expectedName = pytrajResidue.name
        #residue = None
        # In the canonical way this index is equivalent to the prody resiude index
        if index < len(self.residues):
            residue = self.residues[index]
            if (residue.getResnum() == expectedNumber and residue.getResname()[0:3] == expectedName):
                return sourceResidue(residue.getChid(), residue.getResnum(), residue.getIcode())
            
        # Pytraj index may not match the prody index in caotic topologies
        # (i.e. when heavy atoms and hydrogen are listed independently)
        # When this happens we can try to find the residue by comparing resnum and resname
        # WARNING: Note that this alternative method is nos sensitive to chains or icodes
        # WARNING: This is because pytraj does not deal with chains or icodes
        for residue in self.residues:
            #print(str(residue.getResnum()) + ' -> ' + str(expectedNumber) + ' / ' + residue.getResname()[0:3] + ' -> ' + expectedName)
            if (residue.getResnum() == expectedNumber and residue.getResname()[0:3] == expectedName):
                return sourceResidue(residue.getChid(), residue.getResnum(), residue.getIcode())
            
        # Return None if there are no results    
        return None

    # Get the standarized residue array from a prody string selection
    # Residues in the array are in 'sourceResidue' format
    def topology_selection (self, selection : str) -> list:
        
        # Each ATOM in the selection
        sel = self.topology.select(selection)

        if not sel:
            print("WARNING: The selection '" + selection + "' matches no atom in the reference topology")
            return []

        # Getting residue chains
        chains = sel.getChids()

        # Getting residue nums
        residue_numbers = sel.getResnums()

        # Getting icodes
        icodes = sel.getIcodes()

        # Joining chains and nums
        residues = [sourceResidue(i,j,k) for i, j, k in zip(chains, residue_numbers, icodes)]
        
        # Get only the unique residues
        residues = list(set(residues))
        
        # Sort residues by chain letter and residue number
        def by_residue(num):
            return num.residue_number
        def by_chain(num):
            return num.chain_letter
        
        residues.sort(key = by_residue)
        residues.sort(key = by_chain)
        
        return residues

    # Set a function to find the absolute atom index in the corrected topology
    def get_atom_index (self, source, atom_name):
        chain = source.chain_letter
        residue = source.residue_number
        icode = source.icode
        for atom in self.topology.iterAtoms():
            if (atom.getChid() == chain
            and atom.getResnum() == residue
            and atom.getIcode() == icode
            and atom.getName() == atom_name):
                return atom.getIndex()

    # Set a function to find the absolute residue index in the corrected topology
    def get_residue_index (self, source):
        chain = source.chain_letter
        residue = source.residue_number
        icode = source.icode
        for r, res in enumerate(self.topology.iterResidues()):
            if (res.getChid() == chain
            and res.getResnum() == residue
            and res.getIcode() == icode):
                return r

    # Set a function to find the absolute residue index in the corrected topology from an absolute atom index
    def get_atom_residue_index (self, atom_index : int) -> int:
        return self.atoms[atom_index].getResindex()

# Set a pair of independent functions to save and recover chains from a pdb file
# WARNING: These functions must be used only when the pdb has not changed in number of atoms

# Get a list with each atom chain from a pdb
def get_chains (pdb_filename : str) -> list:
    pdb = prody.parsePDB(pdb_filename)
    atoms = list(pdb.iterAtoms())
    chains = [ atom.getChid() for atom in atoms ]
    return chains

# Set each atom chain from a pdb from a list
def set_chains (pdb_filename : str, chains : list):
    pdb = prody.parsePDB(pdb_filename)
    atoms = list(pdb.iterAtoms())
    for a, atom in enumerate(atoms):
        atom.setChid(chains[a])
    prody.writePDB(pdb_filename, pdb)