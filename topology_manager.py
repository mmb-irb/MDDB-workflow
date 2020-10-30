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
import prody
import pytraj as pt

# DANI: Hay que comprobar que no haya residuos de número repetido sin icode
# DANI: Si esto sucede ProDy no verá los residuos y pytraj sí
# DANI: Análisis de pytraj fallarán por este baile de números (e.g. rmsd-perres)

# This is composed of chain, residue number and optionally an icode
class sourceResidue:
    
    def __init__ (self, chainLetter: str, residueNumber: int, icode: str = ''):
        self.chainLetter = chainLetter
        self.residueNumber = residueNumber
        self.icode = self.valid_icode(icode)
    
    def __str__ (self):
        return self.chainLetter + ':' + str(self.residueNumber) + self.icode
    
    def __repr__ (self):
        return self.chainLetter + ':' + str(self.residueNumber) + self.icode
    
    def __eq__ (self, other):
        return (self.chainLetter == other.chainLetter and 
                self.residueNumber == other.residueNumber and
                self.icode == other.icode)
    
    def __hash__(self):
        return hash((self.chainLetter, self.residueNumber, self.icode))
    
    # Replace the default emtpy icode (' ') by an empty string ('')
    def valid_icode (self, icode):
        if icode == ' ':
            return ''
        return icode

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
        # Get pytraj residues
        self.pytraj_residues = list(self.pytraj_topology.residues)
    
    # Transform a topology residue number to the pytraj residue numeration (1, 2, ... n)
    def source2pytraj (self, source):
        residx = None
        # Search all reference residues till we find the one that matches our number
        for index, residue in enumerate(self.residues):
            #print(str(index) + ' -> ' + str(residue.getResnum()) + ':' + residue.getChid())
            if (residue.getResnum() == source.residueNumber and
                residue.getChid() == source.chainLetter and
                residue.getIcode() == source.icode):
                residx = index
        
        # If there is no previous index we stop here
        if (residx == None):
            return None

        # Otherwise, get data from the recently found index residue
        residue = self.residues[residx]
        residueNumber = residue.getResnum()
        residueName = residue.getResname()[0:3]
        
        # And check that this residue data matches the pytraj residues data
        ptres = self.pytraj_residues[residx]
        if (residueNumber == ptres.original_resid and residueName == ptres.name):
            return residx + 1
            
        # If not, we must iterate over all pytraj residues to find a match
        for index, ptres in enumerate(self.pytraj_residues):
            if (residueNumber == ptres.original_resid and residueName == ptres.name):
                return index + 1
        
        # Return None if there is no match
        return None

    # Transform a pytraj residue numeration (1, 2, ... n) to the topology residue number
    def pytraj2source (self, pytraj):
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

    # Get the standarized residue array from a string selection
    def topologySelection (self, selection, form = 'source'):
        
        # Each ATOM in the selection
        sel = self.topology.select(selection)

        # Getting residue chains
        chains = sel.getChids()

        # Getting residue nums
        residueNumbers = sel.getResnums()

        # Getting icodes
        icodes = sel.getIcodes()

        # Joining chains and nums
        residues = [sourceResidue(i,j,k) for i, j, k in zip(chains, residueNumbers, icodes)]
        
        # Get only the unique residues
        residues = list(set(residues))
        
        # Sort residues by chain letter and residue number
        def byResidue(num):
            return num.residueNumber
        def byChain(num):
            return num.chainLetter
        
        residues.sort(key = byResidue)
        residues.sort(key = byChain)
        
        return residues

    # Set a function to find the absolute atom index in the corrected topology
    def getAtomIndex (self, source, atom_name):
        chain = source.chainLetter
        residue = source.residueNumber
        icode = source.icode
        for atom in self.topology.iterAtoms():
            if (atom.getChid() == chain
            and atom.getResnum() == residue
            and atom.getIcode() == icode
            and atom.getName() == atom_name):
                return atom.getIndex()