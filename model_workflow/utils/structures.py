# Main handler of the toolbelt
import os
import math
from scipy.special import comb # DANI: Substituye al math.comb porque fué añadido en python 3.8 y nosotros seguimos en 3.7
from bisect import bisect
from typing import Optional, Union, Tuple, List, Generator, Set
Coords = Tuple[float, float, float]

from model_workflow.utils.file import File
from model_workflow.utils.selections import Selection
from model_workflow.utils.vmd_spells import get_vmd_selection_atom_indices, get_covalent_bonds
from model_workflow.utils.mdt_spells import sort_trajectory_atoms
from model_workflow.utils.auxiliar import InputError, is_imported, residue_name_to_letter, otherwise, warn
from model_workflow.utils.constants import SUPPORTED_POLYMER_ELEMENTS, SUPPORTED_ION_ELEMENTS, SUPPORTED_ELEMENTS
from model_workflow.utils.constants import STANDARD_SOLVENT_RESIDUE_NAMES, STANDARD_COUNTER_ION_ATOM_NAMES

import pytraj
# Import these libraries if they are available
# Otherwise proceed without them
# The error is not raised until a function using the missing module is called
try:
    import prody # type: ignore
except:
    pass
try:
    import MDAnalysis
except:
    pass

# Set all available chains according to pdb standards
available_caps = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K',
    'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
available_lows = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k',
    'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z']
available_chains = available_caps + available_lows

# Set the expected number of bonds for each atom according to its element
coherent_bonds_with_hydrogen = {
    'H': { 'min': 1, 'max': 1 },
    'O': { 'min': 1, 'max': 2 },
    'N': { 'min': 1, 'max': 4 },
    'C': { 'min': 2, 'max': 4 },
    'S': { 'min': 1, 'max': 6 },
    'P': { 'min': 2, 'max': 6 },
}
# WARNING: Not deeply tested
coherent_bonds_without_hydrogen = {
    'O': { 'min': 0, 'max': 2 },
    'N': { 'min': 1, 'max': 3 },
    'C': { 'min': 1, 'max': 3 },
    'S': { 'min': 1, 'max': 4 },
    'P': { 'min': 2, 'max': 4 },
}

# Set the limit of residue numbers according to PDB format (4 characters, starts count at 1)
# This means the last number is 9999 and it is equivalent to index 9998
pdb_last_decimal_residue_index = 9998

# An atom
class Atom:
    def __init__ (self,
        name : Optional[str] = None,
        element : Optional[str] = None,
        coords : Optional[Coords] = None,
        ):
        self.name = name
        self.element = element
        self.coords = tuple(coords) if coords else None
        # Set variables to store references to other related instances
        # These variables will be set further by the structure
        self._structure = None
        self._index = None
        self._residue_index = None

    def __repr__ (self):
        return '<Atom ' + self.name + '>'

    def __eq__ (self, other):
        if type(self) != type(other):
            return False
        return self._residue_index == other._residue_index and self.name == other.name

    def __hash__ (self):
        return hash((self.index))

    # The parent structure (read only)
    # This value is set by the structure itself
    def get_structure (self) -> Optional['Structure']:
        return self._structure
    structure = property(get_structure, None, None, "The parent structure (read only)")

    # The residue index according to parent structure residues (read only)
    # This value is set by the structure itself
    def get_index (self) -> Optional[int]:
        return self._index
    index = property(get_index, None, None, "The residue index according to parent structure residues (read only)")

    # The atom residue index according to parent structure residues
    # If residue index is set then make changes in all the structure to make this change coherent
    def get_residue_index (self) -> int:
        return self._residue_index
    def set_residue_index (self, new_residue_index : int):
        # If the new residue index is the current residue index then do nothing
        # WARNING: It is important to stop this here or it could delete a residue which is not to be deleted
        if new_residue_index == self.residue_index:
            return
        # If there is not strucutre yet it means the residue is beeing set before the structure
        # We just save the residue index and wait for the structure to be set
        if not self.structure:
            self._residue_index = new_residue_index
            return
        # Relational indices are updated through a top-down hierarchy
        # Affected residues are the ones to update this atom internal residue index
        current_residue = self.residue
        current_residue.remove_atom(self)
        new_residue = self.structure.residues[new_residue_index]
        new_residue.add_atom(self)
    residue_index = property(get_residue_index, set_residue_index, None, "The atom residue index according to parent structure residues")

    # The atom residue
    # If residue is set then make changes in all the structure to make this change coherent
    def get_residue (self) -> Optional['Residue']:
        # If there is not strucutre yet it means the atom is beeing set before the structure
        # In this case it is not possible to get related residues in the structure
        # It also may happend that the residue is missing because the atom has been set rawly
        if not self.structure or self.residue_index == None:
            return None
        # Get the residue in the structure according to the residue index
        return self.structure.residues[self.residue_index]
    def set_residue (self, new_residue : 'Residue'):
        # Find the new residue index and set it as the atom residue index
        # Note that the residue must be set in the structure already
        new_residue_index = new_residue.index
        if new_residue_index == None:
            raise ValueError('Residue ' + str(new_residue) + ' is not set in the structure')
        self.set_residue_index(new_residue_index)
    residue = property(get_residue, set_residue, None, "The atom residue")

    # The atom chain index according to parent structure chains (read only)
    # In order to change the chain index it must be changed in the atom residue
    def get_chain_index (self) -> Optional[int]:
        # The residue may be missing if the atom has been set rawly
        if not self.residue:
            return None
        return self.residue.chain_index
    chain_index = property(get_chain_index, None, None, "The atom chain index according to parent structure chains (read only)")

    # The atom chain (read only)
    # In order to change the chain it must be changed in the atom residue
    def get_chain (self) -> Optional['Chain']:
        # If there is not strucutre yet it means the atom is beeing set before the structure
        # In this case it is not possible to get related residues in the structure
        # It also may happend that the chain is missing because the atom has been set rawly
        if not self.structure or self.chain_index == None:
            return None
        # Get the chain in the structure according to the chain index
        return self.structure.chains[self.chain_index]
    chain = property(get_chain, None, None, "The atom chain (read only)")

    # Get indices of other atoms in the structure which are covalently bonded to this atom
    def get_bonds (self, skip_ions : bool = False) -> Optional[ List[int] ]:
        if not self.structure:
            raise ValueError('The atom has not a structure defined')
        if self.index == None:
            raise ValueError('The atom has not an index defined')
        bonds = self.structure.bonds[self.index]
        # If the skip ions argument is passed then remove atom indices of supported ion atoms
        if skip_ions:
            bonds = list(set(bonds) - self.structure.ion_atom_indices)
        return bonds

    # Atoms indices of atoms in the structure which are covalently bonded to this atom
    bonds = property(get_bonds, None, None, 'Atoms indices of atoms in the structure which are covalently bonded to this atom')

    # Get bonded atoms
    def get_bonded_atoms (self) -> List['Atom']:
        return [ self.structure.atoms[atom_index] for atom_index in self.bonds ]

    # Generate a selection for this atom
    def get_selection (self) -> 'Selection':
        return Selection([self.index])

    # Make a copy of the current atom
    def copy (self) -> 'Atom':
        atom_copy = Atom(self.name, self.element, self.coords)
        atom_copy._structure = self._structure
        atom_copy._index = self._index
        atom_copy._residue_index = self._residue_index
        return atom_copy

    # Check if this atom meets specific criteria:
    # 1 - it is a carbon
    # 2 - it is connected only to other carbons and hydrogens
    # 3 - it is connected to 1 or 2 carbons
    def is_fatty_candidate (self) -> bool:
        # Ignore non carbon atoms
        if self.element != 'C':
            return False
        # Get bonded atom elements
        bonded_atoms_elements = [ atom.element for atom in self.get_bonded_atoms() ]
        # Check only carbon and hydrogen atoms to be bonded
        if any((element != 'C' and element != 'H') for element in bonded_atoms_elements):
            return False
        # Check it is connected to 1 or 2 carbons
        connected_carbons_count = bonded_atoms_elements.count('C')
        if connected_carbons_count != 1 and connected_carbons_count != 2:
            return False
        return True

    # Check if this atom meets specific criteria:
    # 1 - it is a carbon
    # 2 - it is connected only to other carbons, hydrogens or oxygens
    # 3 - it is connected to 1 or 2 carbons
    # 4 - It is connected to 1 oxygen
    def is_carbohydrate_candidate (self) -> bool:
        # Ignore non carbon atoms
        if self.element != 'C':
            return False
        # Get bonded atom elements
        bonded_atoms_elements = [ atom.element for atom in self.get_bonded_atoms() ]
        # Check only carbon and hydrogen atoms to be bonded
        if any((element != 'C' and element != 'H' and element != 'O') for element in bonded_atoms_elements):
            return False
        # Check it is connected to 1 or 2 carbons
        connected_carbons_count = bonded_atoms_elements.count('C')
        if connected_carbons_count != 1 and connected_carbons_count != 2:
            return False
        # Check it is connected to 1 oxygen
        connected_oxygens_count = bonded_atoms_elements.count('O')
        if connected_oxygens_count != 1:
            return False
        return True

    # Check if it is an ion by checking if it has no bonds with other atoms
    def is_ion (self) -> bool:
        return len(self.bonds) == 0

    # Guess an atom element from its name and number of bonds
    def guess_element (self) -> str:
        # If the name is SOD and it is a lonly atom then it is clearly sodium, not sulfur
        if self.name.upper() == 'SOD' and self.residue.atom_count == 1:
            return 'Na'
        # Find a obvios element name in the atom name
        element = self.get_name_suggested_element()
        # CA may refer to calcium or alpha carbon, so the number of atoms in the residue is decisive
        if element == 'Ca' and self.residue.atom_count >= 2:
            return 'C'
        # NA may refer to sodium or some ligand nitrogen, so the number of atoms in the residue is decisive
        if element == 'Na' and self.residue.atom_count >= 2:
            return 'N'
        return element

    # Guess an atom element from its name only
    def get_name_suggested_element (self) -> str:
        # Get the atom name and its characters length
        name = self.name
        length = len(name)
        next_character = None
        for i, character in enumerate(name):
            # Get the next character, since element may be formed by 2 letters
            if i < length - 1:
                next_character = name[i+1]
                # If next character is not a string then ignore it
                if not next_character.isalpha():
                    next_character = None
            # Try to get all possible matches between the characters and the supported atoms
            # First letter is always caps
            character = character.upper()
            # First try to match both letters together
            if next_character:
                # Start with the second letter in caps
                next_character = next_character.upper()
                both = character + next_character
                if both in SUPPORTED_ELEMENTS:
                    return both
                # Continue with the second letter in lowers
                next_character = next_character.lower()
                both = character + next_character
                if both in SUPPORTED_ELEMENTS:
                    return both
            # Finally, try with the first character alone
            if character in SUPPORTED_ELEMENTS:
                return character
        raise InputError(f"Not recognized element in '{name}'")

    # Get a standard label
    def get_label (self) -> str:
        return f'{self.residue.label}.{self.name}'
    # The atom standard label (read only)
    label = property(get_label, None, None)

# A residue
class Residue:
    def __init__ (self,
        name : Optional[str] = None,
        number : Optional[int] = None,
        icode : Optional[str] = None,
        ):
        self.name = name
        self.number = number
        self.icode = icode
        # Set variables to store references to other related instaces
        # These variables will be set further by the structure
        self._structure = None
        self._index = None
        self._atom_indices = []
        self._chain_index = None
        self._classification = None

    def __repr__ (self):
        return '<Residue ' + self.name + str(self.number) + (self.icode if self.icode else '') + '>'

    def __eq__ (self, other):
        if type(self) != type(other):
            return False
        return (
            self._chain_index == other._chain_index and
            #self.name == other.name and
            self.number == other.number and
            self.icode == other.icode
        )

    def __hash__ (self):
        # WARNING: This is susceptible to duplicated residues
        return hash((self._chain_index, self.number, self.icode))
        # WARNING: This is not susceptible to duplicated residues
        #return hash(tuple(self._atom_indices))

    # The parent structure (read only)
    # This value is set by the structure itself
    def get_structure (self) -> Optional['Structure']:
        return self._structure
    structure = property(get_structure, None, None, "The parent structure (read only)")

    # The residue index according to parent structure residues (read only)
    # This value is set by the structure itself
    def get_index (self) -> Optional[int]:
        return self._index
    def set_index (self, index):
        # Update residue atoms
        for atom in self.atoms:
            atom._residue_index = index
        # Update residue chain
        chain_residue_index = self.chain._residue_indices.index(self._index)
        self.chain._residue_indices[chain_residue_index] = index
        # Finally update self index
        self._index = index
    index = property(get_index, set_index, None, "The residue index according to parent structure residues (read only)")

    # The atom indices according to parent structure atoms for atoms in this residue
    # If atom indices are set then make changes in all the structure to make this change coherent
    def get_atom_indices (self) -> List[int]:
        return self._atom_indices
    def set_atom_indices (self, new_atom_indices : List[int]):
        # If there is not strucutre yet it means the residue is beeing set before the structure
        # We just save atom indices and wait for the structure to be set
        if not self.structure:
            self._atom_indices = new_atom_indices
            return
        # Update the current atoms
        for atom in self.atoms:
            atom._residue_index = None
        # Update the new atoms
        for index in new_atom_indices:
            atom = self.structure.atoms[index]
            atom._residue_index = self.index
        # Now new indices are coherent and thus we can save them
        self._atom_indices = new_atom_indices
    atom_indices = property(get_atom_indices, set_atom_indices, None, "The atom indices according to parent structure atoms for atoms in this residue")

    # The atoms in this residue
    # If atoms are set then make changes in all the structure to make this change coherent
    def get_atoms (self) -> List['Atom']:
        # If there is not strucutre yet it means the chain is beeing set before the structure
        # In this case it is not possible to get related atoms in the structure
        if not self.structure:
            return []
        # Get atoms in the structure according to atom indices
        atoms = self.structure.atoms
        return [ atoms[atom_index] for atom_index in self.atom_indices ]
    def set_atoms (self, new_atoms : List['Atom']):
        # Find indices for new atoms and set their indices as the new atom indices
        # Note that atoms must be set in the structure already
        new_atom_indices = []
        for new_atom in new_atoms:
            new_atom_index = new_atom.index
            if new_atom_index == None:
                raise ValueError('Atom ' + str(new_atom) + ' is not set in the structure')
            new_atom_indices.append(new_atom_index)
        self.set_atom_indices(new_atom_indices)
    atoms = property(get_atoms, set_atoms, None, "The atoms in this residue")

    # Set the number of atoms in the residue
    def get_atom_count (self) -> int:
        return len(self.atom_indices)
    atom_count = property(get_atom_count, None, None, "The number of atoms in the residue (read only)")

    # Add an atom to the residue
    def add_atom (self, new_atom : 'Atom'):
        # Insert the new atom index in the list of atom indices keeping the order
        new_atom_index = new_atom.index
        sorted_atom_index = bisect(self.atom_indices, new_atom_index)
        self.atom_indices.insert(sorted_atom_index, new_atom_index)
        # Update the atom internal index
        new_atom._residue_index = self.index

    # Remove an atom from the residue
    def remove_atom (self, current_atom : 'Atom'):
        # Remove the current atom index from the atom indices list
        self.atom_indices.remove(current_atom.index) # This index MUST be in the list
        # Update the atom internal index
        current_atom._residue_index = None

    # The residue chain index according to parent structure chains
    # If chain index is set then make changes in all the structure to make this change coherent
    def get_chain_index (self) -> int:
        return self._chain_index
    def set_chain_index (self, new_chain_index : int):
        # If the new chain index is the current chain index do nothing
        # WARNING: It is important to stop this here or it could delete a chain which is not to be deleted
        if new_chain_index == self.chain_index:
            return
        # If there is not strucutre yet it means the chain is beeing set before the structure
        # We just save the chain index and wait for the structure to be set
        if not self.structure:
            self._chain_index = new_chain_index
            return
        # Relational indices are updated through a top-down hierarchy
        # Affected chains are the ones to update this residue internal chain index
        # WARNING: It is critical to find the new chain before removing/adding residues
        # WARNING: It may happend that we remove the last residue in the current chain and the current chain is purged
        # WARNING: Thus the 'new_chain_index' would be obsolete since the structure.chains list would have changed
        new_chain = self.structure.chains[new_chain_index]
        if self.chain:
            self.chain.remove_residue(self)
        new_chain.add_residue(self)
    chain_index = property(get_chain_index, set_chain_index, None, "The residue chain index according to parent structure chains")

    # The residue chain
    # If chain is set then make changes in all the structure to make this change coherent
    def get_chain (self) -> Optional['Chain']:
        # If there is not strucutre yet then it means the residue is set before the structure
        # In this case it is not possible to get related chain in the structure
        if not self.structure or self.chain_index == None:
            return None
        # Get the chain in the structure according to the chain index
        return self.structure.chains[self.chain_index]
    def set_chain (self, new_chain : Union['Chain', str]):
        # In case the chain is just a string we must find/create the corresponding chain
        if type(new_chain) == str:
            letter = new_chain
            # Get the residue structure
            structure = self.structure
            if not structure:
                raise ValueError('Cannot find the corresponding ' + new_chain + ' chain without the structure')
            # Find if the letter belongs to an already existing chain
            new_chain = structure.get_chain_by_name(letter)
            # If the chain does not exist yet then create it
            if not new_chain:
                new_chain = Chain(name=letter)
                structure.set_new_chain(new_chain)
        # Find the new chain index and set it as the residue chain index
        # Note that the chain must be set in the structure already
        new_chain_index = new_chain.index
        if new_chain_index == None:
            raise ValueError('Chain ' + str(new_chain) + ' is not set in the structure')
        self.set_chain_index(new_chain_index)
    chain = property(get_chain, set_chain, None, "The residue chain")

    # Get atom indices from atoms bonded to this residue
    def get_bonded_atom_indices (self) -> List[int]:
        # Get all bonds among residue atoms
        all_bonds = []
        for atom in self.atoms:
            all_bonds += atom.bonds
        # Remove self atom indices from the list
        external_bonds = set(all_bonds) - set(self.atom_indices)
        return list(external_bonds)

    # Get atoms bonded to this residue
    def get_bonded_atoms (self) -> List['Atom']:
        return [ self.structure.atoms[atom_index] for atom_index in self.get_bonded_atom_indices() ]

    # Get residue indices from residues bonded to this residue
    def get_bonded_residue_indices (self) -> List[int]:
        return list(set([ atom.residue_index for atom in self.get_bonded_atoms() ]))

    # Get residues bonded to this residue
    def get_bonded_residues (self) -> List['Residue']:
        return [ self.structure.residues[residue_index] for residue_index in self.get_bonded_residue_indices() ]

    # Get the residue biochemistry classification
    # WARNING: Note that this logic will not work in a structure without hydrogens
    # Available classifications:
    # - protein
    # - nucleic
    # - carbohydrate
    # - fatty
    # - steroid
    # - ion
    # - solvent
    # - acetyl
    # - amide
    # - other
    def get_classification (self) -> str:
        # Return the internal value, if any
        if self._classification:
            return self._classification
        # If we dont have a value then we must classify the residue
        # -------------------------------------------------------------------------------------------------------
        # Ions are single atom residues
        if self.atom_count == 1:
            self._classification = 'ion'
            return self._classification
        # -------------------------------------------------------------------------------------------------------
        # At this point we need atoms to have elements fixed
        # WARNING: missing elements would result in silent failure to recognize some classifications
        if self.structure._fixed_atom_elements == False:
            self.structure.fix_atom_elements()
        # Solvent is water molecules
        if self.atom_count == 3:
            atom_elements = [ atom.element for atom in self.atoms ]
            if atom_elements.count('H') == 2 and atom_elements.count('O') == 1:
                self._classification = 'solvent'
                return self._classification
        # -------------------------------------------------------------------------------------------------------
        # Protein definition according to vmd:
        # a residue with atoms named C, N, CA, and O
        # In our case we accept OC1 and OC2 or OT1 and OT2 instead of O for terminal residues
        atom_names = set([ atom.name for atom in self.atoms ])
        if ( all((name in atom_names) for name in ['C', 'N', 'CA']) and (
            'O' in atom_names or { 'OC1', 'OC2' } <= atom_names or { 'OT1', 'OT2' } <= atom_names
        )):
            self._classification = 'protein'
            return self._classification
        # -------------------------------------------------------------------------------------------------------
        # Nucleic acids definition according to vmd:
        # a residue with atoms named P, O1P, O2P and either O3’, C3’, C4’, C5’, O5’ or O3*, C3*, C4*, C5*, O5*
        # Apparently it has been fixed so now a residue does not need to be phosphorylated to be considered nucleic
        # LORE: This included the condition "all((name in atom_names) for name in ['P', 'OP1', 'OP2'])"
        if (( all((name in atom_names) for name in ["O3'", "C3'", "C4'", "C5'", "O5'"]) or
            all((name in atom_names) for name in ["O3*", "C3*", "C4*", "C5*", "O5*"]) )):
            self._classification = 'nucleic'
            return self._classification
        # -------------------------------------------------------------------------------------------------------
        # To define carbohydrates search for rings made of 1 oxygen and 'n' carbons
        # WARNING: This logic may fail for some very specific molecules such as furine
        # LIMITATION: This logic only aims for cyclical carbohydrates. The linear form of carbohydrates is not yet consider
        rings = list(self.find_rings())
        for ring in rings:
            ring_elements = [ atom.element for atom in ring ]
            # Check the ring has 1 oxygen
            oxygen_count = ring_elements.count('O')
            if oxygen_count != 1:
                continue
            # Check the rest are only carbon
            carbon_count = ring_elements.count('C')
            if carbon_count != len(ring) - 1:
                continue
            self._classification = 'carbohydrate'
            return self._classification
        # -------------------------------------------------------------------------------------------------------
        # To define steroids search for 3 6-carbon rings and 1 5-carbon ring (at least)
        # WARNING: According to this logic some exotical non-steroid molecules may result in false positives
        num_6_carbon_rings = 0
        num_5_carbon_rings = 0
        for ring in rings:
            ring_elements = [ atom.element for atom in ring ]
            if any(element != 'C' for element in ring_elements):
                continue
            ring_lenght = len(ring)
            if ring_lenght == 5:
                num_5_carbon_rings += 1
            if ring_lenght == 6:
                num_6_carbon_rings += 1
        if num_6_carbon_rings >= 3 and num_5_carbon_rings >= 1:
            self._classification = 'steroid'
            return self._classification
        # -------------------------------------------------------------------------------------------------------
        # To define fatty acids search for a series of carbon atoms connected together (3 at least)
        # These carbon atoms must be connected only to hydrogen in addition to 1-2 carbon
        # These carbons must not be bonded in a cycle so check atoms not be repeated in the series
        bonded_fatty_atoms = []
        def get_bonded_fatty_atoms_recuersive (current_atom : Atom, previous_atom : Optional[Atom] = None) -> bool:
            # Iterate over the input atom bonded atoms
            for bonded_atom in current_atom.get_bonded_atoms():
                # Discard the previous atom to avoid an infinite loop
                if bonded_atom == previous_atom:
                    continue
                # If we find an atom which is already in the fatty atoms list then it means we are cycling
                if bonded_atom in bonded_fatty_atoms:
                    return False
                # If this is not a fatty atom then discard it
                if not bonded_atom.is_fatty_candidate():
                    continue
                # Add the current bonded fatty atom to the list
                bonded_fatty_atoms.append(bonded_atom)
                # Find more bonded fatty atoms and add them to list as well
                if get_bonded_fatty_atoms_recuersive(current_atom=bonded_atom, previous_atom=current_atom) == False:
                    return False
            return True
        # Now check all atoms and try to set the series until we find one which works
        already_checked_atoms = []
        for atom in self.atoms:
            # If we already checked this atom trough the recursive logic then skip it
            if atom in already_checked_atoms:
                continue
            # If the atom is not a fatty candidate then skip it
            if not atom.is_fatty_candidate():
                continue
            # Now that we have found a suitable candidate atom we start the series
            bonded_fatty_atoms = [atom]
            if get_bonded_fatty_atoms_recuersive(atom) and len(bonded_fatty_atoms) >= 3:
                self._classification = 'fatty'
                return self._classification
            already_checked_atoms += bonded_fatty_atoms
        # -------------------------------------------------------------------------------------------------------
        # Check if it is an acetylation capping terminal
        if self.name == 'ACE':
            self._classification = 'acetyl'
            return self._classification
        # -------------------------------------------------------------------------------------------------------
        # Check if it is an N-methyl amide capping terminal
        if self.name == 'NME':
            self._classification = 'amide'
            return self._classification
        # -------------------------------------------------------------------------------------------------------
        self._classification = 'other'
        return self._classification


    # The residue biochemistry classification (read only)
    classification = property(get_classification, None, None)

    # Generate a selection for this residue
    def get_selection (self) -> 'Selection':
        return Selection(self.atom_indices)

    # Make a copy of the current residue
    def copy (self) -> 'Residue':
        residue_copy = Residue(self.name, self.number, self.icode)
        residue_copy._structure = self._structure
        residue_copy._index = self._index
        residue_copy._atom_indices = self._atom_indices
        residue_copy._chain_index = self._chain_index
        return residue_copy

    # Get the residue equivalent single letter code
    # Note that this makes sense for aminoacids and nucleotides only
    # Non recognized residue names return 'X'
    def get_letter (self) -> str:
        return residue_name_to_letter(self.name)

    # Get the formula of the residue
    def get_formula (self) -> str:
        elements = [ atom.element for atom in self.atoms ]
        unique_elements = list(set(elements))
        formula = ''
        for unique_element in unique_elements:
            count = elements.count(unique_element)
            number = get_lower_numbers(str(count)) if count > 1 else ''
            formula += unique_element + number
        return formula

    # Find rings in the residue
    def find_rings (self) -> Generator[ List[Atom], None, None ]:
        return self.structure.find_rings(selection=self.get_selection())

    # Split this residue in 2 residues and return them in a tuple
    # Keep things coherent in the structure (renumerate all residues below this one)
    # Note that all residue atoms must be covered by the splits
    def split (self,
        first_residue_atom_indices : List[int],
        second_residue_atom_indices : List[int],
        first_residue_name : Optional[str] = None,
        second_residue_name : Optional[str] = None,
        first_residue_number : Optional[int] = None,
        second_residue_number : Optional[int] = None,
        first_residue_icode : Optional[str] = None,
        second_residue_icode : Optional[str] = None,
    ) -> Tuple['Residue', 'Residue']:
        # This function is expected to be called in a residue with an already set structure
        if not self.structure:
            raise InputError('The split function should be called when the residue has an already defined structure')
        # Make sure all atoms in the residue are covered between both the first and the second residue
        if set(self.atom_indices) != set(first_residue_atom_indices + second_residue_atom_indices):
            print('Residue atoms: ' + ', '.join([ str(v) for v in set(self.atom_indices) ]))
            print('Covered atoms: ' + ', '.join([ str(v) for v in set(first_residue_atom_indices + second_residue_atom_indices) ]))
            raise InputError('All atom indices must be covered between both the first and the second residue')
        # Reuse this first residue instead of creating a new one
        if first_residue_name:
            self.name = first_residue_name
        if first_residue_number:
            self.number = first_residue_number
        if first_residue_icode:
            self.icode = first_residue_icode
        # Set the new second residue
        _second_residue_name = second_residue_name if second_residue_name else self.name
        _second_residue_number = second_residue_number if second_residue_number else self.number
        _second_residue_icode = second_residue_icode if second_residue_icode else get_next_letter(self.icode)
        second_residue = Residue(_second_residue_name, _second_residue_number, _second_residue_icode)
        second_residue._structure = self.structure
        new_residue_index = self.index + 1
        second_residue._index = new_residue_index
        # Insert the second residue in the structure residues list right after this residue
        self.structure.residues.insert(new_residue_index, second_residue)
        # Set the second residue index
        # Update the index of all residues which have been shifted after the insertion
        for residue_index in range(new_residue_index + 1, len(self.structure.residues)):
            residue = self.structure.residues[residue_index]
            residue.index = residue_index
        # Now transfer atoms from residue 1 to residue 2 as it is specified
        # Note that this will automatically update every atom
        second_residue.atom_indices = second_residue_atom_indices
        # Now add the new residue to the chain
        self.chain.add_residue(second_residue)

    # Parse atom names to atom indices and then call the split function
    def split_by_atom_names (self,
        first_residue_atom_names : List[str],
        second_residue_atom_names : List[str],
        first_residue_name : Optional[str] = None,
        second_residue_name : Optional[str] = None,
        first_residue_number : Optional[int] = None,
        second_residue_number : Optional[int] = None,
        first_residue_icode : Optional[str] = None,
        second_residue_icode : Optional[str] = None,
    ) -> Tuple['Residue', 'Residue']:
        # Check all atom names to exist in the residue
        input_atom_names = set(first_residue_atom_names + second_residue_atom_names)
        residue_atom_names = set([ atom.name for atom in self.atoms ])
        if input_atom_names != residue_atom_names:
            print(self)
            print(self.atoms)
            print('Input atom names: ' + ', '.join(input_atom_names))
            print('Residue atom names: ' + ', '.join(residue_atom_names))
            raise InputError('All residue atoms must be covered between both the first and the second residue atom names')
        # Convert atom names to atom indices
        first_residue_atom_indices = [ self.get_atom_by_name(name).index for name in first_residue_atom_names ]
        second_residue_atom_indices = [ self.get_atom_by_name(name).index for name in second_residue_atom_names ]
        # Call the actual split logic
        return self.split(first_residue_atom_indices, second_residue_atom_indices, first_residue_name,
            second_residue_name, first_residue_number, second_residue_number, first_residue_icode, second_residue_icode)

    # Get a residue atom given its name
    def get_atom_by_name (self, atom_name : str) -> 'Atom':
        return next(( atom for atom in self.atoms if atom.name == atom_name ), None)

    # Get a standard label
    def get_label (self) -> str:
        chainname = self.chain.name if self.chain.name.strip() else ''
        return f'{chainname}{self.number}{self.icode}({self.name})'
    # The residue standard label (read only)
    label = property(get_label, None, None)

# A chain
class Chain:
    def __init__ (self,
        name : Optional[str] = None,
        ):
        self.name = name
        # Set variables to store references to other related instaces
        # These variables will be set further by the structure
        self._structure = None
        self._index = None
        self._residue_indices = []

    def __repr__ (self):
        return '<Chain ' + self.name + '>'

    def __eq__ (self, other):
        return self.name == other.name
    
    def __hash__ (self):
        return hash(self.name)

    # The parent structure (read only)
    # This value is set by the structure itself
    def get_structure (self) -> Optional['Structure']:
        return self._structure
    structure = property(get_structure, None, None, "The parent structure (read only)")

    # The residue index according to parent structure residues (read only)
    # This value is set by the structure itself
    def get_index (self) -> Optional[int]:
        return self._index
    # When the index is set all residues are updated with the next chain index
    def set_index (self, index : int):
        for residue in self.residues:
            residue._chain_index = index
        self._index = index
    index = property(get_index, set_index, None, "The residue index according to parent structure residues (read only)")

    # The residue indices according to parent structure residues for residues in this chain
    # If residue indices are set then make changes in all the structure to make this change coherent
    def get_residue_indices (self) -> List[int]:
        return self._residue_indices
    def set_residue_indices (self, new_residue_indices : List[int]):
        # If there is not strucutre yet it means the chain is beeing set before the structure
        # We just save residue indices and wait for the structure to be set
        if not self.structure:
            self._residue_indices = new_residue_indices
            return
        # Update the current residues
        for residue in self.residues:
            residue._chain_index = None
        # Update the new residues
        for index in new_residue_indices:
            residue = self.structure.residues[index]
            residue._chain_index = self.index
        # In case the new residue indices list is empty this chain must be removed from its structure
        if len(new_residue_indices) == 0:
            self.structure.purge_chain(self)
        # Now new indices are coherent and thus we can save them
        self._residue_indices = new_residue_indices
    residue_indices = property(get_residue_indices, set_residue_indices, None, "The residue indices according to parent structure residues for residues in this residue")

    # The residues in this chain
    # If residues are set then make changes in all the structure to make this change coherent
    def get_residues (self) -> List['Residue']:
        # If there is not strucutre yet it means the chain is beeing set before the structure
        # In this case it is not possible to get related residues in the structure
        if not self.structure:
            return []
        # Get residues in the structure according to residue indices
        residues = self.structure.residues
        return [ residues[residue_index] for residue_index in self.residue_indices ]
    def set_residues (self, new_residues : List['Residue']):
        # Find indices for new residues and set their indices as the new residue indices
        # Note that residues must be set in the structure already
        new_residue_indices = []
        for new_residue in new_residues:
            new_residue_index = new_residue.index
            if new_residue_index == None:
                raise ValueError('Residue ' + str(new_residue) + ' is not set in the structure')
            new_residue_indices.append(new_residue_index)
        self.set_residue_indices(new_residue_indices)
    residues = property(get_residues, set_residues, None, "The residues in this chain")

    # Add a residue to the chain
    def add_residue (self, residue : 'Residue'):
        # Insert the new residue index in the list of residue indices keeping the order
        sorted_residue_index = bisect(self.residue_indices, residue.index)
        self.residue_indices.insert(sorted_residue_index, residue.index)
        # Update the residue internal chain index
        residue._chain_index = self.index

    # Remove a residue from the chain
    # WARNING: Note that this function does not trigger the set_residue_indices
    def remove_residue (self, residue : 'Residue'):
        self.residue_indices.remove(residue.index) # This index MUST be in the list
        # If we removed the last index then this chain must be removed from its structure
        if len(self.residue_indices) == 0 and self.structure:
            self.structure.purge_chain(self)
        # Update the residue internal chain index
        residue._chain_index = None

    # Atom indices for all atoms in the chain (read only)
    # In order to change atom indices they must be changed in their corresponding residues
    def get_atom_indices (self) -> List[int]:
        return sum([ residue.atom_indices for residue in self.residues ], [])
    atom_indices = property(get_atom_indices, None, None, "Atom indices for all atoms in the chain (read only)")

    # Atoms in the chain (read only)
    # In order to change atoms they must be changed in their corresponding residues
    def get_atoms (self) -> List[int]:
        return sum([ residue.atoms for residue in self.residues ], [])
    atoms = property(get_atoms, None, None, "Atoms in the chain (read only)")

    # Number of atoms in the chain (read only)
    def get_atom_count (self) -> int:
        return len(self.atom_indices)
    atom_count = property(get_atom_count, None, None, "Number of atoms in the chain (read only)")

    # Get the residues sequence in one-letter code
    def get_sequence (self) -> str:
        return ''.join([ residue_name_to_letter(residue.name) for residue in self.residues ])

    # Generate a selection for this chain
    def get_selection (self) -> 'Selection':
        return Selection(self.atom_indices)

    # Make a copy of the current chain
    def copy (self) -> 'Chain':
        chain_copy = Chain(self.name)
        chain_copy._structure = self._structure
        chain_copy._index = self._index
        chain_copy.residue_indices = self.residue_indices
        return chain_copy

# A structure is a group of atoms organized in chains and residues
class Structure:
    def __init__ (self,
        atoms : List['Atom'] = [],
        residues : List['Residue'] = [],
        chains : List['Chain'] = [],
        ):
        self.atoms = []
        self.residues = []
        self.chains = []
        # Set references between instances
        for atom in atoms:
            self.set_new_atom(atom)
        for residue in residues:
            self.set_new_residue(residue)
        for chain in chains:
            self.set_new_chain(chain)
        # --- Set other internal variables ---
        # Set bonds between atoms
        self._bonds = None
        # Set fragments of bonded atoms
        self._fragments = None
        # --- Set other auxiliar variables ---
        # Trajectory atom sorter is a function used to sort coordinates in a trajectory file
        # This function is generated when sorting indices in the structure
        # Also the new atom order is saved
        self.trajectory_atom_sorter = None
        self.new_atom_order = None
        # Track when atom elements have been fixed
        self._fixed_atom_elements = False
        # Save indices of supported ion atoms
        self._ion_atom_indices = None

    def __repr__ (self):
        return f'<Structure ({self.atom_count} atoms)>'

    # The bonds between atoms
    def get_bonds (self) -> List[ List[int] ]:
        # Return the stored value, if exists
        if self._bonds:
            return self._bonds
        # If not, we must calculate the bonds using vmd
        self._bonds = self.get_covalent_bonds()
        return self._bonds
    # Force specific bonds
    def set_bonds (self, bonds : List[ List[int] ]):
        self._bonds = bonds
        # Reset fragments
        self._fragments = None
    bonds = property(get_bonds, set_bonds, None, "The structure bonds")

    # Get the groups of atoms which are covalently bonded
    def get_fragments (self) -> List['Selection']:
        # Return the stored value, if exists
        if self._fragments != None:
            return self._fragments
        # Otherwise, find fragments in all structure atoms
        self._fragments = list(self.find_fragments())
        return self._fragments
    # Fragments of covalently bonded atoms (read only)
    fragments = property(get_fragments, None, None, "The structure fragments (read only)")

    # Find fragments* in a selection of atoms
    # * A fragment is a selection of colvalently bonded atoms
    # All atoms are searched if no selection is provided
    # WARNING: Note that fragments generated from a specific selection may not match the structure fragments
    # A selection including 2 separated regions of a structure fragment will yield 2 fragments
    def find_fragments (self, selection : Optional['Selection'] = None) -> Generator['Selection', None, None]:
        # If there is no selection we consider all atoms
        if not selection:
            selection = self.select_all()
        # Get/Find covalent bonds between atoms in a new object avoid further corruption (deletion) of the original list
        atom_indexed_covalent_bonds = { atom_index: [ *self.bonds[atom_index] ] for atom_index in selection.atom_indices }
        # Group the connected atoms in "fragments"
        while len(atom_indexed_covalent_bonds) > 0:
            start_atom_index, bonds = next(iter(atom_indexed_covalent_bonds.items()))
            del atom_indexed_covalent_bonds[start_atom_index]
            fragment_atom_indices = [ start_atom_index ]
            while len(bonds) > 0:
                # Get the next bond atom and remove it from the bonds list
                next_atom_index = bonds[0]
                bonds.remove(next_atom_index)
                next_bonds = atom_indexed_covalent_bonds.get(next_atom_index, None)
                # If this atom is out of the selection then skip it
                if next_bonds == None:
                    continue
                next_new_bonds = [ bond for bond in next_bonds if bond not in fragment_atom_indices + bonds ]
                bonds += next_new_bonds
                fragment_atom_indices.append(next_atom_index)
                del atom_indexed_covalent_bonds[next_atom_index]
            yield Selection(fragment_atom_indices)

    # Given a selection of atoms, find all whole structure fragments on them
    def find_whole_fragments (self, selection : 'Selection') -> Generator['Selection', None, None]:
        for fragment in self.fragments:
            if selection & fragment:
                yield fragment

    # Name an atom selection depending on the chains it contains
    # This is used for debug purpouses
    def name_selection (self, selection : 'Selection') -> str:
        atoms = [ self.atoms[index] for index in selection.atom_indices ]
        # Count atoms per chain
        atom_count_per_chain = { chain: 0 for chain in self.chains }
        for atom in atoms:
            atom_count_per_chain[atom.chain] += 1
        # Set labels accoridng to the coverage of every chain
        chain_labels = []
        for chain, atom_count in atom_count_per_chain.items():
            if atom_count == 0: continue
            coverage = atom_count / chain.atom_count
            label = f'chain {chain.name}'
            if coverage < 1:
                percent = round(coverage * 1000) / 10
                label += f' ({percent}%)'
            chain_labels.append(label)
        return ', '.join(chain_labels)

    # Set a new atom in the structure
    def set_new_atom (self, atom : 'Atom'):
        atom._structure = self
        new_atom_index = self.atom_count
        self.atoms.append(atom)
        atom._index = new_atom_index

    # Set a new residue in the structure
    # WARNING: Atoms must be set already before setting residues
    def set_new_residue (self, residue : 'Residue'):
        residue._structure = self
        new_residue_index = len(self.residues)
        self.residues.append(residue)
        residue._index = new_residue_index
        # In case the residue has atom indices, set relational indices on each atom
        for atom_index in residue.atom_indices:
            atom = self.atoms[atom_index]
            atom._residue_index = new_residue_index

    # Purge residue from the structure and its chain
    # This can be done only when the residue has no atoms left in the structure
    # Renumerate all residue indices which have been offsetted as a result of the purge
    def purge_residue (self, residue : 'Residue'):
        # Check the residue can be purged
        if residue not in self.residues:
            raise ValueError(str(residue) + ' is not in the structure already')
        if len(residue.atom_indices) > 0:
            raise ValueError(str(residue) + ' is still having atoms and thus it cannot be purged')
        # Get the current index of the residue to be purged
        purged_index = residue.index
        # Residues and their atoms below this index are not to be modified
        # Residues and their atoms over this index must be renumerated
        for affected_residue in self.residues[purged_index+1:]:
            # Chaging the index automatically changes all residues atoms '_residue_index' values
            # Chaging the index automatically changes its corresponding index in residue chain '_residue_indices'
            affected_residue.index -= 1
        # Finally, remove the current residue from the list of residues in the structure
        del self.residues[purged_index]

    # Set a new chain in the structure
    # WARNING: Residues and atoms must be set already before setting chains
    def set_new_chain (self, chain : 'Chain'):
        chain._structure = self
        new_chain_index = len(self.chains)
        self.chains.append(chain)
        chain._index = new_chain_index
        # In case the chain has residue indices, set relational indices on each residue
        for residue_index in chain.residue_indices:
            residue = self.residues[residue_index]
            residue._chain_index = new_chain_index

    # Purge chain from the structure
    # This can be done only when the chain has no residues left in the structure
    # Renumerate all chain indices which have been offsetted as a result of the purge
    def purge_chain (self, chain : 'Chain'):
        # Check the chain can be purged
        if chain not in self.chains:
            raise ValueError('Chain ' + chain.name + ' is not in the structure already')
        if len(chain.residue_indices) > 0:
            raise ValueError('Chain ' + chain.name + ' is still having residues and thus it cannot be purged')
        # Get the current index of the chain to be purged
        purged_index = chain.index
        # Chains and their residues below this index are not to be modified
        # Chains and their residues over this index must be renumerated
        for affected_chain in self.chains[purged_index+1:]:
            # Chaging the index automatically changes all chain residues '_chain_index' values
            affected_chain.index -= 1
        # Finally, remove the current chain from the list of chains in the structure
        del self.chains[purged_index]

    # Set the structure from a ProDy topology
    @classmethod
    def from_prody (cls, prody_topology):
        parsed_atoms = []
        parsed_residues = []
        parsed_chains = []
        prody_atoms = list(prody_topology.iterAtoms())
        prody_residues = list(prody_topology.iterResidues())
        prody_chains = list(prody_topology.iterChains())
        # Parse atoms
        for prody_atom in prody_atoms:
            name = prody_atom.getName()
            element = prody_atom.getElement()
            coords = tuple(prody_atom.getCoords())
            parsed_atom = Atom(name=name, element=element, coords=coords)
            parsed_atoms.append(parsed_atom)
        # Parse residues
        for prody_residue in prody_residues:
            name = prody_residue.getResname()
            number = int(prody_residue.getResnum())
            icode = prody_residue.getIcode()
            parsed_residue = Residue(name=name, number=number, icode=icode)
            atom_indices = [ int(index) for index in prody_residue.getIndices() ]
            parsed_residue.atom_indices = atom_indices
            parsed_residues.append(parsed_residue)
        # Parse chains
        for prody_chain in prody_chains:
            name = prody_chain.getChid()
            parsed_chain = Chain(name=name)
            residue_indices = [ int(residue.getResindex()) for residue in prody_chain.iterResidues() ]
            parsed_chain.residue_indices = residue_indices
            parsed_chains.append(parsed_chain)
        return cls(atoms=parsed_atoms, residues=parsed_residues, chains=parsed_chains)

    # Set the structure from a pdb file
    @classmethod
    def from_pdb (cls, pdb_content : str):
        # Read the pdb content line by line and set the parsed atoms, residues and chains
        parsed_atoms = []
        parsed_residues = []
        parsed_chains = []
        atom_index = -1
        residue_index = -1
        last_chain_name = None
        last_residue_name = None
        last_residue_number = None
        last_residue_icode = None
        # Set whenever we have reached the limit of residue numbers
        # Thus further residue numbers must be parsed as hexadecimals
        hexadecimal_residue_numbers = False
        pdb_lines = pdb_content.split('\n')
        for line in pdb_lines:
            # Parse atoms only
            start = line[0:6]
            is_atom = start == 'ATOM  ' or start == 'HETATM'
            if not is_atom:
                continue
            # Mine all atom data
            atom_name = line[11:16].strip()
            residue_name = line[17:21].strip()
            chain = line[21:22]
            residue_number = line[22:26]
            icode = line[26:27]
            if icode == ' ':
                icode = ''
            x_coord = float(line[30:38])
            y_coord = float(line[38:46])
            z_coord = float(line[46:54])
            element = line[76:78].strip()
            # Set the parsed atom, residue and chain
            parsed_atom = Atom(name=atom_name, element=element, coords=(x_coord, y_coord, z_coord))
            # Add the parsed atom to the list and update the current atom index
            parsed_atoms.append(parsed_atom)
            atom_index += 1
            # Check if this atom belongs to the same chain than the previous atom
            same_chain = last_chain_name == chain
            # Check if this atom belongs to the same residue than the previous atom
            same_residue = same_chain and last_residue_name == residue_name and last_residue_number == residue_number and last_residue_icode == icode
            # Update last values
            last_chain_name = chain
            last_residue_name = residue_name
            last_residue_number = residue_number
            last_residue_icode = icode
            # Update the residue atom indices
            # If the residue equals the last parsed residue then use the previous instead
            if same_residue:
                parsed_residue = parsed_residues[-1]
                parsed_residue.atom_indices.append(atom_index)
                # If it is the same residue then it will be the same chain as well so we can proceed
                continue
            # Otherwise, include the new residue in the list and update the current residue index
            # When residue index is greater than the limit we start using the hexadecimal parsing
            # Note that residue index 9998 equals residue number 9999 which is the last allowed
            if residue_index == pdb_last_decimal_residue_index:
                hexadecimal_residue_numbers = True
            # Parse the residue number
            parsed_residue_number = int(residue_number, 16) if hexadecimal_residue_numbers else int(residue_number)
            # Now parse the residue and add it to the list
            parsed_residue = Residue(name=residue_name, number=parsed_residue_number, icode=icode)
            parsed_residues.append(parsed_residue)
            residue_index += 1
            # Add current atom to the new residue
            parsed_residue.atom_indices.append(atom_index)
            # If the chain equals the last parsed chain then use the previous instead
            if same_chain:
                parsed_chain = parsed_chains[-1]
                parsed_chain.residue_indices.append(residue_index)
                continue
            # Otherwise, parse the chain and include the new chain in the list
            parsed_chain = Chain(name=chain)
            parsed_chains.append(parsed_chain)
            # Add current atom to the new chain
            parsed_chain.residue_indices.append(residue_index)
        return cls(atoms=parsed_atoms, residues=parsed_residues, chains=parsed_chains)

    # Set the structure from a pdb file
    @classmethod
    def from_pdb_file (cls, pdb_filename : str):
        pdb_file = File(pdb_filename)
        if not pdb_file.exists:
            raise InputError(f'File "{pdb_filename}" not found')
        if not pdb_file.format == 'pdb':
            raise InputError(f'"{pdb_filename}" is not a name for a pdb file')
        # Read the pdb file
        pdb_content = None
        with open(pdb_filename, 'r') as file:
            pdb_content = file.read()
        return cls.from_pdb(pdb_content)

    # Set the structure from an MD analysis object
    @classmethod
    def from_mdanalysis (cls, mdanalysis_universe):
        # Set the final list of atoms to be included in the structure
        parsed_atoms = []
        parsed_residues = []
        parsed_chains = []
        # Setup atoms
        for atom in mdanalysis_universe.atoms:
            name = atom.name
            element = atom.element
            #coords = atom.position # DANI: Esto da error si no hay coordenadas
            parsed_atom = Atom(name=name, element=element)
            parsed_atoms.append(parsed_atom)
        # Setup residues
        for residue in mdanalysis_universe.residues:
            name = residue.resname
            number = residue.resnum
            icode = residue.icode if hasattr(residue, 'icode') else None # DANI: No se ha provado
            parsed_residue = Residue(name=name, number=number, icode=icode)
            atom_indices = [ atom.index for atom in residue.atoms ]
            parsed_residue.atom_indices = atom_indices
            parsed_residues.append(parsed_residue)
        # Setup chains
        for segment in mdanalysis_universe.segments:
            name = segment.segid
            parsed_chain = Chain(name=name)
            residue_indices = [ residue.resindex for residue in segment.residues ]
            parsed_chain.residue_indices = residue_indices
            parsed_chains.append(parsed_chain)
        # Setup the structure
        structure = cls(atoms=parsed_atoms, residues=parsed_residues, chains=parsed_chains)
        # Add bonds and charges which are also available in a topology
        atom_count = len(mdanalysis_universe.atoms)
        atom_bonds = [ [] for i in range(atom_count) ]
        for bond in mdanalysis_universe.bonds:
            a,b = bond.indices
            # Make sure atom indices are regular integers so they are JSON serializables
            atom_bonds[a].append(int(b))
            atom_bonds[b].append(int(a))
        structure.bonds = atom_bonds
        structure.charges = list(mdanalysis_universe._topology.charges.values)
        return structure

    # Set the structure from a prmtop file
    @classmethod
    def from_prmtop_file (cls, prmtop_filepath : str):
        # Make sure the input file exists and has the right format
        prmtop_file = File(prmtop_filepath)
        if not prmtop_file.exists:
            raise InputError(f'File "{prmtop_filepath}" not found')
        if not prmtop_file.format == 'prmtop':
            raise InputError(f'"{prmtop_filepath}" is not a name for a prmtop file')
        # In we do not have mdanalysis in our environment then we cannot proceed
        if not is_imported('MDAnalysis'):
            raise InputError('Missing dependency error: MDAnalysis')
        # Parse the topology using MDanalysis
        parser = MDAnalysis.topology.TOPParser.TOPParser(prmtop_filepath)
        mdanalysis_topology = parser.parse()
        mdanalysis_universe = MDAnalysis.Universe(mdanalysis_topology)
        return cls.from_mdanalysis(mdanalysis_universe)

     # Set the structure from a tpr file
    @classmethod
    def from_tpr_file (cls, tpr_filepath : str):
        # Make sure the input file exists and has the right format
        tpr_file = File(tpr_filepath)
        if not tpr_file.exists:
            raise InputError(f'File "{tpr_filepath}" not found')
        if not tpr_file.format == 'tpr':
            raise InputError(f'"{tpr_filepath}" is not a name for a tpr file')
        # In we do not have mdanalysis in our environment then we cannot proceed
        if not is_imported('MDAnalysis'):
            raise InputError('Missing dependency error: MDAnalysis')
        # Parse the topology using MDanalysis
        parser = MDAnalysis.topology.TPRParser.TPRParser(tpr_filepath)
        mdanalysis_topology = parser.parse()
        mdanalysis_universe = MDAnalysis.Universe(mdanalysis_topology)
        return cls.from_mdanalysis(mdanalysis_universe)

    # Set the structure from a file if the file format is supported
    @classmethod
    def from_file (cls, mysterious_filepath : str):
        mysterious_file = File(mysterious_filepath)
        if mysterious_file.format == 'pdb':
            return cls.from_pdb_file(mysterious_file.path)
        if mysterious_file.format == 'prmtop':
            return cls.from_prmtop_file(mysterious_file.path)
        if mysterious_file.format == 'tpr':
            return cls.from_tpr_file(mysterious_file.path)
        raise InputError(f'Not supported format ({mysterious_file.format}) to setup a Structure')

    # Set the number of atoms in the structure
    def get_atom_count (self) -> int:
        return len(self.atoms)
    atom_count = property(get_atom_count, None, None, "The number of atoms in the structure (read only)")

    # Fix atom elements by guessing them when missing
    # Set all elements with the first letter upper and the second (if any) lower
    # Also check if atom elements are coherent with atom names
    # If 'trust' is set as False then we impose elements according to what we can guess from the atom name
    # Return True if any element was modified or False if not
    def fix_atom_elements (self, trust : bool = True) -> bool:
        modified = False
        # Save the wrong guesses for a final report
        # This way we do not crowd the terminal with logs when a lot of atoms are affected
        reports = {}
        for atom in self.atoms:
            # Make sure elements have the first letter cap and the second letter not cap
            if atom.element:
                new_element = first_cap_only(atom.element)
                if atom.element != new_element:
                    atom.element = new_element
                    modified = True
                # Check the element to match what we would guess from the atom name
                # In case it does not just warn the user
                guess = atom.guess_element()
                if atom.element != guess:
                    report = (atom.name, atom.element, guess)
                    report_count = reports.get(report, 0)
                    reports[report] = report_count + 1
                    if not trust:
                        atom.element = guess
                        modified = True
            # If elements are missing then guess them from atom names
            else:
                atom.element = atom.guess_element()
                modified = True
        # Warn the user about anormalies
        for report, count in reports.items():
            atom_name, atom_element, guess = report
            warn(f"Suspicious element for atom {atom_name}: {atom_element} -> shoudn't it be {guess}? ({count} occurrences)")
        # Warn the user that some elements were modified
        if modified:
            warn('Atom elements have been modified')
        # Set atom elements as fixed in order to avoid repeating this process
        self._fixed_atom_elements = True
        return modified

    # Set new coordinates
    def set_new_coordinates (self, new_coordinates : List[Coords]):
        # Make sure the list of coordinates is as long as the number of atoms
        if len(new_coordinates) != self.atom_count:
            raise ValueError(f'The number of coordinates ({len(new_coordinates)}) does not match the number of atoms ({self.atom_count})')
        # Overwrite current coordinates with the new coordinates
        for i, atom in enumerate(self.atoms):
            atom.coords = tuple(new_coordinates[i])
    
    # Get all supported ion atom indices together in a set
    def get_ion_atom_indices (self) -> Set:
        # If we already did this then return the stored value
        if self._ion_atom_indices != None:
            return self._ion_atom_indices
        # Find ion atom indices
        indices = set()
        for atom in self.atoms:
            if atom.element in SUPPORTED_ION_ELEMENTS:
                indices.add(atom.index)
        self._ion_atom_indices = indices
        return self._ion_atom_indices
    ion_atom_indices = property(get_ion_atom_indices, None, None, "Atom indices for what we consider supported ions")

    # Generate a pdb file with current structure
    def generate_pdb_file (self, pdb_filename : str):
        with open(pdb_filename, "w") as file:
            file.write('REMARK workflow generated pdb file\n')
            for a, atom in enumerate(self.atoms):
                residue = atom.residue
                index = str((a+1) % 100000).rjust(5)
                name = ' ' + atom.name.ljust(3) if len(atom.name) < 4 else atom.name
                residue_name = residue.name.ljust(4) if residue else 'XXX'.ljust(4)
                chain = atom.chain
                chain_name = chain.name if chain.name and len(chain.name) == 1 else 'X'
                residue_number = str(residue.number).rjust(4) if residue else '0'.rjust(4)
                # If residue number is longer than 4 characters then we must parse to hexadecimal
                if len(residue_number) > 4:
                    residue_number = hex(residue.number)[2:].rjust(4)
                icode = residue.icode if residue.icode and len(residue.icode) else ' '
                # Make sure we have atom coordinates
                if atom.coords == None:
                    raise InputError('Trying to write a PDB file from a structure with atoms without coordinates')
                x_coord, y_coord, z_coord = [ "{:.3f}".format(coord).rjust(8) for coord in atom.coords ]
                occupancy = '1.00' # Just a placeholder
                temp_factor = '0.00' # Just a placeholder
                element = atom.element.rjust(2)
                atom_line = ('ATOM  ' + index + ' ' + name + ' ' + residue_name
                    + chain_name + residue_number + icode + '   ' + x_coord + y_coord + z_coord
                    + '  ' + occupancy + '  ' + temp_factor + '          ' + element).ljust(80) + '\n'
                file.write(atom_line)

    # Get the structure equivalent prody topology
    def get_prody_topology (self):
        # In we do not have prody in our environment then we cannot proceed
        if not is_imported('prody'):
            raise InputError('Missing dependency error: prody')
        # Generate the prody topology
        pdb_filename = '.structure.pdb'
        self.generate_pdb_file(pdb_filename)
        prody_topology = prody.parsePDB(pdb_filename)
        os.remove(pdb_filename)
        return prody_topology

    # Get the structure equivalent pytraj topology
    def get_pytraj_topology (self):
        # In we do not have pytraj in our environment then we cannot proceed
        if not is_imported('pytraj'):
            raise InputError('Missing dependency error: pytraj')
        # Generate a pdb file from the current structure to feed pytraj
        pdb_filename = '.structure.pdb'
        self.generate_pdb_file(pdb_filename)
        pytraj_topology = pytraj.load_topology(filename = pdb_filename)
        os.remove(pdb_filename)
        return pytraj_topology

    # Select atoms from the structure thus generating an atom indices list
    # Different tools may be used to make the selection:
    # - vmd (default)
    # - prody
    # - pytraj
    def select (self, selection_string : str, syntax : str = 'vmd') -> Optional['Selection']:
        if syntax == 'vmd':
            # Generate a pdb for vmd to read it
            pdb_filename = '.structure.pdb'
            self.generate_pdb_file(pdb_filename)
            # Use vmd to find atom indices
            atom_indices = get_vmd_selection_atom_indices(pdb_filename, selection_string)
            os.remove(pdb_filename)
            if len(atom_indices) == 0:
                return None
            return Selection(atom_indices)
        if syntax == 'prody':
            # In we do not have prody in our environment then we cannot proceed
            if not is_imported('prody'):
                raise InputError('Missing dependency error: prody')
            prody_topology = self.get_prody_topology()
            prody_selection = prody_topology.select(selection_string)
            if not prody_selection:
                return None
            return Selection.from_prody(prody_selection)
        if syntax == 'pytraj':
            # In we do not have pytraj in our environment then we cannot proceed
            if not is_imported('pytraj'):
                raise InputError('Missing dependency error: pytraj')
            pytraj_topology = self.get_pytraj_topology()
            pytraj_selection = pytraj_topology[selection_string]
            atom_indices = [ atom.index for atom in pytraj_selection.atoms ]
            if len(atom_indices) == 0:
                return None
            return Selection(atom_indices)

        raise InputError('Syntax ' + syntax + ' is not supported. Choose one of the following: vmd, prody, pytraj')

    # Set a function to make selections using atom indices
    def select_atom_indices (self, atom_indices : List[int]) -> 'Selection':
        # Check atom indices to be in the structure
        atom_count = self.atom_count
        for atom_index in atom_indices:
            if atom_index >= atom_count:
                raise InputError(f'Atom index {atom_index} is out of range ({atom_count})')
        return Selection(atom_indices)

    # Set a function to make selections using residue indices
    def select_residue_indices (self, residue_indices : List[int]) -> 'Selection':
        # WARNING: The following line gets stucked sometimes, idk why
        #atom_indices = sum([ self.residues[index].atom_indices for index in residue_indices ], [])
        atom_indices = []
        for i, index in enumerate(residue_indices):
            atom_indices += self.residues[index].atom_indices
        return Selection(atom_indices)

    # Get a selection with all atoms
    def select_all (self) -> 'Selection':
        return Selection(list(range(self.atom_count)))

    # Select water atoms
    # WARNING: This logic is a bit guessy and it may fail for non-standard residue named structures
    def select_water (self) -> 'Selection':
        water_residues_indices = [ residue.index for residue in self.residues if residue.name in STANDARD_SOLVENT_RESIDUE_NAMES ]
        return self.select_residue_indices(water_residues_indices)

    # Select counter ion atoms
    # WARNING: This logic is a bit guessy and it may fail for non-standard atom named structures
    def select_counter_ions (self) -> 'Selection':
        counter_ion_indices = []
        for atom in self.atoms:
            # If the residue has not one single atom then it is not an ion
            if len(atom.residue.atoms) != 1:
                continue
            # Get a simplified version of the atom name
            # Set all letters upper and remove non-letter characters (e.g. '+' and '-' symbols)
            simple_atom_name = ''.join(filter(str.isalpha, atom.name.upper()))
            if simple_atom_name in STANDARD_COUNTER_ION_ATOM_NAMES:
                counter_ion_indices.append(atom.index)
        return Selection(counter_ion_indices)

    # Select heavy atoms
    def select_heavy_atoms (self) -> 'Selection':
        atom_indices = []
        for atom in self.atoms:
            # If the atom is not an hydrogen then add it to the list
            if atom.element != 'H':
                atom_indices.append(atom.index)
        return Selection(atom_indices)

    # Select protein atoms
    # WARNING: Note that there is a small difference between VMD protein and our protein
    # This function is improved to consider terminal residues as proteins
    def select_protein (self) -> 'Selection':
        atom_indices = []
        for residue in self.residues:
            if residue.classification == 'protein':
                atom_indices += residue.atom_indices
        return Selection(atom_indices)

    # Select cartoon representable regions for VMD
    # Rules are:
    # 1. Residues must be protein (i.e. must contain C, CA, N and O atoms) or nucleic (P, OP1, OP2, O3', C3', C4', C5', O5')
    # 2. There must be at least 3 covalently bonded residues
    # It does not matter their chain, numeration or even index order as long as they are bonded
    # * Note that we can represent cartoon while we display one residue alone, but it must be connected anyway
    # Also, we have the option to include terminals in the cartoon selection although they are not representable
    # This is helpful for the screenshot: terminals are better hidden than represented as ligands
    def select_cartoon (self, include_terminals : bool = False) -> 'Selection':
        # Set fragments which are candidate to be cartoon representable
        fragments = []
        # Get protein fragments according to VMD
        protein_selection = self.select_protein() if include_terminals else self.select('protein', syntax='vmd')
        if protein_selection:
            fragments += list(self.find_fragments(protein_selection))
        # Get nucleic fragments according to VMD
        nucleic_selection = self.select('nucleic', syntax='vmd')
        if nucleic_selection:
            fragments += list(self.find_fragments(nucleic_selection))
        # Set the final selection including all valid fragments
        cartoon_selection = Selection()
        # Iterate over every fragment
        for fragment in fragments:
            # Make sure the fragment has at least 3 residues before adding it to the cartoon selection
            fragment_residue_indices = self.get_selection_residue_indices(fragment)
            if len(fragment_residue_indices) >= 3:
                cartoon_selection += fragment
        return cartoon_selection

    # Invert a selection
    def invert_selection (self, selection : 'Selection') -> 'Selection':
        atom_indices = list(range(self.atom_count))
        for atom_index in selection.atom_indices:
            atom_indices[atom_index] = None
        return Selection([ atom_index for atom_index in atom_indices if atom_index != None ])

    # Get protein selection
    # WARNING: Do not rely in VMD's "protein" selection since it misses terminal residues
    # VMD considers protein any resiude including N, C, CA and O while terminals may have OC1 and OC2 instead of O
    def get_protein_selection (self) -> 'Selection':
        atom_indices = []
        for residue in self.residues:
            if residue.classification == 'protein':
                atom_indices += residue.atom_indices
        return Selection(atom_indices)
    
    # Given a selection, get a list of residue indices for residues implicated
    # Note that if a single atom from the residue is in the selection then the residue index is returned
    def get_selection_residue_indices (self, selection : 'Selection') -> List[int]:
        return list(set([ self.atoms[atom_index].residue_index for atom_index in selection.atom_indices ]))
    
    # Given a selection, get a list of chain indices for chains implicated
    # Note that if a single atom from the chain is in the selection then the chain index is returned
    def get_selection_chain_indices (self, selection : 'Selection') -> List[int]:
        return list(set([ self.atoms[atom_index].chain_index for atom_index in selection.atom_indices ]))

    # Create a new structure from the current using a selection to filter atoms
    def filter (self, selection : Union['Selection', str], selection_syntax : str = 'vmd') -> 'Structure':
        if not selection:
            raise InputError('No selection was passed')
        # In case the selection is not an actual Selection, but a string, parse the string into a Selection
        if type(selection) == str:
            selection = self.select(selection, selection_syntax)
        new_atoms = []
        new_residues = []
        new_chains = []
        # Get the selected atoms
        # Generate dictionaries with new indexes as keys and previous indexes as values for atoms, residues and chains
        # This is done with this structure for the residues and chains further to find the new indices fast
        new_atom_indices = {}
        # Collect also original indices to related atom residues and chains
        original_atom_residue_indices = []
        for new_index, original_index in enumerate(selection.atom_indices):
            # Make a copy of the selected atom in order to not modify the original one
            original_atom = self.atoms[original_index]
            new_atom = Atom(
                name=original_atom.name,
                element=original_atom.element,
                coords=original_atom.coords,
            )
            new_atoms.append(new_atom)
            # Save old and new indices in order to reconfigure all indices later
            new_atom_indices[original_index] = new_index
            original_atom_residue_indices.append(original_atom.residue_index)
        # Find the selected residues
        selected_residue_indices = list(set(original_atom_residue_indices))
        # Repeat the original/new indices backup we did before
        new_residue_indices = {}
        original_residue_atom_indices = []
        original_residue_chain_indices = []
        for new_index, original_index in enumerate(selected_residue_indices):
            # Make a copy of the selected residue in order to not modify the original one
            original_residue = self.residues[original_index]
            new_residue = Residue(
                name=original_residue.name,
                number=original_residue.number,
                icode=original_residue.icode,
            )
            new_residues.append(new_residue)
            # Save old and new indices in order to reconfigure all indices later
            new_residue_indices[original_index] = new_index
            original_residue_atom_indices.append(original_residue.atom_indices)
            original_residue_chain_indices.append(original_residue.chain_index)
        # Find the selected chains
        selected_chain_indices = list(set(original_residue_chain_indices))
        # Repeat the original/new indices backup we did before
        new_chain_indices = {}
        original_chain_residue_indices = []
        for new_index, original_index in enumerate(selected_chain_indices):
            # Make a copy of the selected chain in order to not modify the original one
            original_chain = self.chains[original_index]
            new_chain = Chain(
                name=original_chain.name,
            )
            new_chains.append(new_chain)
            # Save old and new indices in order to reconfigure all indices later
            new_chain_indices[original_index] = new_index
            original_chain_residue_indices.append(original_chain.residue_indices)
        # Finally, reset indices in all instances
        for a, atom in enumerate(new_atoms):
            atom.residue_index = new_residue_indices[ original_atom_residue_indices[a] ]
        for r, residue in enumerate(new_residues):
            atom_indices = []
            for original_index in original_residue_atom_indices[r]:
                new_index = new_atom_indices.get(original_index, None)
                if new_index != None:
                    atom_indices.append(new_index)
            residue.atom_indices = atom_indices
            residue.chain_index = new_chain_indices[ original_residue_chain_indices[r] ]
        for c, chain in enumerate(new_chains):
            residue_indices = []
            for original_index in original_chain_residue_indices[c]:
                new_index = new_residue_indices.get(original_index, None)
                if new_index != None:
                    residue_indices.append(new_index)
            chain.residue_indices = residue_indices
        return Structure(atoms=new_atoms, residues=new_residues, chains=new_chains)

    # Set chains on demand
    # If no selection is passed then the whole structure will be affected
    # If no chain is passed then a "chain by fragment" logic will be applied
    def chainer (self, selection : Optional['Selection'] = None, letter : Optional[str] = None, whole_fragments : bool = True):
        # If there is no selection we consider all atoms
        if selection == None:
            selection = self.select_all()
        # If the selection is empty then there is nothing to do here
        if len(selection) == 0:
            return
        # If a letter is specified then the logic is way simpler
        if letter:
            self.set_selection_chain_name(selection, letter)
            return
        # If a letter is not specified we run the "fragments" logic
        fragment_getter = self.find_whole_fragments if whole_fragments else self.find_fragments
        for fragment in fragment_getter(selection):
            chain_name = self.get_next_available_chain_name()
            self.set_selection_chain_name(fragment, chain_name)

    # Smart function to set chains automatically
    # Original chains will be overwritten
    def auto_chainer (self):
        # Reset chains
        self.chainer(letter='X')
        # Set solvent and ions as a unique chain
        ion_residue_indices = [ residue.index for residue in self.residues if residue.classification == 'ion' ]
        ion_selection = self.select_residue_indices(ion_residue_indices)
        solvent_residue_indices = [ residue.index for residue in self.residues if residue.classification == 'solvent' ]
        solvent_selection = self.select_residue_indices(solvent_residue_indices)
        ion_and_indices_selection = ion_selection + solvent_selection
        self.chainer(selection=ion_and_indices_selection, letter='S')
        # Set fatty acids and steroids as a unique chain
        # DANI: Se podrían incluir algunos carbohydrates (glycans)
        # DANI: Se podrían descartarían residuos que no pertenezcan a la membrana por proximidad
        fatty_residue_indices = [ residue.index for residue in self.residues if residue.classification == 'fatty' ]
        fatty_selection = self.select_residue_indices(fatty_residue_indices)
        steroid_residue_indices = [ residue.index for residue in self.residues if residue.classification == 'steroid' ]
        steroid_selection = self.select_residue_indices(steroid_residue_indices)
        membrane_selection = fatty_selection + steroid_selection
        self.chainer(selection=membrane_selection, letter='M')
        # Set carbohydrates as a unique chain as well, just in case we have glycans
        # Note that in case glycan atoms are mixed with protein atoms glycan chains will be overwritten
        # However this is not a problem. It is indeed the best solution if we don't want ro resort atoms
        carbohydrate_residue_indices = [ residue.index for residue in self.residues if residue.classification == 'carbohydrate' ]
        carbohydrate_selection = self.select_residue_indices(carbohydrate_residue_indices)
        self.chainer(selection=carbohydrate_selection, letter='H')
        # Add a chain per fragment for both protein and nucleic acids
        protein_residue_indices = [ residue.index for residue in self.residues if residue.classification == 'protein' ]
        protein_selection = self.select_residue_indices(protein_residue_indices)
        nucleic_residue_indices = [ residue.index for residue in self.residues if residue.classification == 'nucleic' ]
        nucleic_selection = self.select_residue_indices(nucleic_residue_indices)
        protein_and_nucleic_selection = protein_selection + nucleic_selection
        self.chainer(selection=protein_and_nucleic_selection)
        # At this point we should have covered most of the molecules in the structure
        # However, in case there are more molecules, we have already set them all as a single chain ('X')
        # Here we do not apply the chain per fragment logic since it may be dangerous
        # Note that we may have a lot of independent residues (and thus a los of small fragments)
        # This would make us run out of letters in the alphabet and thus there would be no more chains
        # As the last step, fix repeated chains
        self.check_repeated_chains(fix_chains=True)

    # This is an alternative system to find protein chains (anything else is chained as 'X')
    # This system does not depend on VMD
    # It totally overrides previous chains since it is expected to be used only when chains are missing
    def raw_protein_chainer (self):
        current_chain = self.get_next_available_chain_name()
        previous_alpha_carbon = None
        for residue in self.residues:
            alpha_carbon = next((atom for atom in residue.atoms if atom.name == 'CA'), None)
            if not alpha_carbon:
                residue.set_chain('X')
                continue
            # Connected aminoacids have their alpha carbons at a distance of around 3.8 Ångstroms
            residues_are_connected = previous_alpha_carbon and calculate_distance(previous_alpha_carbon, alpha_carbon) < 4
            if not residues_are_connected:
                current_chain = self.get_next_available_chain_name()
            residue.set_chain(current_chain)
            previous_alpha_carbon = alpha_carbon

    # Given an atom selection, set the chain for all these atoms
    # Note that the chain is changed in every whole residue, no matter if only one atom was selected
    def set_selection_chain_name (self, selection : 'Selection', letter : str):
        # Find if the letter belongs to an already existing chain
        chain = self.get_chain_by_name(letter)
        # If the chain does not exist yet then create it
        if not chain:
            chain = Chain(name=letter)
            self.set_new_chain(chain)
        # Convert the selection to residues
        # WARNING: Remember to never use a set when handling residues or you may eliminate 'duplicated' residues
        residue_indices = self.get_selection_residue_indices(selection)
        # Set the chain index of every residue to the new chain
        for residue_index in residue_indices:
            # WARNING: Chain index has to be read every iteration since it may change. Do not save it
            residue = self.residues[residue_index]
            residue.set_chain_index(chain.index)

    # Get the next available chain name
    # Find alphabetically the first letter which is not yet used as a chain name
    # If all letters in the alphabet are used already then return None
    def get_next_available_chain_name (self) -> Optional[str]:
        current_chain_names = [ chain.name for chain in self.chains ]
        next_available_chain_name = next((name for name in available_chains if name not in current_chain_names), None)
        if next_available_chain_name == None:
            raise InputError('There are more chains than available chain letters (' + str(len(available_chains)) + ')')
        return next_available_chain_name

    # Get a chain by its name
    def get_chain_by_name (self, name : str) -> 'Chain':
        return next((c for c in self.chains if c.name == name), None)

    # Get a summary of the structure
    def display_summary (self):
        print('Atoms: ' + str(self.atom_count))
        print('Residues: ' + str(len(self.residues)))
        print('Chains: ' + str(len(self.chains)))
        for chain in self.chains:
            print('Chain ' + chain.name + ' (' + str(len(chain.residue_indices)) + ' residues)')
            print(' -> ' + chain.get_sequence())

    # There may be chains which are equal in the structure (i.e. same chain name)
    # This means we have a duplicated/splitted chain
    # Repeated chains are usual and they are usually supported but with some problems
    # Also, repeated chains ususally come with repeated residues, which means more problems (see explanation below)

    # In the conext of this structure class we may have 2 different problems with a different solution each:
    # 1 - There is more than one chain with the same letter (repeated chain) -> rename the duplicated chains
    # 2 - There is a chain with atom indices which are not consecutive (splitted chain) -> create new chains

    # Rename repeated chains or create new chains if the fix_chains argument is True
    # WARNING: These fixes are possible only if there are less chains than the number of letters in the alphabet
    # Although there is no limitation in this code for chain names, setting long chain names is not compatible with pdb format
    
    # Check repeated chains (two chains with the same name) and return True if there were any repeats
    # Check splitted chains (a chains with non consecuitve residues) and try to fix them if requested
    def check_repeated_chains (self, fix_chains : bool = False, display_summary : bool = False) -> bool:
        # Order chains according to their names
        # Save also those chains which have a previous duplicate
        name_chains = {}
        repeated_chains = []
        for chain in self.chains:
            chain_name = chain.name
            current_name_chains = name_chains.get(chain_name, None)
            if not current_name_chains:
                name_chains[chain_name] = [chain]
            else:
                name_chains[chain_name].append(chain)
                repeated_chains.append(chain)
        # Display the summary of repeated chains if requested
        if display_summary:
            if len(repeated_chains) > 0:
                print('WARNING: There are repeated chains:')
                for chain_name, chains in name_chains.items():
                    chains_count = len(chains)
                    if chains_count > 1:
                        print('- Chain ' + chain_name + ' has ' + str(chains_count) + ' repeats' )
        # Rename repeated chains if requested
        if len(repeated_chains) > 0 and fix_chains:
            n_chains = len(self.chains)
            n_available_chains = len(available_chains)
            if n_chains > n_available_chains:
                # for chain in self.chains:
                #     print(str(chain) + ': ' + str(chain.atom_indices[0]) + ' to ' + str(chain.atom_indices[-1]))
                raise ValueError('There are more chains (' + str(n_chains) + ') than available chain letters (' + str(n_available_chains) + ')')
            current_letters = list(name_chains.keys())
            for repeated_chain in repeated_chains:
                last_chain_letter = repeated_chain.name
                while last_chain_letter in current_letters:
                    last_chain_letter = get_next_letter(last_chain_letter)
                repeated_chain.name = last_chain_letter
                current_letters.append(last_chain_letter)
        # Check if there are splitted chains
        for chain in self.chains:
            residue_indices = sorted(chain.residue_indices)
            # Check if residue indices are consecutive
            if residue_indices[-1] - residue_indices[0] + 1 == len(residue_indices):
                continue
            print('WARNING: splitted chain ' + chain.name)
            # If indices are not consecutive then we must find ranges of consecutive residues and create new chains for them
            previous_residue_index = residue_indices[0]
            consecutive_residues = [previous_residue_index]
            overall_consecutive_residues = []
            for residue_index in residue_indices[1:]:
                # If next residue is consecutive
                if residue_index == previous_residue_index + 1:
                    consecutive_residues.append(residue_index)
                # If next residue is NOT consecutive
                else:
                    overall_consecutive_residues.append(consecutive_residues)
                    consecutive_residues = [residue_index]
                # Update the previous index
                previous_residue_index = residue_index
            # Add the last split
            overall_consecutive_residues.append(consecutive_residues)
            # Now create new chains and reasign residues
            # Skip the first set of consecutive residues since they will stay in the original chain
            for residues_indices in overall_consecutive_residues[1:]:
                chain_name = self.get_next_available_chain_name()
                residues_selection = self.select_residue_indices(residues_indices)
                self.set_selection_chain_name(residues_selection, chain_name)

        # Fix repeated chains if requested
        return len(repeated_chains) > 0

    # There may be residues which are equal in the structure (i.e. same chain, name, number and icode)
    # In case 2 residues in the structure are equal we must check distance between their atoms
    # If atoms are far it means they are different residues with the same notation (duplicated residues)
    # If atoms are close it means they are indeed the same residue (splitted residue)

    # Splitted residues are found in some pdbs and they are supported by some tools
    # These tools consider all atoms with the same 'record' as the same residue
    # However, there are other tools which would consider the splitted residue as two different residues
    # This causes inconsistency along different tools besides a long list of problems
    # The only possible is fix is changing the order of atoms in the topology
    # Note that this is a breaking change for associated trajectories, which must change the order of coordinates
    # However here we provide tools to fix associates trajectories as well

    # Duplicated residues are usual and they are usually supported but with some problems
    # For example, pytraj analysis outputs use to sort results by residues and each residue is tagged
    # If there are duplicated residues with the same tag it may be not possible to know which result belongs to each residue
    # Another example are NGL selections once in the web client
    # If you select residue ':A and 1' and there are multiple residues 1 in chain A all of them will be displayed

    # Check residues to search for duplicated and splitted residues
    # Renumerate repeated residues if the fix_residues argument is True
    # Return True if there were any repeats
    def check_repeated_residues (self, fix_residues : bool = False, display_summary : bool = False) -> bool:
        # Track if residues have to be changed or not
        modified = False
        # Group all residues in the structure according to their chain, number and icode
        grouped_residues = {}
        # Check repeated residues which are one after the other
        # Note that these residues MUST have a different name
        # Otherwise they would have not been considered different residues
        # For these rare cases we use icodes to solve the problem
        non_icoded_residues = []
        last_residue = None
        for residue in self.residues:
            # Check residue to be equal than the previous residue
            if residue == last_residue:
                non_icoded_residues.append(residue)
                last_residue = residue
                continue
            last_residue = residue
            # Add residue to the groups of repeated residues
            current_residue_repeats = grouped_residues.get(residue, None)
            if not current_residue_repeats:
                grouped_residues[residue] = [ residue ]
            else:
                current_residue_repeats.append(residue)
        # In case we have non icoded residues
        if len(non_icoded_residues) > 0:
            if display_summary:
                print('There are non-icoded residues (' + str(len(non_icoded_residues)) + ')')
            # Set new icodes for non icoded residues
            if fix_residues:
                print('    Non icoded residues will recieve an icode')
                for residue in non_icoded_residues:
                    repeated_residues_group = grouped_residues[residue]
                    current_icodes = [ residue.icode for residue in repeated_residues_group if residue.icode ]
                    next_icode = next((cap for cap in available_caps if cap not in current_icodes), None)
                    if not next_icode:
                        raise ValueError('There are no more icodes available')
                    # print('Setting icode ' + next_icode + ' to residue ' + str(residue))
                    residue.icode = next_icode
                modified = True
        # Grouped residues with more than 1 result are considered as repeated
        repeated_residues = [ residues for residues in grouped_residues.values() if len(residues) > 1 ]
        if len(repeated_residues) == 0:
            return modified
        # In case we have repeated residues...
        if display_summary:
            print('WARNING: There are repeated residues (' + str(len(repeated_residues)) + ')')
            print('    e.g. ' + str(repeated_residues[0][0]))
        # Now for each repeated residue, find out which are splitted and which are duplicated
        covalent_bonds = self.bonds
        overall_splitted_residues = []
        overall_duplicated_residues = []
        for residues in repeated_residues:
            # Iterate over repeated residues and check if residues are covalently bonded
            # If any pair of residues are bonded add them both to the splitted residues list
            # At the end, all non-splitted residues will be considered duplicated residues
            # WARNING: Using a set here is not possible since repeated residues have the same hash
            # WARNING: Also comparing residues themselves is not advisable, so we use indices at this point
            splitted_residue_indices = set()
            for residue, other_residues in otherwise(residues):
                if residue.index in splitted_residue_indices:
                    continue
                # Get atom indices for all atoms connected to the current residue
                residue_bonds = sum([ covalent_bonds[index] for index in residue.atom_indices ], [])
                for other_residue in other_residues:
                    # Get all atom indices for each other residue and collate with the current residue bonds
                    if any( index in residue_bonds for index in other_residue.atom_indices ):
                        splitted_residue_indices.add(residue.index)
                        splitted_residue_indices.add(other_residue.index)
            # Finally obtain the splitted residues from its indices
            splitted_residues = [ self.residues[index] for index in splitted_residue_indices ]
            # Repeated residues which are not splitted are thus duplicated
            duplicated_residues = [ residue for residue in residues if residue.index not in splitted_residue_indices ]
            # Update the overall lists
            if len(splitted_residues) > 0:
                overall_splitted_residues.append(splitted_residues)
            if len(duplicated_residues) > 0:
                overall_duplicated_residues += duplicated_residues
        # In case we have splitted residues
        if len(overall_splitted_residues) > 0:
            if display_summary:
                print('    There are splitted residues (' + str(len(overall_splitted_residues)) + ')')
            # Fix splitted residues
            if fix_residues:
                print('        Splitted residues will be merged')
                # Set a function to sort atoms and residues by index
                def by_index (v):
                    return v._index
                # Merge splitted residues
                # WARNING: Merging residues without sorting atoms is possible, but this would be lost after exporting to pdb
                for splitted_residues in overall_splitted_residues:
                    # The first residue (i.e. the residue with the lower index) will be the one which remains
                    # It will absorb other residue atom indices
                    splitted_residues.sort(key=by_index) # Residues are not sorted by default, this is mandatory
                    first_residue = splitted_residues[0]
                    other_residues = splitted_residues[1:]
                    for residue in other_residues:
                        for atom in residue.atoms:
                            atom.residue = first_residue
                        # Then these other residues are eliminated
                        self.purge_residue(residue)
                print('        Atoms will be sorted to be together by residues')
                print('        NEVER FORGET: This will break any associated trajectory if coordinates are not sorted as well')
                # Sort atoms to group residue atoms together
                # Note that each atom index must be updated
                new_atom_indices = sum([ residue.atom_indices for residue in self.residues ], [])
                for i, new_atom_index in enumerate(new_atom_indices):
                    atom = self.atoms[new_atom_index]
                    atom._index = i
                self.atoms.sort(key=by_index)
                # Also residue 'atom_indices' must be updated
                for residue in self.residues:
                    residue._atom_indices = [ new_atom_indices.index(atom_index) for atom_index in residue._atom_indices ]
                # Bonds must be reset since atom indices have changes
                self._bonds = None
                # Prepare the trajectory atom sorter which must be returned
                # Include atom indices already so the user has to provide only the structure and trajectory filenames
                def trajectory_atom_sorter (
                    input_structure_file : 'File',
                    input_trajectory_file : 'File',
                    output_trajectory_file : 'File'
                ):
                    sort_trajectory_atoms(
                        input_structure_file,
                        input_trajectory_file,
                        output_trajectory_file,
                        new_atom_indices
                    )
                self.trajectory_atom_sorter = trajectory_atom_sorter
                self.new_atom_order = new_atom_indices
                modified = True
         # In case we have duplicated residues
        if len(overall_duplicated_residues) > 0:
            if display_summary:
                print('    There are duplicated residues (' + str(len(overall_duplicated_residues)) + ')')
            # Renumerate duplicated residues if requested
            if fix_residues:
                print('        Duplicated residues will be renumerated')
                # Get the next available number in the residue chain
                for duplicated_residue in overall_duplicated_residues:
                    maximum_chain_number = max([ residue.number for residue in duplicated_residue.chain.residues ])
                    duplicated_residue.number = maximum_chain_number + 1
                modified = True
        return modified

    # Atoms with identical chain, residue and name are considered repeated atoms
    # DANI: No recuerdo los problemas que daba tener átomos repetidos

    # Check atoms to search for repeated atoms
    # Rename repeated atoms if the fix_atoms argument is True
    # Return True if there were any repeats
    def check_repeated_atoms (self, fix_atoms : bool = False, display_summary : bool = False) -> bool:
        # Set two trackers for display
        repeated_atoms_count = 0
        example = None
        for residue in self.residues:
            # Iterate over the residue atoms counting how many repeated names we have
            current_names = []
            for atom in residue.atoms:
                # We check if the atom name already exists. If not, go to the next atom
                name = atom.name
                if name not in current_names:
                    current_names.append(name)
                    continue
                # When atom is repeated
                repeated_atoms_count += 1
                # If the fix was not requested we stop here
                if not fix_atoms:
                    continue
                # We set the name of the atom as the element letter + the count of this element
                # If element is missing take the first character of the atom name
                initial = atom.element
                if not initial or initial == ' ':
                    initial = name[0]
                number = 1
                new_name = initial + str(number)
                while new_name in current_names:
                    number += 1
                    new_name = initial + str(number)
                # Save an example for the logs if there is none yet
                if not example:
                    example = atom.name + ' renamed as ' + new_name + ' in residue ' + str(residue)
                atom.name = new_name
                current_names.append(new_name)
        # Display the summary of repeated atoms if requested
        if display_summary:
            if repeated_atoms_count > 0:
                print('WARNING: There are repeated atoms (' + str(repeated_atoms_count) + ')')
                print('    e.g. ' + example)
        return repeated_atoms_count > 0

    # Check bonds to be incoherent
    # i.e. check atoms not to have more or less bonds than expected according to their element
    def check_incoherent_bonds (self) -> bool:
        # Find out if there are hydrogens in the structure
        # It may happen that we encounter an structure without hydrogens
        has_hydrogen = next(( True for atom in self.atoms if atom.element == 'H' ), False)
        coherent_bonds = coherent_bonds_with_hydrogen if has_hydrogen else coherent_bonds_without_hydrogen
        # Cechk coherent bonds atom by atom
        for atom in self.atoms:
            # Do not check this atom bonds if this may be an ion
            # Most authors "force" dummy bonds in these atoms to make them stable
            if atom.element in SUPPORTED_ION_ELEMENTS:
                continue
            # Get actual number of bonds in the current atom both with and without ion bonds
            # LORE: This was a problem since some ions are force-bonded but bonds are actually not real
            # LORE: When an ion is forced we may end up with hyrogens with 2 bonds or carbons with 5 bonds
            # LORE: When an ions is really bonded we can not discard it or we may end up with orphan carbonds (e.g. weird ligands)
            min_nbonds = len(atom.get_bonds(skip_ions = True))
            max_nbonds = len(atom.get_bonds(skip_ions = False))
            # Get the accepted range of number of bonds for the current atom according to its element
            element = atom.element
            element_coherent_bonds = coherent_bonds.get(element, None)
            # If there are no specficiations for the current atom element then skip it
            if not element_coherent_bonds:
                continue
            # Get the range of accepted number of bonds
            min_allowed_bonds = element_coherent_bonds['min']
            max_allowed_bonds = element_coherent_bonds['max']
            # Check the actual number of bonds is insdie the accepted range
            if max_nbonds < min_allowed_bonds or min_nbonds > max_allowed_bonds:
                if min_nbonds == max_nbonds:
                    print(f' An atom with element {element} has {min_nbonds} bonds')
                else:
                    print(f' An atom with element {element} has between {min_nbonds} bonds (withou ions) and {max_nbonds} bonds (with ions)')
                if min_allowed_bonds == max_allowed_bonds:
                    plural_sufix = '' if min_allowed_bonds == 1 else 's'
                    print(f' It should have {min_allowed_bonds} bond{plural_sufix}')
                else:
                    print(f' It should have between {min_allowed_bonds} and {max_allowed_bonds} bonds')
                print(f' -> Atom {atom.label} (index {atom.index})')
                return True
        return False

    # Get all atomic covalent (strong) bonds
    # Bonds are defined as a list of atom indices for each atom in the structure
    # Rely on VMD logic to do so
    def get_covalent_bonds (self, selection : Optional['Selection'] = None) -> List[ List[int] ]:
        # It is important to fix elements before trying to fix bonds, since elements have an impact on bonds
        # VMD logic to find bonds relies in the atom element to set the covalent bond distance cutoff
        self.fix_atom_elements()
        # Generate a pdb strucutre to feed vmd
        auxiliar_pdb_filename = '.structure.pdb'
        self.generate_pdb_file(auxiliar_pdb_filename)
        # Get covalent bonds between both residue atoms
        covalent_bonds = get_covalent_bonds(auxiliar_pdb_filename, selection)
        # Remove the auxiliar pdb file
        os.remove(auxiliar_pdb_filename)
        return covalent_bonds

    # Make a copy of the current structure
    def copy (self) -> 'Structure':
        atom_copies = [ atom.copy() for atom in self.atoms ]
        residue_copies = [ residue.copy() for residue in self.residues ]
        chain_copies = [ chain.copy() for chain in self.chains ]
        structure_copy = Structure(atom_copies, residue_copies, chain_copies)
        structure_copy.bonds = [ [ atom_index for atom_index in atom_bonds ] for atom_bonds in self.bonds ]
        return structure_copy

    # Merge currnet structure with another structure
    # DANI: No lo he testeado en profundidad
    def merge (self, other : 'Structure') -> 'Structure':
        # Copy self atoms, residues and chains
        self_atom_copies = [ atom.copy() for atom in self.atoms ]
        self_residue_copies = [ residue.copy() for residue in self.residues ]
        self_chain_copies = [ chain.copy() for chain in self.chains ]
        # Copy other atoms, residues and chains
        other_atom_copies = [ atom.copy() for atom in other.atoms ]
        other_residue_copies = [ residue.copy() for residue in other.residues ]
        other_chain_copies = [ chain.copy() for chain in other.chains ]
        # Adapt indices in other atoms, residues and chains
        atom_index_offset = len(self_atom_copies)
        residue_index_offset = len(self_residue_copies)
        chain_index_offset = len(self_chain_copies)
        for atom in other_atom_copies:
            atom._index += atom_index_offset
            atom._residue_index += residue_index_offset
        for residue in other_residue_copies:
            residue._index += residue_index_offset
            residue._atom_indices = [ i + atom_index_offset for i in residue._atom_indices ]
            residue._chain_index += chain_index_offset
        for chain in other_chain_copies:
            chain._index += chain_index_offset
            chain._residue_indices = [ i + residue_index_offset for i in chain._residue_indices ]
        # Merge self with other atoms, residues and chains and build the new structure
        merged_atoms = self_atom_copies + other_atom_copies
        merged_residues = self_residue_copies + other_residue_copies
        merged_chains = self_chain_copies + other_chain_copies
        return Structure(merged_atoms, merged_residues, merged_chains)

    # Find rings in the structure and yield them as they are found
    # LIMITATION: Wrong rings may be found including smaller rings
    # LIMITATION: Use this function to see if a specific ring exists, but not to count rings
    def find_rings (self, selection : Optional[Selection] = None) -> Generator[ List[Atom], None, None ]:
        selected_atom_indices = selection.atom_indices if selection else None
        # We will use a logic which finds long paths of bonded heavy atoms
        # These paths have no 'branches' and thus, when there are 2 possible paths then both paths are explored independently
        # Save already searched atoms and the number of times it has been explored already
        # Note that atoms may be visited as many times as their number of heavy atom bonds -1
        already_searched_atoms = {}
        def find_rings_recursive (atom_path : List[Atom]) -> Generator[ List[Atom], None, None ]:
            # Get the current atom to continue the followup: the last atom
            current_atom = atom_path[-1] # Note that the list MUST always have at least 1 atom
            # Get bonded atoms to continue the path
            followup_atoms = current_atom.get_bonded_atoms()
            # Hydrogens are discarded since they have only 1 bond and thus they can only lead to a dead end
            followup_atoms = [ atom for atom in followup_atoms if atom.element != 'H' ]
            # In case there is a selection, check followup atoms not to be in the selection
            if selected_atom_indices:
                followup_atoms = [ atom for atom in followup_atoms if atom.index in selected_atom_indices ]
            # Check if this atom has been visited already too many times and, if so, stop here
            # Allowed search times is the posible number of rings where this atom can belong
            allowed_search_times = comb(len(followup_atoms), 2)
            already_search_times = already_searched_atoms.get(current_atom, 0)
            if already_search_times == allowed_search_times:
                return
            already_searched_atoms[current_atom] = already_search_times + 1
            # Remove the previous atom to avoid going back
            previous_atom = atom_path[-2] if len(atom_path) > 1 else None
            if previous_atom:
                followup_atoms.remove(previous_atom)
            # Iterate over the following bonded atoms
            for followup_atom in followup_atoms:
                # Check if the followup atom is in the path
                path_index = next(( index for index, atom in enumerate(atom_path) if atom == followup_atom ), None)
                # If so, we have found a ring
                if path_index != None:
                    ring = atom_path[path_index:]
                    yield ring
                    continue
                # Otherwise, keep searching
                new_path = atom_path + [followup_atom]
                for ring in find_rings_recursive(atom_path=new_path):
                    yield ring
        # Find a starting atom and run the recursive logic
        candidate_atom_indices = selected_atom_indices if selected_atom_indices else range(self.atom_count)
        for candidate_atom_index in candidate_atom_indices:
            candidate_atom = self.atoms[candidate_atom_index]
            # It must be heavy atom
            if candidate_atom.element == 'H':
                continue
            # It must not be a dead end already
            bonded_atoms = candidate_atom.get_bonded_atoms()
            bonded_atoms = [ atom for atom in bonded_atoms if atom.element != 'H' ]
            # In case there is a selection, check followup atoms not to be in the selection
            if selected_atom_indices:
                bonded_atoms = [ atom for atom in bonded_atoms if atom.index in selected_atom_indices ]
            # Check this atom is not a dead end already
            allowed_search_times = comb(len(bonded_atoms), 2)
            if allowed_search_times < 1:
                continue
            for ring in find_rings_recursive(atom_path=[candidate_atom]):
                yield ring

    # Given an atom selection, get all bonds between these atoms and any other atom in the structure
    # Note that inner bonds between atoms in the selection are discarded
    def get_selection_outer_bonds (self, selection : Selection) -> List[int]:
        # Get bonds from all atoms in the selection
        bonds = set()
        for atom_index in selection.atom_indices:
            atom_bonds = set(self.bonds[atom_index])
            bonds = bonds.union(atom_bonds)
        # Discard selection atoms from the bonds list to discard inner bonds
        bonds -= set(selection.atom_indices)
        return list(bonds)

    # Set the type of PTM according to the classification of the bonded residue
    ptm_options = {
        'protein': ValueError('A PTM residue must never be protein'),
        'nucleic': 'Nucleic acid linkage', # DANI: Esto es posible aunque muy raro y no se como se llama
        'carbohydrate': 'Glycosilation',
        'fatty': 'Lipidation',
        'steroid': 'Steroid linkage', # DANI: Esto es posible aunque muy raro y no se como se llama
        'ion': Warning('Ion is covalently bonded to protein'), # DANI: esto no es "correcto" pero si habitual
        'solvent': Warning('Solvent is covalently bonded to protein'), # DANI: probablemente sea un error
        'acetyl': 'Acetylation', # Typical N-capping terminals
        'amide': 'Amidation', # Typical C-capping terminals
        'other': Warning('Unknow type of PTM'), # Something not supported
    }

    # Find Post Translational Modifications (PTM) in the structure
    def find_ptms (self) -> Generator[dict, None, None]:
        # Find bonds between protein and non-protein atoms
        # First get all protein atoms
        protein_selection = self.get_protein_selection()
        if not protein_selection:
            return
        protein_atom_indices = set(protein_selection.atom_indices) # This is used later
        protein_outer_bonds = set(self.get_selection_outer_bonds(protein_selection))
        non_protein_selection = self.invert_selection(protein_selection)
        if not non_protein_selection:
            return
        non_protein_atom_indices = set(non_protein_selection.atom_indices)
        non_protein_atom_indices_bonded_to_protein = protein_outer_bonds.intersection(non_protein_atom_indices)
        # Get each residue bonded to the protein and based on its 'classification' set the name of the PTM
        for atom_index in non_protein_atom_indices_bonded_to_protein:
            atom = self.atoms[atom_index]
            residue = atom.residue
            residue_classification = residue.get_classification()
            ptm_classification = self.ptm_options[residue_classification]
            # If we found something impossible then raise the error
            if type(ptm_classification) == ValueError:
                print('Problematic residue: ' + str(residue))
                raise ptm_classification
            # If we do not know what it is then do tag it as a PTM but print a warning
            if type(ptm_classification) == Warning:
                print('WARNING: ' + str(ptm_classification) + ' -> ' + str(residue))
                continue
            # At this point we have found a valid PTM
            # Find the protein residue linked to this PTM
            atom_bonds = self.bonds[atom_index]
            protein_bond = next((bond for bond in atom_bonds if bond in protein_atom_indices), None)
            if protein_bond == None:
                raise ValueError('There must be at least one protein bond to the current atom')
            protein_residue_index = self.atoms[protein_bond].residue_index
            # Set the PTM
            yield { 'name': ptm_classification, 'residue_index': protein_residue_index }
            

### Related functions ###

# Calculate the distance between two atoms
def calculate_distance (atom_1 : Atom, atom_2 : Atom) -> float:
    squared_distances_sum = 0
    for i in { 'x': 0, 'y': 1, 'z': 2 }.values():
        squared_distances_sum += (atom_1.coords[i] - atom_2.coords[i])**2
    return math.sqrt(squared_distances_sum)

### Auxiliar functions ###

# Set a function to get the next letter from an input letter in alphabetic order
def get_next_letter (letter : str) -> str:
    if not letter:
        return 'A'
    if letter == 'z' or letter == 'Z':
        raise InputError("Limit of chain letters has been reached")
    next_letter = chr(ord(letter) + 1)
    return next_letter

# Given a string with 1 or 2 characters, return a new string with the first letter cap and the second letter not cap (if any)
def first_cap_only (name : str) -> str:
    if len(name) == 1:
        return name.upper()
    first_character = name[0].upper()
    second_character = name[1].lower()
    return first_character + second_character

# Convert numbers to lower text characters (chr 8020-8029)
lower_numbers = {
    '0': '₀',
    '1': '₁',
    '2': '₂',
    '3': '₃',
    '4': '₄',
    '5': '₅',
    '6': '₆',
    '7': '₇',
    '8': '₈',
    '9': '₉',
}
def get_lower_numbers (numbers_text : str) -> str:
    return ''.join([ lower_numbers[c] for c in numbers_text ])