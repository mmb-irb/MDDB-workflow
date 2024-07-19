from typing import List, Optional

from model_workflow.utils.auxiliar import ranger

# A selection is a list of atom indices from a structure
class Selection:

    def __init__ (self, atom_indices : List[int]):
        self.atom_indices = atom_indices

    def __repr__ (self):
        return f'<Selection ({len(self.atom_indices)} atoms)>'

    def __len__ (self):
        return len(self.atom_indices)

    # Return true if there is at least one atom index on the selection and false otherwise
    def __bool__ (self):
        return len(self.atom_indices) > 0

    # Two selections are equal if they have the same atom indices
    def __eq__ (self, other):
        return set(self.atom_indices) == set(other.atom_indices)

    # Return a new selection with atom indices from both self and the other selection
    def __add__ (self, other):
        return self.merge(other)

    # Return a new selection with self atom indices except for those atom indices in other
    def __sub__ (self, other):
        return self.substract(other)

    # Return a new selection with the intersection of both selections
    def __and__ (self, other):
        return self.intersection(other)

    # Return a new selection with atom indices from both self and the other selection (same as add)
    def __or__ (self, other):
        return self.merge(other)

    @classmethod
    def from_prody (cls, prody_selection):
        indexes = [ atom.getIndex() for atom in prody_selection.iterAtoms() ]
        return cls(indexes)

    # Return a new selection made of self and other selection atom indices
    def merge (self, other : Optional['Selection']) -> 'Selection':
        if not other:
            return self
        unique_atom_indices = list(set( self.atom_indices + other.atom_indices ))
        return Selection(unique_atom_indices)

    # Return a new selection made of self and not other selection atom indices 
    def substract (self, other : Optional['Selection']) -> 'Selection':
        if not other:
            return self
        remaining_atom_indices = [ atom_index for atom_index in self.atom_indices if atom_index not in other.atom_indices ]
        return Selection(remaining_atom_indices)

    # Return a new selection with the intersection of both selections
    def intersection (self, other : Optional['Selection']) -> 'Selection':
        if not other:
            return self
        self_atom_indices = set(self.atom_indices)
        other_atom_indices = set(other.atom_indices)
        intersection_atom_indices = list(self_atom_indices.intersection(other_atom_indices))
        return Selection(intersection_atom_indices)

    def to_prody (self) -> str:
        return 'index ' + ' '.join([ str(index) for index in self.atom_indices ])

    def to_mdanalysis (self) -> str:
        return self.to_prody() # Same notation than prody for atom indices

    def to_pytraj (self) -> str:
        # NEVER FORGET: Pytraj counts atoms starting at 1, not at 0
        indices = [ index + 1 for index in self.atom_indices ]
        # Make ranges for atoms in a row
        return '@' + ranger(indices)

    # Get a string made of all indexes separated by underscores
    # This string can be then passed as a bash argument and easily parsed by other programms
    # Indices can start from 0 or from 1
    def to_bash (self, one_start : bool = False) -> str:
        if one_start:
            return '_'.join([ str(index + 1) for index in self.atom_indices ])
        else:
            return '_'.join([ str(index) for index in self.atom_indices ])

    # Generate a vmd selection in tcl format
    def to_vmd (self) -> str:
        return 'index ' + ' '.join([ str(index) for index in self.atom_indices ])

    # Generate the file content of gromacs ndx file
    def to_ndx (self, selection_name : str = 'Selection') -> str:
        # Add a header
        # WARNING: Sometimes (and only sometimes) if the selection name includes white spaces there is an error
        # WARNING: Gromacs may get the wirst word only as the name of the selection, so we must pruge white spaces
        # WARNING: However, if a call to this function passes a s election names it will expect it to be in the index file
        # WARNING: For this reason we must kill here the proccess and warn the user
        if ' ' in selection_name:
            raise ValueError(f'A Gromacs index file selection name must never include white spaces: {selection_name}')
        content = '[ ' + selection_name + ' ]\n'
        count = 0
        for index in self.atom_indices:
            # Add a breakline each 15 indices
            count += 1
            if count == 15:
                content += '\n'
                count = 0
            # Add a space between indices
            # Atom indices go from 0 to n-1
            # Add +1 to the index since gromacs counts from 1 to n
            content += str(index + 1) + ' '
        content += '\n' 
        return content