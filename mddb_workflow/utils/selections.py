from typing import Optional

from mddb_workflow.utils.auxiliar import ranger


class Selection:
    """List of atom indices from a structure."""

    def __init__(self, atom_indices: Optional[list[int]] = None):
        self.atom_indices = sorted(atom_indices) if atom_indices != None else []

    def __repr__(self):
        return f'<Selection ({len(self.atom_indices)} atoms)>'

    def __hash__(self):
        return hash(tuple(self.atom_indices))

    def __len__(self):
        return len(self.atom_indices)

    def __bool__(self):
        """Return true if there is at least one atom index on the selection and false otherwise."""
        return len(self.atom_indices) > 0

    def __eq__(self, other):
        """Two selections are equal if they have the same atom indices."""
        if not isinstance(other, self.__class__):
            return False
        return set(self.atom_indices) == set(other.atom_indices)

    def __add__(self, other):
        """Return a new selection with atom indices from both self and the other selection."""
        return self.merge(other)

    def __sub__(self, other):
        """Return a new selection with self atom indices except for those atom indices in other."""
        return self.substract(other)

    def __and__(self, other):
        """Return a new selection with the intersection of both selections."""
        return self.intersection(other)

    def __or__(self, other):
        """Return a new selection with atom indices from both self and the other selection (same as add)."""
        return self.merge(other)

    def merge(self, other: Optional['Selection']) -> 'Selection':
        """Return a new selection made of self and other selection atom indices."""
        if not other:
            return self
        unique_atom_indices = list(set(self.atom_indices + other.atom_indices))
        return Selection(unique_atom_indices)

    def substract(self, other: Optional['Selection']) -> 'Selection':
        """Return a new selection with self atom indices except for those atom indices in other."""
        if not other: return self
        remaining_atom_indices = list(set(self.atom_indices) - set(other.atom_indices))
        return Selection(remaining_atom_indices)

    def intersection(self, other: Optional['Selection']) -> 'Selection':
        """Return a new selection with the intersection of both selections."""
        if not other:
            return Selection()
        self_atom_indices = set(self.atom_indices)
        other_atom_indices = set(other.atom_indices)
        intersection_atom_indices = list(self_atom_indices.intersection(other_atom_indices))
        return Selection(intersection_atom_indices)

    def to_mdanalysis(self) -> str:
        """Produce an MDAnalysis selection."""
        # Make sure it is not an empty selection
        if not self: raise ValueError('Trying to get MDAnalysis selection from an empty selection')
        return 'index ' + ' '.join([str(index) for index in self.atom_indices])

    def to_pytraj(self) -> str:
        """Produce a PyTraj selection."""
        # Make sure it is not an empty selection
        if not self: raise ValueError('Trying to get PyTraj selection from an empty selection')
        # NEVER FORGET: Pytraj counts atoms starting at 1, not at 0
        indices = [index + 1 for index in self.atom_indices]
        # Make ranges for atoms in a row
        return '@' + ranger(indices)

    def to_ngl(self) -> str:
        """Produce an NGL selection."""
        return '@' + ','.join([str(index) for index in self.atom_indices])

    def to_list(self) -> list[int]:
        """Return a copy of the atom indices as a list."""
        return self.atom_indices.copy()

    def to_bash(self, one_start: bool = False) -> str:
        """Get a string made of all indexes separated by underscores.
        This string can be then passed as a bash argument and easily parsed by other programms.
        Indices can start from 0 or from 1.
        """
        # Make sure it is not an empty selection
        if not self: raise ValueError('Trying to get Bash selection from an empty selection')
        if one_start:
            return '_'.join([str(index + 1) for index in self.atom_indices])
        else:
            return '_'.join([str(index) for index in self.atom_indices])

    def to_vmd(self) -> str:
        """Produce a vmd selection in tcl format."""
        # Make sure it is not an empty selection
        if not self: raise ValueError('Trying to get VMD selection from an empty selection')
        return 'index ' + ' '.join([str(index) for index in self.atom_indices])

    def to_ndx(self, selection_name: str = 'Selection') -> str:
        """Produce the content of gromacs ndx file."""
        # Make sure it is not an empty selection
        if not self: raise ValueError('Trying to get NDX selection from an empty selection')
        # Add a header
        # WARNING: Sometimes (and only sometimes) if the selection name includes white spaces there is an error
        # WARNING: Gromacs may get the first word only as the name of the selection, so we must purge white spaces
        # WARNING: However, if a call to this function passes a selection name it will expect it to be in the index file
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

    def to_ndx_file(self, selection_name: str = 'Selection', output_filepath: str = 'index.ndx'):
        """Create a gromacs ndx file."""
        index_content = self.to_ndx(selection_name)
        with open(output_filepath, 'w') as file:
            file.write(index_content)
