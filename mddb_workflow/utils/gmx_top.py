"""Read and filter Gromacs topologies (.top) together with their .itp files."""

import os
import re
import shutil
from dataclasses import dataclass, field

from mddb_workflow.utils.constants import GROMACS_EXECUTABLE

# Number of atom columns at the start of each line, in sections which refer to atoms
ATOM_COLUMNS = {
    'bonds': 2, 'pairs': 2, 'pairs_nb': 2, 'angles': 3, 'dihedrals': 4, 'cmap': 5,
    'constraints': 2, 'position_restraints': 1,
    'virtual_sites1': 2, 'virtual_sites2': 3, 'virtual_sites3': 4, 'virtual_sites4': 5,
    'distance_restraints': 2, 'dihedral_restraints': 4, 'orientation_restraints': 2,
    'angle_restraints': 4, 'angle_restraints_z': 2,
    'polarization': 2, 'water_polarization': 1, 'thole_polarization': 4,
}
# Sections where the first atom is a virtual site built from the next atoms
VIRTUAL_SITE_SECTIONS = {'virtual_sites1', 'virtual_sites2', 'virtual_sites3', 'virtual_sites4'}
# Sections which end a molecule definition
NON_MOLECULE_SECTIONS = {'moleculetype', 'system', 'molecules', 'intermolecular_interactions'}

SECTION_PATTERN = re.compile(r'^\s*\[\s*([A-Za-z0-9_]+)\s*\]')
INCLUDE_PATTERN = re.compile(r'^\s*#include\s+["<](.+)[">]')


@dataclass
class External:
    """An include line kept as it is, since the file is not next to the topology (e.g. a force field)."""
    line: str
    # Where the file was found in the Gromacs library, or None if it was not found
    path: str | None
    # Missing files included in the middle of a molecule, such as position restraints
    inside_molecule: bool = False


@dataclass
class MoleculeType:
    name: str | None
    # The whole definition, starting at its [ moleculetype ] line
    lines: list[str] = field(default_factory=list)
    atom_count: int = 0


def get_gromacs_library_dirs() -> list[str]:
    """Get the directories where Gromacs looks for included files, such as force fields."""
    dirs = [d for d in os.environ.get('GMXLIB', '').split(':') if d]
    gmx = shutil.which(GROMACS_EXECUTABLE) if GROMACS_EXECUTABLE else None
    if gmx:
        share = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(gmx))), 'share', 'gromacs', 'top')
        if os.path.isdir(share): dirs.append(share)
    return dirs


def split_comment(line: str) -> tuple[str, str]:
    """Split a line in data and comment."""
    position = line.find(';')
    if position == -1: return line, ''
    return line[:position], line[position:]


def read_flat(path: str, library_dirs: list[str]) -> list[str | External]:
    """Read a topology and put the content of the files next to it in place of their include lines.
    Files found only in the Gromacs library (force fields) or not found are kept as include lines."""
    items = []
    directory = os.path.dirname(os.path.abspath(path))
    with open(path) as file:
        for line in file:
            line = line.rstrip('\n')
            match = INCLUDE_PATTERN.match(line)
            if not match:
                items.append(line)
                continue
            included = match.group(1)
            local_path = os.path.join(directory, included)
            if os.path.isfile(local_path):
                items += read_flat(local_path, library_dirs)
                continue
            library_path = next((os.path.join(d, included) for d in library_dirs
                if os.path.isfile(os.path.join(d, included))), None)
            items.append(External(line, library_path))
    return items


def read_molecule_types(items: list[str | External], library_dirs: list[str]) -> dict[str, MoleculeType]:
    """Find the molecule definitions and count their atoms, also in included library files."""
    molecule_types = {}
    current = None
    section = None
    conditional_depth = 0
    for item in items:
        if isinstance(item, External):
            # A missing file included under #ifdef in the middle of a molecule does not end the molecule
            if current and not item.path and conditional_depth > 0:
                item.inside_molecule = True
                continue
            current = None
            section = None
            if item.path:
                external_items = read_flat(item.path, [os.path.dirname(item.path)] + library_dirs)
                molecule_types.update(read_molecule_types(external_items, library_dirs))
            continue
        match = SECTION_PATTERN.match(item)
        if match:
            section = match.group(1).lower()
            if section in NON_MOLECULE_SECTIONS: current = None
            if section == 'moleculetype':
                current = MoleculeType(name=None)
                current.lines.append(item)
            elif current: current.lines.append(item)
            continue
        if current: current.lines.append(item)
        data = split_comment(item)[0].strip()
        if data.startswith('#if'): conditional_depth += 1
        elif data.startswith('#endif'): conditional_depth -= 1
        if not data or data.startswith('#'): continue
        if section == 'moleculetype' and current and current.name is None:
            current.name = data.split()[0]
            molecule_types[current.name] = current
        elif section == 'atoms' and current:
            current.atom_count += 1
    return molecule_types


def read_molecules(items: list[str | External]) -> tuple[list[tuple[str, int]], list[int]]:
    """Read the [ molecules ] section: name and number of copies of each molecule.
    Also return the positions of its lines in the items list."""
    molecules = []
    positions = []
    section = None
    for i, item in enumerate(items):
        if isinstance(item, External): continue
        match = SECTION_PATTERN.match(item)
        if match:
            section = match.group(1).lower()
            if section == 'intermolecular_interactions':
                raise ValueError('Gromacs topologies with [ intermolecular_interactions ] are not supported')
            continue
        if section != 'molecules': continue
        positions.append(i)
        data = split_comment(item)[0].strip()
        if not data: continue
        if data.startswith('#'):
            raise ValueError(f'Preprocessor lines in the [ molecules ] section are not supported: {item}')
        name, count = data.split()[:2]
        molecules.append((name, int(count)))
    if not molecules:
        raise ValueError('The Gromacs topology has no molecules')
    return molecules, positions


def iterate_molecules(
    molecules: list[tuple[str, int]],
    molecule_types: dict[str, MoleculeType],
    atom_residues: list[int] | None,
):
    """Go through every copy of every molecule and yield its name, first atom and number of atoms.
    When the definition of a molecule is not found (e.g. water in a missing force field)
    each copy is one residue of the structure."""
    offset = 0
    for name, count in molecules:
        molecule = molecule_types.get(name)
        for _ in range(count):
            if molecule:
                size = molecule.atom_count
            else:
                if atom_residues is None:
                    raise ValueError(f'Definition of molecule {name} not found in the Gromacs topology and there is no structure to guess its atoms')
                if offset >= len(atom_residues):
                    raise ValueError(f'The Gromacs topology has more atoms than the structure ({len(atom_residues)})')
                residue = atom_residues[offset]
                size = 1
                while offset + size < len(atom_residues) and atom_residues[offset + size] == residue: size += 1
            yield name, offset, size
            offset += size


def get_gromacs_topology_atom_count(path: str, atom_residues: list[int] | None = None) -> int:
    """Count the atoms in a Gromacs topology.
    Give the residue index of every atom in the structure in case some molecule definitions are missing."""
    library_dirs = get_gromacs_library_dirs()
    items = read_flat(path, library_dirs)
    molecule_types = read_molecule_types(items, library_dirs)
    molecules, _ = read_molecules(items)
    return sum(size for _, _, size in iterate_molecules(molecules, molecule_types, atom_residues))


def filter_molecule_lines(molecule: MoleculeType, kept: list[int], new_name: str) -> list[str]:
    """Write a molecule definition with only some of its atoms (local indices, starting at 1)."""
    new_index = {old: new for new, old in enumerate(kept, 1)}
    lines = []
    section = None
    named = False
    for line in molecule.lines:
        match = SECTION_PATTERN.match(line)
        if match:
            section = match.group(1).lower()
            lines.append(line)
            continue
        data, comment = split_comment(line)
        tokens = data.split()
        if not tokens or tokens[0].startswith('#'):
            lines.append(line)
            continue
        if section == 'moleculetype':
            if named: raise ValueError(f'Unexpected line in molecule {molecule.name}: {line}')
            named = True
            lines.append(f'{new_name} {" ".join(tokens[1:])}')
            continue
        if section == 'atoms':
            # Comments are removed since they show the total charge, which is not right anymore
            atom = int(tokens[0])
            if atom in new_index: lines.append(f'{new_index[atom]:>6} ' + ' '.join(tokens[1:]))
            continue
        if section == 'exclusions':
            atoms = [int(t) for t in tokens]
            if atoms[0] not in new_index: continue
            remaining = [new_index[a] for a in atoms[1:] if a in new_index]
            if remaining: lines.append(' '.join(str(a) for a in [new_index[atoms[0]]] + remaining))
            continue
        # Settles only name the oxygen but they need the 3 atoms of the water
        if section == 'settles':
            raise ValueError(f'Molecule {molecule.name} has settles and only some of its atoms are kept')
        column_count = ATOM_COLUMNS.get(section)
        if column_count is None:
            raise ValueError(f'Section [ {section} ] in molecule {molecule.name} is not supported when removing only some of its atoms')
        atoms = [int(t) for t in tokens[:column_count]]
        kept_atoms = [a in new_index for a in atoms]
        if all(kept_atoms):
            new_tokens = [str(new_index[a]) for a in atoms] + tokens[column_count:]
            lines.append(' '.join(new_tokens) + (f' {comment}' if comment else ''))
            continue
        if section in VIRTUAL_SITE_SECTIONS and kept_atoms[0]:
            raise ValueError(f'Virtual site {atoms[0]} in molecule {molecule.name} is kept but some of the atoms it is built from are removed')
    return lines


def filter_gromacs_topology(
    input_path: str,
    output_path: str,
    kept_atoms: list[int],
    atom_residues: list[int] | None = None,
) -> int:
    """Write a new Gromacs topology with only the kept atoms (global indices, starting at 0).
    The .itp files next to the input topology are written inside the new topology.
    Molecules which keep only some of their atoms get a new definition (e.g. Protein_chain_A_filtered1).
    Give the residue index of every atom in the structure in case some molecule definitions are missing.
    Return the number of atoms in the new topology."""
    library_dirs = get_gromacs_library_dirs()
    items = read_flat(input_path, library_dirs)
    molecule_types = read_molecule_types(items, library_dirs)
    molecules, molecules_positions = read_molecules(items)

    # See which atoms are kept in every copy of every molecule
    kept_set = set(kept_atoms)
    new_molecules = []
    # Molecule name -> { kept atoms: new molecule name }
    variants = {}
    atom_count = 0
    kept_count = 0
    for name, offset, size in iterate_molecules(molecules, molecule_types, atom_residues):
        atom_count = offset + size
        local_kept = tuple(i + 1 for i in range(size) if offset + i in kept_set)
        kept_count += len(local_kept)
        if len(local_kept) == size:
            new_molecules.append(name)
        elif local_kept:
            if name not in molecule_types:
                raise ValueError(f'Only some atoms of molecule {name} are kept but its definition was not found')
            molecule_variants = variants.setdefault(name, {})
            if local_kept not in molecule_variants:
                molecule_variants[local_kept] = f'{name}_filtered{len(molecule_variants) + 1}'
            new_molecules.append(molecule_variants[local_kept])
    if atom_residues is not None and atom_count != len(atom_residues):
        raise ValueError(f'The Gromacs topology has {atom_count} atoms but the structure has {len(atom_residues)}')
    if kept_atoms and max(kept_atoms) >= atom_count:
        raise ValueError(f'Atom {max(kept_atoms)} is kept but the Gromacs topology has only {atom_count} atoms')

    # Join consecutive copies of the same molecule
    grouped = []
    for name in new_molecules:
        if grouped and grouped[-1][0] == name: grouped[-1][1] += 1
        else: grouped.append([name, 1])

    # Write the new topology
    # Molecules defined in local files get their filtered definitions right after the original one
    # Molecules defined in library files get them right before the [ system ] section
    output = []
    defined = set()
    written = set()
    def add_variants(names):
        for name in names:
            if name in written: continue
            for kept, new_name in variants.get(name, {}).items():
                output.append('')
                output.extend(filter_molecule_lines(molecule_types[name], list(kept), new_name))
            written.add(name)
        defined.clear()
    first_molecules_position = molecules_positions[0]
    molecules_positions = set(molecules_positions)
    section = None
    for i, item in enumerate(items):
        if isinstance(item, External):
            if not item.inside_molecule: add_variants(list(defined))
            output.append(item.line)
            continue
        match = SECTION_PATTERN.match(item)
        if match:
            section = match.group(1).lower()
            if section == 'system': add_variants(list(defined) + list(variants))
            elif section in NON_MOLECULE_SECTIONS: add_variants(list(defined))
        if i == first_molecules_position:
            output += [f'{name:<20} {count}' for name, count in grouped]
            continue
        if i in molecules_positions: continue
        output.append(item)
        if section == 'moleculetype' and not match:
            data = split_comment(item)[0].strip()
            if data and not data.startswith('#'): defined.add(data.split()[0])
    with open(output_path, 'w') as file:
        file.write('\n'.join(output) + '\n')
    return kept_count
