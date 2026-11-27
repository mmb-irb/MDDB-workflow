"""Tests for counting and filtering atoms in Gromacs topologies (.top and .itp files)."""
import os
import pytest

from mddb_workflow.utils.gmx_top import (
    get_gromacs_topology_atom_count, filter_gromacs_topology,
    read_flat, read_molecule_types, read_molecules, SECTION_PATTERN, split_comment,
)

# A 4 atoms peptide with position restraints in its own file
PEPTIDE_ITP = """
[ moleculetype ]
; name nrexcl
PEP 3

[ atoms ]
1 N  1 ALA N  1 -0.3 14.0 ; qtot -0.3
2 CT 1 ALA CA 2  0.1 12.0 ; qtot -0.2
3 C  2 GLY C  3  0.5 12.0 ; qtot 0.3
4 O  2 GLY O  4 -0.3 16.0 ; qtot 0.0

[ bonds ]
1 2 1
2 3 1
3 4 1

[ angles ]
1 2 3 1
2 3 4 1

[ exclusions ]
1 3 4

#ifdef POSRES
#include "posre_pep.itp"
#endif
"""
POSRE_ITP = """
[ position_restraints ]
1 1 1000 1000 1000
3 1 1000 1000 1000
"""
# Water is defined here, ions come from a force field which is not available
TOPOLOGY = """
#include "missing.ff/forcefield.itp"
#include "pep.itp"

[ moleculetype ]
SOL 2

[ atoms ]
1 OW 1 SOL OW 1 -0.8 16.0
2 HW 1 SOL HW1 1 0.4 1.0
3 HW 1 SOL HW2 1 0.4 1.0

[ settles ]
1 1 0.1 0.16

#include "missing.ff/ions.itp"

[ system ]
Test

[ molecules ]
PEP 2
SOL 2
NA 1
"""
# Residue of every atom: 2 peptides with 2 residues each, 2 waters and 1 ion
ATOM_RESIDUES = [0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6]


@pytest.fixture
def topology(tmp_path) -> str:
    (tmp_path / 'pep.itp').write_text(PEPTIDE_ITP)
    (tmp_path / 'posre_pep.itp').write_text(POSRE_ITP)
    path = tmp_path / 'topol.top'
    path.write_text(TOPOLOGY)
    return str(path)


def read_output(path: str) -> tuple[list, dict]:
    items = read_flat(path, [])
    return read_molecules(items)[0], read_molecule_types(items, [])


def section_lines(molecule, section_name: str) -> list[list[str]]:
    """Data of the lines in a section of a molecule definition."""
    section = None
    lines = []
    for line in molecule.lines:
        match = SECTION_PATTERN.match(line)
        if match:
            section = match.group(1)
            continue
        data = split_comment(line)[0].split()
        if section == section_name and data and not data[0].startswith('#'): lines.append(data)
    return lines


@pytest.mark.unit_int
def test_atom_count(topology: str):
    """The ion has no definition, so it takes the size of its residue in the structure."""
    assert get_gromacs_topology_atom_count(topology, ATOM_RESIDUES) == 15
    with pytest.raises(ValueError, match='Definition of molecule NA not found'):
        get_gromacs_topology_atom_count(topology)


@pytest.mark.unit_int
def test_remove_whole_molecules(topology: str, tmp_path):
    """Removing the second water and the ion only changes the molecules section."""
    output = str(tmp_path / 'filtered.top')
    assert filter_gromacs_topology(topology, output, list(range(11)), ATOM_RESIDUES) == 11
    molecules, molecule_types = read_output(output)
    assert molecules == [('PEP', 2), ('SOL', 1)]
    assert sorted(molecule_types) == ['PEP', 'SOL']
    # The local files are inside the new topology and the force field ones are kept as includes
    text = open(output).read()
    assert '#include "pep.itp"' not in text
    assert '#include "missing.ff/ions.itp"' in text
    assert section_lines(molecule_types['PEP'], 'position_restraints') == [['1', '1', '1000', '1000', '1000'], ['3', '1', '1000', '1000', '1000']]


@pytest.mark.unit_int
def test_remove_part_of_a_molecule(topology: str, tmp_path):
    """Removing the CA of the second peptide gives it a new definition with renumbered atoms."""
    output = str(tmp_path / 'filtered.top')
    kept = [i for i in range(15) if i != 5]
    assert filter_gromacs_topology(topology, output, kept, ATOM_RESIDUES) == 14
    molecules, molecule_types = read_output(output)
    assert molecules == [('PEP', 1), ('PEP_filtered1', 1), ('SOL', 2), ('NA', 1)]
    filtered = molecule_types['PEP_filtered1']
    assert [line[4] for line in section_lines(filtered, 'atoms')] == ['N', 'C', 'O']
    assert [line[0] for line in section_lines(filtered, 'atoms')] == ['1', '2', '3']
    # Only the bond between C and O remains, and no angle
    assert section_lines(filtered, 'bonds') == [['2', '3', '1']]
    assert section_lines(filtered, 'angles') == []
    assert section_lines(filtered, 'exclusions') == [['1', '2', '3']]
    assert section_lines(filtered, 'position_restraints') == [['1', '1', '1000', '1000', '1000'], ['2', '1', '1000', '1000', '1000']]
    # The original definition is still there for the first peptide
    assert molecule_types['PEP'].atom_count == 4


@pytest.mark.unit_int
def test_same_cut_shares_definition(topology: str, tmp_path):
    """Both peptides lose the same atom, so they share the new definition."""
    output = str(tmp_path / 'filtered.top')
    kept = [i for i in range(15) if i not in (3, 7)]
    filter_gromacs_topology(topology, output, kept, ATOM_RESIDUES)
    molecules, _ = read_output(output)
    assert molecules == [('PEP_filtered1', 2), ('SOL', 2), ('NA', 1)]


@pytest.mark.unit_int
def test_unsupported_cuts(topology: str, tmp_path):
    """Waters with settles can not lose some of their atoms."""
    output = str(tmp_path / 'filtered.top')
    with pytest.raises(ValueError, match='settles'):
        filter_gromacs_topology(topology, output, [i for i in range(15) if i != 9], ATOM_RESIDUES)
    with pytest.raises(ValueError, match='more atoms than the structure'):
        filter_gromacs_topology(topology, output, list(range(10)), ATOM_RESIDUES[:-1])
