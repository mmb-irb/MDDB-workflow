import re

from model_workflow.utils.gmx_spells import get_tpr_content
from model_workflow.utils.type_hints import *

# Set some constants
AMBER_FLAG_PATTERN = r'%FLAG ([A-Z_]*)\s*'
AMBER_FORMAT_PATTERN = r'%FORMAT\(([0-9A-Za-z.]*)\)\s*'
AMBER_FLAG_FORMAT_TYPES = {
    '20a4': str,
    '1a80': str,
    '10I8': int,
    '20I4': int,
    '1I8': int,
    '3I8': int,
    '5E16.8': float,
}

# Unified class to read topologies of multiple formats
# This class relies in other readers and reads what others do not reach
class Topology:
    def __init__ (self, file : 'File'):
        self.file = file
        self.format = file.format
        # Set internal variables
        self._raw_content = None
        self._flag_summary = None

    # Read topology file content
    def read_raw_content (self):
        # If we already read the content then return it
        if self._raw_content: return self._raw_content
        # TPR topologies are not readable as they are, since they are in binary format
        if self.format == 'tpr':
            output, error = get_tpr_content(self.path)
            self._raw_content = output
            return self._raw_content
        # The rest of formats are usually readable ASCII
        with open(self.file.path, 'r') as file:
            self._raw_content = file.readlines()
        return self._raw_content
    raw_content = property(read_raw_content, None, None, "Topology file content (read only)")

    # Amber topologies only
    # Set a summary with all the available falgs and their lines
    def get_flag_summary (self) -> dict:
        # If we already made the summary then return it
        if self._flag_summary: return self._flag_summary
        # Otherwise build the summary
        summary = {}
        last_flag_lines = {}
        # Iterate lines in the topology raw content until we find a target flag
        line_number = 0
        for line in self.raw_content:
            match = re.search(AMBER_FLAG_PATTERN, line)
            if match:
                # Close the previous flag lines dict (first time is redundant)
                last_flag_lines['to'] = line_number-1
                # Get the new flag
                flag = match[1]
                # Get flag format
                flag_format_line = self.raw_content[line_number+1]
                format_match = re.search(AMBER_FORMAT_PATTERN, flag_format_line)
                if not format_match: raise ValueError(f'Flag {flag} has no format line')
                format = format_match[1]
                format_type = AMBER_FLAG_FORMAT_TYPES[format]
                # Set the new flag lines dict
                last_flag_lines = { 'type': format_type, 'from': line_number }
                summary[flag] = last_flag_lines
            line_number += 1
        # Finish the last summary lines
        last_flag_lines['to'] = line_number
        self._flag_summary = summary
        return self._flag_summary
    flag_summary = property(get_flag_summary, None, None, "Amber topologies only. Flags summary and lines (read only)")

    # Amber topologies only
    # Mine values under a specific flag
    # Available flags: TITLE, POINTERS, ATOM_NAME, CHARGE, ATOMIC_NUMBER, MASS, ATOM_TYPE_INDEX,
    #   NUMBER_EXCLUDED_ATOMS, NONBONDED_PARM_INDEX, RESIDUE_LABEL, RESIDUE_POINTER, BOND_FORCE_CONSTANT,                                                      
    #   BOND_EQUIL_VALUE, ANGLE_FORCE_CONSTANT, ANGLE_EQUIL_VALUE, DIHEDRAL_FORCE_CONSTANT,
    #   DIHEDRAL_PERIODICITY, DIHEDRAL_PHASE, SCEE_SCALE_FACTOR, SCNB_SCALE_FACTOR, SOLTY,
    #   LENNARD_JONES_ACOEF, LENNARD_JONES_BCOEF, BONDS_INC_HYDROGEN, BONDS_WITHOUT_HYDROGEN,
    #   ANGLES_INC_HYDROGEN, ANGLES_WITHOUT_HYDROGEN, DIHEDRALS_INC_HYDROGEN, DIHEDRALS_WITHOUT_HYDROGEN,
    #   EXCLUDED_ATOMS_LIST, HBOND_ACOEF, HBOND_BCOEF, HBCUT, AMBER_ATOM_TYPE, TREE_CHAIN_CLASSIFICATION,
    #   JOIN_ARRAY, IROTAT, RADIUS_SET, RADII, SCREEN, IPOL                       
    def mine_amber_flag (self, flag : str) -> list:
        # Find the lines to be read
        summary = self.flag_summary[flag]
        # Get the type of data in the flag content
        flag_type = summary['type']
        # Read the lines
        # Skip the first 2 lines since they are the flag header and the format header
        start_line = summary['from'] + 2
        end_line = summary['to'] + 1
        lines = self.raw_content[start_line:end_line]
        # Parse and add all values within the flag in a list
        values = []
        for line in lines:
            # Remove spaces and breaklines, then split the line by spaces
            string_values = line[0:-1].strip().split()
            # If flag contents are strings then return string values as they are
            if flag_type == str: values += string_values
            # Otherwise parse strings to numeric values if required
            # Note that flag_type may be int or float
            else: values += [ flag_type(value) for value in string_values ]
        return values

    # Get dihedrals data in a standarized format
    def get_dihedrals_data (self) -> List[dict]:
        if self.format == 'prmtop': return self._get_dihedrals_data_amber()
        raise ValueError(f'Reading dihedrals data is not supported in {self.format} topologies')

    # Amber topologies only
    # Get dihedrals data in a standarized format
    # https://ambermd.org/FileFormats.php
    def _get_dihedrals_data_amber (self) -> List[dict]:
        # Get dihedral-related values from the topology
        including_hydrogens = self.mine_amber_flag('DIHEDRALS_INC_HYDROGEN')
        excluding_hydrogens = self.mine_amber_flag('DIHEDRALS_WITHOUT_HYDROGEN')
        force_constants = self.mine_amber_flag('DIHEDRAL_FORCE_CONSTANT')
        preiodicities = self.mine_amber_flag('DIHEDRAL_PERIODICITY')
        phases = self.mine_amber_flag('DIHEDRAL_PHASE')
        ee_scale_factors = self.mine_amber_flag('SCEE_SCALE_FACTOR')
        vdw_scale_factors = self.mine_amber_flag('SCNB_SCALE_FACTOR')
        charges = self.mine_amber_flag('CHARGE')
        atom_types = self.mine_amber_flag('ATOM_TYPE_INDEX')
        ntypes = len(set(atom_types))
        nonbonded_parameter_indices = self.mine_amber_flag('NONBONDED_PARM_INDEX')
        lj_acoefs = self.mine_amber_flag('LENNARD_JONES_ACOEF')
        lj_bcoefs = self.mine_amber_flag('LENNARD_JONES_BCOEF')
        # Start processing dihedrals by finding their implicated atoms and "type"
        # These values are listed in groups of 5 integer numbers (4 atoms + 1 type)
        dihedrals_data = {}
        all_dihedrals = including_hydrogens + excluding_hydrogens
        dihedral_count = int(len(all_dihedrals) / 5)
        for dih in range(dihedral_count):
            start = dih * 5
            dihedral_numbers = all_dihedrals[start:start+5]
            # However indices are a bit "encripted" for performance reasons
            # This is documented somewhere, but you must divide by 3 to get the actual index
            # Some docs say yo also havo to add 1, but this is only to have the 1-based numeration
            # If you want the 0-based index then is just dividing by 3
            # Also the third index may be negative for some dihedrals but this is not literal
            # The negative number means the dihedral was ignored for end interaction calculations
            # This means its 1-4 atom electrostatic and van der waals energies are ignored
            # We make a tuple out of atom indices so we can use them as keys to group dihedrals
            # Note that some dihedrals have multiple terms
            encrypted_atom_indices = dihedral_numbers[0:4]
            atom_indices = tuple([ int(abs(index)/3) for index in encrypted_atom_indices ])
            ignored_1_4_energies = dihedral_numbers[2] < 0
            dihedral_type = dihedral_numbers[4]
            # Get the correspondong dihedral constants
            # Dihedral type numbering starts at 1, not at 0
            dihedral_index = dihedral_type - 1
            force_constant = force_constants[dihedral_index]
            periodicity = preiodicities[dihedral_index]
            if periodicity < 0:
                raise ValueError('Negative preiodicity is documented but not yet supported')
            phase = phases[dihedral_index]
            # Get the current dihedral data
            dihedral_data = dihedrals_data.get(atom_indices, None)
            if dihedral_data == None:
                # Create a new object if this is the first time to access this specific dihedral
                dihedral_data = { 'atom_indices': atom_indices, 'terms': [] }
                # Get end atom indices
                atom1_index = atom_indices[0]
                atom4_index = atom_indices[3]
                # Mine atom charges for the end atoms
                dihedral_data['atom1_charge'] = charges[atom1_index]
                dihedral_data['atom4_charge'] = charges[atom4_index]
                # To get the next two values we need to find a complex index
                # We substract 1 to atom type inidces to get 0-based indices
                atom1_type_index = atom_types[atom1_index] - 1
                atom4_type_index = atom_types[atom4_index] - 1
                # Thus the parameter index becomes also 0-based
                nonbonded_parameter_index_index = (ntypes * atom1_type_index) + atom4_type_index
                nonbonded_parameter_index = nonbonded_parameter_indices[nonbonded_parameter_index_index]
                # And this is the index we use to get the parameters we really want
                dihedral_data['lj_acoef'] = lj_acoefs[nonbonded_parameter_index]
                dihedral_data['lj_bcoef'] = lj_bcoefs[nonbonded_parameter_index]
                # Save the new dihedral dict in the overall dict
                dihedrals_data[atom_indices] = dihedral_data
            # Add current dihedral term
            dihedral_data['terms'].append({
                #'index': dihedral_index,
                'force': force_constant,
                'period': periodicity,
                'phase': phase,
            })
            # If this dihedral term is not the one used to calculate non covalent energies the continue
            if ignored_1_4_energies: continue
            # Mine non covalent energy factors
            # These values are noted for the whole dihedral and not for each term
            # This is because multiterm dihedrals do not have more than one term with end interactions
            ee_scale_factor = ee_scale_factors[dihedral_index]
            vdw_scale_factor = vdw_scale_factors[dihedral_index]
            dihedral_data['ee_scaling'] = ee_scale_factor
            dihedral_data['vdw_scaling'] = vdw_scale_factor
            # Get end atom indices
            atom1_index = atom_indices[0]
            atom4_index = atom_indices[3]

        return list(dihedrals_data.values())