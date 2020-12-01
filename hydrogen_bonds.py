# Hydrogen bonds analysis
# 
# Perform the hydrogen bonds analysis along the different trajectory frames
# Track when an hydrogen bond is expected to be formed between each pair of candidate residues
# The analysis is carried by pytraj

import pytraj as pt
import numpy
import re

import json

# Perform an hydrogen bonds analysis for each interface
# The 'interfaces' input may be an empty list (i.e. there are no interfaces)
# In case there are no interfaces the analysis stops
def hydrogen_bonds (
    pt_trajectory,
    output_analysis_filename : str,
    topology_reference,
    interfaces : list ):

    # Return before doing anything if there are no interfaces
    if len(interfaces) == 0:
        return

    # Save the reference function to get an absolue atom index from a pytraj atom index
    get_atom_index = topology_reference.get_atom_index
    # Save the reference function to get a source residue from a pytraj residue
    pytraj2source = topology_reference.pytraj2source

    # Sabe each analysis to a dict which will be parsed to json
    output_analysis = []

    for interface in interfaces:
        
        # Select all interface residues in pytraj notation
        pt_interface = interface['pt_interface_1'] + interface['pt_interface_2']
        pt_selection = ':' + ','.join(map(str, pt_interface)) + ' @CA'

        # Run the analysis
        hbonds = pt.hbond(pt_trajectory, mask=pt_selection)

        # Get residues in each interface agent
        interface_1 = interface['interface_1']
        interface_2 = interface['interface_2']

        acceptor_atom_index_list = []
        donor_atom_index_list = []
        hydrogen_atom_index_list = []
        hbond_values = []
        
        # Search all predicted hydrogen bonds
        for d0 in hbonds:
            # d0.key example: "ASN15_OD1-LEU373_N-H"
            matchObj = re.match( r'\w{3}(\d*)_(.*)-\w{3}(\d*)_(.*)-(.*)', d0.key, re.M|re.I)
            # Mine all data from the parsed results
            if matchObj is not None:
                acceptor_resnum = matchObj.group(1)
                acceptor_atom = matchObj.group(2)
                donor_resnum = matchObj.group(3)
                donor_atom = matchObj.group(4)
                hydrogen_atom = matchObj.group(5)

                # Get the acceptor and donor residues in source notation
                acceptor = pytraj2source(int(acceptor_resnum))
                donor = pytraj2source(int(donor_resnum))
                
                # WARNING: The analysis may return hydrogen bonds between residues from the same agent
                # Accept the hydrogen bond only if its residues belong to different interface agents
                if ((acceptor in interface_1 and donor in interface_2)
                or (acceptor in interface_2 and donor in interface_1)):
                    
                    # Get the absolute index of each atom
                    acceptor_atom_index_list.append(get_atom_index(acceptor, acceptor_atom))
                    donor_atom_index_list.append(get_atom_index(donor, donor_atom))
                    hydrogen_atom_index_list.append(get_atom_index(donor, hydrogen_atom))
                    hbond_values.append(' '.join(map(str, d0.values)))

        # Write 
        output_analysis.append(
            {
                'name': interface['name'],
                'acceptors': acceptor_atom_index_list,
                'donors': donor_atom_index_list,
                'hydrogens': hydrogen_atom_index_list,
                'hbonds': hbond_values,
            }
        )

    # Export the analysis in json format
    with open(output_analysis_filename, 'w') as file:
        json.dump(output_analysis, file)
