# This script is not part of the workflow, but a tool to set the inputs manually

import json

# INPUTS -----------------------------------------------------------------------

# Set the topology file name
original_topology_filename = 'spike_mutant_prot_glyc_amarolab.psf'

# Set the finally merged output trajectory name or the only input file name
original_trajectory_filename = 'spike_mutant_prot_glyc_amarolab_4.dcd'

# Set how the trajectory must be imaged and fited
# These protocolos may help in some situations, but the imaging step can not be fully automatized
# If protocols do not work the gromacs parameters must be modified manually
# 0 - No imaging, no fitting -> The trajectory is already imaged and fitted
# 1 - No imaging, only fitting -> The trajectory is already imaged but not fitted
# 2 - Single protein
# 3 - Protein in membrane
# 4 - Two interacting proteins
preprocess_protocol = 0

# Set if this is the ACE2, the RBD, both...
unit = 'Spike'

# Chain names
# EXPERIMENTAL INPUT
chainnames = { 'A':'Spike', 'B':'Spike', 'C':'Spike', 'G':'Glycans', 'H':'Glycans', 'I':'Glycans' }

# Set the membrane type as the name of membrane resiudes in the pdb file
# Set the membrane as 'No' if there is no membrane
membrane = 'No'

# Set the source pdb to represent in the overview metadata
# Additional data from the pdb is harvested by the loader while uploading to the database
pdbId = '6vsb'

# Write a breif description or title for this trajectory
name = ('Trajectories of full-length SPIKE protein in the Closed state (replica 1)')

# Additional comments
description = ('All-atom MD simulation of full-length SPIKE protein in the Open state '
               'bearing N165A and N234A mutations, protein + glycans only')

# Authors
authors = ('Lorenzo Casalino, Zied Gaieb, Abigail C. Dommer, Aoife M. Harbison,'
           ' Carl A. Fogarty, Emilia P. Barros, Bryn C. Taylor, Elisa Fadda,'
           ' Rommie E. Amaro')

# Program and version
program = 'NAMD'
version = '2.14'

# License and link to the license web page
license = ("This trajectory dataset is released under a Creative Commons "
           "Attribution 4.0 International Public License")
linkcense = "https://creativecommons.org/licenses/by/4.0/"

# Citation
# To set a citation use the following instructions:
# To add a line break type '(br)' inside the citation string
# To add superior text type '^' before each character
# Do never type '->' inside the string
citation = ('D. E. Shaw Research,(br)"Molecular Dynamics Simulations Related to SARS-CoV-2,"(br)'
            'D. E. Shaw Research Technical Data, 2020.(br)http://www.deshawresearch.com/resources_sarscov2.html/')
        
# Metadata
# These inputs may be automatically mined from the topology and trajectory files
# However they also may be forced here
# DANI: Todo mentira. El minado de metadata no funciona casi nunca para estos valores
# DANI: Dejé de mantenerlo hace tiempo y hay que poner todos los valores a mano

# Length is an important value since it is used in many graph axes in the web client
length = '411.7' # In nanoseconds (ns)
# The rest of values do not affect other outcomes
snapshots = 4117
frequency = '100' # In picoseconds (ps)
temp = 310 # In Kelvin (K)
ensemble = 'NPT'
timestep = '2' # In fs
pcoupling = 'Isotropic'
ff = 'CHARMM36 All-atom additive'
wat = 'TIP3P'
boxType = 'Triclinic'

# ------------------------------------------------------------------------------

# Join all inputs in a single object
inputs = {
    'original_topology_filename': original_topology_filename,
    'original_trajectory_filename': original_trajectory_filename,
    'unit' : unit,
    'chainnames' : chainnames,
    'membrane' : membrane,
    'pdbId' : pdbId,
    'name' : name,
    'description' : description,
    'authors' : authors,
    'program' : program,
    'version' : version,
    'license' : license,
    'linkcense' : linkcense,
    'citation' : citation,
    'length' : length,
    'snapshots' : snapshots,
    'frequency' : frequency,
    'temp' : temp,
    'ensemble' : ensemble,
    'timestep' : timestep,
    'pcoupling' : pcoupling,
    'ff' : ff,
    'wat': wat,
}

# Export it to json
inputs_filename = 'inputs.json'
with open(inputs_filename, 'w') as file:
    json.dump(inputs, file)