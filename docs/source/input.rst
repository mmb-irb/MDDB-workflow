.. _input_file_documentation: generated with generate_input_docs.py

Input File
==========================

.. note::
   This documentation is automatically generated from the YAML template file.

This file defines a set of fields which are then read and used by the workflow.
Every field is explained and example values are provided.

Project metadata
----------------

The following metadata has no effect on the workflow itself, but it will be written to the output metadata file.
These values will be uploaded to the database and then exposed in the project overview.
They may be also useful to search this project using the browser.

``name``
~~~~~~~~

Write a brief description or title for this trajectory for the overview page.
This name may be used by the client to search the trajectory in the database.
The name is displayed in the overview page.

Default: ``My project``

``description``
~~~~~~~~~~~~~~~

Write additional comments about the project.
The description is displayed in the overview page.



::

	description: My project description

``authors``
~~~~~~~~~~~

Write author names.
Authors are displayed in the overview page.

.. warning::
 Check already existing values for this field to avoid making duplicates



::

	authors:
	  - My name
	  - My partner's name

``groups``
~~~~~~~~~~

Write author group name/s.
The group is displayed in the overview page.

.. warning::
 Check already existing values for this field to avoid making duplicates.



::

	groups:
	  - My group
	  - My partner's group

``contact``
~~~~~~~~~~~

How to contact the authors (e.g. a mail).
The contact is displayed in the overview page.

``program``
~~~~~~~~~~~

Program (software) name which carried the trajectory and its version.
Program and version are both displayed in the overview page.

.. warning::
 Check already existing values for this field to avoid making duplicates.



::

	program: GROMACS

``version``
~~~~~~~~~~~



::

	version: 2024.2

``type``
~~~~~~~~

Type of molecular dynamics.
At this moment there are only two options in this field: 'trajectory' and 'ensemble'.
Note that this field has an effect on the client:
Some time-dependent analysis will change the labels of their axes in order to make sense.
For instance RMSD X axis will be labeled as 'frames' instead of 'time'.



::

	type: trajectory

``method``
~~~~~~~~~~

MD method is displayed in the overview.
e.g. Classical MD, Targeted MD, Biased MD (Accelerated Weighted Ensemble), Enhanced sampling (Hamiltonian Replica Exchange) ...



::

	method: Classical MD

``license``
~~~~~~~~~~~

License and link to the license web page.
The license is displayed in the overview page.
Under the license there is a 'More information' button.
The link is used to redirect the user when this button is clicked.



::

	license: This trajectory dataset is released under a Creative Commons Attribution 4.0 International Public License

``linkcense``
~~~~~~~~~~~~~



::

	linkcense: https://creativecommons.org/licenses/by/4.0/

``citation``
~~~~~~~~~~~~

Citation for refering this simulation or related paper.
The citation is displayed in the overview page.
To set a citation use the following instructions:
to add a line break type '(br)' inside the citation string,
to add superior text type '^' before each character.

``thanks``
~~~~~~~~~~

Acknowledgements to be shown in the overview page.
e.g. I would like to thank my funders...

``accession``
~~~~~~~~~~~~~

OPTIONAL: Use this field to force a custom accession.
Accession is a short code to refer this project in your local database.
If no accession is forced then a default formatted accession will be generated when loading the project.

References
----------

References to other databases to enrich our data.

``links``
~~~~~~~~~

Links to somewhere related to the simulation.
These links are displayed in the overview page.
MolSSI uses this field to find simulations in our database and place the embed viewer in their website.
You must fit to the standard when adding a new MolSSI simulation.

.. note::
 This field has no effect in our workflow BUT others may rely on it.



::

	links:
	  - name: First data source
	  url: https://data.source.org/
	  - name: Second data source
	  url: https://mydata.com/

``pdb_ids``
~~~~~~~~~~~

Set the source pdb ids of the trajectory structure
Additional data from the pdb is harvested by the loader while uploading to the database
This data is displayed in the overview page



::

	pdb_ids:
	  - 6ACS
	  - 6M0J

``forced_references``
~~~~~~~~~~~~~~~~~~~~~

Set which reference sequences must be used in order to map residues in the structure of the simulation.
UniProt accession ids are accepted.
Forced references may be not provided or just cover the structure partially.
Then a blast will be run for each orphan chain sequence.
In addition, UniProt accession ids may be guessed from the PDB ids, when provided.

Forced references may be provided as a list.
In this scenario UniProt sequences are aligned to chain sequences to guess which UniProt belongs to each chain.
Forced references may be provided as a dictionary.
Then the user specifies which reference belongs to each chain.
Use the "noref" flag to mark a chain as "no referable" (e.g. antibodies, synthetic constructs).



::

	forced_references:
	  - Q9BYF1
	  - P0DTC2
	forced_references:
	  A: Q9BYF1
	  B: P0DTC2
	  C: noref

``ligands``
~~~~~~~~~~~

Set ligands in the simulation.
The workflow identifies ligands by their pubchem accession.
If a pubchem accession is passed then it is used.
Otherwise, the pubchem accession is found using other database accessions.
Each ligand must have at least one of the following attributes:

- pubchem: the PubChem accession

- drugbank: the DrugBank accession

- chembl: the ChEMBL accession

Optionally, a list of vmd selections may be provided to force the mapping

- vmd_selection: a list of vmd selections (chain D)

Ligands are mapped in the standard topology file
In addition, an RMSD analysis is run for every defined ligand



::

	ligands:
	  - pubchem: 1986
	  - drugbank: DB00945

Simulation metadata
-------------------

Simulation parameters.
DANI: Algún día esto será minado automáticamente

``framestep``
~~~~~~~~~~~~~

Time framestep in nanoseconds (ns).
May be None if this is not a trajectory, but an ensemble.
Framestep is an important value since it is used in many graph axes in the web client.



::

	framestep: 0.01 # ns

``timestep``
~~~~~~~~~~~~

The rest of values are displayed in the web client as trajectory metadata.
These values do not affect other outcomes in the workflow.
Simulation timestep in femtoseconds (fs)



::

	timestep: 2 # fs

``temp``
~~~~~~~~

Temperature in Kelvin (K).



::

	temp: 310 # K

``ensemble``
~~~~~~~~~~~~

Ensemble
e.g. NVT, NPT, etc.

.. warning::
 Check already existing values for this field to avoid making duplicates



::

	ensemble: NPT

``ff``
~~~~~~

Force fields

.. warning::
 Check already existing values for this field to avoid making duplicates



::

	ff:
	  - Amber ff14SB
	  - GLYCAM-06j

``wat``
~~~~~~~

Water force fields.

.. warning::
 Check already existing values for this field to avoid making duplicates



::

	wat: TIP3P

``boxtype``
~~~~~~~~~~~

Boxtype
e.g. Triclinic, Cubic, Dodecahedron.

.. warning::
 Check already existing values for this field to avoid making duplicates



Analysis parameters
-------------------

These fields have an impact in the analysis workflow.

``interactions``
~~~~~~~~~~~~~~~~

Set which are the interesting interactions to be analyzed
A bunch of interaction-specific analyses will be run for each interaction and displayed in the web client

Interactions are defined by the 'agents' which are meant to interact pairwise
An 'agent' may be anything, even a group of unrelated molecules
Atoms of different agents which are close enought will be considered as interface atoms
These atoms will be the ones considered in interface analyses
If no interface atoms are found then the interaction is considered not valid and the user is warned

Interactions are uploaded to the database as part of the project metadata and as an independent analysis
Project metadata includes the interaction name, agents name and agent atom selections (vmd syntax)
Analysis data includes also every agent atom indices (both the whole agent and the interface only)

Each interaction has the following attributes:

- name: a string tag used to relate interaction analyses data with their corresponding atoms. In addition, the name is used to label the corresponding analyses in the web client.

- agent_1: the name of the first agent in the interaction, which is used to label in the client.

- selection_1: the VMD selection of the first agent in the interaction.

- agent_2: the name of the second agent in the interaction, which is used to label in the client.

- selection_2: the VMD selection of the second agent in the interaction.

- distance_cutoff (optional): the distance used to determine which atoms are in the interface (in Å).

The default value is intended for atomistic simulations.
Thus coarse grain interactions may need manual input distance cutoff

VMD atom selection language:
https://www.ks.uiuc.edu/Research/vmd/vmd-1.3/ug/node132.html



::

	interactions:
	  - name: protein-ligand interaction
	  agent_1: protein
	  selection_1: not resname lig
	  agent_2: ligand
	  selection_2: resname lig
	  - name: domain-domain interaction
	  agent_1: domain 1
	  selection_1: resid 2 to 291
	  agent_2: domain 2
	  selection_2: resid 2 to 291
	  - name: dna-dna hybridization
	  agent_1: strain A
	  selection_1: chain A
	  agent_2: strain B
	  selection_2: chain B
	  distance_cutoff: 10

``pbc_selection``
~~~~~~~~~~~~~~~~~

Set those residues which are under periodic boundary conditions (PBC)
These residues are excluded from the imaging centering and fitting
These residues are excluded in the follwoing analyses:

- RMSD: Sudden jumps in PBC residues result in non-sense high peaks

- RMSD per residue: Sudden jumps in PBC residues result in non-sense high peaks

- RMSD pairwise: Sudden jumps in PBC residues result in non-sense high peaks

- TM score: Sudden jumps in PBC residues result in non-sense high peaks

- RGYR: Sudden jumps in PBC residues result in non-sense high changes

- RMSF: Sudden jumps in PBC residues result in non-sense high peaks

- PCA: Sudden jumps make not sense in PCA and they eclipse non-PBC movements

- SASA: Residues close to the boundary will be considered exposed to solvent while they may be not

- Pockets: Residues close to the boundary may be considered to have pockets while they have not <span style="color:red">(
.. tip::
 Esto en realidad no se puede hacer porque fpocket no permite "descartar" átomos de manera inteligente. Si quitas átomos para que no encuentre pockets en ellos entonces pueden aparecer pockets en los sitios que están ocupados por estos átomos. De momento descartamos el análisis entero cuando hay algo en PBC y listo)<span />

- Clusters: Since Clustering is RMSD-based it has the same limitations

If this field is set to 'auto' then PBC residues are set automatically
Solvent, counter ions and membrane lipids are selected in this cases
Note that these are the most tipical residues under periodic boundary conditions

This field is also useful for those scenarions with several protein or nucleic acid molecules floating around
In this situation you can not image and fit all molecules
You must focus in one molecule and let the others stay in periodic boundary conditions

These residues are defined using VMD selection language: https://www.ks.uiuc.edu/Research/vmd/vmd-1.3/ug/node132.html



::

	pbc_selection: water or ions

Default: ``auto``

``cg_selection``
~~~~~~~~~~~~~~~~

Set those atoms which are not actual "atoms" but coarse grained (CG) beads.

.. warning::
 EXPERIMENTAL INPUT



Representation parameters
-------------------------

These fields have an impact in the display of the simulation once in the web client

``chainnames``
~~~~~~~~~~~~~~

Set optional custom chain names which may be longer than a single letter.
This names are used to label chains in the web client.



::

	chainnames:
	  A: Protein
	  B: Ligand

``customs``
~~~~~~~~~~~

The web client sets some default representations (molecular viewer configurations)
They highlight important features in the structure according to the topology reference or interactions
In addition, you may set extra customized representations which are interesting for you
These representations will be available in the web client
WARNING: Make sure whatever you want to represent is not already represented by default or it will be duplicated



::

	customs:
	  - name: A custom view focusing on an interesting residue
	  representations:
	  - name: The very interesing residue
	  selection: VIR
	  type: ball+stick
	  color: element
	  - name: The resting boring molecule
	  selection: not VIR
	  type: cartoon
	  color: chainid

``orientation``
~~~~~~~~~~~~~~~

Set a specific starting orientation for the web client viewer
Normally this is done once the simulation has been uploaded since there is no easy way to get the orientation before.



::

	  [
	  72.05997406618104,
	  21.871748915422142,
	  47.89720038949639,
	  0,
	  34.3234627961572,
	  42.053333152877315,
	  -70.84188126104011,
	  0,
	  -39.93012781662099,
	  75.61943426331311,
	  25.542927052994127,
	  0,
	  -63.015499114990234,
	  -33.07249975204468,
	  -39.439000606536865,
	  1
	  ]

Others
------

Other metadata

``multimeric``
~~~~~~~~~~~~~~

Set if we have any multimeric form
e.g. monomer, dimer, trimer
This field was requested by the referees
Its only use for now is as a parameter in project queries

.. tip::
 Esto es provisional, lo suyo sería automatizarlo



::

	multimeric:
	  - monomer
	  - trimer

Collections
-----------

Set also additional collection related metadata values.

``collections``
~~~~~~~~~~~~~~~

Set to which collection does this simulation belong to.
Currently supported collections:

- cv19

- mcns

- abc

- bigna

- model

``cv19_unit``
~~~~~~~~~~~~~

BioExcel-CV19 specific metadata fields.
Set which family does this trajectory belong to.
Supported units:

- RBD-ACE2

- RBD

- ACE2

- Spike

- 3CLpro

- PLpro

- Polymerase

- E protein

- Exoribonuclease

- Other

``cv19_startconf``
~~~~~~~~~~~~~~~~~~

Set some additional inputs requested by the referees

- Starting conformation of the spike (options: close, 1 open, 2 open, 3 open)

- Are there antibodies? (e.g. true)

- Are there nanobodies? (e.g. false)

``cv19_abs``
~~~~~~~~~~~~

``cv19_nanobs``
~~~~~~~~~~~~~~~

Input files
-----------

Directories for every MD in the project.
Input file paths for topology, trajectory and structure.
Note that all these values may be specified thorugh command line as well.

``mds``
~~~~~~~

Each project may contain several Molecular Dynamics (MD).
Each MD is to be stored in an independent folder when running the workflow.
Each MD must have a different name, which will be used also to assign directories in the workflow.

MDs may include additional metadata to overwrite the project metadata for a specific case.



::

	mds:
	  - name: replica 310 K
	  mdir: replica_1
	  temp: 310
	  - name: replica 311 K
	  mdir: replica_2
	  temp: 311

``mdref``
~~~~~~~~~

Also the reference MD is to be defined by providing the index of the MDs list.
If there is not an MD which is more important than others then simply set the first MD (0) as the reference



::

	mdref: 0

``input_topology_filepath``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Input topology, trajectory and structure files.

LORE
These files have been passed to the workflow through command line traditionally.
Now the inputs file also provides this option.

``input_structure_filepath``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``input_trajectory_filepaths``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

