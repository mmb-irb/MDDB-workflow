.. _task_documentation: generated with generate_task_docs.py

Workflow Tasks
==================

This page documents the available tasks in the Model Workflow system.
These tasks can be specified with the ``-i`` (include) or ``-e`` (exclude) flags.

Project Tasks
---------------

These tasks are executed once per project:

``chains``
~~~~~~~~~~

Get chain references.

``inputs``
~~~~~~~~~~

Set a function to load the inputs file

``itopology``
~~~~~~~~~~~~~

Get the input topology file.
If the file is not found try to download it.

``ligands``
~~~~~~~~~~~

Get the ligand residues mapping.

``mapping``
~~~~~~~~~~~

Map the structure aminoacids sequences against the Uniprot reference sequences.

``membrane``
~~~~~~~~~~~~

Get mapping of residues in the membrane.

``membranes``
~~~~~~~~~~~~~

Get mapping of residues in the membrane.

``pdbs``
~~~~~~~~

Get PDB references.

``pmeta``
~~~~~~~~~

Generate the project metadata file to be upload to the database.

``populations``
~~~~~~~~~~~~~~~

Get the MSM equilibrium populations filename.

``screenshot``
~~~~~~~~~~~~~~

Generate a screenshot of the system.

``stopology``
~~~~~~~~~~~~~

Generate the standardized topology file.

``topology``
~~~~~~~~~~~~

Get the processed topology file.

``transitions``
~~~~~~~~~~~~~~~

Get the MSM transition probabilities filename.

MD Tasks
-----------

These tasks are executed for each MD in the project:

Files
~~~~~~~~

``interactions``
^^^^^^^^^^^^^^^^

The processed interactions.
This is a bit exceptional since it is a value to be used and an analysis file to be generated.

``istructure``
^^^^^^^^^^^^^^

Get the input pdb filename from the inputs.
If the file is not found try to download it.

``itrajectory``
^^^^^^^^^^^^^^^

Get the input trajectory filename(s) from the inputs.
If file(s) are not found try to download it.

``structure``
^^^^^^^^^^^^^



``trajectory``
^^^^^^^^^^^^^^



Analyses
~~~~~~~~~~~~~~

``apl``
^^^^^^^

Area per lipid analysis.

``clusters``
^^^^^^^^^^^^

Run the cluster analysis.

``density``
^^^^^^^^^^^

Membrane density analysis.

``dihedrals``
^^^^^^^^^^^^^

Calculate torsions and then dihedral energies for every dihedral along the trajectory.

``dist``
^^^^^^^^

Calculate the distance mean and standard deviation of each pair of residues*.

``energies``
^^^^^^^^^^^^

Perform the electrostatic and vdw energies analysis for each pair of interaction agents.

``hbonds``
^^^^^^^^^^

Perform an hydrogen bonds analysis for each interaction interface.

``helical``
^^^^^^^^^^^



``linter``
^^^^^^^^^^

Lipid-protein interactions analysis.

``lorder``
^^^^^^^^^^

Calculate lipid order parameters for membranes.

``markov``
^^^^^^^^^^



``mdmeta``
^^^^^^^^^^

Generate the MD metadata file.

``pairwise``
^^^^^^^^^^^^

Perform an analysis for the overall structure and then one more analysis for each interaction.

``pca``
^^^^^^^

PCA, principal component analysis.

``perres``
^^^^^^^^^^

RMSD per residue analysis.

``pockets``
^^^^^^^^^^^

Perform the pockets analysis.

``rgyr``
^^^^^^^^

Radius of gyration analysis.

``rmsds``
^^^^^^^^^

RMSDs analysis.

``rmsf``
^^^^^^^^

RMSF, atom fluctuation analysis.

``sas``
^^^^^^^

Perform the Solvent Accessible Surface Analysis (SASA).

``thickness``
^^^^^^^^^^^^^

Membrane thickness analysis.

``tmscore``
^^^^^^^^^^^

TM scores analysis.

Task Groups
-------------

These are predefined groups of tasks that can be specified with a single flag.

``download``
~~~~~~~~~~~~

Includes the following tasks: ``itopology``, ``inputs``, ``populations``, ``transitions``, ``istructure``, ``itrajectory``

``setup``
~~~~~~~~~

Includes the following tasks: ``topology``, ``structure``, ``trajectory``

``network``
~~~~~~~~~~~

Includes the following tasks: ``mapping``, ``ligands``, ``chains``, ``pdbs``, ``membrane``

``minimal``
~~~~~~~~~~~

Includes the following tasks: ``pmeta``, ``mdmeta``, ``stopology``

``interdeps``
~~~~~~~~~~~~~

Includes the following tasks: ``interactions``, ``pairwise``, ``hbonds``, ``energies``, ``perres``, ``clusters``

