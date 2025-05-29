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



``inputs``
~~~~~~~~~~

Set a function to load the inputs file

``itopology``
~~~~~~~~~~~~~

Get the input topology file.
If the file is not found try to download it.

``ligands``
~~~~~~~~~~~



``mapping``
~~~~~~~~~~~



``membrane``
~~~~~~~~~~~~



``membranes``
~~~~~~~~~~~~~



``pdbs``
~~~~~~~~



``pmeta``
~~~~~~~~~

Project metadata filename

``populations``
~~~~~~~~~~~~~~~



``screenshot``
~~~~~~~~~~~~~~



``stopology``
~~~~~~~~~~~~~



``topology``
~~~~~~~~~~~~



``transitions``
~~~~~~~~~~~~~~~



MD Tasks
-----------

These tasks are executed for each MD in the project:

Files
~~~~~~~~

``interactions``
^^^^^^^^^^^^^^^^



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



``clusters``
^^^^^^^^^^^^



``density``
^^^^^^^^^^^



``dihedrals``
^^^^^^^^^^^^^



``dist``
^^^^^^^^



``energies``
^^^^^^^^^^^^



``hbonds``
^^^^^^^^^^



``helical``
^^^^^^^^^^^



``linter``
^^^^^^^^^^



``lorder``
^^^^^^^^^^



``markov``
^^^^^^^^^^



``mdmeta``
^^^^^^^^^^



``pairwise``
^^^^^^^^^^^^



``pca``
^^^^^^^



``perres``
^^^^^^^^^^



``pockets``
^^^^^^^^^^^



``rgyr``
^^^^^^^^



``rmsds``
^^^^^^^^^



``rmsf``
^^^^^^^^



``sas``
^^^^^^^



``thickness``
^^^^^^^^^^^^^



``tmscore``
^^^^^^^^^^^



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

