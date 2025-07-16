Usage
=====

How to use
----------

.. contents::
   :local:


Folder structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All you need to start processing your files is a topology file (prmtop, tpr, psf, etc.) and any number of trajectory files (nc, xtc, dcd, etc.). Every independent trajectory must be in a different directory, where several output files will be generated for every trajectory. Here is an example of a directory we are about to analyze:

.. code-block:: text

    project/
    ├── raw_topology.prmtop
    ├── raw_trajectory_01.nc
    ├── raw_trajectory_02_part01.nc
    ├── raw_trajectory_02_part02.nc

Note that here we have 2 independent replicas sharing common topology. The second replica is actually splitted in 2 consecutive parts but this is not a problem.

.. warning::
   In order to keep your trajectories safe we recommend making a copy or using symlinks as workflow inputs. e.g. ``ln -s ~/my_data/my_trajectory.nc raw_trajectory.nc``.
   Note that the workflow may generate output files where it is run. Thus it may overwrite a file already present in a directory if its name matches an output file name.

Processing input files
~~~~~~~~~~~~~~~~~~~~~~

The command to start the processing would be as follows:

.. code-block:: bash

   mwf run -top raw_topology.prmtop -md replica_1 raw_trajectory_01.nc -md replica_2 raw_trajectory_02_part01.nc raw_trajectory_02_part02.nc

The ``mwf`` command is the main handler of the workflow and its letters stand for 'MDDB-WorkFlow'. The ``run`` subcommand is the one which triggers the main logic and thus run the workflow. Now lets explain every argument in this command:

* ``-top``: points to the input topology. No big surprises here.
* ``-md``: defines a new "MD" to be processed. You may use this argument multiple times. Then this argument is followed by:

  * MD directory: Directory to be created and filled with all the corresponding MD processed files and outputs.
  * Trajectories: Any number of trajectory files to be merged and then processed.

This command should process all input files and then stop because it is missing the "inputs file".
If everything is fine after the processing, some *new files* may have been generated and now the directory should look like this:

.. code-block:: text
    :emphasize-lines: 6,7,8,9,10,11,12

    project/
    ├── raw_topology.prmtop
    ├── raw_trajectory.nc
    ├── raw_trajectory_part01.nc
    ├── raw_trajectory_part02.nc
    ├── topology.prmtop
    ├── replica_1/
    │   ├── structure.pdb
    │   └── trajectory.xtc
    └── replica_2/
        ├── structure.pdb
        └── trajectory.xtc

Here we must explain a few things more.
First of all, note that trajectory parts have been automatically merged and a structure have been generated for every MD directory.
Note also that the not-raw (or just 'processed') topology is still in amber format (.prmtop) while the processed trajectory is converted to gromacs (.xtc) format. One of the standards in the workflow is that the trajectory is to be called 'trajectory.xtc' as well as the structure is to be called 'structure.pdb'. Thus the trajectory and the structure are to be in 'xtc' and 'pdb' formats respectively.

.. note::

   The workflow has a wide variety of tools and every tool supports a different range of input formats. Both xtc and pdb formats are extensively used and validated formats. Converting all input formats to a unique format to work along the workflow and thus not having to support every possible structure and trajectory format saves a lot of work.

   Conversions between topology formats are difficult however. For this reason the topology is not converted and then used as little as possible. The only processing in topologies is the atom filtering to keep them coherent with both structure and trajectory.

In this example we run the most basic processing, but there are a couple of additional features we may require.

* Filtering: Argument ``-filt`` is used to filter atoms away in all files (topology, structure and trajectory). The ``-filt`` argument alone applies the default filtering: water and counter ions. However the ``-filt`` argument may be followed by some text to apply a custom filtering selection according to |VMD syntax|. Here is an example on how to filter away everything which is not protein or nucleic acids: ``mwf run (...) -filt 'not protein or nucleic'``

.. |VMD syntax| raw:: html

    <a href="https://www.ks.uiuc.edu/Research/vmd/vmd-1.3/ug/node132.html" target="_blank">VMD syntax</a>

* Imaging and fitting: It is not easy to automatize the imaging process so it is recommended that you manually image your trajectories. However, the workflow is provided with a generic imaging protocol which may be useful in some generic cases. Use the ``-img`` argument to image and the ``-fit`` argument to fit the trajectory.

Once this process is over some tests are run.
If they all pass then we can continue with the analyses.

Running the analyses
~~~~~~~~~~~~~~~~~~~~

Before we start, we need the **inputs file**.
This file contains burocratic data, MD parameters and some additional metadata which is used by the workflow to adapt the analyses.

In order to generate this file, a template to build the file explaining every field in detail is provided. You can find it in the workflow repository, at |inputs_file_template.yml| or open it by simply running the following command:

.. |inputs_file_template.yml| raw:: html

    <a href="https://github.com/mmb-irb/MDDB-workflow/blob/release/model_workflow/resources/inputs_file_template.yml" target="_blank">inputs_file_template.yml</a>

.. code-block:: bash

   mwf inputs

Fill every field and then run the whole workflow to run all the analyses by just using the following command:

.. code-block:: bash

   mwf run

Note that no arguments for input files are provided now.
The workflow will remember which are the MD directories and processed files.
In addition, the workflow has some cache in these directories thus remembering some precalculated values and not having to repeat all the process.
Also, the workflow should me smart enough to recalculate any output if any of its implicated inputs change.

Note that this process will also generate some additional files such as 'metadata.json' and 'topology.json' by default. These files are also to be uploaded to the database.

If you want to run only a few specific analyses or exclude some analyses you can use the include (``-i``) and exclude (``-e``) arguments.
To see additional arguments and how to use them you can request help by just running ``mwf run -h``

Once you are done with this process is time to load your files to the database.
To do so, you must head to the |loader|.

.. |loader| raw:: html

    <a href="https://github.com/mmb-irb/MDDB-loader" target="_blank">MDDB-loader</a>


Tests and other checking processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some tests and checkings are run along the workflow to ensure data quality.
By default, a falied test will kill the workflow in the spot.
Tests may be skipped with the ``-t`` or ``--trust`` argument.
Tests may be allowed to fail with the ``-m`` or ``--mercy`` argument.
Either if a test is skipped or failed it will write a warning log in the metadata.
Note that some tests may not be skipped, but only allowed to fail and vice versa.
Both the trust and the mercy arguments have the same behaviour:
They may be followed by the name of the test to be skipped/allowed to fail.
If they are passed alone then they are applied to all available tests.
For example, if you want to allow bonds to be wrong, you may skip both the stable bonds test and the coherent bonds test:

.. code-block:: bash

   mwf run --trust stabonds cohbonds

These are the available tests:

- **stabonds** - Stable bonds test

  Check atom bonds to be stable along the trajectory.
  Bonds are mined from the topology when possible.
  Otherwise they are guessed by atom distance and radii (|VMD logic| under the hood) 10 frames along the trajectory are checked and bonds happening in most of these frames are the ones to be considered real. This way we avoid taking atom clashes as real bonds.
  Then, either if bonds are mined or guessed, we search for the first frame in the trajectory which respects all these bonds using again the atom distance and radii logic. This frame is set as the "reference frame" and it is usually found among the first 10 frames in the trajectory. If no frame matches all bonds after checking the first 100 frames then the test fails.
  If this test is skipped or allowed to fail then the first frame in the trajectory is set as the reference frame.

.. |VMD logic| raw:: html

    <a href="https://www.ks.uiuc.edu/Research/vmd/current/ug/node27.html" target="_blank">VMD logic</a>

- **cohbonds** - Coherent bonds

  Check number of bonds per atom to be coherent according to chemistry knowledge:

  * Hydrogen atoms must have 1 and only 1 bond
  * Oxygen atoms must have between 1 and 2 bonds
  * Nitrogen atoms must have between 1 and 4 bonds
  * Carbon atoms must have between 2 and 4 bonds
  * Sulfur atoms must have between 1 and 6 bonds
  * Phosphorus atoms must have between 2 and 6 bonds

  If any of these rules is not respected then the test fails.

- **intrajrity** - Trajectory integrity

  Make sure there are no sudden jumps in the middle of the trajectory due to imaging problems.
  Compute the RMSD between every pair of consecuitve frames in the trajectory.
  Then calculate the standard deviation among all RMSD differences.
  If there is at least one jump which is greater than 10 times the standard deviation then the test fails.
  First frames are excepcionally allowed to reach this limit since they may be part of the equilibration proccess.

- **elements** - Correct elements (skip only)

  Set wrong or missing atom elements to make them standard.
  This process relies in guessing mostly but it is smart enough to not set alpha carbons as calcium (Ca).
  This process logs warnings for every element which is guessed to be wrong.
  If this process is skipped it raises warning anyway but it keeps the original elements intact.

- **refseq** - Reference sequence match (mercy only)

  Make sure all protein chains are matched with their corresponding UniProt id.
  UniProt ids may be passed through the inputs file.
  Additional UniProt ids may be mined from the PDB in case PDB ids are passed through the inputs file.
  If there is a protein chain which finds no match among the available UniProt ids then a BLAST is run against reviewed UniProt sequences only.
  If the BLAST also fails then the matching process fails.
  If this process is allowed to fail then the unmatched protein region will remain with no UniProt reference.
  Note that proteins which are not to match anything such as antibodies or artifical constructs are to be tagged as 'no referable' in the inputs file.

- **interact** - Stable interactions (mercy only)

  Make sure defined interactions are actually happening and stable enough to be computed.
  In order to find interface atoms, interactions are checked to happend 100 frames along the trajectory.
  Two atoms are considered to be in the interface when they are close enough at least in one of these frames.
  If no interface atoms are found then the process is killed and this can not be allowed.
  This means you defined an interaction which does not exist and thus it must be removed from the inputs file.
  However, it may happen that an interaction has interface atoms but it is almost not happening.
  For example, a ligand which is fitted in a protein pocket but it leaves its place as soon as the trajectroy starts playing.
  For this reason, it is computed the number of frames that the interaction actually happens.
  If the interaction takes place less than the 10% of the total trajectory then the process fails.
  If it is allowed to fail then the wrong interaction is removed from metadata anyway and interaction analyses are not run for this specific interaction.

