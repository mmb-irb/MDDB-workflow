.. _task_documentation: generated with generate_task_docs.py

Workflow Tasks
==================

This page documents the available tasks in the MDDB Workflow system.
These tasks can be specified with the ``-i`` (include) or ``-e`` (exclude) flags.

Project Tasks
---------------

These tasks are executed once per project:

* ``chains``: define the main function that will be called from the main script. This function will get the parsed chains from the structure and request the InterProScan service to obtain the data for each chain.

* ``charges``: extract charges from a source file. 

* ``inputs``: set a function to load the inputs file

* ``intertypes``: given a list of interactions, find their types.

* ``itopology``: get the input topology file. If the file is not found try to download it.

* ``ligmap``: generate a map of residues associated to ligands.

* ``lipmap``: generate the lipid references.

* ``mda_univ``: create a MDAnalysis universe using data in the workflow.

* ``membranes``: generates a list of residue numbers of membrane components from a given structure and topology file.

* ``pdbs``: prepare the PDB references json file to be uploaded to the database.

* ``pmeta``: prepare a JSON file with all project metadata.

* ``populations``: get the MSM equilibrium populations filename.

* ``protmap``: map the structure aminoacids sequences against the Uniprot reference sequences.

* ``refbonds``: find reference safe bonds in the system.

* ``resmap``: build the residue map from both proteins and ligands maps This is formatted as both the standard topology and metadata generators expect them.

* ``screenshot``: obtain a screenshot from the pdb file using VMD. This screenshot of the system is uploaded to the database. Returns the rotation values used to take the photo so they can be saved and reused. 

* ``stopology``: prepare the standard topology file to be uploaded to the database.

* ``topology``: get the processed topology file.

* ``transitions``: get the MSM transition probabilities filename.

MD Tasks
-----------

These tasks are executed for each MD in the project:

Files
~~~~~~~~

* ``inpro``: process input files to generate the processed files. This process corrects and standarizes the topology, the trajectory and the structure.

* ``istructure``: get the input pdb filename from the inputs. If the file is not found try to download it.

* ``itrajectory``: get the input trajectory filename(s) from the inputs. If file(s) are not found try to download it.

* ``structure``: get the processed structure.

* ``trajectory``: get the processed trajectory.

Analyses
~~~~~~~~~~~~~~

* ``apl``: area per lipid analysis.

* ``average``: get an average structure from a trajectory.

* ``clusters``: run the cluster analysis.

* ``density``: membrane density analysis.

* ``dihedrals``: calculate torsions and then dihedral energies for every dihedral along the trajectory.

* ``dist``: calculate the distance mean and standard deviation of each pair of residues of different agents. Note that the distances are calculated for all residues in the agent, not only the interface residues.

* ``energies``: perform the electrostatic and vdw energies analysis for each pair of interaction agents.

* ``firstframe``: get the trajectory first frame in PDB format using Gromacs.

* ``frames``: get the trajectory frames count.

* ``hbonds``: perform an hydrogen bonds analysis for each interaction interface. The 'interactions' input may be an empty list (i.e. there are no interactions). In case there are no interactions the analysis stops. Note that this analysis find hydrogen bonds in a subset of frames along the trajectory. Storing the results for the whole trajectory is not possible due to storage limits.

* ``helical``: helical parameters analysis.

* ``inter``: find the residues of each interacting agent. It can automatically detect interactions based on chain names or ligand information, or use a predefined list of interactions.

* ``linter``: lipid-protein interactions analysis.

* ``lorder``: calculate lipid order parameters for membranes. This function computes the order parameter (S) for lipid acyl chains, defined as: S = 0.5*(3*<cos²θ> - 1) where θ is the angle between the C-H bond and the membrane normal (z-axis). 

* ``markov``: set the data needed to represent a Markov State Model graph in the client. This is finding the most populated frames and calculating an RMSD matrix between these frames.

* ``mdmeta``: produce the MD metadata file to be uploaded to the database.

* ``pairwise``: perform an analysis for the overall structure and then one more analysis for each interaction.

* ``pca``: perform a PCA analysis on the trajectory.

* ``perres``: perform the RMSD analysis for each residue.

* ``pockets``: perform the pockets analysis.

* ``reframe``: return a canonical frame number where all bonds are exactly as they should. This is the frame used when representing the MD.

* ``rgyr``: perform the RMSd analysis. Use the first trajectory frame in .pdb format as a reference.

* ``rmsds``: run multiple RMSD analyses. One with each reference (first frame, average structure)  and each selection (default: protein, nucleic).

* ``rmsf``: perform the fluctuation analysis.

* ``sas``: perform the Solvent Accessible Surface Analysis.

* ``thickness``: membrane thickness analysis.

* ``tmscore``: perform the tm score using the tmscoring package.

Task Groups
-------------

These are predefined groups of tasks that can be specified with a single flag.

* ``download``: ``itopology``, ``inputs``, ``populations``, ``transitions``, ``istructure``, ``itrajectory``.

* ``setup``: ``topology``, ``structure``, ``trajectory``.

* ``meta``: ``pmeta``, ``mdmeta``.

* ``network``: ``mapping``, ``ligands``, ``chains``, ``pdbs``, ``membrane``.

* ``minimal``: ``pmeta``, ``mdmeta``, ``stopology``.

* ``interdeps``: ``interactions``, ``pairwise``, ``hbonds``, ``energies``, ``perres``, ``clusters``, ``dist``.

* ``membs``: ``membranes``, ``density``, ``thickness``, ``apl``, ``lorder``, ``linter``.

