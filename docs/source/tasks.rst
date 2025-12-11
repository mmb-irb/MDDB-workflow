.. _task_documentation: generated with generate_task_docs.py

Workflow Tasks
==================

This page documents the available tasks in the MDDB Workflow system.
These tasks can be specified with the ``-i`` (include) or ``-e`` (exclude) flags.

Project Tasks
---------------

These tasks are executed once per project:

* ``aiidata`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/tools/provenance.py#L8>`__: produce a provenance file containing AiiDA data adapted for our database

* ``chains`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/tools/chains.py#L111>`__: define the main function that will be called from the main script. This function will get the parsed chains from the structure and request the InterProScan service to obtain the data for each chain.

* ``charges`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/tools/get_charges.py#L13>`__: extract charges from a source file. 

* ``inchikeys`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/tools/get_inchi_keys.py#L120>`__: generate a dictionary mapping InChI keys to residue information for non-standard residues.  This function uses MDAnalysis to parse the input structure and topology files and identifies residues that are not classified as 'ion', 'solvent', 'nucleic', or 'protein'. For each identified residue, it converts the structure to RDKit format to obtain the InChI key and InChI string. The resulting data is stored in dictionaries to map InChI keys to residue details and residue names to InChI keys. PDB coordinates are necesary to distinguish stereoisomers. 

* ``inchimap`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/tools/get_inchi_keys.py#L249>`__: generate InChI references for the database.

* ``inputs`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/mwf.py#L1498>`__: set a function to load the inputs yaml file.

* ``itopology`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/mwf.py#L1579>`__: get the input topology file. If the file is not found try to download it.

* ``ligmap`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/tools/get_ligands.py#L40>`__: generate a map of residues associated to ligands.  This function identifies and maps ligands in the molecular structure through a multi-step matching process:  1. **Direct InChIKey matching**: Extracts InChIKeys from structure fragments    (excluding lipids and membrane components) and attempts direct matching with    PubChem database.  2. **Chemical similarity matching**: If direct matching fails, progressively    modifies the molecular structure and calculates Tanimoto coefficient (TC)    for similarity assessment:    1. Neutralize charges    2. Remove stereochemistry information    3. Apply PubChem standardization (tautomer, protonation, etc.)    4. Match against PDB-derived ligands (TC threshold ≥ 0.9)    5. Perform similarity search in PubChem/ChEMBL (TC threshold ≥ 0.9)  3. **Fallback handling**: Unmatched ligands are saved as-is with warnings.  4. **User-forced selections**: Respects user-specified ligand selections from    inputs.yaml, with warnings if TC compared to original fragment is insufficient.

* ``lipmap`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/tools/get_lipids.py#L7>`__: add lipid-specific database information to InChIKeyData objects.  This function queries SwissLipids and LIPID MAPS databases for each InChI key and adds the results directly to the InChIKeyData objects. It also performs quality checks on lipid classifications. 

* ``memmap`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/tools/membrane_mapping.py#L10>`__: generate a list of residue numbers of membrane components from a given structure and topology file. 

* ``pdbs`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/tools/generate_pdb_references.py#L8>`__: prepare the PDB references json file to be uploaded to the database.

* ``pmeta`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/tools/generate_metadata.py#L13>`__: prepare a JSON file with all project metadata.

* ``populations`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/mwf.py#L1642>`__: get the MSM equilibrium populations file.

* ``protmap`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/tools/generate_map.py#L75>`__: map the structure aminoacids sequences against the Uniprot reference sequences.

* ``refbonds`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/tools/get_bonds.py#L251>`__: find reference safe bonds in the system. First try to mine bonds from a topology files. If the mining fails then search for the most stable bonds. If we trust in stable bonds then simply return the structure bonds.

* ``resmap`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/tools/residue_mapping.py#L5>`__: build the residue map from both proteins and ligands maps.  This is formatted as both the standard topology and metadata generators expect them.

* ``screenshot`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/tools/get_screenshot.py#L25>`__: obtain a screenshot from the pdb file using VMD. This screenshot of the system is uploaded to the database. Returns the rotation values used to take the photo so they can be saved and reused.

* ``stopology`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/tools/generate_topology.py#L5>`__: prepare the standard topology file to be uploaded to the database.

* ``topology`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/mwf.py#L1867>`__: get the processed topology from the reference MD.

* ``transitions`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/mwf.py#L1649>`__: get the MSM transition probabilities file.

MD Tasks
-----------

These tasks are executed for each MD in the project:

Files
~~~~~~~~

* ``inpro`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/tools/process_input_files.py#L46>`__: process input files to generate the processed files. This process corrects and standarizes the topology, the trajectory and the structure.

* ``istructure`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/mwf.py#L292>`__: get the input pdb filename from the inputs. If the file is not found try to download it.

* ``itrajectory`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/mwf.py#L438>`__: get the input trajectory filename(s) from the inputs. If file(s) are not found try to download it.

* ``structure`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/utils/tasks.py#L121>`__: 

* ``trajectory`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/utils/tasks.py#L121>`__: 

Analyses
~~~~~~~~~~~~~~

* ``apl`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/analyses/area_per_lipid.py#L12>`__: area per lipid analysis.

* ``average`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/utils/pyt_spells.py#L166>`__: get an average structure from a trajectory.

* ``channels`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/analyses/channels.py#L10>`__: analyze channels in a membrane protein using MDAnalysis' mda_hole.

* ``clusters`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/analyses/clusters.py#L15>`__: run the cluster analysis.

* ``density`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/analyses/density.py#L8>`__: membrane density analysis.

* ``dihedrals`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/analyses/dihedral_energies.py#L10>`__: calculate torsions and then dihedral energies for every dihedral along the trajectory.

* ``dist`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/analyses/distance_per_residue.py#L22>`__: calculate the distance mean and standard deviation of each pair of residues of different agents. Note that the distances are calculated for all residues in the agent, not only the interface residues.

* ``energies`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/analyses/energies.py#L42>`__: perform the electrostatic and vdw energies analysis for each pair of interaction agents.

* ``firstframe`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/tools/get_first_frame.py#L6>`__: get the trajectory first frame in PDB format using Gromacs.

* ``frames`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/utils/pyt_spells.py#L75>`__: get the trajectory frames count.

* ``hbonds`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/analyses/hydrogen_bonds.py#L26>`__: perform an hydrogen bonds analysis for each interaction interface. The 'interactions' input may be an empty list (i.e. there are no interactions). In case there are no interactions the analysis stops. Note that this analysis find hydrogen bonds in a subset of frames along the trajectory. Storing the results for the whole trajectory is not possible due to storage limits.

* ``helical`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/analyses/helical_parameters.py#L115>`__: helical parameters analysis.

* ``inter`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/tools/process_interactions.py#L37>`__: find the residues of each interacting agent. It can automatically detect interactions based on chain names or ligand information, or use a predefined list of interactions.

* ``linter`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/analyses/lipid_interactions.py#L9>`__: lipid-protein interactions analysis.

* ``lorder`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/analyses/lipid_order.py#L9>`__: calculate lipid order parameters for membranes. This function computes the order parameter (S) for lipid acyl chains, defined as: S = 0.5*(3*<cos²θ> - 1) where θ is the angle between the C-H bond and the membrane normal (z-axis). 

* ``markov`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/analyses/markov.py#L9>`__: set the data needed to represent a Markov State Model graph in the client. This is finding the most populated frames and calculating an RMSD matrix between these frames.

* ``mda_univ`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/utils/mda_spells.py#L89>`__: create a MDAnalysis universe using data in the workflow.

* ``mdmeta`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/tools/generate_metadata.py#L268>`__: produce the MD metadata file to be uploaded to the database.

* ``pairwise`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/analyses/rmsd_pairwise.py#L18>`__: perform an analysis for the overall structure and then one more analysis for each interaction.

* ``pca`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/analyses/pca.py#L13>`__: perform a PCA analysis on the trajectory.

* ``perres`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/analyses/rmsd_per_residue.py#L11>`__: perform the RMSD analysis for each residue.

* ``pockets`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/analyses/pockets.py#L45>`__: perform the pockets analysis.

* ``reframe`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/tools/get_bonds.py#L121>`__: return a reference frame number where all bonds are exactly as they should (by VMD standards). This is the frame used when representing the MD.

* ``rgyr`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/analyses/rgyr.py#L15>`__: perform the RMSd analysis. Use the first trajectory frame in .pdb format as a reference.

* ``rmsds`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/analyses/rmsds.py#L11>`__: run multiple RMSD analyses. One with each reference (first frame, average structure)  and each selection (default: protein, nucleic).

* ``rmsf`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/analyses/rmsf.py#L16>`__: perform the fluctuation analysis.

* ``sas`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/analyses/sasa.py#L15>`__: perform the Solvent Accessible Surface Analysis.

* ``thickness`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/analyses/thickness.py#L10>`__: membrane thickness analysis.

* ``tmscore`` `[source] <https://github.com/mmb-irb/MDDB-workflow/blob/master/mddb_workflow/analyses/tmscores.py#L24>`__: perform the tm score using the tmscoring package.

Task Groups
-------------

These are predefined groups of tasks that can be specified with a single flag.

* ``download``: ``itopology``, ``inputs``, ``populations``, ``transitions``, ``istructure``, ``itrajectory``.

* ``setup``: ``topology``, ``structure``, ``trajectory``.

* ``meta``: ``pmeta``, ``mdmeta``.

* ``network``: ``resmap``, ``ligmap``, ``lipmap``, ``chains``, ``pdbs``, ``memmap``.

* ``minimal``: ``pmeta``, ``mdmeta``, ``stopology``.

* ``interdeps``: ``inter``, ``pairwise``, ``hbonds``, ``energies``, ``perres``, ``clusters``, ``dist``.

* ``membs``: ``memmap``, ``density``, ``thickness``, ``apl``, ``lorder``, ``linter``, ``channels``.

