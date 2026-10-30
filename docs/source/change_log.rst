=========
Changelog
=========

All notable changes to this project will be documented in this file.

The format is based on `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_,
and this project adheres to `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_.

[0.1.12] - 2026-10-30
=====================

Added
-----
- Add inputs validation with Pydantic and autofill of the inputs file
- Add the ``mwf upref`` command to update references, using the loader and a cache
- Add the ``mwf uproj`` command to update outdated projects, including their membrane mapping
- Add the ``mwf view`` and ``mwf ds view`` commands to open VMD with the most processed structure and trajectory, with an atom selection option and the PBC box shown
- Add the ``mwf ds logs`` and ``mwf ds load`` commands, the ``-qi`` filter by inputs fields and the ``--get`` option of ``mwf ds inputs``
- Allow ``mwf ds inputs`` to edit or create a single field and to list the existing inputs files
- Reuse the job script from a previous run in ``mwf ds run``
- Add the ``ua_selection`` and ``qm_selection`` inputs and guess coarse grained atoms with ``cg_selection: auto``
- Add the ``forced_lipids`` input and the lipids mercy flag
- Add a flag to skip frames, a local BLAST flag (``-lb``) and BLAST as optional dependency
- Allow some console-only arguments in the inputs file and a manual chainer in ``mwf run``
- Add the simulation box and cross-PBC contacts checks
- Add an early PCA-based clustering prototype
- Add early ComPASS (allostery) and LEGOLAS (chemical shifts) implementations, not yet in conda
- Classify amphiphilic aggregates as bilayers, vesicles or micelles, with midplane leaflets, and map coarse grained lipids
- Recognise lipids from their structure with SMARTS of the LIPID MAPS categories
- Add automatic multimerism, all atom / coarse grained, nucleic and protein-nucleic complex system keywords
- Note protein and nucleic sequence lengths and RMSD time splits in the MD metadata
- Add CHARMM and Amber nucleic acid residue names
- Add descriptions, client data and VMD scripts to warnings
- Get InChI keys from MDDB, run InChI key generation in parallel and version InChI key references
- Export the heaviest topology fields to binary files when they do not fit in Mongo
- Write structures in mmCIF and name chains with no limit (big reform)

Changed
-------
- Replace the Gromacs no-jump with our own no-jump which keeps molecules together
- Add an automatic translation step and smarter defaults to the imaging protocol
- Keep the box rotated after fitting and ignore it when the trajectory was fitted without rotation
- Fit the trajectory integrity bimodality with gammas on jumps scaled per bond, support rotated boxes and detect permanent fragment breaks
- Find leaflets from local lipid normals inside hydrophobic cores, like FATSLiM, and run LeafletFinder with periodic boundaries
- Calculate lipid order only for flat bilayers, reading every coarse grained tail bead from gorder
- Count all lipid contacts with numpy and reduce memory usage
- Bump the versions of the membrane analyses changed by the new membrane mapping
- Rework the interactions and chemical shifts output formats
- Select automatic interaction agents by chain index
- Smarter protein classification, solvent identification and PDB model filtering
- Relax the stable bonds check for big systems and trajectories with few frames
- Move file and CLI inputs handling to its own class
- Create the register only for ``mwf run`` calls
- Show the workflow version also when the closest tag is not annotated
- Require biobb_mem 5.3.1

Fixed
-----
- Fix several issues with missing topologies, missing atoms, dummy atoms and Gromacs masses
- Fix energies with dummy atoms and agents moving away
- Fix trajectories with no box crashing the integrity checks
- Stop before filtering a TPR when Gromacs convert-tpr would split its residues
- Filter the input structure when it has more atoms than the trajectory and there is no topology
- Raise an input error when the filter selection selects no atoms
- Fail when Gromacs exits with error or leaves an empty output file
- Fix reading mmCIF files with missing values, missing elements or repeated residue labels
- Split overflowed '****' residue numbers into separate residues when reading PDBs
- Classify nucleic residues named like both DNA and RNA from their context
- Fix the dataset lock not shared across nodes and retry transient SQLite errors
- Fix the dataset 'waiting' status and the ``ds watch`` options
- Retry remote project requests, downloads and database requests when the network or the server fails
- Cache PubChem requests and stop treating request errors as missing data
- Recover the registered UUID when a cache is lost
- Leave workflow output files out of input glob patterns
- Check if outputs will fit in Mongo


[0.1.11] - 2026-06-04
=====================

Added
-----
- Re-implement the channels analysis with CHAP and add it to the default analyses
- Add the custom metadata additions input
- Add cg, dummy and forced class selections per MD
- Add a global node URL option
- Add Li, W and Ta elements

Changed
-------
- Use PDB numbering and author notation in PDB references
- Write auxiliary files to the specified ``--dir``
- Load the first mmCIF model by default and trust mmCIF elements
- Use strict numeric values for the temperature, frame step and time step inputs
- Skip coarse grained interactions in energies

Fixed
-----
- Restore support for removed MDs
- Fix thickness for trajectories with negative z
- Fix carbon recognized as cobalt and nitrogen as nickel
- Better handle obsolete PDB ids and not found UniProt references
- Throttle PubChem requests and retry failed requests, including SSL handshake timeouts
- Raise an error for trajectories with 0 frames


[0.1.10] - 2026-04-28
=====================

Fixed
-----
- Fix module import issues in CLI and dataset commands


[0.1.9] - 2026-04-27
=====================

Added
-----
- Add dataset summary, status, scan, inputs, and remove workflow commands
- Add project and MD log links in dataset views
- Add waiting state handling for SLURM jobs
- Add membrane analysis support for united-atom and coarse-grained systems

Changed
-------
- Refactor dataset storage around SQLite-backed projects and MD entries
- Update dataset behavior to overwrite entries with the same ``rel_path`` by default and hide UUIDs in ``show``
- Improve topology, ligand, and reference handling, including missing bonds, dummy atoms, and lipid/PubChem lookups
- Improve workflow task naming, task tracking, and error handling across the console
- Switch PDB trajectory handling to PyTraj to reduce memory usage

Fixed
-----
- Fix dataset summaries, empty MD inputs, and inconsistent input handling
- Fix SAS and solvent handling in structure and screenshot generation
- Fix ligand and reference edge cases, including null Mordred values and missing PubChem mappings
- Fix runtime robustness issues around retries, timeouts, multiprocessing, and database locking

[0.1.8] - 2025-12-03
=====================

Added
-----
- Add ``--ssleep`` flag to skip SSL certificate authentication

Changed
-------
- Update mordred dependency
- Update biobb_mem to 5.1.2
- Add ``n_jobs`` option to dataset command
- Add ignore option for dataset config directories

Fixed
-----
 - Fix clustering cutoff when low RMSD difference


[0.1.7] - 2025-11-28
=====================

Added
-----
- Support for URLs as input files
- Test for network flag and top_no to CI
- Ignore bonds flag (``-ib``) for handling systems without bonds
- Support for virtual site "pseudo-bonds"
- Support for command line forced inputs
- RMSD check per fragment functionality
- Enhanced InChIKey handling with PubChem CID retrieval

Changed
-------
- Refactored ligand and lipid references to use InChIKeys (#42)
- Renamed ligands module for better organization
- Changed dataset subcommand from status to groups
- Changed source code license to Apache 2
- Updated HPC instructions (#30)
- Sample trajectory (smp) now picks frames along trajectory instead of at start
- Better function checksums implementation

Fixed
-----
- Fixed PCA memory inefficiency
- Patched MDTraj problem in PCA
- Fixed TM score no anchors problem
- Safe trajectory download implementation
- Fixed issues with missing topology support
- Fixed errors in handling missing bonds throughout codebase
- Fixed handling of missing bonds when finding PTMs
- Fixed CG-related issues
- Improved pdb_to_uniprot function for nucleic acids

Changed
-------
- Use networkx for cycle detection


[0.1.6] - 2025-10-29
=====================

Added
-----
- Dataset documentation and ``generate_inputs_yaml`` method
- Support for missing bonds in molecular systems
- Error log handling and clickable links for logs in Jupyter display
- Dataset subcommand to CLI
- Status dataset command functionality

Changed
-------
- Renamed ``model_workflow`` to ``mddb_workflow`` (#22)
- Updated to Python 3.11 (#41)
- Updated CI Python version
- Added ``cg_test2`` to CI testing
- Fixed recursive transformer functionality
- Improved exception handling with standardized JSON support

Fixed
-----
- Imaging repeated uselessly bug
- Various bugs in molecular system handling
- SLURM job path issues
- Biotite tmscore protein check fixes
- Logging color removal with ``-nc`` flag

[0.1.5] - 2025-10-23
====================

Added
-----
- License file and proper licensing (#16)
- Manual API tools for better automation
- Dummy data for enhanced unit testing
- Channels analysis functionality
- Membrane analysis for coarse-grained (CG) systems
- Plot functions for membrane analysis results
- InChI key support for ligands
- SwissLipids names to lipid references
- Structure classes to documentation
- Parallel processing for InChI key generation

Changed
-------
- Updated minimum Python version requirement to >=3.10.0
- Enhanced type hints and docstring coverage across codebase
- Improved input file processing to be more task-like
- Moved comments to docstrings for better documentation
- Refactored lipid analysis functions to use MDAnalysis Universe directly
- Changed tmscoring implementation to use biotite
- Updated FATSLiM compatibility with numpy 1.26
- Improved structure corrector messages
- Enhanced interactions logic with automatic detection

Fixed
-----
- Critical fixes in Gromacs integration
- Gromacs masses fixer improvements
- Cached exceptions handling
- Automatic no-reference (noref) detection bugs
- Interactions reform with type-critical bug fixes
- Biotite TM score with non-standard amino acids
- Input trajectory overwriting issues
- Topology and trajectory input filepath handling
- UniProt secure connection failures
- Various path handling improvements
- Glucolipids linter fixes

Removed
-------
- Unnecessary dependencies for improved performance
- Path filters from release workflow trigger

Security
--------
- Improved secure connection handling for external APIs

[0.1.4] - 2025-08-11
====================

Changed
-------
- Fixed version name handling
- Added GHCR (GitHub Container Registry) support
- Improved Slack notification system

[0.1.3] - 2025-08-11
====================

Added
-----
- GitHub Container Registry (GHCR) integration

[0.1.2] - 2025-08-11
====================

Added
-----
- GitHub Container Registry (GHCR) support

[0.1.1] - 2025-08-11
====================

Added
-----
- GitHub Container Registry (GHCR) functionality

[0.1.0] - 2025-07-23
====================

This is the first major release of the MDDB workflow package.

Added
-----
- Core molecular dynamics analysis workflow
- Support for multiple MD simulation formats
- Automated topology generation
- Ligand and lipid reference generation
- PDB reference handling
- Metadata generation capabilities
- Structure correction and validation
- Interaction type detection
- Comprehensive testing suite
- Documentation system
- CI/CD pipeline

[0.0.2] - 2025-07-21
====================

Changed
-------
- Simplified version extraction from ``pyproject.toml``

[0.0.1] - 2025-04-14
====================

Added
-----
- Initial release
- Basic project structure
- Core functionality implementation
- Include subfolders support

.. note::
   This changelog was generated from the git commit history on 2025-10-29.
   For detailed information about specific changes, please refer to the git commit history.