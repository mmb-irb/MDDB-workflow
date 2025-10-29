=========
Changelog
=========

All notable changes to this project will be documented in this file.

The format is based on `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_,
and this project adheres to `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_.

[0.1.6] - 2025-10-29
==========

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