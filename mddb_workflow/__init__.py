import os
import sys
import subprocess

# When the mwf command is launched, check the installed dependencies before importing any of them
# This way an old or missing package is reported clearly instead of failing later in the middle of the workflow
# The check can be skipped by setting the MWF_SKIPDEPS environment variable (e.g. MWF_SKIPDEPS=1)
if os.path.basename(sys.argv[0]) == 'mwf' and not os.environ.get('MWF_SKIPDEPS'):
    from mddb_workflow.utils.dependencies import check_dependencies
    check_dependencies()

from rdkit import RDLogger
# Mute RDKit warnings
lg = RDLogger.logger()
lg.setLevel(RDLogger.ERROR)

name = "mddb_workflow"
__all__ = ["analysis", "core", "resources", "tools", "utils"]

try:
    from importlib.metadata import version
    __version__ = version("mddb_workflow")
except Exception:
    __version__ = "unknown"

# Workaround to fix ReadTheDocs build issue with conda environments
# when fix_gromacs_masses and helical_parameters are called.
if "CONDA_PREFIX" not in os.environ:
    try:
        # Get the path of the python executable
        python_path = subprocess.check_output(["which", "python"], text=True).strip()
        # Extract the directory path (one level up from bin/python)
        conda_prefix = os.path.dirname(os.path.dirname(python_path))
        # Set the environment variable
        os.environ["CONDA_PREFIX"] = conda_prefix
    except subprocess.SubprocessError:
        print("Warning: Could not determine CONDA_PREFIX from 'which python'")
