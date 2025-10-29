import os
import subprocess

name = "mddb_workflow"
__all__ = ["resources", "tools", "utils"]

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