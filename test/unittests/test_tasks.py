"""Test for Task class to check for bugs related to tasks that don't generate output files.

Bug scenario: When a task like protmap runs but doesn't generate an output file
(e.g., because there are no protein sequences to map), on the next run the task
logic fails when checking for satisfied output because the expected file doesn't exist.
"""

import os
import shutil
import pathlib
from unittest.mock import patch
from mddb_workflow.mwf import Project

# Set up paths
data_dir = pathlib.Path(__file__).parent.parent / 'data'
dummy_dir = data_dir / 'input/dummy'
test_fld = data_dir / 'output/test_tasks'


def regenerate_test_fld():
    """Remove old test folder and create a fresh one."""
    shutil.rmtree(test_fld, ignore_errors=True)
    test_fld.mkdir(parents=True, exist_ok=True)


def test_task_no_output_file_second_run():
    """Test that tasks which don't generate output files work correctly on subsequent runs."""
    regenerate_test_fld()

    # Copy necessary files
    shutil.copy(dummy_dir / "inputs.yaml", test_fld)
    shutil.copy(dummy_dir / "gromacs/ala_ala.tpr", test_fld)
    (test_fld / "replica_1").mkdir(exist_ok=True)
    shutil.copy(dummy_dir / "gromacs/trajectory.xtc", test_fld / "replica_1")

    cwd = pathlib.Path.cwd()
    os.chdir(test_fld)

    try:
        # Mock generate_protein_mapping to return an empty list (simulating no protein sequences)
        # This causes the task to complete without creating the output file
        # We need to patch the Task's func attribute directly since it stores the function reference at class definition time
        # The mock must have the same signature as the original function
        def mock_generate_protein_mapping(
            structure,
            protein_references_file,
            database,
            cache,
            register,
            mercy=[],
            input_protein_references=[],
            pdb_ids=[],
        ):
            """Mock that simulates no protein sequences found."""
            return []

        with patch.object(Project.get_protein_map, 'func', mock_generate_protein_mapping):
            # First run: Initialize Project and run protmap task
            print("\n ------ First run: protmap with no protein sequences ------ \n")
            project = Project(
                directory=str(test_fld),
                input_topology_filepath=str(test_fld / "ala_ala.tpr"),
                input_trajectory_filepaths=str(test_fld / "replica_1/trajectory.xtc"),
                md_directories=['replica_1']
            )
            # Trigger the protmap task
            result1 = project.protein_map
            print(f"First run result: {result1}")
            assert result1 == [], "First run should return empty list when no proteins"

        # Second run: Same project, protmap task should use cached result
        # This is where the bug manifests - the task tries to check if output file exists
        # but it was never created, causing the satisfied_output check to fail
        with patch.object(Project.get_protein_map, 'func', mock_generate_protein_mapping):
            print("\n ------ Second run: protmap should use cache ------ \n")
            project2 = Project(
                directory=str(test_fld),
                input_topology_filepath=str(test_fld / "ala_ala.tpr"),
                input_trajectory_filepaths=str(test_fld / "replica_1/trajectory.xtc"),
                md_directories=['replica_1']
            )
            # This should not fail - it should use the cached result
            result2 = project2.protein_map
            assert not hasattr(project2.get_protein_map, 'cache_cksums'), "Cache checksums not should exist"

    finally:
        os.chdir(cwd)


if __name__ == "__main__":
    test_task_no_output_file_second_run()
