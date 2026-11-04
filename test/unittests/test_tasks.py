"""Test for Task class to check for bugs related to tasks that don't generate output files.

Bug scenario: When a task like protmap runs but doesn't generate an output file
(e.g., because there are no protein sequences to map), on the next run the task
logic fails when checking for satisfied output because the expected file doesn't exist.
"""

import os
import shutil
import pathlib
import pytest
from unittest.mock import patch, create_autospec
from mddb_workflow.mwf import Project, MD
from mddb_workflow.tools.generate_map import generate_protein_mapping
from mddb_workflow.tools.process_input_files import process_input_files

# Set up paths
data_dir = pathlib.Path(__file__).parent.parent / 'data'
dummy_dir = data_dir / 'input/dummy'
test_fld = data_dir / 'output/test_tasks'
# The box of the dummy Gromacs trajectory is smaller than the peptide, so it touches its own periodic image
DUMMY_GROMACS_TRUST = ['pbcontacts']


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
    shutil.copy(dummy_dir / "gromacs/raw_trajectory.xtc", test_fld / "replica_1")

    cwd = pathlib.Path.cwd()
    os.chdir(test_fld)

    try:
        # Mock generate_protein_mapping to return an empty list (simulating no protein sequences)
        # This causes the task to complete without creating the output file
        # We need to patch the Task's func attribute directly since it stores the function reference at class definition time
        # Use create_autospec to copy the signature from the original function
        mock_generate_protein_mapping = create_autospec(
            generate_protein_mapping, return_value=[]
        )

        with patch.object(Project.get_protein_map, 'func', mock_generate_protein_mapping):
            # First run: Initialize Project and run protmap task
            print("\n ------ First run: protmap with no protein sequences ------ \n")
            project = Project(
                directory=str(test_fld),
                input_topology_filepath=str(test_fld / "ala_ala.tpr"),
                input_trajectory_filepaths=str(test_fld / "replica_1/raw_trajectory.xtc"),
                input_md_directories=['replica_1'],
                trust=DUMMY_GROMACS_TRUST
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
            # Reset the mock call count before second run
            mock_generate_protein_mapping.reset_mock()

            project2 = Project(
                directory=str(test_fld),
                input_topology_filepath=str(test_fld / "ala_ala.tpr"),
                input_trajectory_filepaths=str(test_fld / "replica_1/raw_trajectory.xtc"),
                input_md_directories=['replica_1'],
                trust=DUMMY_GROMACS_TRUST
            )
            # This should not fail - it should use the cached result
            result2 = project2.protein_map
            mock_generate_protein_mapping.assert_not_called()
            assert result2 == [], "Second run should return cached empty list"

    finally:
        os.chdir(cwd)


class SimulatedInterruption(Exception):
    """Exception used to simulate the workflow being killed in the middle of a task."""


def test_md_task_overwrite_keeps_project_output():
    """Test that overwriting an MD task does not delete outputs shared with other MDs.

    Bug scenario: The inpro task of every MD outputs the project topology when all MDs share the same input topology.
    When the inpro of one MD was overwritten (e.g. because its input trajectory changed) the project topology was deleted.
    If the workflow was interrupted before the task finished then the topology was left missing.
    Then, in the next run, the inpro task of the first MD was run again with no input changes.
    """
    regenerate_test_fld()

    shutil.copy(dummy_dir / "inputs.yaml", test_fld)
    shutil.copy(dummy_dir / "gromacs/ala_ala.tpr", test_fld)
    md_directories = ['replica_1', 'replica_2']
    for md_directory in md_directories:
        (test_fld / md_directory).mkdir(exist_ok=True)
        shutil.copy(dummy_dir / "gromacs/raw_trajectory.xtc", test_fld / md_directory)

    def get_project() -> Project:
        return Project(
            directory=str(test_fld),
            input_topology_filepath=str(test_fld / "ala_ala.tpr"),
            input_trajectory_filepaths="raw_trajectory.xtc",
            input_md_directories=md_directories,
            trust=DUMMY_GROMACS_TRUST
        )

    cwd = pathlib.Path.cwd()
    os.chdir(test_fld)

    try:
        # First run: process the input files of both MDs
        project = get_project()
        md1, md2 = project.mds
        md1.topology_file
        md2.topology_file
        project_topology_file = md1.topology_file
        assert md2.topology_file == project_topology_file, "Both MDs should share the project topology"
        assert project_topology_file.exists

        # Second run: overwrite the inpro of the second MD and interrupt it in the middle
        interrupted_processing = create_autospec(process_input_files, side_effect=SimulatedInterruption)
        with patch.object(MD.input_files_processing, 'func', interrupted_processing):
            project = get_project()
            md1, md2 = project.mds
            md2.overwritables.add(MD.input_files_processing.flag)
            with pytest.raises(SimulatedInterruption):
                md2.topology_file
        assert project_topology_file.exists, "The project topology should not be deleted by an MD task"

        # Third run: the inpro of the first MD must not be run again since nothing changed
        unexpected_processing = create_autospec(process_input_files)
        with patch.object(MD.input_files_processing, 'func', unexpected_processing):
            project = get_project()
            md1, md2 = project.mds
            md1.topology_file
            unexpected_processing.assert_not_called()

    finally:
        os.chdir(cwd)


if __name__ == "__main__":
    test_task_no_output_file_second_run()
    test_md_task_overwrite_keeps_project_output()
