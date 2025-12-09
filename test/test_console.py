import subprocess
import pytest
import shutil
import sys
import os
from io import StringIO
from unittest.mock import patch
from mddb_workflow.console import parser, main
from mddb_workflow.utils.type_hints import *


@pytest.mark.CI
@pytest.mark.release
class TestConsoleArgumentParsing:
    """Test the console argument parsing functionality."""

    def test_parser_exists(self):
        """Verify that the parser is properly defined."""
        assert parser is not None
        assert hasattr(parser, 'parse_args')

    def test_no_arguments(self):
        """Test behavior when no arguments are provided."""
        # Set up empty arguments
        sys.argv = ['mddb_workflow']

        # Capture the output of parser.print_help() for comparison
        expected_help_buffer = StringIO()
        with patch('sys.stdout', expected_help_buffer):
            parser.print_help()
            expected_output = expected_help_buffer.getvalue()

        # Run the main function - should print help
        expected_help_buffer = StringIO()
        with patch('sys.stdout', expected_help_buffer):
            main()
            actual_output = expected_help_buffer.getvalue()

        # Compare outputs, ignoring whitespace differences
        assert actual_output.strip() == expected_output.strip(), \
            f"Expected:\n{expected_output}\n\nActual:\n{actual_output}"


@pytest.mark.CI
@pytest.mark.release
class TestConsoleIntegration:
    """Integration tests for console functionality."""

    @pytest.mark.parametrize("subcommand", [
        "convert", "filter", "subset", "chainer", "nassa",
    ])
    def test_subcommand_help(self, subcommand):
        """Test that help text is printed for each subcommand."""
        sys.argv = ['mddb_workflow', subcommand, '-h']

        # Run the main function, which should print help
        with pytest.raises(SystemExit):
            # Capture stdout during execution
            help_buffer = StringIO()
            with patch('sys.stdout', help_buffer):
                main()

        output = help_buffer.getvalue()
        # Check that help text is printed
        assert f"usage: pytest {subcommand}" in output

    @pytest.mark.parametrize("subcommand", ["error", "run error"])
    def test_errors(self, subcommand):
        """Test that errors are raised for invalid commands."""
        sys.argv = ['mddb_workflow', *subcommand.split()]

        # Run the main function, which should raise SystemExit
        with pytest.raises(SystemExit) as exc_info:
            main()

        # Check that the exit code is non-zero
        assert exc_info.value.code != 0


@pytest.mark.CI
@pytest.mark.release
class TestSubcommands:
    """Test the subcommands of the console interface."""

    def test_run(self, test_data_dir: str):
        """Test that the workflow runs without errors, simulating console execution."""
        working_directory = os.path.join(test_data_dir, 'output/test_run')
        # Copy the inputs from raw_project
        shutil.copytree(
            os.path.join(test_data_dir, 'input/raw_project'),
            working_directory,
            dirs_exist_ok=True
        )

        sys.argv = ['mddb_workflow', 'run',
                    '-dir', working_directory,
                    '-stru', 'raw_structure.pdb',
                    '-traj', 'raw_trajectory.xtc',
                    '-top', 'topology.tpr',
                    '-i', 'setup', 'rmsds',
                    '-filt', 'chain A']
        main()
        os.chdir(test_data_dir)

    def test_inputs(self, test_data_dir: str):
        """Test that the inputs subcommand creates an inputs.yaml file."""
        # We change the directory because the inputs command expects to be run in the output directory
        # TODO add a flag to specify the input directory
        cwd = os.getcwd()
        os.chdir(test_data_dir + '/output')
        sys.argv = ['mddb_workflow', 'inputs', '-ed', 'none']
        main()
        os.chdir(cwd)
        assert os.path.exists(f'{test_data_dir}/output/inputs.yaml')

    def clean_run_assert(self, test_data_dir: str, outputs: list = [], args: list = []):
        """Run a command and assert outputs are created."""
        os.makedirs(f'{test_data_dir}/output/subcommand', exist_ok=True)
        # Ensure the outputs are removed before running the command
        for output in outputs:
            if os.path.exists(output):
                os.remove(output)
        # Set the arguments for the command
        sys.argv = args
        main()
        # Check that the outputs are created
        for output in outputs:
            assert os.path.exists(output)

    def test_convert(self, test_data_dir: str):
        """Test that the convert subcommand can be executed without errors."""
        outputs = [f'{test_data_dir}/output/subcommand/converted.gro',
                   f'{test_data_dir}/output/subcommand/converted.xtc']
        args = ['mddb_workflow', 'convert',
                '-is', f'{test_data_dir}/input/raw_project/raw_structure.pdb',
                '-it', f'{test_data_dir}/input/raw_project/2_frames.mdcrd',
                '-os', outputs[0], '-ot', outputs[1]]
        self.clean_run_assert(test_data_dir, outputs, args)

    def test_filter(self, test_data_dir: str):
        """Test that the filter subcommand can be executed without errors."""
        outputs = [f'{test_data_dir}/output/subcommand/filtered.pdb',
                   f'{test_data_dir}/output/subcommand/filtered.xtc']
        args = ['mddb_workflow', 'filter',
                '-is', f'{test_data_dir}/input/raw_project/raw_structure.pdb',
                '-it', f'{test_data_dir}/input/raw_project/raw_trajectory.xtc',
                '-os', outputs[0], '-ot', outputs[1], '-sel', 'chain A']
        self.clean_run_assert(test_data_dir, outputs, args)

    def test_subset(self, test_data_dir: str):
        """Test that the subset subcommand can be executed without errors."""
        outputs = [f'{test_data_dir}/output/subcommand/subset.xtc']
        args = ['mddb_workflow', 'subset',
                '-is', f'{test_data_dir}/input/raw_project/raw_structure.pdb',
                '-it', f'{test_data_dir}/input/raw_project/raw_trajectory.xtc',
                '-ot', outputs[0], '-start', '1', '-end', '8', '-step', '2']
        self.clean_run_assert(test_data_dir, outputs, args)

    def test_chainer(self, test_data_dir: str):
        """Test that the chainer subcommand can be executed without errors."""
        outputs = [f'{test_data_dir}/output/subcommand/chained.pdb']
        args = ['mddb_workflow', 'chainer',
                '-is', f'{test_data_dir}/input/raw_project/raw_structure.pdb',
                '-os', f'{test_data_dir}/output/subcommand/chained.pdb',
                '-sel', '::A',
                '-let', 'Z',
                '-syn', 'pytraj']
        self.clean_run_assert(test_data_dir, outputs, args)


@pytest.mark.release
class TestNassa:
    """Test the nassa subcommand of the console interface."""
    @pytest.fixture(scope="class")
    def test_accession(self):
        """Override the default accession for this test."""
        return "seq001-1"

    @pytest.mark.CI
    def test_helical(self, project: 'Project'):
        """Test the helical analysis functionality."""
        project.mds[0].run_helical_analysis(project)

    def test_nassa_helical(self, test_data_dir: str):
        """Runs helical analysis using nassa subcommand."""
        os.chdir(os.path.join(test_data_dir, 'output'))
        command = "mwf nassa -hp -stru source_topology.prmtop  -pdirs seq001-1 -m -proj seq001-1 -smp"
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        assert result.returncode == 0, f"Error: {result.stderr}"

    def test_nassa_create_config(self, test_data_dir: str):
        """Creates a nassa.yml configuration file."""
        os.chdir(os.path.join(test_data_dir, 'output'))
        command = f"mwf nassa -w seq001-1/replica_1 -seq {test_data_dir}/input"
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        assert result.returncode == 0, f"Error: {result.stderr}"
        shutil.move("nassa.yml", "seq001-1")

    def test_nassa_all(self, test_data_dir: str):
        """Creates nassa_analysis folder with all results."""
        os.chdir(os.path.join(test_data_dir, 'output'))
        command = "mwf nassa -c seq001-1/nassa.yml -all -o seq001-1/nassa_analysis"
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        assert result.returncode == 0, f"Error: {result.stderr}"
