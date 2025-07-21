import os
import sys
import shutil
import pytest
import subprocess
from io import StringIO
from unittest.mock import patch
from model_workflow.console import parser, main
from model_workflow.utils.type_hints import *


@pytest.mark.CI
@pytest.mark.release
class TestConsoleArgumentParsing:
    """Test the console argument parsing functionality"""
    
    def test_parser_exists(self):
        """Verify that the parser is properly defined"""
        assert parser is not None
        assert hasattr(parser, 'parse_args')

    def test_no_arguments(self):
        """Test behavior when no arguments are provided"""
        # Set up empty arguments
        sys.argv = ['model_workflow']
        
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
    """Integration tests for console functionality"""
    
    @pytest.mark.parametrize("subcommand", [
        "convert", 
        "filter",
        "subset",
        "chainer",
        "nassa"
    ])
    def test_subcommand_help(self, subcommand):
        """Test that help text is printed for each subcommand"""
        sys.argv = ['model_workflow', subcommand, '-h']
        
        # Run the main function, which should print help
        with pytest.raises(SystemExit):
            # Capture stdout during execution
            help_buffer = StringIO()
            with patch('sys.stdout', help_buffer):
                main()

        output = help_buffer.getvalue()
        # Check that help text is printed
        assert f"usage: pytest {subcommand}" in output

@pytest.mark.CI
@pytest.mark.release
class TestSubcommands:
    """Test the subcommands of the console interface"""

    def test_subcommand_run(self, test_data_dir: str):
        """Test that the workflow runs without errors, simulating console execution"""

        working_directory = os.path.join(test_data_dir, 'output/test_workflow')
        # Ensure the working directory doesn't exist and it is empty
        if os.path.exists(working_directory):
            # Remove the directory if it exists
            shutil.rmtree(working_directory)
        # Copy the inputs from raw_project
        shutil.copytree(
            os.path.join(test_data_dir, 'input/raw_project'),
            working_directory,
            dirs_exist_ok=True
        )

        sys.argv = ['model_workflow', 'run',
                    '-dir', working_directory,
                    '-i', 'setup',
                    '-filt', 'chain A']
        main()
        os.chdir(test_data_dir)

    @pytest.mark.skip(reason="Test not implemented yet")
    def test_subcommand_inputs(self, test_data_dir: str):
        pass
    
    def test_subcommand_convert(self, test_data_dir: str):
        """Test that the convert subcommand can be executed without errors"""
        os.makedirs(f'{test_data_dir}/output/subcommand', exist_ok=True)
        sys.argv = ['model_workflow', 'convert', 
                    '-is', f'{test_data_dir}/input/raw_project/5ggr-rs1-310k-0ns.pdb',
                    '-os', f'{test_data_dir}/output/subcommand/converted.pdb',
                    '-it', f'{test_data_dir}/input/raw_project/replica_1/trajectory.xtc',
                    '-ot', f'{test_data_dir}/output/subcommand/converted.dcd']
        main()
        assert os.path.exists(f'{test_data_dir}/output/subcommand/converted.pdb')
        assert os.path.exists(f'{test_data_dir}/output/subcommand/converted.dcd')

    def test_subcommand_filter(self, test_data_dir: str):
        """Test that the filter subcommand can be executed without errors"""
        os.makedirs(f'{test_data_dir}/output/subcommand', exist_ok=True)
        sys.argv = ['model_workflow', 'filter', 
                    '-is', f'{test_data_dir}/input/raw_project/5ggr-rs1-310k-0ns.pdb',
                    '-os', f'{test_data_dir}/output/subcommand/filtered.pdb',
                    '-it', f'{test_data_dir}/input/raw_project/replica_1/trajectory.xtc',
                    '-ot', f'{test_data_dir}/output/subcommand/filtered.xtc',
                    '-sel', 'chain A']
        main()
        assert os.path.exists(f'{test_data_dir}/output/subcommand/filtered.pdb')
        assert os.path.exists(f'{test_data_dir}/output/subcommand/filtered.xtc')

    @pytest.mark.skip(reason="Test not implemented yet")
    def test_subcommand_subset(self, test_data_dir: str):
        pass

    @pytest.mark.skip(reason="Test not implemented yet")
    def test_subcommand_chainer(self, test_data_dir: str):
        pass

@pytest.mark.release
class TestNassa:
    @pytest.fixture(scope="class")
    def test_accession(self):
        """Override the default accession for this test"""
        return "seq001-1"
    
    def test_helical(self, test_data_dir, project: 'Project'):
        project.mds[0].run_helical_analysis(project)

    def test_nassa_helical(self, test_data_dir: str):
        os.chdir(os.path.join(test_data_dir, 'output'))
        command = "mwf nassa -hp -stru *.prmtop -traj *.xtc -top *.prmtop  -pdirs seq001-1 -mdir replica_1 -m"
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        assert result.returncode == 0

    def test_nassa_create_config(self, test_data_dir: str):
        os.chdir(os.path.join(test_data_dir, 'output'))
        command = f"mwf nassa -w seq001-1/replica_1 -seq {test_data_dir}/input"
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        assert result.returncode == 0

    def test_nassa_all(self, test_data_dir: str):
        os.chdir(os.path.join(test_data_dir, 'output'))
        command = "mwf nassa -c nassa.yml -all"
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        assert result.returncode == 0