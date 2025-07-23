import subprocess
import os, sys, shutil, pytest
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
        "convert", "filter", "subset", "chainer", "nassa",
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

    @pytest.mark.parametrize("subcommand", [ "error", "run error"])
    def test_errors(self, subcommand):
        """Test that errors are raised for invalid commands"""
        sys.argv = ['model_workflow', *subcommand.split()]

        # Run the main function, which should raise SystemExit
        with pytest.raises(SystemExit) as exc_info:
            main()
        
        # Check that the exit code is non-zero
        assert exc_info.value.code != 0 

@pytest.mark.CI
@pytest.mark.release
class TestSubcommands:
    """Test the subcommands of the console interface"""

    def test_run(self, test_data_dir: str):
        """Test that the workflow runs without errors, simulating console execution"""

        working_directory = os.path.join(test_data_dir, 'output/test_run')
        # Copy the inputs from raw_project
        shutil.copytree(
            os.path.join(test_data_dir, 'input/raw_project'),
            working_directory,
            dirs_exist_ok=True
        )

        sys.argv = ['model_workflow', 'run',
                    '-dir', working_directory,
                    '-stru', 'raw_structure.pdb',
                    '-traj', 'raw_trajectory.xtc',
                    '-top', 'topology.tpr',
                    '-i', 'setup', 'rmsds',
                    '-filt', 'chain A']
        main()
        os.chdir(test_data_dir)

    def test_inputs(self, test_data_dir: str):
        # We change the directpry because the inputs command expects to be run in the output directory
        # TODO add a flag to specify the input directory
        cwd = os.getcwd()
        os.chdir(test_data_dir + '/output')
        sys.argv = ['model_workflow', 'inputs', '-ed', 'none']
        main()
        os.chdir(cwd)
        assert os.path.exists(f'{test_data_dir}/output/inputs.yaml')

    def test_convert(self, test_data_dir: str):
        """Test that the convert subcommand can be executed without errors"""
        os.makedirs(f'{test_data_dir}/output/subcommand', exist_ok=True)
        sys.argv = ['model_workflow', 'convert', 
                    '-is', f'{test_data_dir}/input/raw_project/raw_structure.pdb', #mdcrd xtc
                    '-os', f'{test_data_dir}/output/subcommand/converted.gro',
                    '-it', f'{test_data_dir}/input/raw_project/2_frames.mdcrd',
                    '-ot', f'{test_data_dir}/output/subcommand/converted.xtc']
        main()
        assert os.path.exists(f'{test_data_dir}/output/subcommand/converted.gro')
        assert os.path.exists(f'{test_data_dir}/output/subcommand/converted.xtc')

    def test_filter(self, test_data_dir: str):
        """Test that the filter subcommand can be executed without errors"""
        os.makedirs(f'{test_data_dir}/output/subcommand', exist_ok=True)
        sys.argv = ['model_workflow', 'filter', 
                    '-is', f'{test_data_dir}/input/raw_project/raw_structure.pdb',
                    '-os', f'{test_data_dir}/output/subcommand/filtered.pdb',
                    '-it', f'{test_data_dir}/input/raw_project/raw_trajectory.xtc',
                    '-ot', f'{test_data_dir}/output/subcommand/filtered.xtc',
                    '-sel', 'chain A']
        main()
        assert os.path.exists(f'{test_data_dir}/output/subcommand/filtered.pdb')
        assert os.path.exists(f'{test_data_dir}/output/subcommand/filtered.xtc')

    def test_subset(self, test_data_dir: str):
        os.makedirs(f'{test_data_dir}/output/subcommand', exist_ok=True)
        sys.argv = ['model_workflow', 'subset', 
                '-is', f'{test_data_dir}/input/raw_project/raw_structure.pdb',
                '-it', f'{test_data_dir}/input/raw_project/raw_trajectory.xtc',
                '-ot', f'{test_data_dir}/output/subcommand/subset.xtc',
                '-start', '1',
                '-end', '8',
                '-step', '2']
        main()
        assert os.path.exists(f'{test_data_dir}/output/subcommand/subset.xtc')

    def test_chainer(self, test_data_dir: str):
        os.makedirs(f'{test_data_dir}/output/subcommand', exist_ok=True)
        sys.argv = ['model_workflow', 'chainer', 
                '-is', f'{test_data_dir}/input/raw_project/raw_structure.pdb',
                '-os', f'{test_data_dir}/output/subcommand/chained.pdb',
                '-sel', '::A',
                '-let', 'Z',
                '-syn', 'pytraj']
        main()
        assert os.path.exists(f'{test_data_dir}/output/subcommand/chained.pdb')

@pytest.mark.release
class TestNassa:
    @pytest.fixture(scope="class")
    def test_accession(self):
        """Override the default accession for this test"""
        return "seq001-1"
    
    def test_helical(self, project: 'Project'):
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
        