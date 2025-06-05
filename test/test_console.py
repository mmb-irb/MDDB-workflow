import sys
import pytest
from io import StringIO
from unittest.mock import patch
from model_workflow.console import parser, main

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

class TestConsoleIntegration:
    """Integration tests for console functionality"""
    
    @pytest.mark.parametrize("subcommand", [
        "convert", 
        "filter",
        "subset",
        "chainer",
        "nassa"
    ])
    def test_subcommand_help(self, capture_stdout, subcommand):
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