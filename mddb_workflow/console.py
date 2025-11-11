from os.path import exists
from shutil import copyfile
from subprocess import call
from argparse import ArgumentParser, RawTextHelpFormatter, Action, _SubParsersAction
from textwrap import wrap, dedent
import re

from mddb_workflow.mwf import workflow, Project, requestables, DEPENDENCY_FLAGS
from mddb_workflow.utils.structures import Structure
from mddb_workflow.utils.file import File
from mddb_workflow.utils.filters import filter_atoms
from mddb_workflow.utils.subsets import get_trajectory_subset
from mddb_workflow.utils.constants import *
from mddb_workflow.utils.auxiliar import InputError
from mddb_workflow.utils.nassa_file import generate_nassa_config
from mddb_workflow.tools.conversions import convert
from mddb_workflow.analyses.nassa import workflow_nassa
from mddb_workflow.core.dataset import Dataset

expected_project_args = set(Project.__init__.__code__.co_varnames)

test_docs_url = 'https://mddb-workflow.readthedocs.io/en/latest/usage.html#tests-and-other-checking-processes'
task_docs_url = 'https://mddb-workflow.readthedocs.io/en/latest/tasks.html'


class CustomHelpFormatter(RawTextHelpFormatter):
    """Custom formatter for argparse help text with better organization and spacing."""

    def __init__(self, prog, indent_increment=2, max_help_position=6, width=None):
        super().__init__(prog, indent_increment, max_help_position, width)

    def _get_help_string(self, action):
        import argparse
        help = action.help
        if '%(default)' not in action.help:
            if action.default is not argparse.SUPPRESS and \
                action.default and action.default != '.':
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    help += '\nDefault: %(default)s'
        return help

    def _split_lines(self, text, width):
        lines = []
        for line in text.splitlines():
            if line.strip() != '':
                if line.startswith('https'):
                    lines.append(line)
                else:
                    lines.extend(wrap(line, width, break_long_words=False, replace_whitespace=False))
        return lines

    def _format_usage(self, usage, actions, groups, prefix):
        essential_usage = super()._format_usage(usage, actions, groups, prefix)
        # Only for mwf run
        if 'run' in self._prog:
            # Combine the aguments for -i, -e, -ow
            lines = essential_usage.split('\n')
            filtered_lines = []
            for line in lines:
                if line.strip().startswith("[-i "):
                    line = line.replace("[-i", "[-i/-e/-ow")
                    filtered_lines.append(line)
                elif line.strip().startswith("[-e") or line.strip().startswith("[-ow"):
                    continue
                else:
                    filtered_lines.append(line)
            essential_usage = '\n'.join(filtered_lines)
        return essential_usage

    def _format_action_invocation(self, action):
        """Format the display of options with choices more cleanly."""
        if not action.option_strings:
            # This is a positional argument
            return super()._format_action_invocation(action)

        # For options with choices, format them nicely
        opts = ', '.join(action.option_strings)

        # Special case for include, exclude, and overwrite
        if action.dest in ['include', 'exclude', 'overwrite']:
            opts = ', '.join(action.option_strings)
            metavar = 'TASKS'
            return f"{opts} {metavar}"
        if action.nargs == 0:
            # Boolean flag
            return opts
        else:
            # Format with metavar or choices
            metavar = self._format_args(action, action.dest.upper())
            if action.choices:
                choice_str = '{' + ','.join(str(c) for c in action.choices) + '}'
                # if action.nargs is not None and action.nargs != 1:
                #     choice_str += ' ...'
                return f"{opts} [{choice_str}]"
            else:
                return f"{opts} {metavar}"


class CustomArgumentParser(ArgumentParser):
    """Extends the ArgumentParser to handle subcommands and errors more gracefully."""

    def error(self, message):
        # Check for subcommand in sys.argv
        import sys
        # Extract subcommand from command line if it exists
        if hasattr(self, '_subparsers') and self._subparsers is not None:
            subcommands = [choice for action in self._subparsers._actions
                          if isinstance(action, _SubParsersAction)
                          for choice in action.choices]
            if len(sys.argv) > 1 and sys.argv[1] in subcommands:
                self.subcommand = sys.argv[1]

        # Now continue with your existing logic
        if hasattr(self, 'subcommand') and self.subcommand:
            self._print_message(f"{self.prog} {self.subcommand}: error: {message}\n", sys.stderr)
            # Show help for the specific subparser
            for action in self._subparsers._actions:
                if isinstance(action, _SubParsersAction):
                    for choice, subparser in action.choices.items():
                        if choice == self.subcommand:
                            subparser.print_usage()
                            break
        else:
            # Default error behavior for main parser
            self.print_usage(sys.stderr)
            self._print_message(f"{self.prog}: error: {message}\n", sys.stderr)
        sys.exit(2)


def parse_docstring_for_help(docstring):
    """Parse a docstring to extract help for arguments."""
    if not docstring:
        return {}

    docstring = dedent(docstring)
    arg_section_match = re.search(r'Args:\n(.*?)(?=\n\n|\Z)', docstring, re.S)
    if not arg_section_match:
        return {}

    arg_section = arg_section_match.group(1)
    help_dict = {}
    # Regex to find argument name and its full description, including newlines
    arg_blocks = re.findall(r'^\s*([a-zA-Z0-9_]+)\s*\(.*?\):\s*(.*?)(?=\n\s*[a-zA-Z0-9_]+\s*\(|\Z)', arg_section, re.S | re.M)

    for arg_name, help_text in arg_blocks:
        # Clean up the help text: remove leading/trailing whitespace from each line and join
        lines = []
        for line in help_text.strip().split('\n'):
            # For lines that are part of a list, just rstrip to keep indentation
            if line.lstrip().startswith('-'):
                lines.append('  ' + line.strip())
            else:
                lines.append(line.strip())
        help_dict[arg_name] = '\n'.join(lines)

    return help_dict


def pretty_list(availables: list[str]) -> str:
    """Pretty print a list of available checkings / failures."""
    final_line = '\n'
    for available in availables:
        nice_name = NICE_NAMES.get(available, None)
        if not nice_name:
            raise Exception('Flag "' + available + '" has not a defined nice name')
        final_line += '\n  - ' + available + ': ' + nice_name
    final_line += f'\nTo know more about each test please visit:\n{test_docs_url}'
    return final_line


class custom (Action):
    """Custom argparse action to handle the following 2 arguments.

    This is done becuase it is not possible to combine nargs='*' with const
    https://stackoverflow.com/questions/72803090/argparse-how-to-create-the-equivalent-of-const-with-nargs
    """

    def __call__(self, parser, namespace, values, option_string=None):
        """Handle optional argument values with smart defaults.

        If argument is not passed -> default
        If argument is passed empty -> const
        If argument is passed with values -> values
        """
        if values:
            setattr(namespace, self.dest, values)
        else:
            setattr(namespace, self.dest, self.const)


def main():
    """Parse the arguments and calls the workflow accordingly."""
    # Parse input arguments from the console
    # The vars function converts the args object to a dictionary
    args = parser.parse_args()
    if hasattr(args, 'subcommand') and args.subcommand:
        parser.subcommand = args.subcommand
    # Apply common arguments as necessary
    if hasattr(args, 'no_symlinks') and args.no_symlinks:
        GLOBALS['no_symlinks'] = True
    if hasattr(args, 'no_colors') and args.no_colors:
        GLOBALS['no_colors'] = True
    # Find which subcommand was called
    subcommand = args.subcommand
    # If there is not subcommand then print help
    if not subcommand:
        parser.print_help()
    # If user wants to run the workflow
    elif subcommand == "run":
        # Ger all parsed arguments
        dict_args = vars(args)
        # Remove arguments not related to this subcommand
        del dict_args['subcommand']
        # Remove common arguments from the dict as well
        common_args = [action.dest for action in common_parser._actions]
        for arg in common_args:
            del dict_args[arg]
        # Find out which arguments are for the Project class and which ones are for the workflow
        project_args = {}
        workflow_args = {}
        for k, v in dict_args.items():
            if k in expected_project_args:
                project_args[k] = v
            else:
                workflow_args[k] = v
        # Call the actual main function
        workflow(project_parameters=project_args, **workflow_args)
    # If user wants to setup the inputs
    elif subcommand == "inputs":
        # Make a copy of the template in the local directory if there is not an inputs file yet
        if exists(DEFAULT_INPUTS_FILENAME):
            print(f"File {DEFAULT_INPUTS_FILENAME} already exists")
        else:
            copyfile(INPUTS_TEMPLATE_FILEPATH, DEFAULT_INPUTS_FILENAME)
            print(f"File {DEFAULT_INPUTS_FILENAME} has been generated")
        # Set the editor to be used to modify the inputs file
        editor_command = args.editor
        if editor_command:
            if editor_command == 'none': return
            return call([editor_command, DEFAULT_INPUTS_FILENAME])
        # If no editor argument is passed then ask the user for one
        print("Choose your preferred editor:")
        available_editors = list(AVAILABLE_TEXT_EDITORS.keys())
        for i, editor_name in enumerate(available_editors, 1):
            print(f"{i}. {editor_name}")
        print("*. exit")
        try:
            choice = int(input("Number: ").strip())
            if not (1 <= choice <= len(available_editors)): raise ValueError
            editor_name = available_editors[choice - 1]
            editor_command = AVAILABLE_TEXT_EDITORS[editor_name]
            # Open a text editor for the user
            print(f"{editor_name} was selected")
            call([editor_command, DEFAULT_INPUTS_FILENAME])
        except ValueError:
            print("No editor was selected")

    # In case the convert tool was called
    elif subcommand == 'convert':
        # If no input arguments are passed print help
        if args.input_structure is None and args.input_trajectories is None:
            convert_parser.print_help()
            return
        if args.input_trajectories is None:
            args.input_trajectories = []
        # Run the convert command
        convert(
            input_structure_filepath=args.input_structure,
            output_structure_filepath=args.output_structure,
            input_trajectory_filepaths=args.input_trajectories,
            output_trajectory_filepath=args.output_trajectory,
        )

    # In case the filter tool was called
    elif subcommand == 'filter':
        # Run the convert command
        filter_atoms(
            input_structure_file=File(args.input_structure),
            output_structure_file=File(args.output_structure),
            input_trajectory_file=File(args.input_trajectory),
            output_trajectory_file=File(args.output_trajectory),
            selection_string=args.selection_string,
            selection_syntax=args.selection_syntax
        )
        print('There you have it :)')

    # In case the subset tool was called
    elif subcommand == 'subset':
        output_trajectory = args.output_trajectory if args.output_trajectory else args.input_trajectory
        get_trajectory_subset(
            input_structure_file=File(args.input_structure),
            input_trajectory_file=File(args.input_trajectory),
            output_trajectory_file=File(output_trajectory),
            start=args.start,
            end=args.end,
            step=args.step,
            skip=args.skip,
            frames=args.frames
        )
        print('All done :)')

    # In case the chainer tool was called
    elif subcommand == 'chainer':
        # Parse the structure
        structure = Structure.from_pdb_file(args.input_structure)
        # Select atom accoridng to inputs
        selection = structure.select(args.selection_string, args.selection_syntax) if args.selection_string else structure.select_all()
        if not selection: raise InputError(f'Empty selection {selection}')
        # Run the chainer logic
        structure.chainer(selection, args.letter, args.whole_fragments)
        # Generate the output file from the modified structure
        structure.generate_pdb_file(args.output_structure)
        print(f'Changes written to {args.output_structure}')
    elif subcommand == 'dataset':
        if not hasattr(args, 'dataset_subcommand') or not args.dataset_subcommand:
            dataset_parser.print_help()
            return

        dataset = Dataset(dataset_yaml_path=args.dataset_yaml)

        if args.dataset_subcommand == 'run':
            dataset.launch_workflow(
                include_groups=args.include_groups,
                exclude_groups=args.exclude_groups,
                slurm=args.slurm,
                job_template=args.job_template
            )
        elif args.dataset_subcommand == 'status':
            dataset.show_groups(cmd=True)
    # If user wants to run the NASSA analysis
    elif subcommand == "nassa":
        # If no input arguments are passed print help
        if args.config is None and args.make_config is False:
            nassa_parser.print_help()
            print('Please provide a configuration file or make one with the -m flag')
            return
        # If the user wants to make a configuration file
        if args.make_config:
            # print('args.make_config: ', args.make_config)
            if args.make_config is True or args.make_config == []:
                # Make a copy of the template in the local directory if there is not an inputs file yet
                if not exists(DEFAULT_NASSA_CONFIG_FILENAME):
                    copyfile(NASSA_TEMPLATE_FILEPATH, DEFAULT_NASSA_CONFIG_FILENAME)
                # Open a text editor for the user
                call(["vim", DEFAULT_NASSA_CONFIG_FILENAME])
                print('Configuration file created as nassa.json\nNow you can run the analysis with the -c flag.')
                return
            # If the user provides a path to the files
            else:
                generate_nassa_config(
                    args.make_config,
                    args.seq_path,
                    args.output,
                    args.unit_len,
                    args.n_sequences
                    )
            print('Configuration file created as nassa.json\nNow you can run the analysis with the -c flag.')
            return
        # If the user wants to run the analysis. With the config file an analysis name must be provided, or the all flag must be set
        if args.config and args.analysis_names is None and args.all is False:
            nassa_parser.print_help()
            print('Please provide an analysis name to run:', ', '.join(NASSA_ANALYSES_LIST))
            return
        # If the user wants to run the helical parameters analysis we must check if the necessary files are provided (structure, topology and trajectory)
        if args.helical_parameters:
            # Also, it is necesary to provide the project directories. Each of the project directories must contain an independent MD
            if args.proj_directories is None:
                nassa_parser.print_help()
                print('Please provide a project directory to run the helical parameters analysis with the -pdirs flag')
                return
            if args.input_structure_filepath is None:
                raise InputError('Please provide a structure file to run the helical parameters analysis with the -stru flag')
            elif args.input_trajectory_filepath is None:
                raise InputError('Please provide a trajectory file to run the helical parameters analysis with the -traj flag')
            elif args.input_topology_filepath is None:
                raise InputError('Please provide a topology file to run the helical parameters analysis with the -top flag')
            # If the all flag is set, the user must provide the path to the sequences because it is necessary to create the nassa.yml and run the NASSA analysis
            if args.all:
                if not args.seq_path:
                    raise InputError('Please, if all option is selected provide the path to the sequences (--seq_path)')
            # If all the flags are correctly set, we can run the analysis
            workflow_nassa(
                config_file_path=None,  # The configuration file is not needed in this case because we are going to run the helical parameters analysis so it will be created then
                analysis_names=args.analysis_names,
                overwrite=args.overwrite,
                overwrite_nassa=args.overwrite_nassa,
                helical_par=args.helical_parameters,
                proj_dirs=args.proj_directories,
                input_structure_file=args.input_structure_filepath,
                input_trajectory_file=args.input_trajectory_filepath,
                input_top_file=args.input_topology_filepath,
                all=args.all,
                unit_len=args.unit_len,
                n_sequences=args.n_sequences,
                seq_path=args.seq_path,
                md_directories=args.md_directories,
                trust=args.trust,
                mercy=args.mercy
            )
        # If the user wants to run the NASSA analysis with the config file already created and the analysis name provided
        else:
            dict_args = vars(args)
            del dict_args['subcommand']  # preguntar Dani ¿?
            # Call the actual main function
            workflow_nassa(
                    config_file_path=args.config,
                    analysis_names=args.analysis_names,
                    make_config=args.make_config,
                    output=args.output,
                    working_directory=args.working_directory,
                    overwrite=args.overwrite,
                    overwrite_nassa=args.overwrite_nassa,
                    n_sequences=args.n_sequences,
                    unit_len=args.unit_len,
                    all=args.all,
                    md_directories=args.md_directories,
                    trust=args.trust,
                    mercy=args.mercy
            )


# Define a common parser running in top of all others
# This arguments declared here are available among all subparsers
common_parser = ArgumentParser(add_help=False)

# If this argument is passed then no symlinks will be used anywhere
# Files will be copied instead thus taking more time and disk
# However symlinks are not always allowed in all file systems so this is sometimes necessary
common_parser.add_argument("-ns", "--no_symlinks", default=False, action='store_true', help="Do not use symlinks internally")
common_parser.add_argument("-nc", "--no_colors", default=False, action='store_true', help="Do not use colors for logging")

# Define console arguments to call the workflow
parser = CustomArgumentParser(description="MDDB Workflow")
subparsers = parser.add_subparsers(help='Name of the subcommand to be used', dest="subcommand")

project_init_help = parse_docstring_for_help(Project.__init__.__doc__)
workflow_help = parse_docstring_for_help(workflow.__doc__)

# Set the run subcommand
run_parser = subparsers.add_parser("run",
    help="Run the workflow",
    formatter_class=CustomHelpFormatter,
    parents=[common_parser]
)

# Set optional arguments
run_parser_input_group = run_parser.add_argument_group('INPUT OPTIONS')
run_parser_input_args = [
    # There is no default since many formats may be possible
    (['-top', '--input_topology_filepath'], {'default': None, 'help': project_init_help['input_topology_filepath']}),
    (['-stru', '--input_structure_filepath'], {'default': None, 'help': project_init_help['input_structure_filepath']}),
    (['-traj', '--input_trajectory_filepaths'], {'default': None, 'nargs': '*', 'help': project_init_help['input_trajectory_filepaths']}),
    (['-dir', '--working_directory'], {'default': '.', 'help': "Directory where the whole workflow is run."}),
    (['-mdir', '--md_directories'], {'default': None, 'nargs': '*', 'help': project_init_help['md_directories']}),
    (['-md', '--md_config'], {'action': 'append', 'default': None, 'nargs': '*', 'help': project_init_help['md_config']}),
    (['-proj', '--accession'], {'default': None, 'help': project_init_help['accession']}),
    (['-url', '--database_url'], {'default': DEFAULT_API_URL, 'help': project_init_help['database_url']}),
    (['-inp', '--inputs_filepath'], {'default': None, 'help': "Path to inputs file"}),
    (['-fin', '--forced_inputs'], {'action': 'append', 'nargs': '*', 'default': None, 'help': project_init_help['forced_inputs']}),
    (['-pop', '--populations_filepath'], {'default': DEFAULT_POPULATIONS_FILENAME, 'help': project_init_help['populations_filepath']}),
    (['-tpro', '--transitions_filepath'], {'default': DEFAULT_TRANSITIONS_FILENAME, 'help': project_init_help['transitions_filepath']}),
    (['-ad', '--aiida_data_filepath'], {'default': None, 'help': project_init_help['aiida_data_filepath']}),
]
for flags, kwargs in run_parser_input_args:
    run_parser_input_group.add_argument(*flags, **kwargs)

# Set a group for the workflow control options
run_parser_workflow_group = run_parser.add_argument_group('WORKFLOW CONTROL OPTIONS')
run_parser_workflow_args = [
    (['-img', '--image'], {'action': 'store_true', 'help': project_init_help['image']}),
    (['-fit', '--fit'], {'action': 'store_true', 'help': project_init_help['fit']}),
    (['-trans', '--translation'], {'nargs': '*', 'default': [0, 0, 0], 'help': project_init_help['translation']}),
    (['-d', '--download'], {'action': 'store_true', 'help': workflow_help['download']}),
    (['-s', '--setup'], {'action': 'store_true', 'help': workflow_help['setup']}),
    (['-smp', '--sample_trajectory'], {'type': int, 'nargs': '?', 'default': None, 'const': 10, 'metavar': 'N_FRAMES', 'help': project_init_help['sample_trajectory']}),
    (['-rcut', '--rmsd_cutoff'], {'type': float, 'default': DEFAULT_RMSD_CUTOFF, 'help': project_init_help['rmsd_cutoff']}),
    (['-icut', '--interaction_cutoff'], {'type': float, 'default': DEFAULT_INTERACTION_CUTOFF, 'help': project_init_help['interaction_cutoff']}),
    (['-iauto', '--interactions_auto'], {'type': str, 'nargs': '?', 'const': True, 'help': project_init_help['interactions_auto']}),
    (['-gb', '--guess_bonds'], {'action': 'store_true', 'help': project_init_help['guess_bonds']}),
    (['-ib', '--ignore_bonds'], {'action': 'store_true', 'help': project_init_help['ignore_bonds']}),
]
for flags, kwargs in run_parser_workflow_args:
    run_parser_workflow_group.add_argument(*flags, **kwargs)

# Set a group for the selection options
run_parser_selection_group = run_parser.add_argument_group('SELECTION OPTIONS')
run_parser_selection_args = [
    (['-filt', '--filter_selection'], {'nargs': '?', 'default': False, 'const': True, 'help': project_init_help['filter_selection']}),
    (['-pbc', '--pbc_selection'], {'default': None, 'help': project_init_help['pbc_selection']}),
    (['-cg', '--cg_selection'], {'default': None, 'help': project_init_help['cg_selection']}),
    (['-pcafit', '--pca_fit_selection'], {'default': PROTEIN_AND_NUCLEIC_BACKBONE, 'help': project_init_help['pca_fit_selection']}),
    (['-pcana', '--pca_analysis_selection'], {'default': PROTEIN_AND_NUCLEIC_BACKBONE, 'help': project_init_help['pca_analysis_selection']}),
]
for flags, kwargs in run_parser_selection_args:
    run_parser_selection_group.add_argument(*flags, **kwargs)

# Set a group with all input checking options
run_parser_checks_group = run_parser.add_argument_group('INPUT CHECKS OPTIONS', description=f"For more information about each check please visit:\n{test_docs_url}")
run_parser_checks_args = [
    (['-t', '--trust'], {'default': [], 'nargs': '*', 'action': custom, 'const': AVAILABLE_CHECKINGS, 'choices': AVAILABLE_CHECKINGS,
      'help': ("If passed, do not run the specified checking. Note that all checkings are skipped if passed alone. "
               "Available checkings:" + pretty_list(AVAILABLE_CHECKINGS))}),
    (['-m', '--mercy'], {'default': [], 'nargs': '*', 'action': custom, 'const': AVAILABLE_FAILURES, 'choices': AVAILABLE_FAILURES,
     'help': ("If passed, do not kill the process when any of the specfied checkings fail and proceed with the workflow. "
              "Note that all checkings are allowed to fail if the argument is passed alone. "
              "Available checkings:" + pretty_list(AVAILABLE_FAILURES))}),
    (['-f', '--faith'], {'default': False, 'action': 'store_true',
     'help': ("Use this flag to force-skip all data processing thus asuming inputs are already processed.\n"
              "WARNING: Do not use this flag if you don't know what you are doing.\n"
              "This may lead to several silent errors.")}),
]
for flags, kwargs in run_parser_checks_args:
    run_parser_checks_group.add_argument(*flags, **kwargs)

# Set a list with the alias of all requestable dependencies
choices = sorted(list(requestables.keys()) + list(DEPENDENCY_FLAGS.keys()))
task_groups = [
  "download: Check/download input files (already ran with analyses)",
  "setup: Process and test input files (already ran with analyses)",
  "meta: Run project and MD metadata analyses",
  "network: Run dependencies which require internet connection",
  "minimal: Run dependencies required by the web client to work",
  "interdeps: Run interactions and all its dependent analyses",
  "membs: Run all membrane-related analyses",
]
assert len(DEPENDENCY_FLAGS.keys()) == len(task_groups), "The number of dependency flags and task groups must be the same"

run_parser_analysis_group = run_parser.add_argument_group('TASKS OPTIONS',
    description=f"Available tasks: {choices}\nFor more information about each task, please visit:\n{task_docs_url}")
run_parser_analysis_args = [
    (['-i', '--include'], {'nargs': '*', 'choices': choices,
      'help': ("Set the unique analyses or tools to be run. All other steps will be skipped.\n"
               "There are also some additional flags to define a preconfigured group of dependencies:"
               + '\n  - ' + '\n  - '.join(task_groups))}),
    (['-e', '--exclude'], {'nargs': '*', 'choices': choices, 'help': workflow_help['exclude']}),
    (['-ow', '--overwrite'], {'type': str, 'nargs': '*', 'default': [], 'action': custom, 'const': True, 'choices': choices, 'help': workflow_help['overwrite']}),
]
for flags, kwargs in run_parser_analysis_args:
    run_parser_analysis_group.add_argument(*flags, **kwargs)


# Add a new to command to aid in the inputs file setup
inputs_parser = subparsers.add_parser("inputs",
    help="Set the inputs file",
    formatter_class=CustomHelpFormatter,
    parents=[common_parser]
)
# Chose the editor in advance
inputs_parser.add_argument(
    "-ed", "--editor",
    choices=[*AVAILABLE_TEXT_EDITORS.values(), 'none'],
    help="Set the text editor to modify the inputs file")


# The convert command
convert_parser = subparsers.add_parser("convert",
    help="Convert a structure and/or several trajectories to other formats\n" +
        "If several input trajectories are passed they will be merged previously",
    formatter_class=CustomHelpFormatter,
    parents=[common_parser])
convert_parser_args = [
    (['-is', '--input_structure'], {'help': "Path to input structure file"}),
    (['-os', '--output_structure'], {'help': "Path to output structure file"}),
    (['-it', '--input_trajectories'], {'nargs': '*', 'help': "Path to input trajectory file or same format files."}),
    (['-ot', '--output_trajectory'], {'help': "Path to output trajectory file"}),
]
for flags, kwargs in convert_parser_args:
    convert_parser.add_argument(*flags, **kwargs)


# The filter command
filter_parser = subparsers.add_parser("filter",
    help="Filter atoms in a structure and/or a trajectory\n",
    formatter_class=CustomHelpFormatter,
    parents=[common_parser])
filter_parser_args = [
    (["-is", "--input_structure"], {'required': True, 'help': "Path to input structure file"}),
    (["-os", "--output_structure"], {'help': "Path to output structure file"}),
    (["-it", "--input_trajectory"], {'help': "Path to input trajectory file"}),
    (["-ot", "--output_trajectory"], {'help': "Path to output trajectory file"}),
    (["-sel", "--selection_string"], {'help': "Atom selection"}),
    (["-syn", "--selection_syntax"], {'default': 'vmd', 'help': "Atom selection syntax (vmd by default)"}),
]
for flags, kwargs in filter_parser_args:
    filter_parser.add_argument(*flags, **kwargs)


# The subset command
subset_parser = subparsers.add_parser("subset",
    help="Get a subset of frames from the current trajectory",
    formatter_class=CustomHelpFormatter,
    parents=[common_parser])
subset_parser_args = [
    (["-is", "--input_structure"], {'required': True, 'help': "Path to input structure file"}),
    (["-it", "--input_trajectory"], {'help': "Path to input trajectory file"}),
    (["-ot", "--output_trajectory"], {'help': "Path to output trajectory file"}),
    (["-start", "--start"], {'type': int, 'default': 0, 'help': "Start frame (0-based)"}),
    (["-end", "--end"], {'type': int, 'default': None, 'help': "End frame (0-based)"}),
    (["-step", "--step"], {'type': int, 'default': 1, 'help': "Frame step"}),
    (["-skip", "--skip"], {'nargs': '*', 'type': int, 'default': [], 'help': "Frames to be skipped (0-based)"}),
    (["-fr", "--frames"], {'nargs': '*', 'type': int, 'default': [], 'help': "Frames to be returned (0-based). Input frame order is ignored as original frame order is conserved."}),
]
for flags, kwargs in subset_parser_args:
    subset_parser.add_argument(*flags, **kwargs)


# The chainer command
chainer_parser = subparsers.add_parser("chainer",
    help="Edit structure (pdb) chains",
    formatter_class=CustomHelpFormatter,
    parents=[common_parser])
chainer_parser_args = [
    (["-is", "--input_structure"], {'required': True, 'help': "Path to input structure file"}),
    (["-os", "--output_structure"], {'default': 'chained.pdb', 'help': "Path to output structure file"}),
    (["-sel", "--selection_string"], {'help': "Atom selection (the whole structure by default)"}),
    (["-syn", "--selection_syntax"], {'default': 'vmd', 'choices': Structure.SUPPORTED_SELECTION_SYNTAXES, 'help': "Atom selection syntax (VMD syntax by default)"}),
    (["-let", "--letter"], {'help': "New chain letter (one letter per fragment by default)"}),
    (["-whfr", "--whole_fragments"], {'type': bool, 'default': False, 'help': "Consider fragments beyond the atom selection. Otherwise a fragment could end up having multiple chains."}),
]
for flags, kwargs in chainer_parser_args:
    chainer_parser.add_argument(*flags, **kwargs)


# The NASSA commands
nassa_parser = subparsers.add_parser("nassa", formatter_class=CustomHelpFormatter,
    help="Run and set the configuration of the NASSA analysis",
    parents=[common_parser])
nassa_parser.add_argument(
    "-c", "--config",
    help="Configuration file for the NASSA analysis")
nassa_parser.add_argument(
    "-n", "--analysis_names",
    nargs='*',
    default=None,
    help="Name of the analysis to be run. It can be: " + ', '.join(NASSA_ANALYSES_LIST))
nassa_parser.add_argument(
    "-w", "--make_config",
    # type=str,
    nargs='*',
    default=None,
    # const=True,
    # action=custom,
    help="Make a configuration file for the NASSA analysis: makecfg.\nThe base path could be given as an argument. If not, an example of configuration file is created.")
nassa_parser.add_argument(
    "-seq", "--seq_path",
    type=str,
    const=False,
    action=custom,
    help="Set the base path of the sequences. If not given, the sequences are searched in the current directory.")
nassa_parser.add_argument(
    "-o", "--output",
    help="Output path for the NASSA analysis")
nassa_parser.add_argument(
    "-dir", "--working_directory",
    default='.',
    help="Directory where the whole workflow is run. Current directory by default.")
nassa_parser.add_argument(
    "-ow", "--overwrite",
    type=str,
    nargs='*',
    default=[],
    action=custom,
    const=True,
    help="Set the output files to be overwritten thus re-runing its corresponding analysis or tool")
nassa_parser.add_argument(
    "-own", "--overwrite_nassa",
    type=str,
    nargs='*',
    default=[],
    action=custom,
    const=True,
    help="Set the output files to be overwritten thus re-runing its corresponding analysis or tool for the NASSA analysis")
nassa_parser.add_argument(
    "-nseq", "--n_sequences",
    type=int,
    help="Number of sequences to be analyzed")
nassa_parser.add_argument(
    "-i", "--unit_len",
    type=int,
    default=6,
    help="Number of base pairs to be analyzed")
nassa_parser.add_argument(
    "-hp", "--helical_parameters",
    action='store_true',
    default=False,
    help="Run the helical parameters analysis")
nassa_parser.add_argument(
    "-pdirs", "--proj_directories",
    nargs='*',
    default=None,
    help=("Path to the different project directories. Each directory is to contain an independent project.\n"
        "Several output files will be generated in the same folder directory"))
nassa_parser.add_argument(
    "-all", "--all",
    action='store_true',
    default=False,
    help="Run all the helical parameters and NASSA analyses")
nassa_parser.add_argument(
    "-stru", "--input_structure_filepath",
    default=None,
    help=("Path to input structure file. It may be relative to the project or to each MD directory.\n"
        "If this value is not passed then the standard structure file is used as input by default"))
nassa_parser.add_argument(
    "-traj", "--input_trajectory_filepath",
    nargs='*',
    default=None,
    help=("Path to input trajectory file. It is relative to each MD directory.\n"
        "If this value is not passed then the standard trajectory file path is used as input by default"))
nassa_parser.add_argument(
    "-top", "--input_topology_filepath",
    default=None,  # There is no default since many formats may be possible
    help="Path to input topology file. It is relative to the project directory.")
nassa_parser.add_argument(
    "-mdir", "--md_directories",
    nargs='*',
    default=None,
    help=("Path to the different MD directories. Each directory is to contain an independent trajectory and structure.\n"
        "Several output files will be generated in every MD directory"))
nassa_parser.add_argument(
    "-t", "--trust",
    type=str,
    nargs='*',
    default=[],
    action=custom,
    const=AVAILABLE_CHECKINGS,
    choices=AVAILABLE_CHECKINGS,
    help="If passed, do not run the specified checking. Note that all checkings are skipped if passed alone.  Available checkings:" + pretty_list(AVAILABLE_CHECKINGS)
)
nassa_parser.add_argument(
    "-m", "--mercy",
    type=str,
    nargs='*',
    default=[],
    action=custom,
    const=AVAILABLE_FAILURES,
    choices=AVAILABLE_FAILURES,
    help=("If passed, do not kill the process when any of the specfied checkings fail and proceed with the workflow.\n"
        "Note that all checkings are allowed to fail if the argument is passed alone. Available checkings:" + pretty_list(AVAILABLE_FAILURES))
)
nassa_parser.add_argument(
    "-dup", "--duplicates",
    default=False,
    action='store_true',
    help="If passed, merge duplicate subunits if there is more than one, in the sequences. if not only the last will be selected"
)

# Dataset subcommand
dataset_parser = subparsers.add_parser("dataset", formatter_class=CustomHelpFormatter,
    help="Manage and process a dataset of MDDB projects.")
dataset_subparsers = dataset_parser.add_subparsers(dest='dataset_subcommand', help='Dataset subcommands')

# Dataset run subcommand
dataset_run_parser = dataset_subparsers.add_parser("run", formatter_class=CustomHelpFormatter,
help="Run the workflow for a dataset of MDDB projects.",
    parents=[common_parser])
dataset_run_parser.add_argument("dataset_yaml", help="Path to the dataset YAML file.")
dataset_run_parser.add_argument("--slurm", action="store_true", help="Submit the workflow to SLURM.")
dataset_run_parser.add_argument("-jt", "--job-template", help="Path to the SLURM job template file. Required if --slurm is used.")
dataset_run_parser.add_argument("-ig", "--include-groups", nargs='*', type=int, default=[], help="List of group IDs to be run.")
dataset_run_parser.add_argument("-eg", "--exclude-groups", nargs='*', type=int, default=[], help="List of group IDs to be excluded.")

# Dataset status subcommand
dataset_status_parser = dataset_subparsers.add_parser("groups", formatter_class=CustomHelpFormatter,
    help="Show the status of projects in a dataset, grouped by their last log message.",
    parents=[common_parser])
dataset_status_parser.add_argument("dataset_yaml", help="Path to the dataset YAML file.")
