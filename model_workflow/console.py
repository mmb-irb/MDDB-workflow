from pathlib import Path
from os.path import exists
from shutil import copyfile
from subprocess import call
from typing import List
from argparse import ArgumentParser, RawTextHelpFormatter, Action
from textwrap import wrap

from model_workflow.mwf import workflow, Project, requestables, DEPENDENCY_FLAGS
from model_workflow.utils.structures import Structure
from model_workflow.utils.file import File
from model_workflow.utils.conversions import convert
from model_workflow.utils.filters import filter_atoms
from model_workflow.utils.subsets import get_trajectory_subset
from model_workflow.utils.constants import *
from model_workflow.utils.auxiliar import InputError
from model_workflow.utils.nassa_file import generate_nassa_config
from model_workflow.analyses.nassa import workflow_nassa

# Set the path to the input setter jupyter notebook
inputs_template = str(Path(__file__).parent / "resources" / "inputs_file_template.yml")
nassa_template = str(Path(__file__).parent / "resources" / "nassa_template.yml")

expected_project_args = set(Project.__init__.__code__.co_varnames)

test_docs_url = 'https://mddb-workflow.readthedocs.io/en/latest/usage.html#tests-and-other-checking-processes'
task_docs_url = 'https://mddb-workflow.readthedocs.io/en/latest/tasks.html'

class CustomHelpFormatter(RawTextHelpFormatter):
    """Custom formatter for argparse help text with better organization and spacing"""
    
    def __init__(self, prog, indent_increment=2, max_help_position=6, width=None):
        super().__init__(prog, indent_increment, max_help_position, width)
        
    def _split_lines(self, text, width):
        lines = []
        for line in text.splitlines():
            if line.strip() != '':
                lines.extend(wrap(line, width, break_long_words=False, replace_whitespace=False))
        return lines

    def _format_usage(self, usage, actions, groups, prefix):
        essential_usage = super()._format_usage(usage, actions, groups, prefix)
        # Remove lines about arguments -i, -e, -ow from usage
        lines = essential_usage.split('\n')
        filtered_lines = []
        for line in lines:
            if line.strip().startswith("[-i"):
                line = line.replace("[-i", "[-i/-e/-ow")
                filtered_lines.append(line)
            elif line.strip().startswith("[-e") or line.strip().startswith("[-ow"):
                continue
            else:
                filtered_lines.append(line)
        essential_usage = '\n'.join(filtered_lines)
        return essential_usage
    
    def _format_action_invocation(self, action):
        """Format the display of options with choices more cleanly"""
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
                choice_str = '{' + ','.join(sorted(str(c) for c in action.choices)) + '}'
                # if action.nargs is not None and action.nargs != 1:
                #     choice_str += ' ...'
                return f"{opts} [{choice_str}]"
            else:
                return f"{opts} {metavar}"

    
# Main ---------------------------------------------------------------------------------            

# Function called through argparse
def main ():
    # Parse input arguments from the console
    # The vars function converts the args object to a dictionary
    args = parser.parse_args()
    # Apply common arguments as necessary
    if hasattr(args, 'no_symlinks') and args.no_symlinks:
        GLOBALS['no_symlinks'] = True
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
        common_args = [ action.dest for action in common_parser._actions ]
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
        workflow(project_parameters = project_args, **workflow_args)
    # If user wants to setup the inputs
    elif subcommand == "inputs":
        # Make a copy of the template in the local directory if there is not an inputs file yet
        if exists(DEFAULT_INPUTS_FILENAME):
            print(f"File {DEFAULT_INPUTS_FILENAME} already exists")
        else:
            copyfile(inputs_template, DEFAULT_INPUTS_FILENAME)
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
            print(f"No editor was selected")
        

    # In case the convert tool was called
    elif subcommand == 'convert':
        # If no input arguments are passed print help
        if args.input_structure == None and args.input_trajectories == None:
            convert_parser.print_help()
            return
        if args.input_trajectories == None:
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
            input_structure_file = File(args.input_structure),
            output_structure_file = File(args.output_structure),
            input_trajectory_file = File(args.input_trajectory),
            output_trajectory_file = File(args.output_trajectory),
            selection_string = args.selection_string,
            selection_syntax = args.selection_syntax
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
    # If user wants to run the NASSA analysis
    elif subcommand == "nassa":
        # If no input arguments are passed print help
        if args.config == None and args.make_config == False:
            nassa_parser.print_help()
            print('Please provide a configuration file or make one with the -m flag')
            return
        # If the user wants to make a configuration file
        if args.make_config:
            #print('args.make_config: ', args.make_config)
            if args.make_config == True or args.make_config == []:
            # Make a copy of the template in the local directory if there is not an inputs file yet
                if not exists(DEFAULT_NASSA_CONFIG_FILENAME):
                    copyfile(nassa_template, DEFAULT_NASSA_CONFIG_FILENAME)
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
        if args.config and args.analysis_names == None and args.all == False:
            nassa_parser.print_help()
            print('Please provide an analysis name to run:', ', '.join(NASSA_ANALYSES_LIST))
            return
        # If the user wants to run the helical parameters analysis we must check if the necessary files are provided (structure, topology and trajectory)
        if args.helical_parameters:
            # Also, it is necesary to provide the project directories. Each of the project directories must contain an independent MD
            if args.proj_directories == None:
                nassa_parser.print_help()
                print('Please provide a project directory to run the helical parameters analysis')
                return
            if args.input_structure_filepath == None:
                raise InputError('Please provide a structure file to run the helical parameters analysis')
            elif args.input_trajectory_filepath == None:
                raise InputError('Please provide a trajectory file to run the helical parameters analysis')
            elif args.input_topology_filepath == None:
                raise InputError('Please provide a topology file to run the helical parameters analysis')
            # If the all flag is set, the user must provide the path to the sequences because it is necessary to create the nassa.yml and run the NASSA analysis
            if args.all:
                if not args.seq_path:
                    raise InputError('Please, if all option is selected provide the path to the sequences (--seq_path)')
            # If all the flags are correctly set, we can run the analysis
            workflow_nassa(
                config_file_path=None, # The configuration file is not needed in this case because we are going to run the helical parameters analysis so it will be created then
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
            del dict_args['subcommand'] # preguntar Dani ¿?
            # Call the actual main function
            workflow_nassa(
                    config_file_path = args.config, 
                    analysis_names = args.analysis_names, 
                    make_config = args.make_config, 
                    output =  args.output,
                    working_directory = args.working_directory,
                    overwrite = args.overwrite, 
                    overwrite_nassa = args.overwrite_nassa,
                    n_sequences = args.n_sequences,
                    unit_len = args.unit_len,
                    all= args.all,
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

# Define console arguments to call the workflow
parser = ArgumentParser(description="MDDB Workflow", formatter_class=RawTextHelpFormatter)
subparsers = parser.add_subparsers(help='Name of the subcommand to be used', dest="subcommand")

# Set the run subcommand
run_parser = subparsers.add_parser("run",
    help="Run the workflow",
    formatter_class=CustomHelpFormatter,
    parents=[common_parser]
)

# Set optional arguments
run_parser_input_group = run_parser.add_argument_group('INPUT OPTIONS')
run_parser_input_group.add_argument(
    "-top", "--input_topology_filepath",
    default=None, # There is no default since many formats may be possible
    help="Path to input topology file. It is relative to the project directory.")
run_parser_input_group.add_argument(
    "-stru", "--input_structure_filepath",
    default=None,
    help=("Path to input structure file. It may be relative to the project or to each MD directory.\n"
        "If this value is not passed then the standard structure file is used as input by default"))
run_parser_input_group.add_argument(
    "-traj", "--input_trajectory_filepaths",
    #type=argparse.FileType('r'),
    nargs='*',
    default=None,
    help=("Path to input trajectory file. It is relative to each MD directory.\n"
        "If this value is not passed then the standard trajectory file path is used as input by default"))
run_parser_input_group.add_argument(
    "-dir", "--working_directory",
    default='.',
    help="Directory where the whole workflow is run. Current directory by default.")
run_parser_input_group.add_argument(
    "-mdir", "--md_directories",
    nargs='*',
    default=None,
    help=("Path to the different MD directories. Each directory is to contain an independent trajectory and structure.\n"
        "Several output files will be generated in every MD directory")
)
run_parser_input_group.add_argument(
    "-md", "--md_config",
    action='append',
    nargs='*',
    default=None,
    help=("Configuration of a specific MD. You may declare as many as you want."
          "Every MD requires a directory name, a structure path and at least one trajectory path."
          "The structure is -md <directory> <structure> <trajectory 1> <trajectory 2>... "
          "Note that all trajectories from the same MD will be merged.")
)
run_parser_input_group.add_argument(
    "-proj", "--accession",
    default=None,
    help="Project accession to download missing input files from the database.")
run_parser_input_group.add_argument(
    "-url", "--database_url",
    default=DEFAULT_API_URL,
    help=f"API URL to download missing data. Default value is {DEFAULT_API_URL}")
run_parser_input_group.add_argument(
    "-inp", "--inputs_filepath",
    default=None,
    help="Path to inputs file")
run_parser_input_group.add_argument(
    "-pop", "--populations_filepath",
    default=DEFAULT_POPULATIONS_FILENAME,
    help="Path to equilibrium populations file (Markov State Model only)")
run_parser_input_group.add_argument(
    "-tpro", "--transitions_filepath",
    default=DEFAULT_TRANSITIONS_FILENAME,
    help="Path to transition probabilities file (Markov State Model only)")
# Set a group for the workflow control options
run_parser_workflow_group = run_parser.add_argument_group('WORKFLOW CONTROL OPTIONS')
run_parser_workflow_group.add_argument(
    "-img", "--image",
    action='store_true',
    help="Set if the trajectory is to be imaged")
run_parser_workflow_group.add_argument(
    "-fit", "--fit",
    action='store_true',
    help="Set if the trajectory is to be fitted (both rotation and translation)")
run_parser_workflow_group.add_argument(
    "-trans", "--translation",
    nargs='*',
    default=[0,0,0],
    help=("Set the x y z translation for the imaging process\n"
        "e.g. -trans 0.5 -1 0"))
run_parser_workflow_group.add_argument(
    "-d", "--download",
    action='store_true',
    help="If passed, only download required files. Then exits.")
run_parser_workflow_group.add_argument(
    "-s", "--setup",
    action='store_true',
    help="If passed, only download required files and run mandatory dependencies. Then exits.")
run_parser_workflow_group.add_argument(
    "-smp", "--sample_trajectory",
    type=int,
    nargs='?',
    default=None,
    const=10,
    help="If passed, download just a few frames (10 by default) from the trajectory instead of it all")
run_parser_workflow_group.add_argument(
    "-rcut", "--rmsd_cutoff",
    type=float,
    default=DEFAULT_RMSD_CUTOFF,
    help=(f"Set the cutoff for the RMSD sudden jumps analysis to fail (default {DEFAULT_RMSD_CUTOFF}).\n"
        "This cutoff stands for the number of standard deviations away from the mean an RMSD value is to be.\n"))
run_parser_workflow_group.add_argument(
    "-icut", "--interaction_cutoff",
    type=float,
    default=DEFAULT_INTERACTION_CUTOFF,
    help=(f"Set the cutoff for the interactions analysis to fail (default {(DEFAULT_INTERACTION_CUTOFF)}).\n"
        "This cutoff stands for percent of the trajectory where the interaction happens (from 0 to 1).\n"))
run_parser_workflow_group.add_argument(
    "-iauto", "--interactions_auto",
    type=str,
    nargs='?',
    const=True,
    help=("""Guess input interactions automatically. Available options:
  - greedy (default): All chains against all chains
  - humble: If there are only two chains then select the interaction between them
  - <chain letter>: All chains against this specific chain
  - ligands: All chains against every ligand""")
)

# Set a group for the selection options
run_parser_selection_group = run_parser.add_argument_group('SELECTION OPTIONS')
run_parser_selection_group.add_argument(
    "-filt", "--filter_selection",
    nargs='?',
    default=False,
    const=True,
    help=("Atoms selection to be filtered in VMD format. "
        "If the argument is passed alone (i.e. with no selection) then water and counter ions are filtered"))
run_parser_selection_group.add_argument(
    "-pbc", "--pbc_selection",
    default=None,
    help=("Selection of atoms which stay in Periodic Boundary Conditions even after imaging the trajectory. "
        "e.g. remaining solvent and counter ion molecules, membrane lipids, etc."))
run_parser_selection_group.add_argument(
    "-cg", "--cg_selection",
    default=None,
    help="Selection of atoms which are not actual atoms but Coarse Grained beads.")
run_parser_selection_group.add_argument(
    "-pcafit", "--pca_fit_selection",
    default=PROTEIN_AND_NUCLEIC_BACKBONE,
    help="Atom selection for the pca fitting in vmd syntax")
run_parser_selection_group.add_argument(
    "-pcasel", "--pca_selection",
    default=PROTEIN_AND_NUCLEIC_BACKBONE,
    help="Atom selection for pca analysis in vmd syntax")


# Set a custom argparse action to handle the following 2 arguments
# This is done becuase it is not possible to combine nargs='*' with const
# https://stackoverflow.com/questions/72803090/argparse-how-to-create-the-equivalent-of-const-with-nargs
class custom (Action):
    # If argument is not passed -> default
    # If argument is passed empty -> const
    # If argument is passed with values -> values
    def __call__(self, parser, namespace, values, option_string=None):
        if values:
            setattr(namespace, self.dest, values)
        else:
            setattr(namespace, self.dest, self.const)

# Set a function to pretty print a list of available checkings / failures
def pretty_list (availables : List[str]) -> str:
    final_line = 'Available protocols:'
    for available in availables:
        nice_name = NICE_NAMES.get(available, None)
        if not nice_name:
            raise Exception('Flag "' + available + '" has not a defined nice name')
        final_line += '\n  - ' + available + ' -> ' +  nice_name
    final_line += f'\nTo know more about each test please visit:\n{test_docs_url}'
    return final_line

run_parser_checks_group = run_parser.add_argument_group('INPUT CHECKS OPTIONS', description=f"For more information about each check please visit:\n{test_docs_url}")
run_parser_checks_group.add_argument(
    "-t", "--trust",
    type=str,
    nargs='*',
    default=[],
    action=custom,
    const=AVAILABLE_CHECKINGS,
    choices=AVAILABLE_CHECKINGS,
    help="If passed, do not run the specified checking. Note that all checkings are skipped if passed alone. " + pretty_list(AVAILABLE_CHECKINGS)
)
run_parser_checks_group.add_argument(
    "-m", "--mercy",
    type=str,
    nargs='*',
    default=[],
    action=custom,
    const=AVAILABLE_FAILURES,
    choices=AVAILABLE_FAILURES,
    help=("If passed, do not kill the process when any of the specfied checkings fail and proceed with the workflow. "
        "Note that all checkings are allowed to fail if the argument is passed alone. " + pretty_list(AVAILABLE_FAILURES))
)

# Set a list with the alias of all requestable dependencies
choices = sorted(list(requestables.keys()) + list(DEPENDENCY_FLAGS.keys()))

run_parser_analysis_group = run_parser.add_argument_group('TASKS OPTIONS', description=f"Available tasks: {choices}\nFor more information about each task please visit:\n{task_docs_url}")
run_parser_analysis_group.add_argument(
    "-i", "--include",
    nargs='*',
    choices=choices,
    help="""Set the unique analyses or tools to be run. All other steps will be skipped. There are also some additional flags to define a preconfigured group of dependencies:
  - download: Check/download input files (already ran with analyses)
  - setup: Process and test input files (already ran with analyses)
  - network: Run dependencies which require internet connection
  - minimal: Run dependencies required by the web client to work
  - interdeps: Run interactions and all its dependent analyses""")
run_parser_analysis_group.add_argument(
    "-e", "--exclude",
    nargs='*',
    choices=choices,
    help=("Set the analyses or tools to be skipped. All other steps will be run. Note that if we exclude a dependency of something else then it will be run anyway."))
run_parser_analysis_group.add_argument(
    "-ow", "--overwrite",
    type=str,
    nargs='*',
    default=[],
    action=custom,
    const=True,
    choices=choices,
    help=("Set the output files to be overwritten thus re-runing its corresponding analysis or tool. Use this flag alone to overwrite everything."))


# Add a new to command to aid in the inputs file setup
inputs_parser = subparsers.add_parser("inputs",
    help="Set the inputs file",
    formatter_class=RawTextHelpFormatter,
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
        "If several input trajectories are passed they will be merged previously.",
    parents=[common_parser])
convert_parser.add_argument(
    "-is", "--input_structure",
    help="Path to input structure file")
convert_parser.add_argument(
    "-os", "--output_structure",
    help="Path to output structure file")
convert_parser.add_argument(
    "-it", "--input_trajectories", nargs='*',
    help="Path to input trajectory file(s)")
convert_parser.add_argument(
    "-ot", "--output_trajectory",
    help="Path to output trajectory file")

# The filter command
filter_parser = subparsers.add_parser("filter",
    help="Filter atoms in a structure and/or a trajectory\n",
    parents=[common_parser])
filter_parser.add_argument(
    "-is", "--input_structure", required=True,
    help="Path to input structure file")
filter_parser.add_argument(
    "-os", "--output_structure",
    help="Path to output structure file")
filter_parser.add_argument(
    "-it", "--input_trajectory",
    help="Path to input trajectory file")
filter_parser.add_argument(
    "-ot", "--output_trajectory",
    help="Path to output trajectory file")
filter_parser.add_argument(
    "-sel", "--selection_string",
    help="Atom selection")
filter_parser.add_argument(
    "-syn", "--selection_syntax", default='vmd',
    help="Atom selection syntax (vmd by default)")

# The subset command
subset_parser = subparsers.add_parser("subset",
    help="Get a subset of frames from the current trajectory.",
    parents=[common_parser])
subset_parser.add_argument(
    "-is", "--input_structure", required=True,
    help="Path to input structure file")
subset_parser.add_argument(
    "-it", "--input_trajectory",
    help="Path to input trajectory file")
subset_parser.add_argument(
    "-ot", "--output_trajectory",
    help="Path to output trajectory file")
subset_parser.add_argument(
    "-start", "--start", type=int, default=0,
    help="Start frame")
subset_parser.add_argument(
    "-end", "--end", type=int, default=None,
    help="End frame")
subset_parser.add_argument(
    "-step", "--step", type=int, default=1,
    help="Frame step")
subset_parser.add_argument(
    "-skip", "--skip", nargs='*', type=int, default=[],
    help="Frames to be skipped")
subset_parser.add_argument(
    "-fr", "--frames", nargs='*', type=int, default=[],
    help="Frames to be returned. Input frame order is ignored as original frame order is conserved.")

# The chainer command
chainer_parser = subparsers.add_parser("chainer",
    help="Edit structure (pdb) chains",
    parents=[common_parser])
chainer_parser.add_argument(
    "-is", "--input_structure", required=True,
    help="Path to input structure file")
chainer_parser.add_argument(
    "-os", "--output_structure", default='chained.pdb',
    help="Path to output structure file")
chainer_parser.add_argument(
    "-sel", "--selection_string",
    help="Atom selection (the whole structure by default)")
chainer_parser.add_argument(
    "-syn", "--selection_syntax", default='vmd',
    choices=Structure.SUPPORTED_SELECTION_SYNTAXES,
    help="Atom selection syntax (VMD syntax by default)")
chainer_parser.add_argument(
    "-let", "--letter",
    help="New chain letter (one letter per fragment by default)")
chainer_parser.add_argument(
    "-whfr", "--whole_fragments", type=bool, default=True,
    help="Consider fragments beyond the atom selection. Otherwise a fragment could end up having multiple chains.")

# The NASSA commands 
nassa_parser = subparsers.add_parser("nassa", formatter_class=RawTextHelpFormatter,
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
    #type=str,
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
    default=None, # There is no default since many formats may be possible
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
    help="If passed, do not run the specified checking. Note that all checkings are skipped if passed alone.\n" + pretty_list(AVAILABLE_CHECKINGS)
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
        "Note that all checkings are allowed to fail if the argument is passed alone.\n" + pretty_list(AVAILABLE_FAILURES))
)