from subprocess import Popen
from pathlib import Path

from argparse import ArgumentParser, RawTextHelpFormatter, Action
from model_workflow.mwf import workflow, Project, requestables
from model_workflow.constants import *

# Set the path to the input setter jupyter notebook
input_setter = str(Path(__file__).parent / "utils" / "input_setter.ipynb")

expected_project_args = set(Project.__init__.__code__.co_varnames)

# Main ---------------------------------------------------------------------------------            

# Function called through argparse
def main ():
    # Parse input arguments from the console
    # The vars function converts the args object to a dictionary
    args = parser.parse_args()
    # Find which subcommand was called
    subcommand = args.subcommand
    # If there is not subcommand then print help
    if not subcommand:
        parser.print_help()
        return
    # If user wants to run the workflow
    if subcommand == "run":
        dict_args = vars(args)
        del dict_args['subcommand']
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
        command = "jupyter-notebook " + input_setter
        Popen(command, shell=True)


# Define console arguments to call the workflow
parser = ArgumentParser(description="MoDEL Workflow", formatter_class=RawTextHelpFormatter)
subparsers = parser.add_subparsers(help='Name of the subcommand to be used', dest="subcommand")

# Set the run subcommand
run_parser = subparsers.add_parser("run", help="Run the workflow")

# Set optional arguments
run_parser.add_argument(
    "-dir", "--working_directory",
    default='.',
    help="Directory where the whole workflow is run. Current directory by default.")

run_parser.add_argument(
    "-mdir", "--md_directories",
    nargs='*',
    default=None,
    help="Path to the different MD directories. Each directory is to contain an independent trajectory and structure."
)

run_parser.add_argument(
    "-proj", "--accession",
    default=None,
    help="Project accession to download missing input files from the database.")

run_parser.add_argument(
    "-url", "--database_url",
    default=DEFAULT_API_URL,
    help="API URL to download missing data")

run_parser.add_argument(
    "-stru", "--input_structure_filepath",
    default=STRUCTURE_FILENAME,
    help="Path to input structure file. It may be relative to the project or to each MD directory.")

run_parser.add_argument(
    "-traj", "--input_trajectory_filepaths",
    #type=argparse.FileType('r'),
    nargs='*',
    default=TRAJECTORY_FILENAME,
    help="Path to input trajectory file. It is relative to each MD directory.")

run_parser.add_argument(
    "-top", "--input_topology_filepath",
    default=None, # There is no default since many formats may be possible
    help="Path to input topology file. It is relative to the project directory.")

run_parser.add_argument(
    "-inp", "--inputs_filepath",
    default=DEFAULT_INPUTS_FILENAME,
    help="Path to inputs file")

run_parser.add_argument(
    "-pop", "--populations_filepath",
    default=DEFAULT_POPULATIONS_FILENAME,
    help="Path to equilibrium populations file (Markov State Model only)")

run_parser.add_argument(
    "-tpro", "--transitions_filepath",
    default=DEFAULT_TRANSITIONS_FILENAME,
    help="Path to transition probabilities file (Markov State Model only)")

run_parser.add_argument(
    "-img", "--image",
    action='store_true',
    help="Set if the trajectory is to be imaged")

run_parser.add_argument(
    "-fit", "--fit",
    action='store_true',
    help="Set if the trajectory is to be fitted (both rotation and translation)")

run_parser.add_argument(
    "-trans", "--translation",
    nargs='*',
    default=[0,0,0],
    help=("Set the x y z translation for the imaging process\n"
        "e.g. -trans 0.5 -1 0"))

run_parser.add_argument(
    "-filt", "--filter_selection",
    nargs='?',
    default=False,
    const=True,
    help=("Atoms selection to be filtered in VMD format\n"
        "If the argument is passed alone (i.e. with no selection) then water and counter ions are filtered"))

run_parser.add_argument(
    "-pcafit", "--pca_fit_selection",
    default=PROTEIN_AND_NUCLEIC_BACKBONE,
    help="Atom selection for the pca fitting in vmd syntax")

run_parser.add_argument(
    "-pcasel", "--pca_selection",
    default=PROTEIN_AND_NUCLEIC_BACKBONE,
    help="Atom selection for pca analysis in vmd syntax")

run_parser.add_argument(
    "-d", "--download",
    action='store_true',
    help="If passed, only download required files. Then exits.")

run_parser.add_argument(
    "-s", "--setup",
    action='store_true',
    help="If passed, only download required files and run mandatory dependencies. Then exits.")

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

run_parser.add_argument(
    "-t", "--trust",
    type=str,
    nargs='*',
    default=[],
    action=custom,
    const=AVAILABLE_CHECKINGS,
    choices=AVAILABLE_CHECKINGS,
    help=("If passed, do not run the specified checking. Note that all checkings are skipped if passed alone.\n"
        "Available protocols:\n"
        "- stabonds - Stable bonds\n"
        "- cohbonds - Coherent bonds\n"
        "- intrajrity - Trajectory integrity")
)

run_parser.add_argument(
    "-m", "--mercy",
    type=str,
    nargs='*',
    default=[],
    action=custom,
    const=AVAILABLE_FAILURES,
    choices=AVAILABLE_FAILURES,
    help=("If passed, do not kill the process when any of the specfied checkings fail and proceed with the workflow.\n"
        "Note that all checkings are allowed to fail if the argument is passed alone.\n"
        "Available protocols:\n"
        "- stabonds - Stable bonds\n"
        "- cohbonds - Coherent bonds\n"
        "- intrajrity - Trajectory integrity\n"
        "- refseq - Reference sequence")
)

run_parser.add_argument(
    "-smp", "--sample_trajectory",
    action='store_true',
    help="If passed, download just the 10 first frames of the trajectory instead of it all")

# Set a list with the alias of all requestable dependencies
choices = list(requestables.keys())

run_parser.add_argument(
    "-i", "--include",
    nargs='*',
    choices=choices,
    help="Set the unique analyses or tools to be run. All other steps will be skipped")

run_parser.add_argument(
    "-e", "--exclude",
    nargs='*',
    choices=choices,
    help="Set the analyses or tools to be skipped. All other steps will be run")


run_parser.add_argument(
    "-rcut", "--rmsd_cutoff",
    type=float,
    default=DEFAULT_RMSD_CUTOFF,
    help=("Set the cutoff for the RMSD sudden jumps analysis to fail (default " + str(DEFAULT_RMSD_CUTOFF) + ").\n"
        "This cutoff stands for the number of standard deviations away from the mean an RMSD value is to be.\n"))

run_parser.add_argument(
    "-icut", "--interaction_cutoff",
    type=float,
    default=DEFAULT_INTERACTION_CUTOFF,
    help=("Set the cutoff for the interactions analysis to fail (default " + str(DEFAULT_INTERACTION_CUTOFF) + ").\n"
        "This cutoff stands for percent of the trajectory where the interaction happens (from 0 to 1).\n"))

# Add a new to command to aid in the inputs.json setup
inputs_parser = subparsers.add_parser("inputs", help="Set the inputs.json file")
