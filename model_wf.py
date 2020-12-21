#!/usr/bin/env python

from os import chdir

import argparse
import json

import pytraj as pt

from topology_manager import TopologyReference
import index as model_index

from generic_analyses import rmsd, rmsf, rgyr
from pca import pca
from pca_contacts import pca_contacts
from rmsd_per_residue import rmsd_per_residue
from rmsd_pairwise import rmsd_pairwise
from distance_per_residue import distance_per_residue
from hydrogen_bonds import hydrogen_bonds
from energies import energies
from pockets import pockets


def required_length(nmin, nmax):
    """Custom Action to require minimum and maximum number of values for argument."""
    class RequiredLength(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not nmin <= len(values) <= nmax:
                msg = 'argument "{f}" requires between {nmin} and {nmax} arguments'.format(
                    f=self.dest, nmin=nmin, nmax=nmax)
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return RequiredLength


# define arguments
parser = argparse.ArgumentParser(description="MoDEL Workflow")
action = parser.add_mutually_exclusive_group(required=True)

# action arguments
action.add_argument(
    "-analyze_directory",
    help="Indicate directory to analyze, which is meant to "
    "include the topology and trajectory files. "
    "Topology and trajectory filenames can be indicated using "
    "-topology_filename and -trajectory_filename args.")

action.add_argument(
    "-analyze_project",
    nargs='+',
    action=required_length(2, 3),
    help="Indicate project name, directory and (optionally) URL from which to "
    "download topology and trajectory files. "
    "Topology and trajectory filenames can be indicated using "
    "-topology_filename and -trajectory_filename args.")

action.add_argument(
    "-single_analysis",
    choices=[
        "rmsd",
        "rmsf",
        "rgyr",
        "pca",
        "pca_contacts",
        "rmsd_per_residue",
        "rmsd_pairwise",
        "distance_per_residue",
        "hydrogen_bonds",
        "energies",
        "pockets"
    ],
    help="Perform single analysis.")

# optional arguments
parser.add_argument(
    "-output_analysis_filename",
    help="output filename for single analysis.")

parser.add_argument(
    "-pca_eigenval_file",
    default="pca.eigenval.xvg",
    help="output filename for eigenvalues in PCA analysis.")

parser.add_argument(
    "-pca_eigenvec_file",
    default="eigenvec.trr",
    help="output filename for eigenvectors in PCA analysis.")

parser.add_argument(
    "-topology_filename",
    default="md.imaged.rot.dry.pdb",
    help="path to topology filename")

parser.add_argument(
    "-trajectory_filename",
    default="md.imaged.rot.xtc",
    help="path to trajectory filename")

parser.add_argument(
    "-inputs_filename",
    default="inputs.json",
    help="path to inputs filename")

parser.add_argument(
    "-metadata_filename",
    default="metadata.json",
    help="path to metadata filename")

parser.add_argument(
    "-working_dir",
    help="directory where to perform single analysis. "
    "If empty, will use current directory.")

parser.add_argument(
    "-interface_cutoff_distance",
    type=int,
    default=5,
    help="cut-off distance to use when computing interfaces")

# flags
parser.add_argument(
    "--analysis_prep",
    action="store_true",
    help="execute analyses preparation.")

args = parser.parse_args()

# execution flow
if args.analyze_directory:
    model_index.analyze_directory(
        directory=args.analyze_directory,
        topology_filename=args.topology_filename,
        trajectory_filename=args.trajectory_filename)
elif args.analyze_project:
    model_index.analyze_project(
        *args.analyze_project,
        topology_filename=args.topology_filename,
        trajectory_filename=args.trajectory_filename)
elif args.single_analysis:
    chdir(args.working_dir)

    if args.analysis_prep:
        # compute from scratch
        print("Preparing analysis...")
        metadata, topology_reference, interfaces = model_index.analysis_prep(
            trajectory_filename=args.trajectory_filename,
            topology_filename=args.topology_filename,
            inputs_filename=args.inputs_filename,
            interface_cutoff_distance=args.interface_cutoff_distance)
    else:
        # load from files
        with open(args.metadata_filename, "r") as f:
            metadata = json.load(f)
        with open(args.inputs_filename, "r") as f:
            inputs = json.load(f)
        interfaces = inputs['interfaces']
        topology_reference = TopologyReference(args.topology_filename)

    pt_trajectory = pt.iterload(
        args.trajectory_filename, args.topology_filename)[0:2000:10]

    analysis_functions = {
        # "key": [func, args],
        "rmsd": (rmsd, {
            "input_first_frame_filename": "firstFrame.pdb",
            "input_trajectory_filename": args.trajectory_filename,
            "output_analysis_filename": args.output_analysis_filename or "md.rmsd.xvg"
        }),
        "rmsf": (rmsf, {
            "input_topology_filename": args.topology_filename,
            "input_trajectory_filename": args.trajectory_filename,
            "output_analysis_filename": args.output_analysis_filename or "md.rmsf.xvg"
        }),
        "rgyr": (rgyr, {
            "input_topology_filename": args.topology_filename,
            "input_trajectory_filename": args.trajectory_filename,
            "output_analysis_filename": args.output_analysis_filename or "md.rgyr.xvg"
        }),
        "pca": (pca, {
            "input_topology_filename": args.topology_filename,
            "input_trajectory_filename": args.trajectory_filename,
            "output_eigenvalues_filename": args.pca_eigenval_file,
            "output_eigenvectors_filename": args.pca_eigenvec_file,
            "snapshots": metadata["SNAPSHOTS"]
        }),
        "pca_contacts": (pca_contacts, {
            "trajectory": args.trajectory_filename,
            "topology": args.topology_filename,
            "interfaces": metadata["INTERFACES"],
            "output_analysis_filename": args.output_analysis_filename or "contacts_PCA.json"
        }),
        "rmsd_per_residue": (rmsd_per_residue, {
            "pt_trajectory": pt_trajectory,
            "output_analysis_filename": args.output_analysis_filename or "md.rmsd.perres.json",
            "topology_reference": topology_reference
        }),
        "rmsd_pairwise": (rmsd_pairwise, {
            "pt_trajectory": pt_trajectory,
            "output_analysis_filename": args.output_analysis_filename or "md.rmsd.pairwise.json",
            "interfaces": interfaces
        }),
        "distance_per_residue": (distance_per_residue, {
            "pt_trajectory": pt_trajectory,
            "output_analysis_filename": args.output_analysis_filename or "md.dist.perres.json",
            "interfaces": interfaces
        }),
        "hydrogen_bonds": (hydrogen_bonds, {
            "pt_trajectory": pt_trajectory,
            "output_analysis_filename": args.output_analysis_filename or "md.hbonds.json",
            "topology_reference": topology_reference,
            "interfaces": interfaces
        }),
        "energies": (energies, {
            "input_topology_filename": args.topology_filename,
            "input_trajectory_filename": args.trajectory_filename,
            "output_analysis_filename": args.output_analysis_filename or "md.energies.json",
            "reference": topology_reference,
            "snapshots": metadata["SNAPSHOTS"],
            "ligands": metadata["LIGANDS"]
        }),
        "pockets": (pockets, {
            "input_topology_filename": args.topology_filename,
            "input_trajectory_filename": args.trajectory_filename,
            "output_analysis_filename": args.output_analysis_filename or "md.pockets.json",
            "topology_reference": topology_reference,
            "snapshots": metadata["SNAPSHOTS"]
        })
    }

    # execute analysis_func
    print(f"\nExecuting analys function {args.single_analysis}...")
    analysis_function, analysis_args = analysis_functions[args.single_analysis]
    analysis_function(**analysis_args)

print("\nDone!")
