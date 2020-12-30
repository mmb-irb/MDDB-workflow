#!/usr/bin/env python

from pathlib import Path

import math
import argparse

import pytraj as pt

import model_workflow.analyses.index as model_index

from model_workflow.analyses.generic_analyses import rmsd, rmsf, rgyr
from model_workflow.analyses.pca import pca
from model_workflow.analyses.pca_contacts import pca_contacts
from model_workflow.analyses.rmsd_per_residue import rmsd_per_residue
from model_workflow.analyses.rmsd_pairwise import rmsd_pairwise
from model_workflow.analyses.distance_per_residue import distance_per_residue
from model_workflow.analyses.hydrogen_bonds import hydrogen_bonds
from model_workflow.analyses.energies import energies
from model_workflow.analyses.pockets import pockets


def required_length(nmin, nmax):
    """Require minimum and maximum number of values for argument."""
    class RequiredLength(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not nmin <= len(values) <= nmax:
                msg = (
                    f"argument \"{self.dest}\" requires between "
                    f"{nmin} and {nmax} arguments")
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return RequiredLength


def execute_workflow(args):
    start_args = {
        "directory": Path(args.working_dir).resolve(),
        "project": args.project,
        "url": args.url,
        "inputs_filename": args.inputs_filename}
    model_index.start(**start_args)

    trajectory_filename = "md.imaged.rot.xtc"
    # trajectory_filename = "md.trr"
    topology_filename = "md.imaged.rot.dry.pdb"
    if not args.single_analysis:
        # run all analyses
        model_index.run_analyses(
            topology_filename, trajectory_filename)
    else:
        # run analysis prep
        topology_reference, interactions, ligands, snapshots = model_index.analysis_prep(
            inputs_filename=args.inputs_filename)

        # get reduced trajectory
        pt_trajectory = pt.iterload(
            trajectory_filename,
            topology_filename)
        trajectory_frames = pt_trajectory.n_frames
        # Set a reduced trajectory used for heavy analyses
        reduced_pt_trajectory = None
        # First, set the maximum number of frames for the reduced trajectory
        reduced_trajectory_frames = 200
        # If the current trajectory has already less frames
        # than the maximum then use it also as reduced
        if trajectory_frames < reduced_trajectory_frames:
            reduced_pt_trajectory = pt_trajectory
            # Add a step value which will be required later
            reduced_pt_trajectory.step = 1
        # Otherwise, create a reduced trajectory with
        # as much frames as specified above
        # These frames are picked along the trajectory
        else:
            # Calculate how many frames we must jump between each
            # reduced frame to never exceed the limit
            # The '- 1' is because the first frame is 0
            # (you have to do the math to understand)
            step = math.floor(trajectory_frames /
                              (reduced_trajectory_frames - 1))
            reduced_pt_trajectory = pt_trajectory[0:trajectory_frames:step]
            # Add the step value to the reduced trajectory explicitly.
            # It will be required later
            reduced_pt_trajectory.step = step

        analysis_functions = {
            # "key": [func, args],
            "rmsd": (rmsd, {
                "input_first_frame_filename": "firstFrame.pdb",
                "input_trajectory_filename": trajectory_filename,
                "output_analysis_filename": "md.rmsd.xvg"
            }),
            "rmsf": (rmsf, {
                "input_topology_filename": topology_filename,
                "input_trajectory_filename": trajectory_filename,
                "output_analysis_filename": "md.rmsf.xvg"
            }),
            "rgyr": (rgyr, {
                "input_topology_filename": topology_filename,
                "input_trajectory_filename": trajectory_filename,
                "output_analysis_filename": "md.rgyr.xvg"
            }),
            "pca": (pca, {
                "input_topology_filename": topology_filename,
                "input_trajectory_filename": trajectory_filename,
                "output_eigenvalues_filename": "pca.eigenval.xvg",
                "output_eigenvectors_filename": "eigenvec.trr",
                "snapshots": snapshots
            }),
            "pca_contacts": (pca_contacts, {
                "trajectory": trajectory_filename,
                "topology": topology_filename,
                "interactions": interactions,
                "output_analysis_filename": "contacts_PCA.json"
            }),
            "rmsd_per_residue": (rmsd_per_residue, {
                "pt_trajectory": reduced_pt_trajectory,
                "output_analysis_filename": "md.rmsd.perres.json",
                "topology_reference": topology_reference
            }),
            "rmsd_pairwise": (rmsd_pairwise, {
                "pt_trajectory": reduced_pt_trajectory,
                "output_analysis_filename": "md.rmsd.pairwise.json",
                "interactions": interactions
            }),
            "distance_per_residue": (distance_per_residue, {
                "pt_trajectory": reduced_pt_trajectory,
                "output_analysis_filename": "md.dist.perres.json",
                "interactions": interactions
            }),
            "hydrogen_bonds": (hydrogen_bonds, {
                "pt_trajectory": reduced_pt_trajectory,
                "output_analysis_filename": "md.hbonds.json",
                "topology_reference": topology_reference,
                "interactions": interactions
            }),
            "energies": (energies, {
                "input_topology_filename": topology_filename,
                "input_trajectory_filename": trajectory_filename,
                "output_analysis_filename": "md.energies.json",
                "reference": topology_reference,
                "snapshots": snapshots,
                "ligands": ligands
            }),
            "pockets": (pockets, {
                "input_topology_filename": topology_filename,
                "input_trajectory_filename": trajectory_filename,
                "output_analysis_filename": "md.pockets.json",
                "topology_reference": topology_reference,
                "snapshots": snapshots,
            })
        }

        # execute single analysis function
        print(f"\nExecuting analysis function {args.single_analysis}...")
        analysis_func, analysis_args = analysis_functions[args.single_analysis]
        analysis_func(**analysis_args)

    print("\nDone!")


# define arguments
parser = argparse.ArgumentParser(description="MoDEL Workflow")

# optional args
parser.add_argument(
    "-dir", "--working_dir",
    default=Path.cwd(),
    help="directory where to perform analysis. "
    "If empty, will use current directory.")

parser.add_argument(
    "-p", "--project",
    default=None,
    help="If given a project name, trajectory and "
    "topology files will be downloaded from remote server.")

parser.add_argument(
    "-url",
    default="https://bioexcel-cv19-dev.bsc.es",
    help="URL from where to download project")

parser.add_argument(
    "-in", "--inputs_filename",
    default="inputs.json",
    help="path to inputs filename")

parser.add_argument(
    "-a", "--single_analysis",
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
        "pockets"],
    help="Indicate single analysis to perform.")


def main():
    args = parser.parse_args()
    execute_workflow(args)
