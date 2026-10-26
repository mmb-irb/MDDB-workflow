import mdtraj as mdt
import pytraj as pt
import numpy as np
import math
import sys
from os.path import abspath

# Visual output tools
from tqdm import tqdm
import plotext as plt
from sklearn.mixture import GaussianMixture

from mddb_workflow.utils.auxiliar import delete_previous_log, reprint, TestFailure, warn, fail
from mddb_workflow.utils.auxiliar import round_to_hundredths
from mddb_workflow.utils.constants import TRAJECTORY_INTEGRITY_FLAG, BIMODALITY_BIC_CUTOFF, BIMODALITY_SEPARATION_CUTOFF
from mddb_workflow.utils.constants import RED_HEADER, COLOR_END
from mddb_workflow.utils.pyt_spells import get_pytraj_trajectory
from mddb_workflow.utils.mdt_spells import iterate_xtc_frames, unrotate_box
from mddb_workflow.utils.selections import Selection
from mddb_workflow.analyses.mindist import remove_view_script, MIN_FOCUS_SIZE
from mddb_workflow.tools.generate_topology import format_selection
from mddb_workflow.utils.type_hints import *

# Number of bonds which change the most in the biggest jump to be highlighted in the VMD script
HIGHLIGHTED_BONDS = 5

# Number of failed fragments to be shown in the logs, since there may be thousands
SHOWN_FAILED_FRAGMENTS = 3

# The meta-fragment is not a real molecule, so it must be explained when it fails
META_FRAGMENT_NAME = 'Inter-fragment distances'
META_FRAGMENT_EXPLANATION = ('Note that "Inter-fragment distances" is not a real molecule, but the first atom of every molecule not under PBC.'
    ' Its distances are between the first atom of every molecule and the next one, which should not change suddenly.'
    ' If they do, then whole molecules jumped at once (e.g. across the box boundaries).')

# Check if the output is going to a terminal or not
is_terminal = sys.stdout.isatty()

# LORE
# This test was originaly intended to use a RMSD jump cutoff based on number of atoms and timestep
# However, after a deep study, it was observed that simulations with similar features may show very different RMSD jumps
# For this reason now we comptue RMSD jumps along the whole trajectory and check that the biggest jump is not an outlier
# The outlier is defined according to how many times the standard deviation far from the mean is a value

# Look for sudden raises of RMSd values from one frame to another
# To do so, we check the RMSD of every frame using its previous frame as reference
# This was allowed thanks to people from MDtraj https://github.com/mdtraj/mdtraj/issues/1966
def check_trajectory_integrity (
    input_structure_filename : str,
    input_trajectory_filename : str,
    view_script_filename : str,
    structure : 'Structure',
    pbc_selection : 'Selection',
    mercy : list[str],
    trust: list[str],
    register : 'Register',
    #time_length : float,
    check_selection : str,
    # DANI: He visto saltos 'correctos' pasar de 6
    # DANI: He visto saltos 'incorrectos' no bajar de 10
    standard_deviations_cutoff : float,
    snapshots : int,
    # Used only to suggest how to skip non-equilibrated frames
    md_index : int,
    skipped_frames : int) -> bool:

    # Skip the test if we trust
    if TRAJECTORY_INTEGRITY_FLAG in trust:
        return True

    # Skip the test if it is already passed according to the register
    if register.tests.get(TRAJECTORY_INTEGRITY_FLAG, None):
        return True

    # If the trajectory has only 1 or 2 frames then there is no test to do
    if snapshots < 3:
        register.update_test(TRAJECTORY_INTEGRITY_FLAG, 'na')
        return True

    # Remove old warnings
    register.remove_warnings(TRAJECTORY_INTEGRITY_FLAG)
    remove_view_script(view_script_filename)

    # Parse the selection in VMD selection syntax
    parsed_selection = structure.select(check_selection, syntax='vmd')

    # If there is nothing to check then warn the user and stop here
    if not parsed_selection:
        raise Exception('WARNING: There are not atoms to be analyzed for the RMSD analysis')

    # # Discard PBC residues from the selection to be checked
    # parsed_selection -= pbc_selection

    # # If there is nothing to check then warn the user and stop here
    # if not parsed_selection:
    #     warn('There are no atoms to be analyzed for the RMSD checking after PBC substraction')
    #     register.update_test(TRAJECTORY_INTEGRITY_FLAG, 'na')
    #     return True

    # Get fragments out of the parsed selection
    # Fragments will be analyzed independently
    # A sudden jump of a small fragment may cause a small RMSD perturbation when
    #  it is part of a large fragment which is not jumping
    # Jumps or artifacts of partial fragments along boundaries are rare imaging problems, althoug possible
    # For this reason splitting the test in fragments is a good meassure
    # Although big fragments could be splitted in even smaller parts in the future
    # In the other hand a very small fragment may overcome the cutoff with small RMSD jumps, so beware
    bonded_fragments = {}
    unbonded_fragments = {}

    # If the structure is missing bonds then there are no "fragments" to select
    # In this case we just split the structure in chains
    missing_bonds_selection = structure.select_missing_any_bonds()
    if missing_bonds_selection:
        # Parse the selection to chain indices
        missing_bonds_chain_indices = structure.get_selection_chain_indices(missing_bonds_selection)
        # Make sure the selection wraps whole chains, for later
        missing_bonds_selection = structure.select_chain_indices(missing_bonds_chain_indices)
        # Now split the definitive selection in chains
        for chain_index in missing_bonds_chain_indices:
            chain = structure.chains[chain_index]
            chain_fragment = chain.get_selection()
            chain_fragment_name = structure.name_selection(chain_fragment)
            unbonded_fragments[chain_fragment_name] = chain_fragment

    # We also split PBC regions in chains
    # Otherwise they would have a lot of fragments (e.g. every water or lipid molecule)
    bonded_pbc_selection = pbc_selection - missing_bonds_selection
    if bonded_pbc_selection:
        # Parse the selection to chain indices
        pbc_chain_indices = structure.get_selection_chain_indices(bonded_pbc_selection)
        # Make sure the chain split selection selects whole chains, for later
        bonded_pbc_selection = structure.select_chain_indices(pbc_chain_indices)
        # Now split the definitive selection in chains
        for chain_index in pbc_chain_indices:
            chain = structure.chains[chain_index]
            chain_fragment = chain.get_selection()
            # Single-atom fragments are excluded
            # Note that single-atoms are not checked by atom distance pairwsie but normal RMSD
            # However if the atoms are under PBC they may jump freely, so we do not check them
            if len(chain_fragment) == 1: continue
            chain_fragment_name = structure.name_selection(chain_fragment)
            bonded_fragments[chain_fragment_name] = chain_fragment

    # Now split the rest of the structure in fragments
    chain_split_selection = missing_bonds_selection + bonded_pbc_selection
    fragment_split_selection = structure.invert_selection(chain_split_selection)
    if fragment_split_selection:
        for fragment in structure.find_fragments(fragment_split_selection):
            # Many fragments may cover the same part of a chain (e.g. lipids) so we add their first atom to the name
            # Otherwise they would share the name and overwrite each other
            first_atom = structure.atoms[fragment.atom_indices[0]]
            fragment_name = f'{structure.name_selection(fragment)} starting at {first_atom.label}'
            bonded_fragments[fragment_name] = fragment

    # Merge all fragments together
    all_fragments = bonded_fragments | unbonded_fragments

    # Meta-fragment: one representative atom (first) per real non-PBC fragment
    # Catches whole-fragment jumps that intra-fragment bond distances cannot see
    non_pbc_fragments = [ fragment for fragment in all_fragments.values() if not fragment.intersection(pbc_selection) ]
    if len(non_pbc_fragments) >= 2:
        header_atoms_indices = [ fragment.atom_indices[0] for fragment in non_pbc_fragments ]
        meta_fragment = structure.select_atom_indices(header_atoms_indices)
        unbonded_fragments[META_FRAGMENT_NAME] = meta_fragment
        all_fragments[META_FRAGMENT_NAME] = meta_fragment

    print(f'Checking trajectory integrity ({len(all_fragments)} fragments)')

    # Read the trajectory frame by frame with boxes as they are in the file
    # MDtraj would rebuild them not rotated, which is wrong for fitted trajectories
    frames = iterate_xtc_frames(input_trajectory_filename)

    # Save the previous frame every time
    previous_positions, previous_box = next(frames)

    # Build atom pairs per fragment to track via distance RMSD
    # Bonded fragments use the safe bonds
    fragment_bond_pairs = {}
    # We make a list from bonded_fragments.items since some fragments may be moved to unbonded while iterating
    for fragment_name, fragment in list(bonded_fragments.items()):
        # If it has one atom only then skip it
        if len(fragment) == 1:
            fragment_bond_pairs[fragment_name] = None
            continue 
        atom_set = set(fragment.atom_indices)
        # Now get bonded atom pairs inside the fragment
        pairs = []
        # Iterate fragment atom indices
        for atom_index in fragment.atom_indices:
            bonded_atom_indices = structure.bonds[atom_index]
            # Iterate bonded atom indices
            for bonded_atom_index in bonded_atom_indices:
                # In order to avoid having duplicates (a-b and b-a) we only record the lower-higher bond
                if bonded_atom_index < atom_index: continue
                # Skip atoms which do not belong to the fragment
                if bonded_atom_index not in atom_set: continue
                # Add this pair to the list
                pairs.append((atom_index, bonded_atom_index))
        # If there are no pairs then it means all atoms inside this fragment are not connected
        # This may happen, for instance, with a chain made of counter ions, which are under PBC
        # If this happens then move the fragment to the unbonded fragments dict
        if not pairs:
            unbonded_fragments[fragment_name] = fragment
            del bonded_fragments[fragment_name]
            continue
        # Assign these pairs to this specific fragment
        fragment_bond_pairs[fragment_name] = np.array(pairs)

    # Unbonded fragments (with missing bonds) use sequential pairs as a fake backbone
    for fragment_name, fragment in unbonded_fragments.items():
        # If it has one atom only then skip it
        if len(fragment) == 1:
            fragment_bond_pairs[fragment_name] = None
            continue
        indices = fragment.atom_indices
        fragment_bond_pairs[fragment_name] = np.array([(indices[i], indices[i + 1]) for i in range(len(indices) - 1)])

    # Join the pairs of all fragments so distances are computed with a single call per frame
    # Then every fragment takes its own slice of distances
    paired_fragment_names = [ name for name, pairs in fragment_bond_pairs.items() if pairs is not None ]
    single_atom_fragment_names = [ name for name, pairs in fragment_bond_pairs.items() if pairs is None ]
    pairs_per_fragment = [ len(fragment_bond_pairs[name]) for name in paired_fragment_names ]
    fragment_starts = np.cumsum([0] + pairs_per_fragment[:-1])
    all_pairs = np.concatenate([ fragment_bond_pairs[name] for name in paired_fragment_names ])

    # Compute distances for the first frame so the loop always has a previous reference
    previous_distances = compute_periodic_distances(previous_positions, previous_box, all_pairs)

    # Save all RMSD jumps
    # Jumps of fragments with pairs are saved as one array per frame, with one value per fragment
    paired_rmsd_jumps = []
    single_atom_rmsd_jumps = { fragment_name: [] for fragment_name in single_atom_fragment_names }

    # Initialize progress bar first
    pbar = tqdm(frames, total=snapshots, desc=' Frame', unit='frame', initial=1)

    # Iterate trajectory frames
    for positions, box_vectors in pbar:

        # Compute bond distance RMSD between consecutive frames for all fragments with pairs at once
        current_distances = compute_periodic_distances(positions, box_vectors, all_pairs)
        squared_differences = (current_distances - previous_distances) ** 2
        squared_sums = np.add.reduceat(squared_differences, fragment_starts)
        paired_rmsd_jumps.append(np.sqrt(squared_sums / pairs_per_fragment))
        previous_distances = current_distances

        # If the fragment has a single atom then its RMSD is just how much it moved
        # Boxes are not used here, so jumps across the boundaries are also caught
        for fragment_name in single_atom_fragment_names:
            atom_index = all_fragments[fragment_name].atom_indices[0]
            rmsd_value = np.linalg.norm(positions[atom_index] - previous_positions[atom_index])
            single_atom_rmsd_jumps[fragment_name].append(rmsd_value)

        # Update the previous frame as the current one
        previous_positions = positions

    # Gather the RMSD jumps of every fragment
    paired_rmsd_jumps = np.array(paired_rmsd_jumps)
    fragment_rmsd_jumps = { name: paired_rmsd_jumps[:, i] for i, name in enumerate(paired_fragment_names) }
    fragment_rmsd_jumps.update(single_atom_rmsd_jumps)

    # Check all fragments and store their test result here before issuing the final report
    fragment_reports = []

    # Count failed fragments, since only the first ones are shown in the logs
    failed_fragments = 0

    # Iterate over the different fragments
    for fragment_name, fragment in all_fragments.items():
        show_logs = failed_fragments < SHOWN_FAILED_FRAGMENTS

        # get the results
        rmsd_jumps = fragment_rmsd_jumps[fragment_name]

        # Skip fragments for which no metric could be computed (single atoms with no pairs)
        # A multi-atom fragment with no intra-fragment bonds means all its bonds point outside
        # the fragment — likely a cofactor or ligand fully bonded to another fragment in the topology
        if len(rmsd_jumps) == 0:
            raise RuntimeError(f'Fragment "{fragment_name}" could not be analyzed')

        # Capture outliers
        # If we capture more than 5 we stop searching
        outliers_count = 0
        max_z_score = 0
        max_z_score_frame = 0

        # Get the maximum RMSD value and check it is a reasonable deviation from the average values
        # Otherwise, if it is an outlier, the test fails
        mean_rmsd_jump = np.mean(rmsd_jumps)
        stdv_rmsd_jump = np.std(rmsd_jumps)

        # First frames may not be perfectly equilibrated and thus have stronger RMSD jumps
        # For this reason we allow the first frames to bypass the check
        # As soon as one frame is below the cutoff the bypass is finished for the following frames
        # Count the number of bypassed frames for this fragment
        bypassed_frames = 0

        for i, rmsd_jump in enumerate(rmsd_jumps, 1):
            z_score = abs( (rmsd_jump - mean_rmsd_jump) / stdv_rmsd_jump )
            # Keep track of the maximum z score
            if z_score > max_z_score:
                max_z_score = z_score
                max_z_score_frame = i
            # If z score bypassed the limit then report it
            if z_score > standard_deviations_cutoff:
                # If there are as many bypassed frames as the index then it means no frame has passed the cutoff yet
                if i - 1 == bypassed_frames:
                    bypassed_frames += 1
                    continue
                # If we are no long bypassing frist frames then outliers mean the test has failed
                # But don't kill the check just yet
                if show_logs and outliers_count < 4:
                    print(f' FAIL: Sudden RMSD jump in fragment "{fragment_name}" between frames {i} and {i+1} (1-based numeration) with a z-score of {round_to_hundredths(z_score)}')
                if show_logs and outliers_count == 4:
                    print(' etc...')
                outliers_count += 1

        # A jump which breaks the fragment for good leaves big jumps in all the following frames
        # These jumps make the standard deviation bigger and then the first jump does not stand out in the z score
        # However, the RMSD jumps are then split in two groups: before and after the break
        # The same happens when a region jumps and comes back repeatedly along the trajectory
        # So we also check if two gaussians fit the RMSD jumps much better than one
        # Jumps are divided by their median so the fit does not depend on their scale
        scaled_jumps = (rmsd_jumps / np.median(rmsd_jumps)).reshape(-1, 1)
        one_gaussian = GaussianMixture(1, random_state=0).fit(scaled_jumps)
        two_gaussians = GaussianMixture(2, random_state=0).fit(scaled_jumps)
        bimodality_score = one_gaussian.bic(scaled_jumps) - two_gaussians.bic(scaled_jumps)
        # Two gaussians may also fit better when there is a single group of jumps which is skewed
        # e.g. a rigid molecule with a single flexible bond has a minimum jump and then a tail of bigger jumps
        # So the groups must be also far from each other: bigger jumps must be several times bigger on average
        group_means = two_gaussians.means_[:, 0]
        separation = group_means.max() / group_means.min()
        bimodal_frames = 0
        if bimodality_score > BIMODALITY_BIC_CUTOFF and separation > BIMODALITY_SEPARATION_CUTOFF:
            # Find the frames in the group of bigger jumps
            bigger_group = np.argmax(group_means)
            bigger_group_frames = np.where(two_gaussians.predict(scaled_jumps) == bigger_group)[0] + 1
            # If bigger jumps are all in the first frames then these frames are not equilibrated
            if bigger_group_frames[-1] == len(bigger_group_frames):
                bypassed_frames = max(bypassed_frames, len(bigger_group_frames))
            else:
                bimodal_frames = len(bigger_group_frames)
                first_frame = bigger_group_frames[0]
                if show_logs:
                    print(f' FAIL: RMSD jumps in fragment "{fragment_name}" are split in two groups'
                        f' (bimodality score = {round_to_hundredths(bimodality_score)}, separation = {round_to_hundredths(separation)}).'
                        f' Bigger jumps start between frames {first_frame} and {first_frame + 1}')

        if outliers_count > 0 or bimodal_frames > 0: failed_fragments += 1

        # Save the report for this fragment
        fragment_reports.append({
            'name': fragment_name,
            'score': max_z_score,
            'frame': max_z_score_frame,
            'outliers': outliers_count,
            'bypass': bypassed_frames,
            'bimodality': bimodality_score,
            'separation': separation,
            'bimodal_frames': bimodal_frames,
            'jumps' : rmsd_jumps
        })

    hidden_fragments = failed_fragments - SHOWN_FAILED_FRAGMENTS
    if hidden_fragments > 0: print(f' (and {hidden_fragments} fragments more)')

    # If there was any outliers or any bypassed frames then the test has failed
    # LORE: Back in the day having bypassed frames was not a test failure, but now it is
    any_outliers = any(report['outliers'] > 0 or report['bimodal_frames'] > 0 for report in fragment_reports)
    any_bypassed_frames = any(report['bypass'] > 0 for report in fragment_reports)

    # If the test has failed then display a full report
    # Only failed fragments are shown since there may be thousands of fragments
    if any_outliers or any_bypassed_frames:
        print('-- RMSD Check final report --')
        passed_count = 0
        shown_count = 0
        for report in fragment_reports:
            name = report['name']
            max_z_score = round_to_hundredths(report['score'])
            max_z_score_frame = report['frame']
            outliers = report['outliers']
            bypassed_frames = report['bypass']
            bimodal_frames = report['bimodal_frames']
            passed = outliers == 0 and bypassed_frames == 0 and bimodal_frames == 0
            # If the fragment has no problems
            if passed:
                passed_count += 1
                continue
            # If the fragment failed to pass the test
            if shown_count == SHOWN_FAILED_FRAGMENTS: continue
            shown_count += 1
            fail(f'Fragment "{name}" FAILED with a maximum z-score of {max_z_score}'
                 f' reported between frames {max_z_score_frame} and {max_z_score_frame + 1}'
                 f', a bimodality score of {round_to_hundredths(report["bimodality"])}'
                 f' and a separation of {round_to_hundredths(report["separation"])}{COLOR_END}')
            if outliers > 0: print(f' Outliers: {outliers}')
            if bimodal_frames > 0: print(f' Frames in the group of bigger jumps: {bimodal_frames}')
            if bypassed_frames > 0: print(f' Bypassed frames: {bypassed_frames}')
        hidden_fragments = len(fragment_reports) - passed_count - shown_count
        if hidden_fragments > 0: print(f'(and {hidden_fragments} failed fragments more)')
        print(f'{passed_count} other fragments PASSED')
        print('*The z-score of a value means how many times the standard deviation away it is from the average')
        print('*The bimodality score means how much better two gaussians fit the jumps than one (BIC difference)')
        print('*The separation means how many times bigger the jumps of the bigger gaussian are on average')
        print('*The bimodality only fails when both the bimodality score and the separation are above their cutoffs')

        # Find the worst fragment
        # Show the most bimodal fragment if any, since its jumps may not stand out in the z score
        bimodal_reports = [ report for report in fragment_reports if report['bimodal_frames'] > 0 ]
        if bimodal_reports:
            worst_report = max(bimodal_reports, key=lambda report: report['bimodality'])
        elif any_outliers:
            worst_report = max(fragment_reports, key=lambda report: report['score'])
        else:
            worst_report = max(fragment_reports, key=lambda report: report['bypass'])

        # Write a VMD script for the user to see the biggest jump in the worst fragment
        # This is not part of the warning since it is only useful for the workflow user
        worst_fragment = all_fragments[worst_report['name']]
        # Frame after the biggest jump, starting at 0
        jump_frame = int(np.argmax(worst_report['jumps'])) + 1
        # Show the two frames before the jump and the two frames after it, as long as they exist
        last_frame = len(worst_report['jumps'])
        shown_frames = list(range(max(jump_frame - 2, 0), min(jump_frame + 1, last_frame) + 1))
        jumping_atoms = find_jumping_atoms(input_trajectory_filename,
            worst_fragment, fragment_bond_pairs[worst_report['name']], jump_frame)
        is_meta_fragment = worst_report['name'] == META_FRAGMENT_NAME
        write_view_script(view_script_filename, input_structure_filename, input_trajectory_filename,
            jump_frame, shown_frames, worst_fragment, jumping_atoms, is_meta_fragment)
        view_hint = f'To see the clearest jump run:\n  vmd -e {view_script_filename}'

        # Describe the worst fragment and its biggest jump in the warning
        failed_count = len(fragment_reports) - passed_count
        description = (f'{failed_count} of {len(fragment_reports)} fragments failed. The worst one is "{worst_report["name"]}"'
            f', with a maximum z-score of {round_to_hundredths(worst_report["score"])} (cutoff {standard_deviations_cutoff})'
            f', a bimodality score of {round_to_hundredths(worst_report["bimodality"])} (cutoff {BIMODALITY_BIC_CUTOFF})'
            f' and a separation of {round_to_hundredths(worst_report["separation"])} (cutoff {BIMODALITY_SEPARATION_CUTOFF})')
        if worst_report['bypass'] > 0:
            description += f', and {worst_report["bypass"]} bypassed frames at the beginning'
        description += f'.\nIts biggest jump is between frames {jump_frame - 1} and {jump_frame}.'
        if is_meta_fragment: description += f'\n{META_FRAGMENT_EXPLANATION}'
        # Frames and atoms to be shown in the client
        data = {
            'frames': shown_frames,
            'selection': format_selection(structure, worst_fragment),
        }

    # If there were any outlier then the check has failed
    if any_outliers > 0:
        # Add a warning an return True since the test failed in case we have mercy
        message = 'RMSD check has failed: there may be sudden jumps along the trajectory'
        if TRAJECTORY_INTEGRITY_FLAG in mercy:
            register.add_warning(TRAJECTORY_INTEGRITY_FLAG, message, description, data)
            print(f' {view_hint}')
            register.update_test(TRAJECTORY_INTEGRITY_FLAG, False)
            return False
        # Otherwise kill the process right away, after displaying a graph
        title = f'RMSD jumps along the trajectory in fragment "{worst_report["name"]}"'
        report_data = worst_report['jumps']
        display_rmsd_jumps_graph(report_data, title)
        raise TestFailure(f'{message}\n{description}\n{view_hint}')

    # Warn the user if we had bypassed frames
    if any_bypassed_frames > 0:
        # Set the error message, which is to be shown in the web client when mercy is asked
        max_bypassed_frames = max(report['bypass'] for report in fragment_reports)
        message = f'First {max_bypassed_frames} frames may be not equilibrated.'
        # Add a warning an return True since the test failed in case we have mercy
        if TRAJECTORY_INTEGRITY_FLAG in mercy:
            register.add_warning(TRAJECTORY_INTEGRITY_FLAG, message, description, data)
            print(f' {view_hint}')
            register.update_test(TRAJECTORY_INTEGRITY_FLAG, False)
            return False
        # Otherwise kill the process, after displaying a graph
        graph_frames_limit = max_bypassed_frames + 100
        report_data = worst_report['jumps'][0:graph_frames_limit]
        title = f'RMSD jumps along the trajectory in fragment "{worst_report["name"]}" (first {graph_frames_limit} frames)'
        display_rmsd_jumps_graph(report_data, title)
        # Note that frames already skipped must be added to the suggested number of frames to skip
        message += f' Use the "--mercy {TRAJECTORY_INTEGRITY_FLAG}" flag to skip this error' \
            f' or exclude these frames from the trajectory with "-skip {md_index} {skipped_frames + max_bypassed_frames}"' \
            ' (or the "skip_frames" input).'
        raise TestFailure(f'{message}\n{description}\n{view_hint}')

    print(' Test has passed successfully')
    register.update_test(TRAJECTORY_INTEGRITY_FLAG, True)
    return True


def compute_periodic_distances(
    positions : np.ndarray,
    box_vectors : np.ndarray,
    pairs : np.ndarray,
) -> np.ndarray:
    """Get the distance between the atoms of every pair, using the closest periodic image.
    Boxes may be rotated (e.g. in fitted trajectories), which MDtraj does not support.
    For this reason the whole frame is rotated back before calculating distances.
    A box full of zeros means there are no periodic boundaries.
    """
    has_box = box_vectors.any()
    if has_box: positions, box_vectors = unrotate_box(positions, box_vectors)
    frame = mdt.Trajectory(positions[np.newaxis], None)
    if has_box: frame.unitcell_vectors = box_vectors[np.newaxis]
    return mdt.compute_distances(frame, pairs, periodic=True)[0]


def find_jumping_atoms(
    trajectory_filename : str,
    fragment : 'Selection',
    pairs : Optional[np.ndarray],
    jump_frame : int,
) -> 'Selection':
    """Find the atoms of the bonds which change the most in a jump, from the previous frame to the jump frame.
    Fragments with a single atom have no bonds, so the atom itself is returned.
    """
    if pairs is None: return fragment
    with mdt.formats.XTCTrajectoryFile(trajectory_filename) as reader:
        reader.seek(jump_frame - 1)
        xyz, _, _, boxes = reader.read(n_frames=2)
    # Trajectories with no box have no periodic boundaries, as a box full of zeros
    if boxes is None: boxes = np.zeros((2, 3, 3), dtype=np.float32)
    previous_distances = compute_periodic_distances(xyz[0], boxes[0], pairs)
    distances = compute_periodic_distances(xyz[1], boxes[1], pairs)
    changes = np.abs(distances - previous_distances)
    biggest_changes = np.argsort(-changes)[:HIGHLIGHTED_BONDS]
    # Keep only bonds which change at least half as much as the bond which changes the most
    # Otherwise a single broken bond would be shown along with bonds which barely change
    biggest_changes = biggest_changes[changes[biggest_changes] >= changes[biggest_changes[0]] / 2]
    return Selection([ int(index) for index in np.unique(pairs[biggest_changes]) ])


def write_view_script(
    script_filename : str,
    structure_filename : str,
    trajectory_filename : str,
    jump_frame : int,
    shown_frames : list[int],
    fragment : 'Selection',
    jumping_atoms : 'Selection',
    is_meta_fragment : bool,
):
    """Write a VMD script to see a jump, from the previous frame to the jump frame.

    Only the shown frames are loaded.
    The fragment is shown in yellow licorice and the atoms of the bonds which change the most in CPK.
    """
    previous_frame = jump_frame - 1
    commands = [
        f'# Biggest RMSD jump, between frames {previous_frame} and {jump_frame}',
        f'mol new {{{abspath(structure_filename)}}} type pdb waitfor all',
        # Frames from the structure file are removed so frame numbers match the trajectory
        'animate delete all',
        f'mol addfile {{{abspath(trajectory_filename)}}} type xtc first {shown_frames[0]} last {shown_frames[-1]} waitfor all',
        'set main [molinfo top]',
        # Fragment
        'mol representation Licorice',
        'mol color ColorID 4',
        f'mol selection {{{fragment.to_vmd()}}}',
        'mol material Opaque',
        'mol addrep $main',
        # Atoms of the bonds which change the most
        'mol representation CPK',
        'mol color Name',
        f'mol selection {{{jumping_atoms.to_vmd()}}}',
        'mol addrep $main',
        # Center the view in the fragment before the jump and zoom in
        f'animate goto {shown_frames.index(previous_frame)}',
        f'set fragment [atomselect $main {{{fragment.to_vmd()}}}]',
        'proc get_size {atoms} {',
        '    lassign [measure minmax $atoms] lowest highest',
        '    return [tcl::mathfunc::max {*}[vecsub $highest $lowest]]',
        '}',
        'display resetview',
        'molinfo $main set center_matrix [list [transoffset [vecscale -1 [measure center $fragment]]]]',
        f'set zoom [expr {{[get_size [atomselect $main all]] / max([get_size $fragment], {MIN_FOCUS_SIZE})}}]',
        'if {$zoom > 1} { scale by $zoom }',
        f'puts "Frames {shown_frames[0]} to {shown_frames[-1]} are loaded as frames 0 to {len(shown_frames) - 1}"',
        f'puts "The jump is between frames {shown_frames.index(previous_frame)} and {shown_frames.index(jump_frame)}"',
        'puts "The fragment is shown in yellow licorice and the atoms of the bonds which change the most in CPK"',
    ]
    if is_meta_fragment:
        commands.append(f'puts {{{META_FRAGMENT_EXPLANATION}}}')
    with open(script_filename, 'w') as script_file:
        script_file.write('\n'.join(commands) + '\n')


def display_rmsd_jumps_graph (data : list, title : str):
    """Display a graph to show the distribution of sudden jumps along the trajectory in the terminal itself."""
    if not is_terminal: return
    # Display a graph to show the distribution of sudden jumps along the trajectory
    plt.scatter(data)
    plt.title(title)
    plt.xlabel('Frame')
    plt.ylabel('RMSD jump')
    n_ticks = 5
    n_jumps = len(data)
    tickstep = math.ceil(n_jumps / n_ticks)
    xticks = [ t+1 for t in range(0, n_jumps, tickstep) ]
    xlabels = [ str(t) for t in xticks ]
    plt.xticks(xticks, xlabels)
    plt.show()

# Compute every residue RMSD to check if there are sudden jumps along the trajectory
# HARDCODE: This function is not fully implemented but enabled manually for specific cases
def check_trajectory_integrity_per_residue (
    input_structure_filename : str,
    input_trajectory_filename : str,
    structure : 'Structure',
    pbc_selection : 'Selection',
    mercy : list[str],
    trust: list[str],
    register : 'Register',
    #time_length : float,
    check_selection : str,
    # DANI: He visto saltos 'correctos' pasar de 11
    # DANI: He visto saltos 'incorrectos' no bajar de 14
    standard_deviations_cutoff : float):

    # HARDCODE: The default value does not work for a single residue
    standard_deviations_cutoff = 12

    # Skip the test if we trust
    if TRAJECTORY_INTEGRITY_FLAG in trust:
        return True

    # Skip the test if it is already passed according to the register
    if register.tests.get(TRAJECTORY_INTEGRITY_FLAG, None):
        return True

    # Remove old warnings
    register.remove_warnings(TRAJECTORY_INTEGRITY_FLAG)

    # Parse the selection in VMD selection syntax
    parsed_selection = structure.select(check_selection, syntax='vmd')

    # If there is nothing to check then warn the user and stop here
    if not parsed_selection:
        raise Exception('WARNING: There are not atoms to be analyzed for the RMSD analysis')

    # Discard PBC residues from the selection to be checked
    parsed_selection -= pbc_selection

    # If there is nothing to check then warn the user and stop here
    if not parsed_selection:
        warn('There are no atoms to be analyzed for the RMSD checking after PBC substraction')
        register.update_test(TRAJECTORY_INTEGRITY_FLAG, 'na')
        return True

    # We must filter out residues which only have 1 atom (e.g. ions)
    # This is because sometimes pytraj does not return results for them and then the number of results and residues does not match
    # More info: https://github.com/Amber-MD/pytraj/issues/1580
    ion_atom_indices = []
    for residue in structure.residues:
        if len(residue.atom_indices) == 1:
            ion_atom_indices += residue.atom_indices
    ions_selection = structure.select_atom_indices(ion_atom_indices)
    parsed_selection -= ions_selection

    # Parse the selection to a pytraj mask
    pytraj_selection = parsed_selection.to_pytraj()

    # Calculate the residue indices of the overall structure remaining in the filtered trajectory
    residue_indices = structure.get_selection_residue_indices(parsed_selection)
    n_residues = len(residue_indices)

    print('Checking trajectory integrity per residue')

    # Parse the trajectory into pytraj and apply the mask
    # NEVER FORGET: The pytraj iterload does not accept a mask, but we apply the mask later in the analysis
    pt_trajectory = get_pytraj_trajectory(input_structure_filename, input_trajectory_filename, atom_selection = parsed_selection)

    # Make sure the expected output number of residues to match with the pytraj results
    # These numbers may not match when ions are included so we better check
    # NEVER FORGET: The pytraj TrajectoryIterator is not an iterator
    first_frame = pt_trajectory[0:1]
    # DANI: When the 'resname' argument is missing it prints "Error: Range::SetRange(None): Range is -1 for None"
    # DANI: However there is no problem and the analysis runs flawlessly
    # DANI: For this reason we call this function with no resname and then we remove the log
    data_sample = pt.rmsd_perres(first_frame, ref=first_frame, perres_mask=pytraj_selection)
    # We remove the previous error log
    delete_previous_log()
    # We remove the first result, which is meant to be the whole rmsd and whose key is 'RMSD_00001'
    del data_sample[0]
    if n_residues != len(data_sample):
        raise ValueError(f'Number of target residues ({n_residues}) does not match number of residues in data ({len(data_sample)})')

    # Saving all RMSD jumps may take a lot of memory
    # Instead we will store the sum of values and the maximum
    # This way we can caluclate the average value at the end and check if the maximum is too far from it
    rmsd_per_residue_per_frame = []

    # Add an extra breakline before the first log
    print()

    # Iterate trajectory frames
    previous_frame_trajectory = first_frame
    frame_number = 1
    for frame in pt_trajectory:
        # Update the current frame log
        reprint(f' Frame {frame_number}')
        # Set a pytraj trajectory out of a pytraj frame
        frame_trajectory = pt.Trajectory(top=pt_trajectory.top)
        frame_trajectory.append(frame)
        # Run the analysis in pytraj
        # The result data is a custom pytraj class: pytraj.datasets.datasetlist.DatasetList
        # This class has keys but its attributes can not be accessed through the key
        # They must be accessed thorugh the index
        # DANI: When the 'resname' argument is missing it prints "Error: Range::SetRange(None): Range is -1 for None"
        # DANI: However there is no problem and the analysis runs flawlessly
        # DANI: Adding resrage as a list/range was tried and did not work, only string works
        # DANI: Adding a string resrange however strongly impacts the speed when this function is called repeatedly
        # DANI: For this reason we call this function with no resname and then we remove the log
        rmsd_per_residue = pt.rmsd_perres(frame_trajectory, ref=previous_frame_trajectory, mask=pytraj_selection)
        # We remove the previous error log
        delete_previous_log()
        # We remove the first result, which is meant to be the whole rmsd and whose key is 'RMSD_00001'
        del rmsd_per_residue[0]
        # Check we have no NaNs
        if np.isnan(rmsd_per_residue[0][0]):
            raise ValueError(f'We are having NaNs at frame {frame_number}')
        # Add last values to the list
        rmsd_per_residue_per_frame.append(rmsd_per_residue)
        # rmsd_per_residue_per_frame[:, frame]
        # Now update data for every residue
        # for index, residue_rmsd in enumerate(rmsd_per_residue):
        #     current_rmsd = residue_rmsd[0]
        #     # Get the current residue rmsd data
        #     total_rmsd = rmsd_per_residue_per_frame[index]
        #     total_rmsd['accumulated'] += current_rmsd
        #     total_rmsd['maximum'] = max(total_rmsd['maximum'], current_rmsd)
        # Update previous coordinates
        previous_frame_trajectory = frame_trajectory
        # Update the frame_number
        frame_number += 1

    # If the trajectory has only 1 or 2 frames then there is no test to do
    n_jumps = len(rmsd_per_residue_per_frame)
    if n_jumps <= 1:
        register.update_test(TRAJECTORY_INTEGRITY_FLAG, True)
        return True

    # Keep the overall maximum z score, and its residue and frame for the logs
    overall_max_z_score = 0
    overall_max_z_score_frame = None
    overall_max_z_score_residue = None

    # Keep the overall maximum bypassed frames number
    overall_bypassed_frames = 0

    # Keep the overall count of residues with outliers
    overall_outliered_residues = 0

    # Add an extra breakline before the next log
    print()

    # Now check there are not sudden jumps for each residue separattely
    for residue_number in range(n_residues):
        reprint(f' Residue {residue_number+1}')
        # Get the rmsd jumps for each frame for this specific residue
        rmsd_jumps = [ frame[residue_number] for frame in rmsd_per_residue_per_frame ]

        # Get the maximum RMSD value and check it is a reasonable deviation from the average values
        # Otherwise, if it is an outlier, the test fails
        mean_rmsd_jump = np.mean(rmsd_jumps)
        stdv_rmsd_jump = np.std(rmsd_jumps)

        # First frames may not be perfectly equilibrated and thus have stronger RMSD jumps
        # For this reason we allow the first frames to bypass the check
        # As soon as one frame is below the cutoff the bypass is finished for the following frames
        # Count the number of bypassed frames and warn the user in case there are any
        bypassed_frames = 0

        # Capture outliers
        # If we capture more than 5 we stop searching
        outliers_count = 0
        max_z_score = 0
        max_z_score_frame = 0
        for i, rmsd_jump in enumerate(rmsd_jumps, 1):
            z_score = abs( (rmsd_jump - mean_rmsd_jump) / stdv_rmsd_jump )
            # Keep track of the maixmum z score
            if z_score > max_z_score:
                max_z_score = z_score
                max_z_score_frame = i
            # If z score bypassed the limit then report it
            if z_score > standard_deviations_cutoff:
                # If there are as many bypassed frames as the index then it means no frame has passed the cutoff yet
                if i - 1 == bypassed_frames:
                    bypassed_frames += 1
                    continue
                # Otherwise we consider this as an outlier and thus the test has failed
                # However we keep checking just to find and report the highest outlier
                outliers_count += 1

        # Update the overall bypassed frames if we overcomed it
        if overall_bypassed_frames < bypassed_frames:
            overall_bypassed_frames = bypassed_frames

        # Update the overall max z score if we overcome it
        if max_z_score > overall_max_z_score:
            overall_max_z_score = max_z_score
            overall_max_z_score_frame = max_z_score_frame
            overall_max_z_score_residue = residue_number

        # If there were any outlier then add one to the overall count
        if outliers_count > 0:
            overall_outliered_residues += 1

    # Always print the overall maximum z score and its frames and residue
    overall_max_z_score_residue_label = pt_trajectory.top.residue(overall_max_z_score_residue)
    print(f' Maximum z score {overall_max_z_score} reported for residue {overall_max_z_score_residue_label} between frames {overall_max_z_score_frame} and {overall_max_z_score_frame + 1}')

    # If there were any outlier then the check has failed
    if overall_outliered_residues > 0:
        # Add a warning an return True since the test failed in case we have mercy
        message = 'RMSD check has failed: there may be sudden jumps along the trajectory'
        if TRAJECTORY_INTEGRITY_FLAG in mercy:
            register.add_warning(TRAJECTORY_INTEGRITY_FLAG, message)
            register.update_test(TRAJECTORY_INTEGRITY_FLAG, False)
            return False
        # Otherwise kill the process right away
        raise TestFailure(message)

    # Warn the user if we had bypassed frames
    if overall_bypassed_frames > 0:
        register.add_warning(TRAJECTORY_INTEGRITY_FLAG, f'First {overall_bypassed_frames} frames may be not equilibrated')

    print(' Test has passed successfully')
    register.update_test(TRAJECTORY_INTEGRITY_FLAG, True)
    return True
