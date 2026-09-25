from os.path import exists
import numpy as np
import mdtraj as mdt
from sklearn.decomposition import PCA
from sklearn.cluster import HDBSCAN
from scipy.linalg import eigh

from mddb_workflow.utils.auxiliar import save_json
from mddb_workflow.utils.auxiliar import numerate_filename, get_analysis_name
from mddb_workflow.utils.auxiliar import reprint
from mddb_workflow.utils.constants import OUTPUT_CLUSTERS_FILENAME, OUTPUT_CLUSTER_SCREENSHOT_FILENAMES
from mddb_workflow.utils.file import File
from mddb_workflow.tools.get_screenshot import get_screenshot
from mddb_workflow.tools.get_reduced_trajectory import get_reduced_trajectory
from mddb_workflow.utils.type_hints import *


def clusters_analysis(
    structure_file: 'File',
    trajectory_file: 'File',
    interactions: list,
    structure: 'Structure',
    pbc_selection: 'Selection',
    cg_selection : 'Selection',
    output_directory: str,
    # False for ensembles (frames without temporal order), where TICA is not applicable
    is_time_dependent: bool,
    # Time between frames in ns; without it frames are not known to be evenly spaced in time
    input_framestep: Optional[float],
    # Atom selection for the overall clustering
    overall_selection: str = "name CA or name C5'",
    # PCA: fraction of variance to retain (e.g. 0.95 keeps enough components for 95%)
    pca_variance_threshold: float = 0.95,
    # Minimum frames a conformational state must occupy to be considered a real cluster.
    # None = 1% of sampled frames (max(3, n_frames // 100)).
    min_cluster_size: int = None,
    # How many close neighbors a frame needs to count as "core" density, rather than a thin bridge.
    # Low values (e.g. 1) make HDBSCAN prone to chaining: elongated clusters linked by a path of
    # nearby frames even when their extremes are structurally very different.
    # None = half of min_cluster_size (at least 3).
    min_samples: int = None,
    # TICA lag time, in ns: the minimum time something must last to count as a conformational state.
    # It is converted to frames with the framestep, and it is never shorter than one frame.
    tica_lag_time: float = 1,
    # The trajectory must last at least this many times the slowest TICA relaxation time,
    # otherwise transitions are not observed repeatedly and TICA is not reliable.
    tica_min_time_ratio: float = 10,
    # TICA: fraction of kinetic variance to retain
    tica_kinetic_variance_threshold: float = 0.95,
    # Maximum number of PCA/TICA components passed to the clustering.
    # Density-based clustering degrades quickly with dimensionality: points become almost
    # equidistant and no dense region stands out. It is also much slower.
    max_clustering_components: int = 10,
):
    """Run the cluster analysis using HDBSCAN on Cartesian coordinates reduced with PCA,
    followed by TICA when frames are a time series."""
    runs = []

    parsed_overall_selection = structure.select(overall_selection)
    if not parsed_overall_selection:
        parsed_overall_selection = structure.select_heavy_atoms()
    parsed_overall_selection -= pbc_selection
    if parsed_overall_selection:
        runs.append({'name': 'Overall', 'selection': parsed_overall_selection})

    for interaction in interactions:
        interface_residue_indices = interaction['interface_residue_indices_1'] \
            + interaction['interface_residue_indices_2']
        interface_selection = structure.select_residue_indices(interface_residue_indices)
        heavy_atoms_selection = structure.select_heavy_atoms()
        final_selection = interface_selection & heavy_atoms_selection
        final_selection -= pbc_selection
        if final_selection:
            runs.append({'name': interaction['name'], 'selection': final_selection})

    if len(runs) == 0:
        print(' No clusters to analyze')
        return

    output_analysis_filepath = f'{output_directory}/{OUTPUT_CLUSTERS_FILENAME}'
    output_screenshot_filepath = f'{output_directory}/{OUTPUT_CLUSTER_SCREENSHOT_FILENAMES}'

    traj = mdt.load(trajectory_file.path, top=structure_file.path)

    auxiliar_structure = structure.copy()
    output_summary = []

    for r, run in enumerate(runs):
        numbered_output_analysis_filepath = numerate_filename(output_analysis_filepath, r)
        name = run['name']
        analysis_name = get_analysis_name(numbered_output_analysis_filepath)
        output_summary.append({'name': name, 'analysis': analysis_name})

        if exists(numbered_output_analysis_filepath):
            continue

        atom_indices = run['selection'].atom_indices
        n_atoms = len(atom_indices)
        print(f'Extracting features for {name} -> {analysis_name} ({n_atoms} atoms, {traj.n_frames} frames)')

        # Superpose all frames to frame 0 on the selected atoms to remove rigid-body motion
        traj.superpose(traj, 0, atom_indices=atom_indices)

        # Build feature matrix: (n_frames, n_atoms * 3) Cartesian coordinates
        features = traj.xyz[:, atom_indices, :].reshape(traj.n_frames, -1)

        # Reduce dimensionality with PCA
        # Coordinates are not rescaled: they are all displacements in the same units, so atoms that
        # move more must weigh more, instead of amplifying the fluctuations of barely moving atoms
        max_components = min(features.shape[0] - 1, features.shape[1])
        pca = PCA(n_components=pca_variance_threshold, svd_solver='full')
        try:
            features_pca = pca.fit_transform(features)
        except ValueError:
            pca = PCA(n_components=min(50, max_components), svd_solver='full')
            features_pca = pca.fit_transform(features)

        n_components_kept = features_pca.shape[1]
        variance_explained = float(pca.explained_variance_ratio_.sum())
        print(f' PCA: kept {n_components_kept} components ({variance_explained*100:.1f}% variance)')

        # TICA needs frames that are a time series evenly spaced in time, and long enough to observe
        # repeated transitions between states. Otherwise fall back to PCA.
        tica_lag = None
        tica_skip_reason = None
        if not is_time_dependent:
            tica_skip_reason = 'frames are not a time series (ensemble)'
        elif not input_framestep:
            tica_skip_reason = 'the time between frames (framestep) is unknown'
        else:
            tica_lag = max(1, round(tica_lag_time / input_framestep))
            if traj.n_frames <= tica_lag:
                tica_skip_reason = f'not enough simulated time: the trajectory is shorter than the TICA lag ({tica_lag_time:g} ns)'
        if tica_skip_reason is None:
            tica_coordinates, tica_eigenvalues = tica(features_pca, tica_lag)
            slowest_relaxation_time = relaxation_time(tica_eigenvalues[0], tica_lag)
            if traj.n_frames < tica_min_time_ratio * slowest_relaxation_time:
                tica_skip_reason = (f'not enough simulated time: the slowest motion takes'
                    f' ~{slowest_relaxation_time * input_framestep:.3g} ns to relax, which is too long compared'
                    f' to the trajectory length ({traj.n_frames * input_framestep:.3g} ns) to observe repeated transitions')
        use_tica = tica_skip_reason is None
        if use_tica:
            kinetic_variance = tica_eigenvalues ** 2 / np.sum(tica_eigenvalues ** 2)
            cumulative_kinetic_variance = np.cumsum(kinetic_variance)
            n_tica_components = min(len(kinetic_variance), max_clustering_components,
                int(np.searchsorted(cumulative_kinetic_variance, tica_kinetic_variance_threshold)) + 1)
            # TICA components come out with equal scale, so minor fast modes would weigh as much
            # as the slowest one and dilute the density; weight each by its share of kinetic variance
            clustering_features = tica_coordinates[:, :n_tica_components] * kinetic_variance[:n_tica_components]
            kept_kinetic_variance = float(cumulative_kinetic_variance[n_tica_components - 1])
            print(f' TICA (lag {tica_lag} frames = {tica_lag * input_framestep:g} ns): kept {n_tica_components} components'
                f' ({kept_kinetic_variance*100:.1f}% kinetic variance)')
        else:
            clustering_features = features_pca[:, :max_clustering_components]
            print(f' TICA not applied: {tica_skip_reason}')
            print(f' Using {clustering_features.shape[1]} PCA components for clustering')

        # Cluster with HDBSCAN — the algorithm determines the number of clusters.
        # min_cluster_size is the minimum density threshold: a state needs this many frames
        # to be considered real rather than transient noise.
        effective_min_cluster_size = min_cluster_size if min_cluster_size is not None \
            else max(3, traj.n_frames // 100)
        # min_samples controls how much a frame's neighborhood must be genuinely dense to
        # count as "core", rather than just a thin bridge between two otherwise separate states
        effective_min_samples = min_samples if min_samples is not None \
            else max(3, effective_min_cluster_size // 2)

        hdb = HDBSCAN(min_cluster_size=effective_min_cluster_size, min_samples=effective_min_samples, copy=True)
        labels = hdb.fit_predict(clustering_features)

        # If no region is denser than the rest, there are no separate states: the whole system is one.
        # HDBSCAN's allow_single_cluster is not used for this since, on a single blob whose density
        # keeps growing towards its center, it selects only the dense core and leaves the rest as noise.
        unique_labels = sorted(set(labels) - {-1})
        if len(unique_labels) == 0:
            print(' No density-separated states found: all frames are treated as a single cluster')
            labels = np.zeros(traj.n_frames, dtype=int)
            unique_labels = [0]

        # Build cluster lists (noise frames are excluded and handled separately below)
        first_frame_per_label = {lbl: int(np.where(labels == lbl)[0][0]) for lbl in unique_labels}
        ordered_labels = sorted(unique_labels, key=lambda lbl: first_frame_per_label[lbl])
        label_to_idx = {lbl: idx for idx, lbl in enumerate(ordered_labels)}
        clusters = [[] for _ in ordered_labels]
        for frame, lbl in enumerate(labels):
            if lbl != -1:
                clusters[label_to_idx[lbl]].append(frame)

        n_clusters = len(clusters)
        n_noise = int((labels == -1).sum())
        print(f' Found {n_clusters} clusters, {n_noise} noise frames')

        # Walk through noise segments and classify each one by its temporal context:
        # – noise between cluster A and cluster B (A≠B): transition frames A→B
        # – noise between cluster A and cluster A:       excursion frames from A
        # – leading / trailing noise:                    excursion of the bordering cluster
        transition_frame_map = {}   # (c_from, c_to) -> [[frame_indices], ...] one sublist per event
        excursion_frame_map  = {}   # c_cluster      -> [[frame_indices], ...] one sublist per event

        i = 0
        while i < traj.n_frames:
            if labels[i] != -1:
                i += 1
                continue
            # Collect the full contiguous noise segment
            noise_segment = []
            while i < traj.n_frames and labels[i] == -1:
                noise_segment.append(i)
                i += 1
            # Nearest non-noise cluster before and after the segment
            seg_start = noise_segment[0]
            before_lbl = next((labels[j] for j in range(seg_start - 1, -1, -1) if labels[j] != -1), None)
            after_lbl  = next((labels[j] for j in range(i, traj.n_frames)       if labels[j] != -1), None)
            # Resolve to cluster indices (fall back to whichever side exists)
            before_c = label_to_idx[before_lbl] if before_lbl is not None else None
            after_c  = label_to_idx[after_lbl]  if after_lbl  is not None else None
            anchor_c = before_c if before_c is not None else after_c
            if before_c is not None and after_c is not None and before_c != after_c:
                key = (before_c, after_c)
                transition_frame_map.setdefault(key, []).append(noise_segment)
            else:
                excursion_frame_map.setdefault(anchor_c, []).append(noise_segment)

        # Count transitions from the non-noise sequence only
        non_noise_clusters = [label_to_idx[labels[i]] for i in range(traj.n_frames) if labels[i] != -1]
        transition_counts = {}
        prev = non_noise_clusters[0]
        for curr in non_noise_clusters[1:]:
            if prev != curr:
                transition_counts[(prev, curr)] = transition_counts.get((prev, curr), 0) + 1
            prev = curr
        print(f' Found {sum(transition_counts.values())} transitions')

        # For each cluster find the most representative frame (closest to the centroid in PCA space)
        # and the least representative frame (farthest from the centroid, i.e. the cluster's own outlier)
        representative_frames = []
        least_representative_frames = []
        screenshot_parameters = None
        print(' Generating cluster screenshots')
        for c, cluster in enumerate(clusters):
            cluster_arr = np.array(cluster)
            centroid = clustering_features[cluster_arr].mean(axis=0)
            dists_to_centroid = np.sum((clustering_features[cluster_arr] - centroid) ** 2, axis=1)
            most_representative_frame = int(cluster_arr[np.argmin(dists_to_centroid)])
            least_representative_frame = int(cluster_arr[np.argmax(dists_to_centroid)])
            representative_frames.append(most_representative_frame)
            least_representative_frames.append(least_representative_frame)

            # Build one structure per frame so the screenshot can overlay both
            most_structure = auxiliar_structure.copy()
            most_structure.set_new_coordinates(traj[most_representative_frame].xyz[0] * 10)
            # Skip the overlay entirely for singleton clusters, where both frames are the same
            least_structure = None
            if least_representative_frame != most_representative_frame:
                least_structure = auxiliar_structure.copy()
                least_structure.set_new_coordinates(traj[least_representative_frame].xyz[0] * 10)

            screenshot_filepath = output_screenshot_filepath.replace('*', str(r).zfill(2)).replace('??', str(c).zfill(2))
            screenshot_file = File(screenshot_filepath)
            reprint(f' Generating cluster screenshot {c+1}/{n_clusters}')
            screenshot_parameters = get_screenshot(most_structure, screenshot_file, cg_selection,
                parameters=screenshot_parameters, reference_structure=least_structure)

        output_clusters = []
        for cluster_frames, most_representative_frame, least_representative_frame in zip(
            clusters, representative_frames, least_representative_frames
        ):
            output_clusters.append({
                'frames': cluster_frames,
                'main': most_representative_frame,
                'less': least_representative_frame,
            })

        # A transition may have a count but no frames (direct jump, no noise in between)
        # or frames but no separate count entry — union both sources by (from, to)
        all_transition_keys = set(transition_counts.keys()) | set(transition_frame_map.keys())
        output_transitions = [
            {
                'from': int(key[0]),
                'to': int(key[1]),
                'count': transition_counts.get(key, 0),
                'frames': transition_frame_map.get(key, []),
            }
            for key in all_transition_keys
        ]
        output_excursions = [
            {'cluster': int(c), 'frames': frames}
            for c, frames in excursion_frame_map.items()
        ]

        output_analysis = {
            'name': name,
            'cutoff': None,
            'clusters': output_clusters,
            'transitions': output_transitions,
            'excursions': output_excursions,
            'reduction': 'tica' if use_tica else 'pca',
            'version': '0.2.0',
        }

        save_json(output_analysis, numbered_output_analysis_filepath)

    save_json(output_summary, output_analysis_filepath)


def tica(features: np.ndarray, lag: int) -> tuple[np.ndarray, np.ndarray]:
    """Project frames onto their slowest collective motions (TICA).
    Frames must be consecutive in time. Returns all components sorted from slowest to fastest,
    together with their eigenvalues (autocorrelation at the given lag)."""
    n_frames, n_features = features.shape
    centered = features - features.mean(axis=0)
    instant, lagged = centered[:-lag], centered[lag:]
    # Symmetrized estimators, so the problem stays a real symmetric generalized eigenproblem
    covariance = (instant.T @ instant + lagged.T @ lagged) / (2 * (n_frames - lag))
    covariance += np.eye(n_features) * 1e-10
    lagged_covariance = (instant.T @ lagged) / (n_frames - lag)
    lagged_covariance = (lagged_covariance + lagged_covariance.T) / 2
    eigenvalues, eigenvectors = eigh(lagged_covariance, covariance)
    order = np.argsort(-np.abs(eigenvalues))
    return centered @ eigenvectors[:, order], eigenvalues[order]


def relaxation_time(eigenvalue: float, lag: int) -> float:
    """Time, in frames, a TICA component takes to lose its memory, from its autocorrelation at the lag."""
    autocorrelation = np.clip(abs(eigenvalue), 1e-12, 1 - 1e-12)
    return -lag / np.log(autocorrelation)
