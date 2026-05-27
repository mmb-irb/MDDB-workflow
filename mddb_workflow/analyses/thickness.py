from mddb_workflow.tools.get_reduced_trajectory import calculate_frame_step
from mddb_workflow.utils.auxiliar import save_json, load_json
from mddb_workflow.utils.constants import OUTPUT_THICKNESS_FILENAME
from mddb_workflow.utils.type_hints import *
from lipyphilic.analysis.z_positions import ZPositions
import numpy as np
from MDAnalysis.transformations import translate, set_dimensions


def set_box(u: 'Universe') -> np.ndarray:
    """Set the box dimensions of the universe based on the positions of the atoms."""
    # Initialize min and max positions with extreme values
    min_pos = np.full(3, np.inf)
    max_pos = np.full(3, -np.inf)

    # Iterate over all frames to find the overall min and max positions
    for ts in u.trajectory:
        positions = u.atoms.positions
        min_pos = np.minimum(min_pos, positions.min(axis=0))
        max_pos = np.maximum(max_pos, positions.max(axis=0))

    # Calculate the dimensions of the box
    box_dimensions = max_pos - min_pos
    zshift = -min_pos[2]  # Shift to ensure the minimum z is at 0
    transformations = [
        set_dimensions([*box_dimensions, 90, 90, 90]),
        translate(np.array([0.0, 0.0, zshift]))
    ]
    u.trajectory.add_transformations(*transformations)
    return zshift


def thickness(
    membrane_map: dict,
    universe: 'Universe',
    output_directory: str,
    snapshots: int,
    frames_limit: int = 100):
    """Membrane thickness analysis."""
    if membrane_map is None or membrane_map['n_mems'] == 0:
        print('-> No membranes found, skipping thickness analysis')
        return

    frame_step, frame_count = calculate_frame_step(snapshots, frames_limit)
    head_sel = []
    for n in range(membrane_map['n_mems']):
        head_sel.extend(membrane_map['mems'][str(n)]['polar_atoms']['top'])
        head_sel.extend(membrane_map['mems'][str(n)]['polar_atoms']['bot'])
    head_sel_mda = 'index ' + " ".join(map(str, (head_sel)))
    # Run the analysis on the whole membrane
    print(' Running LiPyphilic ZPositions')
    u_copy = universe.copy()
    zshift = set_box(u_copy)
    positions = ZPositions(u_copy, head_sel_mda, head_sel_mda, return_midpoint=True)
    positions.run(step=frame_step)

    # Calculate the mean z position of midplane wrt the box axis
    midplane_z = positions.memb_midpoint.flatten() - zshift
    zpos = positions.z_positions
    # All postive/negative of each frame
    zpos_per_frame = [zpos[:, fr][zpos[:, fr] > 0] for fr in range(frame_count)]
    zneg_per_frame = [zpos[:, fr][zpos[:, fr] < 0] for fr in range(frame_count)]
    # Mean and std per frame
    zpos_means = [zpos_frame.mean() for zpos_frame in zpos_per_frame]
    zneg_means = [zneg_frame.mean() for zneg_frame in zneg_per_frame]
    zpos_stds = [zpos_frame.std() for zpos_frame in zpos_per_frame]
    zneg_stds = [zneg_frame.std() for zneg_frame in zneg_per_frame]

    # Calculate the mean and std of the z positions of the head groups
    data = {'data': {
        'mean_positive': zpos_means,
        'mean_negative': zneg_means,
        'std_positive': zpos_stds,
        'std_negative': zneg_stds,
        'thickness': [zpos_mean - zneg_mean for zpos_mean, zneg_mean in zip(zpos_means, zneg_means)],
        # Sum of standard deviations
        'std_thickness': np.sqrt([zpos_std**2 for zpos_std in zpos_stds] + [zneg_std**2 for zneg_std in zneg_stds]),
        'midplane_z': midplane_z,
        },
        'version': '0.1.0',
    }
    # Convert any numpy arrays to lists for JSON serialization
    for key in data['data']:
        data['data'][key] = data['data'][key].tolist() if isinstance(data['data'][key], np.ndarray) else data['data'][key]
    # Save the results to a JSON file
    output_analysis_filepath = f'{output_directory}/{OUTPUT_THICKNESS_FILENAME}'
    save_json(data, output_analysis_filepath)
    return data


def plot_thickness(output_analysis_filepath):
    """Plot membrane thickness analysis results."""
    import matplotlib.pyplot as plt

    data = load_json(output_analysis_filepath)["data"]

    # Extract data for plotting
    frames = data['frame']
    thickness = data['thickness']
    std_thickness = data['std_thickness']
    midplane_z = data['midplane_z']
    mean_positive = data['mean_positive']
    mean_negative = data['mean_negative']
    std_positive = data['std_positive']
    std_negative = data['std_negative']

    # Create the plot
    plt.figure(figsize=(12, 8))

    # Plot thickness with error bars
    plt.errorbar(frames, thickness, yerr=std_thickness, label='Thickness', fmt='-o', capsize=3)

    # Plot midplane z position
    plt.plot(frames, midplane_z, label='Midplane Z', linestyle='--', color='orange')

    # Plot mean positive and mean negative with error bars
    eb1 = plt.errorbar(frames, mean_positive, yerr=std_positive, label='Mean Positive', fmt='-o', color='green', capsize=3)
    eb1[-1][0].set_linestyle('--')
    eb2 = plt.errorbar(frames, mean_negative, yerr=std_negative, label='Mean Negative', fmt='-o', color='red', capsize=3)
    eb2[-1][0].set_linestyle('--')

    # Add labels, legend, and title
    plt.xlabel('Frame'); plt.ylabel('Value')
    plt.title('Membrane Thickness, Midplane Z Position, and Mean Positive/Negative with Errors')
    plt.legend(); plt.grid(); plt.show()
