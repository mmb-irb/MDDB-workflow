# Numba kernel of the grid-free energies engine
#
# Per-atom sums of the pair energies of each atom with every atom of the other group: Coulomb with the Mehler-Solmajer
# sigmoidal dielectric (no cutoff) and CMIP's van der Waals within a cutoff. energies.get_pair_energy_sums imports this
# module on its first call, so numba is only loaded when the grid-free engine runs.
#
# The pairs are the same as in the numpy version this kernel replaced (all pairs, in double precision); only the order
# of the additions differs, so per-atom sums agree to float rounding.

import os
import math
import numba
import numpy as np

# OpenMP threads spin-wait between parallel regions and slow everything down when the machine is busy (measured
# 2026-11-12: 0.45 s instead of 5 ms per call), so use numba's own thread pool, which does not spin, unless the user
# chose a layer. This must be set before the first parallel call in the process.
if 'NUMBA_THREADING_LAYER' not in os.environ:
    numba.config.THREADING_LAYER = 'workqueue'

# Rows (atoms of the larger group) are split in this many chunks, whatever the number of threads, and the column sums
# of each chunk are added in chunk order: results do not depend on the number of threads.
ROW_CHUNKS = 64


@numba.njit(parallel=True, cache=True)
def _pair_energy_sums(
    coords_1: np.ndarray,
    coords_2: np.ndarray,
    coulomb_charges_1: np.ndarray,
    charges_2: np.ndarray,
    epsilons_1: np.ndarray,
    epsilons_2: np.ndarray,
    sigmas_1: np.ndarray,
    sigmas_2: np.ndarray,
    bulk_dielectric: float,
    radii_scale: float,
    vdw_cutoff: float,
    minimum_distance: float,
):
    count_1, count_2 = len(coords_1), len(coords_2)
    vdw_1 = np.zeros(count_1)
    es_1 = np.zeros(count_1)
    chunk_size = (count_1 + ROW_CHUNKS - 1) // ROW_CHUNKS
    chunks = (count_1 + chunk_size - 1) // chunk_size
    vdw_columns = np.zeros((chunks, count_2))
    es_columns = np.zeros((chunks, count_2))
    a = -8.5525
    b = bulk_dielectric - a
    exponent_factor = -0.003627 * b
    for chunk in numba.prange(chunks):
        vdw_column = vdw_columns[chunk]
        es_column = es_columns[chunk]
        for i in range(chunk * chunk_size, min(count_1, (chunk + 1) * chunk_size)):
            x, y, z = coords_1[i, 0], coords_1[i, 1], coords_1[i, 2]
            charge, epsilon, sigma = coulomb_charges_1[i], epsilons_1[i], sigmas_1[i]
            vdw_sum = 0.0
            es_sum = 0.0
            for j in range(count_2):
                dx = x - coords_2[j, 0]
                dy = y - coords_2[j, 1]
                dz = z - coords_2[j, 2]
                distance = max(math.sqrt(dx * dx + dy * dy + dz * dz), minimum_distance)
                dielectric = a + b / (1 + 7.7839 * math.exp(exponent_factor * distance))
                es = charge * charges_2[j] / (dielectric * distance)
                es_sum += es
                es_column[j] += es
                if distance <= vdw_cutoff:
                    ratio = radii_scale * (sigma + sigmas_2[j]) / 2 / distance
                    ratio_3 = ratio * ratio * ratio
                    ratio_6 = ratio_3 * ratio_3
                    vdw = 4 * math.sqrt(epsilon * epsilons_2[j]) * (ratio_6 * ratio_6 - ratio_6)
                    vdw_sum += vdw
                    vdw_column[j] += vdw
            vdw_1[i] = vdw_sum
            es_1[i] = es_sum
    vdw_2 = np.zeros(count_2)
    es_2 = np.zeros(count_2)
    for chunk in range(chunks):
        vdw_2 += vdw_columns[chunk]
        es_2 += es_columns[chunk]
    return vdw_1, es_1, vdw_2, es_2


def get_thread_count() -> int:
    """Get the threads for the kernel: OMP_NUM_THREADS if set (the variable CMIP follows), else all numba threads."""
    maximum = numba.config.NUMBA_NUM_THREADS
    requested = os.environ.get('OMP_NUM_THREADS', '')
    if requested.isdigit() and int(requested) > 0:
        return min(int(requested), maximum)
    return maximum


def pair_energy_sums(
    coords_1: np.ndarray,
    coords_2: np.ndarray,
    charges_1: np.ndarray,
    charges_2: np.ndarray,
    epsilons_1: np.ndarray,
    epsilons_2: np.ndarray,
    sigmas_1: np.ndarray,
    sigmas_2: np.ndarray,
    coulomb_constant: float,
    bulk_dielectric: float,
    radii_scale: float,
    vdw_cutoff: float,
    minimum_distance: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Sum the van der Waals and electrostatic energies of the pairs of each atom with the atoms of the other group.

    Return the vdw and es sums of the atoms of group 1 and then those of group 2.
    The larger group is split in rows, so the per-chunk column sums take the least memory.
    """
    groups = [
        (np.ascontiguousarray(coords, dtype=np.float64), np.asarray(charges, dtype=np.float64),
         np.asarray(epsilons, dtype=np.float64), np.asarray(sigmas, dtype=np.float64))
        for coords, charges, epsilons, sigmas in [
            (coords_1, charges_1, epsilons_1, sigmas_1),
            (coords_2, charges_2, epsilons_2, sigmas_2),
        ]
    ]
    swap = len(groups[1][0]) > len(groups[0][0])
    rows, columns = (groups[1], groups[0]) if swap else (groups[0], groups[1])
    numba.set_num_threads(get_thread_count())
    vdw_rows, es_rows, vdw_columns, es_columns = _pair_energy_sums(
        rows[0], columns[0], coulomb_constant * rows[1], columns[1], rows[2], columns[2], rows[3], columns[3],
        bulk_dielectric, radii_scale, vdw_cutoff, minimum_distance,
    )
    if swap:
        return vdw_columns, es_columns, vdw_rows, es_rows
    return vdw_rows, es_rows, vdw_columns, es_columns
