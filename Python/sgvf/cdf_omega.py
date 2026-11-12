# SPDX-License-Identifier: CC-BY-NC-SA-4.0
# Copyright (c) 2026 ehsan010@umn.edu

import numpy as np
from scipy.ndimage import find_objects

from ._utils import ecdf, label_columnwise, mround, run_frames, tictoc


def components(labels, n):
    """Yield ``(rows, cols)``, the 0-based pixel indices of every labelled component."""
    for j, slc in enumerate(find_objects(labels, max_label=n), start=1):
        if slc is None:
            continue
        local = labels[slc] == j
        rr, cc = np.nonzero(local)
        yield rr + slc[0].start, cc + slc[1].start


def _frame_values(lambda_ci_frame, omega_frame, threshold, omega_rms):
    """Normalized mean vorticity of every strong shear layer in one frame."""
    sz1, sz2 = omega_frame.shape
    mask = np.abs(lambda_ci_frame) >= threshold
    Matrix = omega_frame * mask
    binaryMatrix = Matrix != 0
    labeledMatrix, numComponents = label_columnwise(binaryMatrix)

    frame_values = []
    for rows, cols in components(labeledMatrix, numComponents):
        values = Matrix[rows, cols]

        sum_val = values.sum()
        if sum_val == 0:
            continue  # Avoid division by zero

        # 1-based centroids, as in MATLAB
        row_centroid = int(mround(np.sum((rows + 1) * values) / sum_val))
        col_centroid = int(mround(np.sum((cols + 1) * values) / sum_val))

        # Bounds and validity check
        if (row_centroid <= 0 or row_centroid > sz1
                or col_centroid <= 0 or col_centroid > sz2
                or Matrix[row_centroid - 1, col_centroid - 1] == 0):
            continue

        # Normalize and store
        frame_values.append(values.mean() / omega_rms[row_centroid - 1])
    return frame_values


def cdf_omega(obj, gamma_filter=0.35, parallel=False, n_workers=None):
    """C.d.f. of the mean vorticity of the shear layers (``cdf_omega.m``).

    Regions with ``|Lambda_ci| >= gamma_filter * Lambda_cirms`` are kept; the
    mean vorticity of every connected region, normalized by ``omega_rms`` at
    the elevation of its centroid, is collected and its empirical c.d.f. is
    stored in ``obj.cdf_omega`` and ``obj.mean_omega``. The result is used in
    the percentile matching of the vortex insertion.

    ``parallel=True`` processes the frames in parallel on ``n_workers``
    processes (``None``: all CPU cores), like ``parfor`` in the MATLAB
    version; ``parallel=False`` (default) processes them in series. Both give
    identical results. Processes are used because the per-region loop is
    pure Python; scripts that use ``parallel=True`` need an
    ``if __name__ == "__main__":`` guard.
    """
    with tictoc("cdf_omega"):
        print("Finding percentile of the omega in the field...")
        # Filtering (keep strong shear layers)
        threshold = gamma_filter * np.asarray(obj.Lambda_cirms)[:, None]
        omega_rms = np.asarray(obj.omega_rms)
        num_frames = obj.omega.shape[2]

        args = ((np.ascontiguousarray(obj.Lambda_ci[:, :, S]),
                 np.ascontiguousarray(obj.omega[:, :, S]), threshold, omega_rms)
                for S in range(num_frames))
        per_frame = run_frames(_frame_values, args, num_frames, "cdf of omega",
                               parallel=parallel, n_workers=n_workers, processes=True)
        mean_omega_y_norm = [v for frame_values in per_frame for v in frame_values]

        cdf_, mean_omega_ = ecdf(np.abs(np.array(mean_omega_y_norm)))
        obj.cdf_omega = cdf_
        obj.mean_omega = mean_omega_
