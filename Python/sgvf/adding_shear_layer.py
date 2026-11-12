# SPDX-License-Identifier: CC-BY-NC-SA-4.0
# Copyright (c) 2026 ehsan010@umn.edu

import numpy as np

from ._utils import progressbar


def adding_shear_layer(obj):
    """Smooth the sharp shear layer between the UMZs (``adding_shear_layer.m``).

    Details are in section 2.2 of the paper:
    https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/
    stochastic-modal-velocity-field-in-roughwall-turbulence/800D97E60609F4FE4D92847CCDB1C7A0

    ``obj.Gen_u_prof`` and ``obj.Gen_w_prof`` are modified in place.
    """
    delta_omega = 0.4 * obj.lambda_T
    theta = np.linspace(0, np.pi / 2, 25)
    # Harmonic shape of the shear layer (the same for every jump).
    shape_uw = np.concatenate((1 - np.cos(theta), 1 + np.sin(theta)))
    shape_z = np.concatenate((np.sin(theta) - 1, 1 - np.cos(theta)))
    z = np.asarray(obj.z)
    u_prof = obj.Gen_u_prof
    w_prof = obj.Gen_w_prof
    # One sample of step-like velocity profile. You can plot it before and
    # after smoothing shear-layer:
    # plt.plot(obj.Gen_u_prof[:, 0], obj.z)
    for i in progressbar(range(u_prof.shape[1]), "Adding Shear Layer"):
        u_col = u_prof[:, i]
        w_col = w_prof[:, i]
        Delta_u_m = np.diff(u_col)
        Delta_w_m = np.diff(w_col)

        change_indices = np.flatnonzero(Delta_u_m != 0.0)

        for c_idx in change_indices:
            Delta_u_mi = Delta_u_m[c_idx]
            Delta_w_mi = Delta_w_m[c_idx]

            u = u_col[c_idx] + Delta_u_mi / 2 * shape_uw
            w = w_col[c_idx] + Delta_w_mi / 2 * shape_uw
            zc = (z[c_idx] + z[c_idx + 1]) / 2 + delta_omega / 2 * shape_z

            # grid points with zc[0] <= z <= zc[-1]
            lo = np.searchsorted(z, zc[0], side="left")
            hi = np.searchsorted(z, zc[-1], side="right")
            if hi <= lo:
                continue
            # first point of the smoothed curve with zc >= z(k)
            cc = np.searchsorted(zc, z[lo:hi], side="left")

            u_col[lo:hi] = u[cc]
            w_col[lo:hi] = w[cc]
    return obj
