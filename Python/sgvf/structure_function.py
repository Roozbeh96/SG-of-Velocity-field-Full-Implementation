# SPDX-License-Identifier: CC-BY-NC-SA-4.0
# Copyright (c) 2026 ehsan010@umn.edu

"""Second-order structure functions of the generated field.

Python versions of ``str_func_D11.m`` and ``str_func_D22.m``. More info about
how to compute the structure function in:
https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/
article/stochastic-modal-velocity-field-in-roughwall-turbulence/
800D97E60609F4FE4D92847CCDB1C7A0
"""

import math

import numpy as np

from ._utils import concat_fluctuation, first_index, frame_mean, mround, progressbar

#: Elevations z/delta at which the structure functions are computed.
Z_OVER_DELTA = (0.015, 0.05, 0.0919, 0.14, 0.24)


def _setup(obj):
    dx = obj.Gen_x_HRVF[1] - obj.Gen_x_HRVF[0]
    win15del = math.ceil(15 * obj.delta / dx)
    n_concat = obj.Gen_x_HRVF.size * obj.Gen_u_HRVF.shape[2]
    Numb_frame = n_concat // win15del

    # --- LOGARITHMIC SKIP SETUP ---
    # Define how many points you want in your log scale (e.g., 50 or 100 points)
    num_log_points = 50
    # Generate unique, rounded integer steps logarithmically from 0 to win15del
    r_indices = np.unique(np.concatenate(
        ([0], mround(np.logspace(0, np.log10(win15del), num_log_points))))).astype(int)
    # Ensure it doesn't accidentally exceed bounds
    r_indices = r_indices[r_indices <= win15del]

    z_over_delta = np.asarray(obj.z) / obj.delta
    Z = [i for i in (first_index(z_over_delta >= v) for v in Z_OVER_DELTA) if i is not None]
    return dx, win15del, Numb_frame, r_indices, Z


def _structure_function(obj, field, desc):
    dx, win15del, Numb_frame, r_indices, Z = _setup(obj)
    num_r = r_indices.size
    # Preallocate only for the skipped points
    D = np.zeros((field.shape[0], num_r, Numb_frame))
    mean_rc = frame_mean(field)
    concat = {elev: concat_fluctuation(field, elev, mean_rc) for elev in Z}

    for S in progressbar(range(Numb_frame), desc):
        frame_start = S * win15del
        frame_end = (S + 1) * win15del  # inclusive
        for elev in Z:
            seg = concat[elev][frame_start:frame_end + 1]
            # Loop through our log-spaced indices instead of every integer
            for ri, r in enumerate(r_indices):
                D[elev, ri, S] = np.mean((seg[:seg.size - r] - seg[r:]) ** 2)
    return D, dx * r_indices


def str_func_D11(obj):
    """Streamwise structure function D11(r); also sets epsilon_D11 and eta_D11."""
    obj.D11, obj.r11 = _structure_function(obj, obj.Gen_u_HRVF, "Computing Structure Function D11")
    # Avoid division by zero at r=0 for the epsilon calculation: ignore the
    # first point (r=0).
    with np.errstate(divide="ignore", invalid="ignore"):
        max_val = np.max(obj.D11[:, 1:, :] * (obj.r11[1:] ** (-2 / 3))[None, :, None], axis=1)
        obj.epsilon_D11 = (2.3 / max_val) ** (-3 / 2)
        obj.eta_D11 = (obj.nu ** 3 / obj.epsilon_D11) ** 0.25


def str_func_D22(obj):
    """Vertical structure function D22(r)."""
    obj.D22, obj.r11 = _structure_function(obj, obj.Gen_w_HRVF, "Computing Structure Function D22")
