# SPDX-License-Identifier: CC-BY-NC-SA-4.0
# Copyright (c) 2026 ehsan010@umn.edu

import math

import numpy as np

from ._utils import first_index, frame_mean, label_columnwise, mround, progressbar, tictoc
from .cdf_omega import components

#: Column names of ``obj.near_wall_VorX_data``.
NEAR_WALL_COLUMNS = ("frame", "xc_VorX[m]", "zc_VorX[m]", "u_omega[m/s]", "r_omega[m]",
                     "Dist_to_prev_VorX[m]")
#: Column names of ``obj.prograde_VorX_data`` and ``obj.retrograde_VorX_data``.
LAMBDACI_COLUMNS = ("frame", "xc_VorX[m]", "zc_VorX[m]", "u_omega[m/s]", "r_omega[m]",
                    "rotation[{-1,1}]")


def _transform(u_VorX, w_VorX, rho_uw_VorX):
    """``[u w] * TVorX`` with ``TVorX = [1, rho; 0, sqrt(1-rho^2)]``."""
    return u_VorX, u_VorX * rho_uw_VorX + w_VorX * math.sqrt(1 - rho_uw_VorX ** 2)


def lambdaci_vorx(obj, i, S, xcs, zcs, rot, r_wOlambda_T, u_wOu_tau, rho_uw_VorX, ind_,
                  N, cdf_, mean_omega_, omega_, omega_deciles, r_wOlambda_T_deciles,
                  u_wOu_tau_deciles, rng, prograde, retrograde):
    """Insert one vortex at the centroid of a shear layer (``Lambdaci_VorX.m``).

    The vortex is selected from the decile of generated vortices that matches
    the percentile of the shear-layer vorticity. ``r_wOlambda_T[ind_]`` and
    ``u_wOu_tau[ind_]`` are overwritten with the selected vortex (in place).
    Prograde/retrograde vortex records are appended to the given lists.
    """
    # Matching the percentile of the vorticity of the shear layer with the
    # vorticity of the generated vortex.
    abs_omega = abs(omega_)
    idx = np.searchsorted(mean_omega_, abs_omega, side="right") - 1  # last mean_omega_ <= |omega_|
    if np.isnan(abs_omega) or idx < 0:
        percentile_rank = 1
    else:
        percentile_rank = min(max(int(math.ceil(cdf_[idx] * 100 / 10)), 1), 10)
    d = percentile_rank - 1
    omega_samples = omega_deciles[d]
    randomIndex = rng.integers(omega_samples.shape[0])

    r_omega = r_wOlambda_T_deciles[d][randomIndex] * obj.lambda_T
    u_omega = u_wOu_tau_deciles[d][randomIndex] * obj.u_tau
    Gama = 2 * math.pi * u_omega * r_omega / (1 - math.exp(-1))

    r_wOlambda_T[ind_] = r_wOlambda_T_deciles[d][randomIndex]
    u_wOu_tau[ind_] = u_wOu_tau_deciles[d][randomIndex]

    z = obj.z
    x = obj.Gen_x_HRVF
    indz = np.flatnonzero((z >= zcs[i] - N * r_omega) & (z <= zcs[i] + N * r_omega))
    indx = np.flatnonzero((x >= xcs[i] - N * r_omega) & (x <= xcs[i] + N * r_omega))
    X, Z = np.meshgrid(x[indx], z[indz])
    r = np.sqrt((X - xcs[i]) ** 2 + (Z - zcs[i]) ** 2)
    # Making the vortex circular shape.
    r[r > N * r_omega] = 0

    if r.shape[0] == 1:
        return

    with np.errstate(divide="ignore", invalid="ignore"):
        uazi = Gama / (2 * math.pi * r) * (1 - np.exp(-r ** 2 / r_omega ** 2))
    uazi[np.isnan(uazi)] = 0

    Theta = np.arctan2((Z - zcs[i]), (X - xcs[i]))
    Theta[np.isnan(Theta)] = 0

    record = (S, xcs[i], zcs[i], u_omega, r_omega, rot[i])
    if rot[i] < 0:
        prograde.append(record)
        u_VorX = uazi * np.sin(Theta)
        w_VorX = -uazi * np.cos(Theta)
    else:
        retrograde.append(record)
        u_VorX = -uazi * np.sin(Theta)
        w_VorX = uazi * np.cos(Theta)

    if u_VorX.size == 0:
        return

    u_TVorX, w_TVorX = _transform(u_VorX, w_VorX, rho_uw_VorX)

    # Adding to the velocity field
    iz, ix = np.ix_(indz, indx)
    obj.Gen_u_HRVF[iz, ix, S] += u_TVorX
    obj.Gen_w_HRVF[iz, ix, S] += w_TVorX


def sg_vorx(obj, u_wOu_tau_near, r_wOlambda_T_near, rho_uw_VorX_near,
            u_wOu_tau_far, r_wOlambda_T_far, rho_uw_VorX_far,
            omega_near_deciles, u_wOu_tau_near_deciles, r_wOlambda_T_near_deciles,
            omega_far_deciles, u_wOu_tau_far_deciles, r_wOlambda_T_far_deciles,
            gamma_filter=0.35, seed=None):
    """Add near-wall vortices and shear-layer vortices to the HRVF (``SG_VorX.m``).

    Near : z/k_s < 2
    Far : z/k_s >= 2

    ``obj.Gen_u_HRVF`` and ``obj.Gen_w_HRVF`` are modified in place. The
    inserted vortices are recorded in ``obj.near_wall_VorX_data``,
    ``obj.prograde_VorX_data`` and ``obj.retrograde_VorX_data`` (arrays with
    columns ``NEAR_WALL_COLUMNS`` / ``LAMBDACI_COLUMNS``, 0-based frame index).
    """
    with tictoc("SGVorX"):
        rng = np.random.default_rng(seed)
        # The selected vortices overwrite entries of these arrays (as in the
        # MATLAB version); work on copies so the caller's arrays are unchanged.
        r_wOlambda_T_near = np.array(r_wOlambda_T_near, dtype=float)
        u_wOu_tau_near = np.array(u_wOu_tau_near, dtype=float)
        r_wOlambda_T_far = np.array(r_wOlambda_T_far, dtype=float)
        u_wOu_tau_far = np.array(u_wOu_tau_far, dtype=float)

        z = np.asarray(obj.z)
        x = np.asarray(obj.Gen_x_HRVF)
        obj.z = z
        obj.Gen_x_HRVF = x
        nz, nx = z.size, x.size
        x_max = x.max()

        # some counter parameters
        ii = 0  # counter for r_omega and u_omega of z/ks<2 vortices
        # (could be near-wall or lambda_ci vortices)
        iii = 0  # counter for z/ks>=2 vortices (lambda_ci vortices)

        # Filter low intense lambda_ci regions.
        threshold = gamma_filter * np.asarray(obj.Lambda_cirms)[:, None]

        cdf_omega_ = obj.cdf_omega
        mean_omega_ = obj.mean_omega

        # fluctuating velocity uprime = u - mean(u, 3) of the field before
        # adding vortices
        mean_u = frame_mean(obj.Gen_u_HRVF)

        near_wall, prograde, retrograde = [], [], []

        for S in progressbar(range(obj.Gen_u_HRVF.shape[2]), "Adding VorX to the Velocity Field"):
            # S is iterating over the Velocity field with a length of delta
            uprime = obj.Gen_u_HRVF[:, :, S] - mean_u
            x_r_margin_near_wall = 0.0
            temp_near_wall_VorX_counter = 1
            prev_pos_xc_near_wall_VorX = 0.0
            # Adding near-wall vortices (square shape). As long as we are in
            # the velocity field region, we keep add near-wall vortices.
            while x_r_margin_near_wall < x_max:
                N = 1
                r_omega = r_wOlambda_T_near[ii] * obj.lambda_T
                u_omega = u_wOu_tau_near[ii] * obj.u_tau
                rho = rho_uw_VorX_near[ii]
                Gama = 2 * math.pi * u_omega * r_omega / (1 - math.exp(-1))
                ii += 1
                # Find center position of the generated vortex.
                zc_near_wall = N * r_omega + z[0]
                xc_near_wall = N * r_omega + x_r_margin_near_wall
                x_r_margin_near_wall = xc_near_wall + N * r_omega
                indzmid = first_index(z >= zc_near_wall)
                indxmid = first_index(x >= xc_near_wall)

                # Two conditions for putting near-wall vortices:
                # 1. uprime at center of the vortex should be positive.
                # 2. randomly generated binary number should be one.
                if indzmid is not None and indxmid is not None and uprime[indzmid, indxmid] < 0:
                    continue
                if rng.integers(0, 2) == 0:
                    continue
                # Find distance between the near-wall vortices.
                if temp_near_wall_VorX_counter > 1:
                    dist = xc_near_wall - prev_pos_xc_near_wall_VorX
                else:
                    dist = 0.0
                near_wall.append((S, xc_near_wall, zc_near_wall, u_omega, r_omega, dist))
                prev_pos_xc_near_wall_VorX = xc_near_wall
                temp_near_wall_VorX_counter += 1

                # Generating vortex field. Explained in methodology of the paper.
                indz = np.flatnonzero((z >= zc_near_wall - N * r_omega) & (z <= zc_near_wall + N * r_omega))
                indx = np.flatnonzero((x >= xc_near_wall - N * r_omega) & (x <= xc_near_wall + N * r_omega))
                X, Z = np.meshgrid(x[indx], z[indz])
                # square shape vortex
                r = np.sqrt((X - xc_near_wall) ** 2 + (Z - zc_near_wall) ** 2)

                # Generating the azimuthal velocity using equation 3.1.
                with np.errstate(divide="ignore", invalid="ignore"):
                    uazi = Gama / (2 * math.pi * r) * (1 - np.exp(-r ** 2 / r_omega ** 2))
                uazi[np.isnan(uazi)] = 0

                Theta = np.arctan2((Z - zc_near_wall), (X - xc_near_wall))
                Theta[np.isnan(Theta)] = 0

                u_VorX = uazi * np.sin(Theta)
                # Check whether generated vortex can be resolved in our mesh
                # domain or not.
                if u_VorX.size == 0:
                    continue
                # Subtracting the center velocity from u_VorX
                nz_mask = uazi != 0
                u_VorX[nz_mask] = u_VorX[nz_mask] - Gama / (2 * math.pi * r_omega) * (1 - math.exp(-1))
                w_VorX = -uazi * np.cos(Theta)

                u_TVorX, w_TVorX = _transform(u_VorX, w_VorX, rho)

                # Adding to the velocity field.
                iz, ix = np.ix_(indz, indx)
                obj.Gen_u_HRVF[iz, ix, S] += u_TVorX
                obj.Gen_w_HRVF[iz, ix, S] += w_TVorX

            # Adding vortices into shear layers (Primary vortices)
            mask = np.abs(obj.Lambda_ci[:, :, S]) >= threshold
            Matrix_lambdaci = obj.Lambda_ci[:, :, S] * mask
            Matrix_omega = obj.omega[:, :, S] * mask
            binaryMatrix = Matrix_lambdaci != 0

            # Label connected components (4-connectivity)
            labeledMatrix, numb_clusters = label_columnwise(binaryMatrix)

            # Calculate the mass center (centroid) for each connected component
            centroids = []
            for rows, cols in components(labeledMatrix, numb_clusters):
                values = Matrix_lambdaci[rows, cols]
                sign_region_mean = np.sign(values.mean())
                omega_region_mean = Matrix_omega[rows, cols].mean()
                # Calculate the centroid (1-based, as in MATLAB) and round it
                with np.errstate(divide="ignore", invalid="ignore"):
                    row_centroid = mround(np.sum((rows + 1) * values) / values.sum())
                    col_centroid = mround(np.sum((cols + 1) * values) / values.sum())
                if not (np.isfinite(row_centroid) and np.isfinite(col_centroid)):
                    continue
                row_centroid = int(row_centroid)
                col_centroid = int(col_centroid)
                if row_centroid <= 0 or row_centroid > nz:
                    continue
                if col_centroid <= 0 or col_centroid > nx:
                    continue
                # Store the centroid (0-based indices)
                centroids.append((row_centroid - 1, col_centroid - 1, sign_region_mean,
                                  cols.min(), cols.max(),
                                  omega_region_mean / obj.omega_rms[row_centroid - 1]))

            if not centroids:
                continue
            centroids = np.array(centroids)
            row_c = centroids[:, 0].astype(int)
            col_c = centroids[:, 1].astype(int)
            col_min = centroids[:, 3].astype(int)
            col_max = centroids[:, 4].astype(int)
            zcs = z[row_c]
            xcs = x[col_c]
            rot = centroids[:, 2]
            omega_ = centroids[:, 5]
            Numb_vortices = centroids.shape[0]
            N = 2
            lam = obj.lambda_T
            for i in range(Numb_vortices):
                if zcs[i] < 2 * obj.ks:  # close to the wall
                    r_arr, u_arr, rho_arr = r_wOlambda_T_near, u_wOu_tau_near, rho_uw_VorX_near
                    dec = (omega_near_deciles, r_wOlambda_T_near_deciles, u_wOu_tau_near_deciles)
                    k = ii
                else:  # Far from the wall (z/ks >= 2)
                    r_arr, u_arr, rho_arr = r_wOlambda_T_far, u_wOu_tau_far, rho_uw_VorX_far
                    dec = (omega_far_deciles, r_wOlambda_T_far_deciles, u_wOu_tau_far_deciles)
                    k = iii

                def insert(ind_):
                    lambdaci_vorx(obj, i, S, xcs, zcs, rot, r_arr, u_arr, rho_arr[ind_], ind_, N,
                                  cdf_omega_, mean_omega_, omega_[i], dec[0], dec[1], dec[2],
                                  rng, prograde, retrograde)

                # Primary vortex
                ind_Pri_VorX = k
                insert(k)
                x_ = xcs[i] + N * r_arr[ind_Pri_VorX] * lam
                k += 1
                # Filling the right side till reaching to right margin (Secondary vortices)
                while x_ + N * r_arr[k] * lam < x[col_max[i]]:
                    insert(k)
                    x_ = x_ + 2 * N * r_arr[k] * lam
                    k += 1
                x_ = xcs[i] - N * r_arr[ind_Pri_VorX] * lam
                k += 1
                # Filling the left side till reaching to left margin (Secondary vortices)
                while x_ - N * r_arr[k] * lam > x[col_min[i]]:
                    insert(k)
                    x_ = x_ - 2 * N * r_arr[k] * lam
                    k += 1
                k += 1

                if zcs[i] < 2 * obj.ks:
                    ii = k
                else:
                    iii = k

        obj.near_wall_VorX_data = np.array(near_wall, dtype=float).reshape(-1, len(NEAR_WALL_COLUMNS))
        obj.prograde_VorX_data = np.array(prograde, dtype=float).reshape(-1, len(LAMBDACI_COLUMNS))
        obj.retrograde_VorX_data = np.array(retrograde, dtype=float).reshape(-1, len(LAMBDACI_COLUMNS))
