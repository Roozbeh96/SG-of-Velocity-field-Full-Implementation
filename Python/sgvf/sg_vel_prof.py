# SPDX-License-Identifier: CC-BY-NC-SA-4.0
# Copyright (c) 2026 ehsan010@umn.edu

import math
from bisect import bisect_right

import numpy as np
from scipy.special import ndtri

from ._utils import progressbar, tictoc
from .gen_gauss_cop import gen_gauss_cop

#: Column names of ``obj.log_data_SGVP``.
LOG_COLUMNS = ("prof_num", "hm_i[m]", "zm_i[m]", "um_i[m/s]", "wm_i[m/s]")


def sg_vel_prof(obj):
    """Generate step-like velocity profiles (Python version of ``SG_VelProf.m``).

    Returns
    -------
    Gen_u_prof, Gen_w_prof : ndarray, shape (len(z), N_prof)
    log_data : ndarray, shape (n_steps, 5)
        One row per uniform momentum zone with columns ``LOG_COLUMNS``
        (the MATLAB version stores the same values in nested cells of
        dictionaries). ``prof_num`` is 0-based.
    """
    with tictoc("SGVP"):
        N_rand = int(1e6)
        repro_numb = 137  # for reproducibility, we generate one set of random numbers
        # gen_gauss_cop will generate random numbers with a prescribed correlation
        # coefficient out of the Gaussian copula.
        rand_u, rand_w = gen_gauss_cop(obj.ro_uw, N_rand, repro_numb)
        N_prof = int(obj.N_prof)
        z = np.asarray(obj.z)
        nz = z.size
        # Column-major storage: every profile (column) is contiguous in memory.
        Gen_u_prof = np.zeros((nz, N_prof), order="F")
        Gen_w_prof = np.zeros((nz, N_prof), order="F")
        kappa = 0.39
        rng = np.random.default_rng(37)
        rand_h = rng.random(N_rand)

        # sqrt(2)*erfinv(2*rand-1) is the standard normal quantile; evaluate it
        # once for all pre-generated random numbers (lists are faster to index
        # in the loop below).
        g_h = ndtri(rand_h).tolist()
        g_u = ndtri(rand_u).tolist()
        g_w = ndtri(rand_w).tolist()

        z_list = z.tolist()
        z_start = z_list[0]
        z_end = z_list[-1]
        delta = obj.delta
        z_0 = obj.z_0
        u_tau = obj.u_tau
        log_data = []
        counter_1 = 0
        # Townsend-Perry constants
        A = 1.26
        B = 1.00
        for prof_num in progressbar(range(N_prof), "Profile Generation"):
            z_i = z_start
            z_ind = 0
            prof_log_start = len(log_data)
            while z_i < z_end:
                # Fitted model is available in section 5.2:
                # https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/
                # stochastic-modelling-of-the-instantaneous-velocity-profile-in-roughwall-
                # turbulent-boundary-layers/492F3CD03C8C3E7ED306E9117B848B5E
                mu_Hm_zi = -3.59 * (z_i / delta) ** 0.91
                std_Hm_zi = 1.0

                hm_i = z_i * math.exp(mu_Hm_zi + std_Hm_zi * g_h[counter_1])

                # In case if you want to remove profiles with just one-step
                # if z_i + hm_i > z_end and counter_2 == 1:
                #     counter_1 += 1
                #     continue

                var_uu_plus = -A * math.log((z_i + hm_i / 2) / delta) + B
                mean_u_plus = 1 / kappa * math.log((z_i + hm_i / 2) / z_0)

                if var_uu_plus <= 0:
                    # Restart the profile from z_min.
                    z_i = z_start
                    z_ind = 0
                    del log_data[prof_log_start:]
                    counter_1 += 1
                    if counter_1 >= rand_u.size:
                        raise RuntimeError(
                            f"Not enough pre-generated random numbers for profile #{prof_num} "
                            "after restarting. Please increase N_rand.")
                    continue
                std_uu_plus = math.sqrt(var_uu_plus)

                um_i = u_tau * (mean_u_plus + std_uu_plus * g_u[counter_1])

                std_ww_plus = 0.85
                mean_w_plus = 0.0

                wm_i = u_tau * (mean_w_plus + std_ww_plus * g_w[counter_1])

                # Check if we are beyond the z_{max}
                row = bisect_right(z_list, z_i + hm_i)  # first index with z > z_i + hm_i

                if row < nz:
                    log_data.append((prof_num, hm_i, z_i + hm_i / 2, um_i, wm_i))
                    Gen_u_prof[z_ind:row, prof_num] = um_i
                    Gen_w_prof[z_ind:row, prof_num] = wm_i
                    z_i = z_i + hm_i
                else:
                    log_data.append((prof_num, z_end - z_i, 0.5 * (z_i + z_end), um_i, wm_i))
                    Gen_u_prof[z_ind:, prof_num] = um_i
                    Gen_w_prof[z_ind:, prof_num] = wm_i
                    z_i = z_end
                z_ind = row
                counter_1 += 1
                if counter_1 >= rand_u.size:
                    raise RuntimeError(
                        f"Not enough pre-generated random numbers for profile #{prof_num}. "
                        "Please increase the number of generated random samples. "
                        f"Current limit: {rand_u.size}, required index: {counter_1 + 1}.")

    return Gen_u_prof, Gen_w_prof, np.array(log_data, dtype=float).reshape(-1, len(LOG_COLUMNS))
