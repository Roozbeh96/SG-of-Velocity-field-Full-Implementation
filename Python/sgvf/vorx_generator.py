# SPDX-License-Identifier: CC-BY-NC-SA-4.0
# Copyright (c) 2026 ehsan010@umn.edu

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.special import ndtr, ndtri

from ._utils import tictoc
from .gen_gauss_cop import gen_gauss_cop


def _norm_pdf(x):
    return np.exp(-0.5 * x ** 2) / np.sqrt(2 * np.pi)


def _skewed_gauss_samples(C1, omega, kesi, beta, n, seed):
    """Samples of rho^{VorX}_{u,w} by inverting the c.d.f. of a skewed Gaussian."""
    def pdf_skewedGauss(x):
        return C1 * (1 - x ** 2) * (2 / omega * _norm_pdf((x - kesi) / omega)
                                    * ndtr(beta * (x - kesi) / omega))

    x_vals = np.linspace(-1, 1, 1000)
    # Step 1: Calculate the p.d.f. over the range
    pdf_vals = pdf_skewedGauss(x_vals)
    # Step 2: Compute the c.d.f. using cumulative trapezoidal integration
    cdf_vals = cumulative_trapezoid(pdf_vals, x_vals, initial=0)
    cdf_vals = cdf_vals / cdf_vals.max()  # Normalize to ensure c.d.f. ends at 1
    # Step 4: Generate uniform random numbers in [0, 1]
    rand_rho_uw_VorX = np.random.default_rng(seed).random(n)
    # Step 3 & 5: Invert the c.d.f. using interpolation
    return np.interp(rand_rho_uw_VorX, cdf_vals, x_vals)


def _radius_samples(rand_r_omega, mu, sigma, x_t, rand_x_t, alpha):
    """Log-normal core with a power-law tail for r_omega/lambda_T."""
    r = np.zeros_like(rand_r_omega)
    # rand_r_omega <= rand_x_t (Log-normal region)
    idx_log_normal = rand_r_omega <= rand_x_t
    r[idx_log_normal] = np.exp(mu + sigma * ndtri(rand_r_omega[idx_log_normal]))
    # rand_r_omega > rand_x_t (Power-law region)
    idx_power_law = rand_r_omega > rand_x_t
    r[idx_power_law] = x_t * ((1 - rand_r_omega[idx_power_law]) / (1 - rand_x_t)) ** (-1 / alpha)
    return r


def _deciles(omega_norm, u_w, r_w):
    """Split the generated vortices into 10 groups of the normalized vorticity."""
    omega_deciles, u_deciles, r_deciles = [], [], []
    for i in range(1, 11):
        L = (i - 1) * 10
        U = i * 10
        # MATLAB prctile corresponds to the 'hazen' method
        valL = np.percentile(omega_norm, L, method="hazen")
        valU = np.percentile(omega_norm, U, method="hazen")
        sel = (omega_norm >= valL) & (omega_norm <= valU)
        # Store the samples for this specific bracket
        omega_deciles.append(omega_norm[sel])
        u_deciles.append(u_w[sel])
        r_deciles.append(r_w[sel])
    return omega_deciles, u_deciles, r_deciles


def vorx_generator(obj):
    """Generate radius, maximum azimuthal velocity and corrcoef(u,w) of vortices (``VorX_generator.m``).

    Parameters can be found in Table 4 of the paper:
    https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/
    stochastic-generation-of-velocity-fields-to-reproduce-the-energy-spectrum-in-wall-turbulence/
    1D7BDC29D1D577DAE7672E86E21641B7

    Returns a dictionary with the 12 outputs of the MATLAB function (same names).
    ``r_wOlambda_T_*`` are normalized by lambda_T and ``u_wOu_tau_*`` by u_tau.
    The ``*_deciles`` entries are lists of 10 arrays.
    """
    with tictoc("VorX generation"):
        print("Generation of Vortices...")
        # ---------------- near the wall ----------------
        obj.ro_r_u_omega_near_wall = 0.4
        n = int(1e7)  # Number of samples
        rand_u_omega, rand_r_omega = gen_gauss_cop(obj.ro_r_u_omega_near_wall, n, 97)

        # Generating u_omega (truncated exponential)
        a, b, m = 0.65, 4.5, 0.55
        u_wOu_tau_near = -m * np.log(-rand_u_omega * (np.exp(-a / m) - np.exp(-b / m)) + np.exp(-a / m))

        # Generating r_omega
        r_wOlambda_T_near = _radius_samples(rand_r_omega, mu=-1.55, sigma=0.36, x_t=0.36,
                                            rand_x_t=0.93, alpha=5.0)
        # Check the p.d.f.:
        # plt.hist(r_wOlambda_T_near, bins=500, density=True)

        # Generating rho^{VorX}_{u,w}
        rho_uw_VorX_near = _skewed_gauss_samples(C1=1.23, omega=0.58, kesi=-0.61, beta=2.37,
                                                 n=n, seed=18)

        # ---------------- far from the wall ----------------
        obj.ro_r_u_omega_far_wall = 0.45
        rand_u_omega, rand_r_omega = gen_gauss_cop(obj.ro_r_u_omega_far_wall, n, 143)

        a, b, m = 0.4, 4.5, 0.44
        u_wOu_tau_far = -m * np.log(-rand_u_omega * (np.exp(-a / m) - np.exp(-b / m)) + np.exp(-a / m))

        r_wOlambda_T_far = _radius_samples(rand_r_omega, mu=-1.94, sigma=0.36, x_t=0.25,
                                           rand_x_t=0.93, alpha=4.5)

        rho_uw_VorX_far = _skewed_gauss_samples(C1=1.28, omega=0.66, kesi=-0.54, beta=1.80,
                                                n=n, seed=61)

        # Find the omega of the generated vortices
        omega_near = u_wOu_tau_near / r_wOlambda_T_near * obj.u_tau / obj.lambda_T
        omega_near_norm = omega_near / np.std(omega_near, ddof=1)

        omega_far = u_wOu_tau_far / r_wOlambda_T_far * obj.u_tau / obj.lambda_T
        omega_far_norm = omega_far / np.std(omega_far, ddof=1)

        # Dividing the generated vortices into deciles for the matter of matching
        omega_near_deciles, u_wOu_tau_near_deciles, r_wOlambda_T_near_deciles = _deciles(
            omega_near_norm, u_wOu_tau_near, r_wOlambda_T_near)
        omega_far_deciles, u_wOu_tau_far_deciles, r_wOlambda_T_far_deciles = _deciles(
            omega_far_norm, u_wOu_tau_far, r_wOlambda_T_far)

    return dict(
        u_wOu_tau_near=u_wOu_tau_near, r_wOlambda_T_near=r_wOlambda_T_near,
        rho_uw_VorX_near=rho_uw_VorX_near,
        u_wOu_tau_far=u_wOu_tau_far, r_wOlambda_T_far=r_wOlambda_T_far,
        rho_uw_VorX_far=rho_uw_VorX_far,
        omega_near_deciles=omega_near_deciles, u_wOu_tau_near_deciles=u_wOu_tau_near_deciles,
        r_wOlambda_T_near_deciles=r_wOlambda_T_near_deciles,
        omega_far_deciles=omega_far_deciles, u_wOu_tau_far_deciles=u_wOu_tau_far_deciles,
        r_wOlambda_T_far_deciles=r_wOlambda_T_far_deciles,
    )
