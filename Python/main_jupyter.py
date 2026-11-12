# SPDX-License-Identifier: CC-BY-NC-SA-4.0
# Copyright (c) 2026 ehsan010@umn.edu

"""
To Stochastically Generate Velocity Profiles (SGVP), reorganize, and add vortex
cores, you need to provide initial physical parameters such as
u_tau [m/s] (friction velocity), z_0 [m] (aerodynamic roughness length),
delta [m] (boundary layer thickness), lambda [m] (Taylor micro-scale). The
details about the rest of the parameters are provided in the following lines.
If you need details about the statistics, refer to the following papers:

**For SG of Velocity Profiles:
https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/
stochastic-modelling-of-the-instantaneous-velocity-profile-in-roughwall-
turbulent-boundary-layers/492F3CD03C8C3E7ED306E9117B848B5E

**For SG of Velocity field:
https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/
stochastic-modal-velocity-field-in-roughwall-turbulence/
800D97E60609F4FE4D92847CCDB1C7A0

**For SG of Vortex core:
https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/
stochastic-generation-of-velocity-fields-to-reproduce-the-energy-spectrum-
in-wall-turbulence/1D7BDC29D1D577DAE7672E86E21641B7

Two ways to use this file:

1. As a script: ``poetry run python main_jupyter.py`` (``--help`` lists the
   options) or "Run Python File" goes through the whole pipeline from zero to
   the end and saves the object, as main.py does.

2. Interactive (like MATLAB sections): in VS Code click "Run Cell" above a
   ``# %%`` line. Run the first two cells once: the first defines the
   functions, the second loads the saved object ``Gen_sample`` (skipped if it
   is already in memory, e.g. loaded by run_plots.py in the same interactive
   window). Then run any step cell on its
   own, e.g. "Increase the resolution" followed by "Analysis" gives the
   spectra and structure functions of the HRVF without vortices, and
   "Add vortices", "Gaussian filter" and "Analysis" give them for the final
   field again.

The steps overwrite the attributes of the object in place (e.g.
``Gen_u_HRVF``, ``E11``, ``D11``). Results of later steps are not cleared, so
after going back, rerun the steps whose results you want for the new field.
Steps from "Increase the resolution" on use only the saved LRVF and are
deterministic, except the vortex positions in "Add vortices" (``SEED``).
"""

# %% Pipeline functions and setup (run this cell first)
import argparse
import math
from pathlib import Path

from IPython import get_ipython

from sgvf import (StochasticGeneration, cdf_omega, increasing_res, lambda_ci,
                  production_comp, small_scale_sigma, vorx_generator)


def setup(u_tau, z_0, ks, delta, lambda_T, nu):
    """Create the object for the given physical parameters.

    u_tau [m/s] (friction velocity), z_0 [m] (aerodynamic roughness length),
    ks [m] (roughness height), delta [m] (boundary layer thickness),
    lambda_T [m] (Taylor micro-scale), nu [m^2/s] (kinematic viscosity).
    """
    # Define the wall-normal extent of the generated profiles which covers the
    # logarithmic layer. The generated profiles start from z = 1.5*ks and
    # continue till z/delta = 0.25. The values here are conservatively
    # considered. The wall-normal resolution of profiles is 0.4*lambda/6.
    # Since the average thickness of a shear layer is 0.4*lambda, it is
    # divided into 6 to resolve the shear layer.
    z_min = 1.5 * ks  # [m]
    z_max = 0.25 * delta  # [m]
    res_z = 0.4 * lambda_T / 6  # [m]

    return StochasticGeneration(u_tau, z_0, ks, delta, lambda_T, nu, z_min, z_max, res_z)


def sgvp(Gen_sample, N_prof=15e4, ro_uw=-0.33):
    """Stochastic Generation of Velocity Profiles (SGVP).

    Sets ``Gen_u_prof``, ``Gen_w_prof`` and ``log_data_SGVP``.
    """
    # For the stochastic generation of velocity profiles, we need to start
    # from z_min, then stochastically generate thickness (step), streamwise
    # (or modal) velocity, and vertical velocity. Repeat the process of step
    # generation until we reach z_max, and then stop. We can independently
    # generate as many profiles as we want. The minimum number of profiles
    # for convergence in statistics (stationary dataset) is ~100 profiles.
    # However, since we want to generate a long velocity field in the next
    # step, we need more profiles.
    # The parameter we need here is the Pearson correlation coefficient
    # between streamwise and vertical velocity components (ro_uw). This
    # parameter is necessary for the Reynolds shear stress (covariance)
    # reproduction. This value in the wind tunnel (m1) dataset is -0.4. You
    # can tweak it w.r.t. your dataset. N_prof is the number of profiles we
    # want to generate. The parameters for the Townsend-Perry constant (A, B)
    # can be set inside the function.
    # Can also use Priyadarshana and Klewicki (2004) equation for the
    # correlation coefficient:
    # U_inf = 50
    # n = 4
    # Theta = Gen_sample.delta * n / ((n + 1) * (n + 2))
    # Re_tau = Gen_sample.u_tau * Gen_sample.delta / Gen_sample.nu
    # ro_uw = -0.03*math.log(U_inf*Theta/(Gen_sample.u_tau*Gen_sample.delta)*Re_tau) + 0.63
    Gen_sample.SGVP(ro_uw, N_prof)


def sgvf(Gen_sample, seed=None):
    """Stochastic Generation of Velocity Field (SGVF).

    Sets ``Gen_u_LRVF``, ``Gen_w_LRVF`` and ``hist_corr``.
    """
    # For the generated velocity field, you need to prescribe the distance
    # between the velocity profiles which here is lambda_T. The storage size M
    # is formulated through a Re dependent equation.
    # In the function, first, the shear layers between the adjacent UMZs are
    # smoothed using harmonic functions.
    Delta_x = Gen_sample.lambda_T
    M = math.ceil(150 * math.sqrt(Gen_sample.u_tau * Gen_sample.delta / Gen_sample.nu))
    Gen_sample.SGVF(Delta_x, M, seed=seed)


def increase_resolution(Gen_sample, L=1, Kr=10, parallel=False):
    """Increase the resolution of the LRVF (HRVF without vortices).

    Overwrites ``Gen_u_HRVF``, ``Gen_w_HRVF``, ``Gen_x_HRVF``, ``Gen_x_LRVF``
    and ``Kr`` from ``Gen_u_LRVF`` and ``Gen_w_LRVF``.
    """
    # At the first step, we need to divide the VF into L*delta intervals, and
    # increase the resolution. For increasing the resolution, the factor is
    # another input parameter. Here Kr is the factor parameter which is
    # consistent with the parameter in the paper. All these steps are done in
    # the increasing_res function. The band-stop filter is implemented inside
    # the function.
    # With parallel=True, increasing_res, lambda_ci and cdf_omega process the
    # frames in parallel (like parfor in MATLAB); otherwise in series.
    increasing_res(Gen_sample, L, Kr, parallel=parallel)


def add_vortices(Gen_sample, seed=None, parallel=False):
    """Add the stochastically generated vortex cores to the HRVF.

    Computes the Lambda_ci field and the c.d.f. of omega of the current HRVF,
    generates the vortices and adds them to ``Gen_u_HRVF`` and ``Gen_w_HRVF``
    in place. Run it once on a HRVF without vortices (after
    ``increase_resolution``).
    """
    # Lambda_ci (swirling strength) field
    # In the next step, we need to find the Lambda_ci field (swirling
    # strength) to insert the vortices in the shear layer.
    lambda_ci(Gen_sample, parallel=parallel)

    # Generating radius, maximum azimuthal velocity, and corrcoeff(u,w)
    # At this step, we generate radius (r_omega), maximum azimuthal velocity
    # (u_omega), and the correlation coefficient between u and w component of
    # the generated vortex rho^{VorX}_{u,w}. Generated values for r_omega are
    # normalized by lambda_T, for u_omega are normalized by u_tau.
    # For the decile matching algorithm, the vorticity of the generated
    # vortices is also computed.
    vorx = vorx_generator(Gen_sample)

    # Find the c.d.f. of the omega distribution
    # At this step, we find the c.d.f. of the omega distribution in the
    # non-zero lambda_ci regions. The output of this step is used in the
    # decile matching process for adaptive selection of vortices.
    cdf_omega(Gen_sample, parallel=parallel)

    # Adding Stochastically Generated Vortex Cores to the HRVF (both near
    # wall and shear layers).
    Gen_sample.SGVorX(**vorx, seed=seed)


def gaussian_filter(Gen_sample, kern_xsize=0.4, kern_ysize=0.4):
    """Convolve the HRVF with the small-scale 2-D Gaussian filter (in place)."""
    # computing std of the small-scale 2-D Gaussian filter
    # s0 is the initial guess (change it if it is necessary)
    Sigma_small, sol = small_scale_sigma(Gen_sample.lambda_T, w0=1e-5, li_factor=0.4, s0=1e-7)
    print(sol)

    # Convolving the small-scale 2-D Gaussian filter
    Gen_sample.Gauss_kernel(kern_xsize, kern_ysize, Sigma_small)


def analysis(Gen_sample):
    """Spectra, structure functions and production of the current HRVF."""
    # The spatial and spectral analysis can be done on different stages of
    # the generated velocity field: LRVF, HRVF (w/ or w/o band-stop filter),
    # HRVF (w/ or w/o inserting VorX), HRVF (w/ or w/o passing 2-D filter).
    # If you want to compute it on LRVF, run increase_resolution(Gen_sample,
    # Kr=1) first, then compute the spectra or structure function.
    Gen_sample.Spectral_analysis_E11()
    Gen_sample.Spectral_analysis_E22()
    Gen_sample.Spectral_analysis_E12()
    Gen_sample.Structure_function_D11()
    Gen_sample.Structure_function_D22()

    # Computing production
    production_comp(Gen_sample)


def save(Gen_sample, output="Gen_sample.pkl"):
    """Save the object with pickle."""
    Gen_sample.save(output)
    print(f"Saved to {output}")


def main(u_tau, z_0, ks, delta, lambda_T, nu, N_prof=15e4, output="Gen_sample.pkl", seed=None,
         parallel=False):
    """Run the whole pipeline for the given physical parameters."""
    # Let us generate an object.
    Gen_sample = setup(u_tau, z_0, ks, delta, lambda_T, nu)

    sgvp(Gen_sample, N_prof)
    sgvf(Gen_sample, seed=seed)
    increase_resolution(Gen_sample, parallel=parallel)
    add_vortices(Gen_sample, seed=seed, parallel=parallel)
    gaussian_filter(Gen_sample)
    analysis(Gen_sample)

    if output:
        save(Gen_sample, output)
    return Gen_sample


# get_ipython() returns the running IPython/Jupyter session (VS Code's
# Interactive Window), or None when the file is run as a normal script.
_ipython = get_ipython()
INTERACTIVE = _ipython is not None

if __name__ == "__main__" and not INTERACTIVE:
    # Running as a script: the whole pipeline from zero to the end.
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--n-prof", type=float, default=16e4,
                        help="number of generated profiles (default: %(default)g)")
    parser.add_argument("--output", default="Gen_sample.pkl",
                        help="output pickle file; pass '' to skip saving")
    parser.add_argument("--seed", type=int, default=None,
                        help="seed of the field reorganization and vortex insertion "
                             "(default: different every run)")
    parser.add_argument("--parallel", action="store_true",
                        help="run increasing_res, lambda_ci and cdf_omega in parallel "
                             "(default: in series)")
    args = parser.parse_args()

    # The physical properties of the Saddoughi High-velocity (50 [m/s]) are
    # used. You can find details about this dataset in the papers.
    u_tau = 2.325  # [m/s]
    z_0 = 3.8e-4  # [m]
    kappa = 0.39
    AFR = 8.5
    ks = z_0 / math.exp(-kappa * AFR)
    delta = 1.09  # [m]
    lambda_T = 0.007  # [m]
    nu = 1.57e-5  # [m^2/s]

    main(u_tau, z_0, ks, delta, lambda_T, nu,
         N_prof=args.n_prof, output=args.output, seed=args.seed,
         parallel=args.parallel)
    raise SystemExit

# %% Load the saved object (interactive use; run once per session)
if INTERACTIVE and not any(e.endswith("autoreload") for e in _ipython.extension_manager.loaded):
    # Reload sgvf automatically when its files are saved.
    _ipython.run_line_magic("load_ext", "autoreload")
    _ipython.run_line_magic("autoreload", "2")

try:
    HERE = Path(__file__).resolve().parent
except NameError:  # the interactive window runs in the folder of this file
    HERE = Path.cwd()

# Saved object to work with (written by the script mode above or main.py)
PKL_FILE = HERE / "Gen_sample_before_VorX.pkl"
# Options of the interactive steps
PARALLEL = False
SEED = None

# Reading the file is skipped if it is already in memory, so rerunning this
# cell is cheap. Delete Gen_sample (del Gen_sample) to force a new read.
if "Gen_sample" not in globals():
    print(f"Loading {PKL_FILE} ...")
    Gen_sample = StochasticGeneration.load(PKL_FILE)

# %% Increase the resolution (HRVF without vortices, from the saved LRVF)
increase_resolution(Gen_sample, L=1, Kr=10, parallel=PARALLEL)

# %% Add vortices
add_vortices(Gen_sample, seed=SEED, parallel=PARALLEL)

# %% Gaussian filter
gaussian_filter(Gen_sample, kern_xsize=0.4, kern_ysize=0.4)

# %% Analysis (spectra, structure functions and production of the current HRVF)
analysis(Gen_sample)

# %% Save (use a new name to keep the original file)
save(Gen_sample, HERE / "Gen_sample_VorX_before_filter.pkl")

