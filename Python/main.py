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

Run with ``poetry run python main.py`` (``--help`` lists the options).
"""

import argparse
import math


from sgvf import (StochasticGeneration, cdf_omega, increasing_res, lambda_ci,
                  production_comp, small_scale_sigma, vorx_generator)


def main(u_tau, z_0, ks, delta, lambda_T, nu, N_prof=15e4, output="Gen_sample.pkl", seed=None,
         parallel=False):
    """Run the whole pipeline for the given physical parameters.

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

    # Let us generate an object.
    Gen_sample = StochasticGeneration(u_tau, z_0, ks, delta, lambda_T, nu, z_min, z_max, res_z)

    # %% Stochastic Generation of Velocity Profiles (SGVP)
    # For the stochastic generation of velocity profiles, we need to start
    # from z_min, then stochastically generate thickness (step), streamwise
    # (or modal) velocity, and vertical velocity. Repeat the process of step
    # generation until we reach z_max, and then stop. We can independently
    # generate as many profiles as we want. The minimum number of profiles
    # for convergence in statistics (stationary dataset) is ~100 profiles.
    # However, since we want to generate a long velocity field in the next
    # step, we need more profiles.
    # The parameter we need here is the Pearson correlation coefficient
    # between streamwise and vertical velocity components. This parameter is
    # necessary for the Reynolds shear stress (covariance) reproduction. This
    # value in the wind tunnel (m1) dataset is -0.4. You can tweak it w.r.t.
    # your dataset. N_prof is the number of profiles we want to generate.
    # The parameters for the Townsend-Perry constant (A, B) can be set inside
    # the function.
    ro_uw = -0.33
    # Can also use Priyadarshana and Klewicki (2004) equation for the
    # correlation coefficient:
    U_inf = 50
    n = 4
    Theta = delta * n / ((n + 1) * (n + 2))
    Re_tau = u_tau * delta / nu
    # ro_uw = -0.03*math.log(U_inf*Theta/(u_tau*delta)*Re_tau) + 0.63

    Gen_sample.SGVP(ro_uw, N_prof)

    # %% Stochastic Generation of Velocity Field (SGVF)
    # For the generated velocity field, you need to prescribe the distance
    # between the velocity profiles which here is lambda_T. The storage size M
    # is formulated through a Re dependent equation.
    # In the function, first, the shear layers between the adjacent UMZs are
    # smoothed using harmonic functions.
    Delta_x = lambda_T
    M = math.ceil(150 * math.sqrt(Gen_sample.u_tau * Gen_sample.delta / Gen_sample.nu))
    Gen_sample.SGVF(Delta_x, M, seed=seed)

    # %% Increasing the resolution of the SGVF and compute the lambda_ci field
    # At the first step, we need to divide the VF into 1*delta (L = 1)
    # intervals, and increase the resolution. For increasing the resolution,
    # the factor is another input parameter. Here Kr is the factor parameter
    # which is consistent with the parameter in the paper. All these steps
    # are done in the increasing_res function. The band-stop filter is
    # implemented inside the function.
    # With parallel=True, increasing_res, lambda_ci and cdf_omega process the
    # frames in parallel (like parfor in MATLAB); otherwise in series.
    L = 1
    Kr = 10
    increasing_res(Gen_sample, L, Kr, parallel=parallel)

    # %% Lambda_ci (swirling strength) field
    # In the next step, we need to find the Lambda_ci field (swirling
    # strength) to insert the vortices in the shear layer.
    lambda_ci(Gen_sample, parallel=parallel)

    # %% Generating radius, maximum azimuthal velocity, and corrcoeff(u,w)
    # At this step, we generate radius (r_omega), maximum azimuthal velocity
    # (u_omega), and the correlation coefficient between u and w component of
    # the generated vortex rho^{VorX}_{u,w}. Generated values for r_omega are
    # normalized by lambda_T, for u_omega are normalized by u_tau.
    # For the decile matching algorithm, the vorticity of the generated
    # vortices is also computed.
    vorx = vorx_generator(Gen_sample)

    # %% Find the c.d.f. of the omega distribution
    # At this step, we find the c.d.f. of the omega distribution in the
    # non-zero lambda_ci regions. The output of this step is used in the
    # decile matching process for adaptive selection of vortices.
    cdf_omega(Gen_sample, parallel=parallel)

    # %% Adding Stochastically Generated Vortex Cores to the HRVF (both near
    # wall and shear layers).
    Gen_sample.SGVorX(**vorx, seed=seed)

    # %% computing std of the small-scale 2-D Gaussian filter
    # s0 is the initial guess (change it if it is necessary)
    Sigma_small, sol = small_scale_sigma(Gen_sample.lambda_T, w0=1e-5, li_factor=0.4, s0=1e-7)
    print(sol)

    # %% Convolving the small-scale 2-D Gaussian filter
    kern_xsize = 0.4
    kern_ysize = 0.4
    Gen_sample.Gauss_kernel(kern_xsize, kern_ysize, Sigma_small)

    # %% Spectral and structure function analysis
    # The spatial and spectral analysis can be done on different stages of
    # the generated velocity field: LRVF, HRVF (w/ or w/o band-stop filter),
    # HRVF (w/ or w/o inserting VorX), HRVF (w/ or w/o passing 2-D filter).
    # If you want to compute it on LRVF, run increasing_res(Gen_sample, L, Kr)
    # with Kr = 1 first, then compute the spectra or structure function.
    Gen_sample.Spectral_analysis_E11()
    Gen_sample.Spectral_analysis_E22()
    Gen_sample.Spectral_analysis_E12()
    Gen_sample.Structure_function_D11()
    Gen_sample.Structure_function_D22()

    # %% Computing production
    production_comp(Gen_sample)

    # %% Saving
    if output:
        Gen_sample.save(output)
        print(f"Saved to {output}")
    return Gen_sample


if __name__ == "__main__":
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
