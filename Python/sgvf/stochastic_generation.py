# SPDX-License-Identifier: CC-BY-NC-SA-4.0
# Copyright (c) 2026 ehsan010@umn.edu

import pickle
from pathlib import Path

import numpy as np

from ._utils import colon
from .mvg_kernel import mvg_kernel
from .sg_vel_field import sg_vel_field
from .sg_vel_prof import sg_vel_prof
from .sg_vorx import sg_vorx
from .spectral import spec_analysis_E11, spec_analysis_E12, spec_analysis_E22
from .structure_function import str_func_D11, str_func_D22


class StochasticGeneration:
    """Methods for SGVP, SGVF, SG of vortex cores and their addition.

    Python version of the ``stochastic_generation`` MATLAB class. Array
    properties keep the MATLAB layout: profiles are ``(len(z), N)`` and the
    high-resolution fields are ``(len(z), len(Gen_x_HRVF), n_frames)``.
    Indices stored in the object (frames, profile numbers) are 0-based.
    """

    def __init__(self, u_tau, z_0, ks, delta, lambda_T, nu, z_min, z_max, res_z):
        self.u_tau = u_tau
        self.z_0 = z_0
        self.ks = ks
        self.delta = delta
        self.lambda_T = lambda_T
        self.nu = nu
        self.z = colon(z_min, res_z, z_max)

        # Properties filled by the different steps
        self.ro_uw = self.N_prof = self.Delta_x = None
        self.Gen_u_prof = self.Gen_w_prof = self.log_data_SGVP = None
        self.Gen_u_LRVF = self.Gen_w_LRVF = self.Gen_x_LRVF = self.hist_corr = None
        self.Kr = self.Gen_u_HRVF = self.Gen_w_HRVF = self.Gen_x_HRVF = None
        self.Lambda_ci = self.Lambda_cirms = self.omega = self.omega_rms = None
        self.ro_r_u_omega_near_wall = self.ro_r_u_omega_far_wall = None
        self.near_wall_VorX_data = self.prograde_VorX_data = self.retrograde_VorX_data = None
        self.E11 = self.k1 = self.E22 = self.E12 = self.epsilon_E11 = self.eta_E11 = None
        self.D11 = self.r11 = self.epsilon_D11 = self.eta_D11 = self.D22 = None
        self.cdf_omega = self.mean_omega = None
        self.production = self.eta_prod = None

    def SGVP(self, ro_uw, N):
        """Generate N velocity profiles stochastically."""
        self.ro_uw = ro_uw
        self.N_prof = int(N)
        self.Gen_u_prof, self.Gen_w_prof, self.log_data_SGVP = sg_vel_prof(self)

    def SGVF(self, Delta_x, M, seed=None):
        """Reorganize the profiles generated in SGVP into a correlated velocity field."""
        self.Delta_x = Delta_x
        self.Gen_u_LRVF, self.Gen_w_LRVF, self.hist_corr = sg_vel_field(self, M, seed=seed)

    def SGVorX(self, u_wOu_tau_near, r_wOlambda_T_near, rho_uw_VorX_near,
               u_wOu_tau_far, r_wOlambda_T_far, rho_uw_VorX_far,
               omega_near_deciles, u_wOu_tau_near_deciles, r_wOlambda_T_near_deciles,
               omega_far_deciles, u_wOu_tau_far_deciles, r_wOlambda_T_far_deciles, seed=None):
        """Add the generated vortices to the SGVF."""
        sg_vorx(self, u_wOu_tau_near, r_wOlambda_T_near, rho_uw_VorX_near,
                u_wOu_tau_far, r_wOlambda_T_far, rho_uw_VorX_far,
                omega_near_deciles, u_wOu_tau_near_deciles, r_wOlambda_T_near_deciles,
                omega_far_deciles, u_wOu_tau_far_deciles, r_wOlambda_T_far_deciles, seed=seed)

    def Structure_function_D11(self):
        str_func_D11(self)

    def Structure_function_D22(self):
        str_func_D22(self)

    def Spectral_analysis_E11(self):
        spec_analysis_E11(self)

    def Spectral_analysis_E22(self):
        spec_analysis_E22(self)

    def Spectral_analysis_E12(self):
        spec_analysis_E12(self)

    def Gauss_kernel(self, kern_xsize, kern_ysize, Sigma):
        """Apply the multivariate Gaussian kernel."""
        mvg_kernel(self, kern_xsize, kern_ysize, Sigma)

    # ------------------------------------------------------------------
    def save(self, path="Gen_sample.pkl", exclude=()):
        """Save the object with pickle.

        ``exclude`` lists attribute names that are not saved (e.g.
        ``("Lambda_ci", "omega")`` to reduce the file size).
        """
        state = {k: v for k, v in self.__dict__.items() if k not in set(exclude)}
        with open(Path(path), "wb") as f:
            pickle.dump(state, f, protocol=pickle.HIGHEST_PROTOCOL)

    @classmethod
    def load(cls, path="Gen_sample.pkl"):
        """Load an object written by :meth:`save`."""
        with open(Path(path), "rb") as f:
            state = pickle.load(f)
        obj = cls.__new__(cls)
        obj.__dict__.update(state)
        return obj

    def __repr__(self):
        filled = [k for k, v in self.__dict__.items() if isinstance(v, np.ndarray)]
        return f"StochasticGeneration(u_tau={self.u_tau}, delta={self.delta}, arrays={filled})"
