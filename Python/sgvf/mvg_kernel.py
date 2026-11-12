# SPDX-License-Identifier: CC-BY-NC-SA-4.0
# Copyright (c) 2026 ehsan010@umn.edu

import math

import numpy as np
from scipy.optimize import brentq

from ._utils import conv2_same, progressbar, tictoc


def small_scale_sigma(lambda_T, w0=1e-5, li_factor=0.4, s0=1e-7):
    """Covariance of the small-scale 2-D Gaussian filter (section in ``main.m``).

    Solves ``2*sqrt(-2*s*log(2*w0*pi*s)) = li`` with ``li = li_factor*lambda_T``
    starting the search from ``s0`` (like ``fzero(f, s0)``), and returns the
    2x2 covariance matrix ``[[s, 0], [0, s]]`` together with ``s``.
    """
    li = li_factor * lambda_T

    def f(s):
        return 2 * (-2 * s * math.log(2 * w0 * math.pi * s)) ** 0.5 - li

    # Search for a sign change around s0 (as fzero does), staying inside the
    # domain s > 0, then use Brent's method.
    dx = s0 / 50
    a = b = s0
    fa = fb = f(s0)
    while fa * fb > 0:
        dx *= math.sqrt(2)
        a = max(s0 - dx, s0 * 1e-12)
        b = s0 + dx
        fa, fb = f(a), f(b)
        if dx > 1e6 * s0:
            raise RuntimeError("Could not bracket the root; change the initial guess s0.")
    sol = brentq(f, a, b, xtol=1e-30, rtol=4 * np.finfo(float).eps)
    # Check: li should be equal to li_factor*lambda_T
    # li = 2*(-2*sol*log(2*pi*w0*sol))**0.5
    return np.array([[sol, 0.0], [0.0, sol]]), sol


def mvg_kernel(obj, kern_xsize, kern_ysize, Sigma):
    """Apply the multivariate Gaussian kernel to the HRVF (``MVG_kernel.m``)."""
    with tictoc("Gaussian kernel"):
        mu = np.zeros(2)
        Sigma = np.asarray(Sigma, dtype=float)
        dx = obj.Gen_x_HRVF[1] - obj.Gen_x_HRVF[0]
        dz = obj.z[1] - obj.z[0]

        X, Y = np.meshgrid(
            np.linspace(-kern_xsize / 2 * obj.lambda_T, kern_xsize / 2 * obj.lambda_T,
                        math.ceil(kern_xsize * obj.lambda_T / dx) + 1),
            np.linspace(-kern_ysize / 2 * obj.lambda_T, kern_ysize / 2 * obj.lambda_T,
                        math.ceil(kern_ysize * obj.lambda_T / dz) + 1))
        # mvnpdf
        pos = np.stack((X.ravel(), Y.ravel()), axis=1) - mu
        inv = np.linalg.inv(Sigma)
        expo = -0.5 * np.einsum("ij,jk,ik->i", pos, inv, pos)
        kernel = np.exp(expo) / (2 * np.pi * np.sqrt(np.linalg.det(Sigma)))
        kernel = kernel.reshape(X.shape)
        # If you want to see the distribution of the Kernel:
        # plt.contourf(X/obj.lambda_T, Y/obj.lambda_T, kernel, 100); plt.axis('equal')

        nrows, ncols, nframes = obj.Gen_u_HRVF.shape
        normalization_matrix = conv2_same(np.ones((nrows, ncols)), kernel)
        for S in progressbar(range(nframes), "Passing Gaussian Kernel"):
            obj.Gen_u_HRVF[:, :, S] = conv2_same(obj.Gen_u_HRVF[:, :, S], kernel) / normalization_matrix
            obj.Gen_w_HRVF[:, :, S] = conv2_same(obj.Gen_w_HRVF[:, :, S], kernel) / normalization_matrix
