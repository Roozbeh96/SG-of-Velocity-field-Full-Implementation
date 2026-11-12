# SPDX-License-Identifier: CC-BY-NC-SA-4.0
# Copyright (c) 2026 ehsan010@umn.edu

import numpy as np

from ._utils import run_frames, std_matlab, tictoc


def lambda_ci(obj, parallel=False, n_workers=None):
    """Swirling strength field of the HRVF (``Lambda_ci.m``).

    ``parallel=True`` processes the frames in parallel on ``n_workers``
    threads (``None``: all CPU cores), like ``parfor`` in the MATLAB version;
    ``parallel=False`` (default) processes them in series. Both give
    identical results.

    Sets ``obj.Lambda_ci`` (signed with the vorticity), ``obj.omega`` (the
    vorticity where swirling exists, zero elsewhere), and the per-elevation
    rms of their non-zero values ``obj.Lambda_cirms`` and ``obj.omega_rms``.
    """
    with tictoc("Lambda_ci"):
        Gen_u_HRVF = obj.Gen_u_HRVF
        Gen_w_HRVF = obj.Gen_w_HRVF
        z = np.asarray(obj.z)

        nRows, nCols, nFrames = Gen_u_HRVF.shape

        Lambda_ci_local = np.zeros((nRows, nCols, nFrames))
        omega_local_out = np.zeros((nRows, nCols, nFrames))

        DelxGen = obj.Gen_x_HRVF[1] - obj.Gen_x_HRVF[0]
        DelzGen = z[1] - z[0]

        print("Computing Lambda_{ci} Field...")

        def process_frame(S):
            # Compute gradients for the current frame (rows: z, columns: x)
            dudz, dudx = np.gradient(Gen_u_HRVF[:, :, S], DelzGen, DelxGen)
            dwdz, dwdx = np.gradient(Gen_w_HRVF[:, :, S], DelzGen, DelxGen)

            omega = dwdx - dudz

            # Vectorized eigenvalue calculation for 2x2 matrix
            # The imaginary part of eigenvalues for [a b; c d] is
            # sqrt(abs(min(0, (a+d)^2 - 4*(a*d - b*c)))) / 2
            tr = dudx + dwdz
            detJ = (dudx * dwdz) - (dudz * dwdx)
            discriminant = (tr ** 2) - 4 * detJ

            # Calculate swirling strength
            swirl = np.zeros((nRows, nCols))
            mask = discriminant < 0
            swirl[mask] = np.sqrt(np.abs(discriminant[mask])) / 2

            Lambda_ci_local[:, :, S] = swirl * np.sign(omega)

            # Store omega where swirling exists
            temp_omega = np.zeros((nRows, nCols))
            temp_omega[mask] = omega[mask]
            omega_local_out[:, :, S] = temp_omega

        # numpy releases the GIL, so in parallel the frames run on threads.
        run_frames(process_frame, ((S,) for S in range(nFrames)), nFrames,
                   "Computing Lambda_ci Field", parallel=parallel, n_workers=n_workers)

        # Statistics Calculation (RMS)
        Lambda_cirms_ = np.zeros(nRows)
        omega_rms_ = np.zeros(nRows)

        for z_ in range(nRows):
            temp_ = Lambda_ci_local[z_]
            Lambda_cirms_[z_] = std_matlab(temp_[temp_ != 0])

            temp__ = omega_local_out[z_]
            omega_rms_[z_] = std_matlab(temp__[temp__ != 0])

        obj.Lambda_ci = Lambda_ci_local
        obj.omega = omega_local_out
        obj.Lambda_cirms = Lambda_cirms_
        obj.omega_rms = omega_rms_
