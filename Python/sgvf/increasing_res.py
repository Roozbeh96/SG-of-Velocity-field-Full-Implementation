# SPDX-License-Identifier: CC-BY-NC-SA-4.0
# Copyright (c) 2026 ehsan010@umn.edu

import math

import numpy as np
from scipy.interpolate import Akima1DInterpolator

from ._utils import mround, run_frames, tictoc


def increasing_res(obj, L0, factor, parallel=False, n_workers=None):
    """Increase the streamwise resolution of the LRVF (``increasing_res.m``).

    The LRVF is cut into frames of length ``L0*delta``. Each row of a frame is
    interpolated with modified Akima (``makima``) interpolation onto a grid
    that is ``factor`` times finer, and a band-stop filter is applied in
    Fourier space. The first 20% of the frames are discarded.

    Parameters
    ----------
    parallel : bool
        ``False`` (default) processes the frames one after the other.
        ``True`` processes them in parallel with a pool of threads (like
        ``parfor`` in the MATLAB version). Both give identical results.
    n_workers : int or None
        Number of threads when ``parallel`` is ``True``. ``None`` uses all
        CPU cores.

    Sets ``obj.Kr``, ``obj.Gen_x_HRVF``, ``obj.Gen_x_LRVF``, ``obj.Gen_u_HRVF``
    and ``obj.Gen_w_HRVF`` (shape ``(len(z), len(Gen_x_HRVF), n_frames)``).
    """
    with tictoc("Increasing resolution"):
        L = math.ceil(L0 * obj.delta / obj.Delta_x)
        Numrator = (obj.Gen_u_LRVF.shape[1] - 1) // L
        obj.Kr = factor

        # Pre-calculate coordinate vectors
        Gen_x_HRVF = np.arange(L * obj.Kr + 1) * (obj.Delta_x / obj.Kr)
        Gen_x_LRVF = np.arange(L + 1) * obj.Delta_x
        obj.Gen_x_HRVF = Gen_x_HRVF
        obj.Gen_x_LRVF = Gen_x_LRVF

        LRVF_u = obj.Gen_u_LRVF
        LRVF_w = obj.Gen_w_LRVF
        num_rows = np.asarray(obj.z).size
        num_cols_HR = Gen_x_HRVF.size

        print("Increasing Resolution...")

        k_start = mround(2 * np.pi / (10 * obj.lambda_T))
        k_mid = mround(2 * np.pi / (5 * obj.lambda_T))
        k_end = mround(2 * np.pi / (2 * obj.lambda_T))
        min_val = 0.80

        N = Gen_x_HRVF.size
        dx = Gen_x_HRVF[1] - Gen_x_HRVF[0]
        k = (2 * np.pi / (N * dx)) * np.concatenate(
            (np.arange(0, N // 2 + 1), np.arange(-((N - 1) // 2), 0)))
        k_abs = np.abs(k)

        control_k = np.array([0, k_start, k_mid, k_end, k_abs.max()])
        control_v = np.array([1, 1, min_val, 1, 1])
        order = np.argsort(control_k, kind="stable")
        mask = np.interp(k_abs, control_k[order], control_v[order])

        # To see the filter:
        # eta_Sadd = 0.055e-3  # eta_sadd at z/delta = 0.09
        # plt.semilogx(k_abs*eta_Sadd, mask, 'k-')

        n_drop = int(mround(0.2 * Numrator))
        n_keep = Numrator - n_drop
        HRVF_u = np.zeros((num_rows, num_cols_HR, n_keep))
        HRVF_w = np.zeros((num_rows, num_cols_HR, n_keep))

        def refine(slice_):
            # makima interpolation of every row at once (rows are along axis 1)
            interp = Akima1DInterpolator(Gen_x_LRVF, slice_.T, method="makima",
                                         extrapolate=True)
            fine = interp(Gen_x_HRVF).T
            return np.real(np.fft.ifft(np.fft.fft(fine, axis=1) * mask, axis=1))

        def process_frame(S):
            cols = slice(S * L, (S + 1) * L + 1)
            HRVF_u[:, :, S - n_drop] = refine(LRVF_u[:, cols])
            HRVF_w[:, :, S - n_drop] = refine(LRVF_w[:, cols])

        # The first 20% of the frames are discarded at the end of the MATLAB
        # code; they are skipped here to save memory and time.
        frames = range(n_drop, Numrator)
        # numpy's FFT and scipy's interpolation release the GIL, so in
        # parallel the frames run concurrently on threads. Every frame is
        # written to its own slice, so no locking is needed.
        run_frames(process_frame, ((S,) for S in frames), len(frames), "Increasing Resolution",
                   parallel=parallel, n_workers=n_workers)

        obj.Gen_u_HRVF = HRVF_u
        obj.Gen_w_HRVF = HRVF_w

        print("Done.")
