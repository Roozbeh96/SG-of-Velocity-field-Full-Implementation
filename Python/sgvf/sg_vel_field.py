# SPDX-License-Identifier: CC-BY-NC-SA-4.0
# Copyright (c) 2026 ehsan010@umn.edu

import numpy as np

from ._utils import progressbar, tictoc
from .adding_shear_layer import adding_shear_layer


def sg_vel_field(obj, M, seed=None):
    """Reorganize the generated profiles into a correlated velocity field (``SG_VelField.m``).

    Starting from a random profile, the next profile of the field is the one
    in the storage whose fluctuation u' = u - mean(u) has the largest
    zero-lag normalized cross-correlation with the previous one. The storage
    is refilled with a randomly selected profile from the repository after
    every selection (sampling without replacement).

    Parameters
    ----------
    obj : StochasticGeneration
    M : int
        Size of the storage (Reynolds number dependent).
    seed : int or None
        Seed of the random number generator. ``None`` gives a different field
        every run, like ``rng('shuffle')`` in the MATLAB version.

    Returns
    -------
    Gen_u_prof_reorg, Gen_w_prof_reorg : ndarray, shape (len(z), N_prof - M)
    hist_corr : ndarray, shape (N_prof - M - 1,)
        Cross-correlation of each profile with the next profile.
    """
    # If you do not want to smooth the shear layer in the generated
    # profiles, please comment this function out.
    adding_shear_layer(obj)

    with tictoc("SGVF"):
        u = obj.Gen_u_prof
        w = obj.Gen_w_prof
        # Size of storage is Reynold's number dependent.
        NP_storage = int(M)
        n_prof = u.shape[1]
        size_reorg_VField = n_prof - NP_storage
        if size_reorg_VField < 2:
            raise ValueError("N_prof must be larger than the storage size M + 1.")

        mean_u = u.mean(axis=1)

        def unit_prime(cols):
            """Normalized u' of the given profiles, one profile per row."""
            up = u[:, cols] - mean_u[:, None]
            norms = np.sqrt(np.sum(up ** 2, axis=0))
            norms[norms == 0] = np.inf
            return (up / norms).T

        rng = np.random.default_rng(seed)

        # Since we want to sample W/O replacement, the profiles that are not in
        # the storage form the repository. Profiles are tracked by their column
        # index so the (large) profile matrices are not copied.
        perm = rng.permutation(n_prof)
        storage_idx = perm[:NP_storage].copy()
        repo = perm[NP_storage:].copy()
        n_repo = repo.size
        # Sampled W/O replacement from repo, transferred to storage.
        storage = np.ascontiguousarray(unit_prime(storage_idx))

        order = np.empty(size_reorg_VField, dtype=np.int64)

        def refill(slot):
            """Replace the selected storage profile with a random profile from the repo."""
            nonlocal n_repo
            j = rng.integers(n_repo)
            p = repo[j]
            repo[j] = repo[n_repo - 1]
            n_repo -= 1
            storage_idx[slot] = p
            storage[slot] = unit_prime(np.array([p]))[0]

        # sample first profile for the velocity field.
        r = rng.integers(NP_storage)
        order[0] = storage_idx[r]
        ref_signal = storage[r].copy()
        # Remove the selected profile from the storage and refill it with a
        # randomly selected profile from repo.
        refill(r)

        hist_corr = np.empty(size_reorg_VField - 1)

        for prof_num_Vfield in progressbar(range(1, size_reorg_VField), "Generating Velocity Field"):
            # Normalized correlation (lag 0) of the previous profile with all
            # profiles in the storage at once.
            C = storage @ ref_signal

            # Find the best match index
            best_idx = int(np.argmax(C))

            # We can save the history of the correlation coefficient.
            hist_corr[prof_num_Vfield - 1] = C[best_idx]

            # Filling the VField w/ the most correlated profile.
            order[prof_num_Vfield] = storage_idx[best_idx]
            ref_signal = storage[best_idx].copy()

            # Remove the selected profile from the storage and refill the
            # storage with a randomly selected profile from repo.
            refill(best_idx)

        Gen_u_prof_reorg = np.ascontiguousarray(u[:, order])
        Gen_w_prof_reorg = np.ascontiguousarray(w[:, order])

    return Gen_u_prof_reorg, Gen_w_prof_reorg, hist_corr
