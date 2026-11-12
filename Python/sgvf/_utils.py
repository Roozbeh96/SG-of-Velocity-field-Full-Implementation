# SPDX-License-Identifier: CC-BY-NC-SA-4.0
# Copyright (c) 2026 ehsan010@umn.edu

"""Small helpers that reproduce MATLAB behaviour used throughout the package."""

import time
from contextlib import contextmanager

import numpy as np
from scipy.signal import convolve2d
from tqdm import tqdm


def mround(x):
    """MATLAB ``round``: halves are rounded away from zero (numpy rounds to even)."""
    x = np.asarray(x, dtype=float)
    out = np.sign(x) * np.floor(np.abs(x) + 0.5)
    return out if out.ndim else float(out)


def colon(start, step, stop):
    """MATLAB ``start:step:stop`` for floating point values."""
    n = int(np.floor((stop - start) / step + 1e-10)) + 1
    return start + step * np.arange(n)


def first_index(mask):
    """MATLAB ``find(mask, 1, 'first')`` (0-based). Returns ``None`` if empty."""
    idx = np.flatnonzero(mask)
    return int(idx[0]) if idx.size else None


def log_indices(n, num_points, start=0.0):
    """0-based version of ``unique(round(logspace(start, log10(n), num_points)))``."""
    return np.unique(mround(np.logspace(start, np.log10(n), num_points))).astype(int) - 1


def movmean(x, kb, kf):
    """MATLAB ``movmean(x, [kb kf])`` (window shrinks at the end points)."""
    x = np.asarray(x, dtype=float)
    c = np.concatenate(([0.0], np.cumsum(x)))
    i = np.arange(x.size)
    lo = np.maximum(i - kb, 0)
    hi = np.minimum(i + kf + 1, x.size)
    return (c[hi] - c[lo]) / (hi - lo)


def conv2_same(a, kernel):
    """MATLAB ``conv2(a, kernel, 'same')``.

    scipy centres even-sized kernels differently, so the central part of the
    full convolution is cut out the way MATLAB does it.
    """
    full = convolve2d(a, kernel, mode="full")
    r0 = kernel.shape[0] // 2
    c0 = kernel.shape[1] // 2
    return full[r0:r0 + a.shape[0], c0:c0 + a.shape[1]]


def ecdf(y):
    """MATLAB ``[f, x] = ecdf(y)``: the first point is (min(y), 0)."""
    y = np.sort(np.asarray(y, dtype=float).ravel())
    vals, counts = np.unique(y, return_counts=True)
    f = np.concatenate(([0.0], np.cumsum(counts) / y.size))
    x = np.concatenate(([vals[0]], vals))
    return f, x


def label_columnwise(binary):
    """``bwlabel(binary, 4)``: 4-connected labels numbered in MATLAB (column-major) order."""
    from scipy.ndimage import label

    labels_t, n = label(binary.T)
    return labels_t.T, n


def std_matlab(x):
    """MATLAB ``std(x, 0)``: N-1 normalisation, and 0 for a single sample."""
    x = np.asarray(x)
    if x.size == 0:
        return 0.0
    if x.size == 1:
        return 0.0
    return float(np.std(x, ddof=1))


def progressbar(iterable, desc, **kwargs):
    """Text replacement for MATLAB's ``progressbar`` figure."""
    return tqdm(iterable, desc=desc, leave=False, **kwargs)


@contextmanager
def tictoc(label=None):
    """MATLAB ``tic``/``toc``."""
    t0 = time.perf_counter()
    yield
    dt = time.perf_counter() - t0
    print(f"{label + ': ' if label else ''}Elapsed time is {dt:.6f} seconds.")


def frame_mean(field):
    """``mean(field, 3)`` computed one row at a time to limit memory use."""
    out = np.empty(field.shape[:2])
    for r in range(field.shape[0]):
        out[r] = field[r].mean(axis=1)
    return out


def frame_var(field):
    """``var(field, 0, 3)`` computed one row at a time to limit memory use."""
    out = np.empty(field.shape[:2])
    for r in range(field.shape[0]):
        out[r] = field[r].var(axis=1, ddof=1)
    return out


def frame_cov(field_a, field_b):
    """``sum(a'.*b', 3)/(n-1)`` with primes being fluctuations over the 3rd dimension."""
    n = field_a.shape[2]
    out = np.empty(field_a.shape[:2])
    for r in range(field_a.shape[0]):
        a = field_a[r] - field_a[r].mean(axis=1, keepdims=True)
        b = field_b[r] - field_b[r].mean(axis=1, keepdims=True)
        out[r] = (a * b).sum(axis=1) / (n - 1)
    return out


def concat_fluctuation(field, r, mean_rc=None):
    """Row ``r`` of ``reshape(field - mean(field,3), nz, [])`` (MATLAB column-major order)."""
    row = field[r]
    m = row.mean(axis=1, keepdims=True) if mean_rc is None else mean_rc[r][:, None]
    return (row - m).ravel(order="F")


def run_frames(func, args_iter, n_tasks, desc, parallel=False, n_workers=None,
               processes=False):
    """Apply ``func(*args)`` to every item of ``args_iter`` and return the results in order.

    ``parallel=False`` runs the items one after the other. ``parallel=True``
    runs them on a pool of threads, or on a pool of processes when
    ``processes=True`` (for work that holds Python's GIL; ``func`` must then
    be a module-level function). At most ``2*n_workers`` items are in flight
    at a time so that the frames are not all copied at once.
    """
    import os
    from collections import deque
    from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

    if not parallel:
        return [func(*args) for args in progressbar(args_iter, desc, total=n_tasks)]

    n_workers = n_workers or os.cpu_count()
    print(f"Running on {n_workers} {'processes' if processes else 'threads'}")
    pool_cls = ProcessPoolExecutor if processes else ThreadPoolExecutor
    results = []
    pending = deque()
    with pool_cls(max_workers=n_workers) as pool, tqdm(total=n_tasks, desc=desc, leave=False) as bar:
        for args in args_iter:
            pending.append(pool.submit(func, *args))
            if len(pending) >= 2 * n_workers:
                results.append(pending.popleft().result())
                bar.update()
        while pending:
            results.append(pending.popleft().result())
            bar.update()
    return results
