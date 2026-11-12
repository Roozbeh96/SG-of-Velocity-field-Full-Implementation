# SPDX-License-Identifier: CC-BY-NC-SA-4.0
# Copyright (c) 2026 ehsan010@umn.edu

import numpy as np
from scipy.special import ndtr


def gen_gauss_cop(rho, n, repro_numb=None):
    """Uniform random numbers with a prescribed Pearson correlation (Gaussian copula).

    Python version of ``genGaussCop.m``.

    Parameters
    ----------
    rho : float
        Target Pearson correlation coefficient between the two uniform samples.
    n : int
        Number of samples.
    repro_numb : int or None
        Seed for reproducibility. Pass ``None`` to generate different numbers
        every time.
    """
    rng = np.random.default_rng(repro_numb)
    # Convert target Pearson correlation for Uniforms to underlying Gaussian
    # correlation
    r = 2 * np.sin((np.pi / 6) * rho)
    P = np.array([[1.0, r], [r, 1.0]])  # Correlation matrix
    d = P.shape[0]  # Dimension
    # Generate sample
    rands = rng.multivariate_normal(np.zeros(d), P, size=int(n))
    # rands have Gaussian distribution. We need to uniform it.
    rands = ndtr(rands)
    # Now random numbers are uniformly distributed. You can double check with
    # a histogram of rands[:, 0] and rands[:, 1].
    return rands[:, 0], rands[:, 1]
