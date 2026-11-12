# SPDX-License-Identifier: CC-BY-NC-SA-4.0
# Copyright (c) 2026 ehsan010@umn.edu

"""Energy spectra and co-spectra of the generated field.

Python versions of ``spec_analysis_E11.m``, ``spec_analysis_E22.m`` and
``spec_analysis_E12.m``. More info about how to compute the spectra in:
https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/
article/stochastic-modal-velocity-field-in-roughwall-turbulence/
800D97E60609F4FE4D92847CCDB1C7A0
"""

import math

import numpy as np
from scipy.integrate import trapezoid
from scipy.signal import csd, welch
from scipy.signal.windows import hann

from ._utils import concat_fluctuation, frame_cov, frame_mean, frame_var, progressbar


def _welch_settings(obj):
    dx = obj.Gen_x_HRVF[1] - obj.Gen_x_HRVF[0]
    ksmax = math.floor(2 * np.pi / dx)
    win15del = math.floor(15 * obj.delta / dx)
    window = hann(win15del, sym=True)
    NFFT = 2 ** math.ceil(math.log2(win15del))
    return dict(fs=ksmax, window=window, nperseg=win15del, noverlap=NFFT // 2,
                nfft=win15del, detrend=False, return_onesided=True, scaling="density")


def _pwelch_rows(obj, field, desc):
    """pwelch of every row of ``reshape(field - mean(field,3), nz, [])``."""
    settings = _welch_settings(obj)
    mean_rc = frame_mean(field)
    psd = None
    k = None
    for r in progressbar(range(field.shape[0]), desc):
        k, p = welch(concat_fluctuation(field, r, mean_rc), **settings)
        if psd is None:
            psd = np.empty((field.shape[0], p.size))
        psd[r] = p
    return k, psd


def spec_analysis_E11(obj):
    """Streamwise energy spectrum E11(k1); also sets epsilon_E11 and eta_E11."""
    print("Spectral Analysis E11...")
    k, PSD_k = _pwelch_rows(obj, obj.Gen_u_HRVF, "Spectral Analysis E11")  # PSD_k: (nz, nk)
    Apk = trapezoid(PSD_k, k, axis=1)
    Q = frame_var(obj.Gen_u_HRVF).mean(axis=1) / Apk  # [m^3/s^2]=[m^2/s]*[m^2/s^2]*1/[m/s]
    obj.E11 = PSD_k * Q[:, None]

    obj.k1 = k
    obj.epsilon_E11 = 15 * obj.nu * trapezoid(obj.E11 * k ** 2, k, axis=1)
    obj.eta_E11 = (obj.nu ** 3 / obj.epsilon_E11) ** 0.25


def spec_analysis_E22(obj):
    """Vertical energy spectrum E22(k1)."""
    print("Spectral Analysis E22...")
    k, PSD_k = _pwelch_rows(obj, obj.Gen_w_HRVF, "Spectral Analysis E22")
    Apk = trapezoid(PSD_k, k, axis=1)
    Q = frame_var(obj.Gen_w_HRVF).mean(axis=1) / Apk
    obj.E22 = PSD_k * Q[:, None]


def spec_analysis_E12(obj):
    """Co-spectrum E12(k1) of u and w."""
    print("Spectral Analysis E12...")
    settings = _welch_settings(obj)
    mean_u = frame_mean(obj.Gen_u_HRVF)
    mean_w = frame_mean(obj.Gen_w_HRVF)
    nz = obj.Gen_u_HRVF.shape[0]
    CPSD_k = None
    k = None
    for r in progressbar(range(nz), "Spectral Analysis E12"):
        k, p = csd(concat_fluctuation(obj.Gen_u_HRVF, r, mean_u),
                   concat_fluctuation(obj.Gen_w_HRVF, r, mean_w), **settings)
        if CPSD_k is None:
            CPSD_k = np.empty((nz, p.size))
        CPSD_k[r] = np.real(p)
    Apk = trapezoid(CPSD_k, k, axis=1)
    Cov = frame_cov(obj.Gen_u_HRVF, obj.Gen_w_HRVF)
    Q = Cov.mean(axis=1) / Apk
    obj.E12 = CPSD_k * Q[:, None]
