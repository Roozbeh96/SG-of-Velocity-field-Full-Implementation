# SPDX-License-Identifier: CC-BY-NC-SA-4.0
# Copyright (c) 2026 ehsan010@umn.edu

import numpy as np

from ._utils import frame_cov, frame_mean


def production_comp(obj):
    """Turbulent kinetic energy production -<u'w'> dU/dz (``production_comp.m``)."""
    print("Computing Production...")
    dUdz = np.gradient(frame_mean(obj.Gen_u_HRVF).mean(axis=1), obj.z[1] - obj.z[0])
    obj.production = (-frame_cov(obj.Gen_u_HRVF, obj.Gen_w_HRVF)).mean(axis=1) * dUdz

    with np.errstate(invalid="ignore"):
        obj.eta_prod = (obj.nu ** 3 / obj.production) ** 0.25
