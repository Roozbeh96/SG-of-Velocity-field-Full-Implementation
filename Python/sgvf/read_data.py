# SPDX-License-Identifier: CC-BY-NC-SA-4.0
# Copyright (c) 2026 ehsan010@umn.edu

"""Read the Saddoughi wind-tunnel data (``read_data.m``).

The JSON files are shared with the MATLAB package and live in
``Matlab/src/Data/Saddoughi`` of the repository.
"""

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

#: Default location of the Saddoughi data (shared with the MATLAB package).
DATA_DIR = Path(__file__).resolve().parents[2] / "Matlab" / "src" / "Data" / "Saddoughi"


def read_json_xy(path):
    """Return the ``x`` and ``y`` columns of a digitized-data JSON file as float arrays."""
    with open(path) as f:
        data = json.load(f)
    x = np.array([float(d["x"]) for d in data])
    y = np.array([float(d["y"]) for d in data])
    return x, y


def read_saddoughi(data_dir=DATA_DIR):
    """Read the Saddoughi data for both U_e = 10 and 50 [m/s].

    Returns a namespace with the same variable names as ``read_data.m``
    (e.g. ``sadd.k1_eta_10_Sadd``, ``sadd.u_tau_50_Sadd``).
    """
    data_dir = Path(data_dir)
    d10 = data_dir / "Ue=10[ms]"
    d50 = data_dir / "Ue=50[ms]"
    s = SimpleNamespace()

    s.k1_eta_10_Sadd, s.energy_10_Sadd = read_json_xy(d10 / "mid-layer-u-spectra(Ue_10).json")
    s.k1_eta_compen_inner_10_Sadd, s.compen_energy_inner_10_Sadd = read_json_xy(
        d10 / "inner-layer-compen-u-spectra(Ue_10).json")
    s.k1_eta_compen_mid_10_Sadd, s.compen_energy_mid_10_Sadd = read_json_xy(
        d10 / "mid-layer-compen-u-spectra(Ue_10).json")
    s.z_U_o_Ue_10_Sadd, s.U_o_Ue_10_Sadd = read_json_xy(d10 / "profile_meanvelocity.json")
    s.z_uu_o_Ue2_10_Sadd, s.uu_o_Ue2_10_Sadd = read_json_xy(d10 / "profile_uvar.json")
    s.z_ww_o_Ue2_10_Sadd, s.ww_o_Ue2_10_Sadd = read_json_xy(d10 / "profile_wvar.json")
    s.z_uw_o_Ue2_10_Sadd, s.uw_o_Ue2_10_Sadd = read_json_xy(d10 / "profile_uwcovar.json")
    s.r_o_eta_10_Sadd, s.compen_D11_10_Sadd = read_json_xy(d10 / "mid-layer-D11(Ue_10).json")

    s.Ue_10_Sadd = 10
    s.u_tau_10_Sadd = 0.465
    s.z_0_10_Sadd = 3.8e-4
    s.z_mid_o_delta_10_Sadd = 0.515 / 1.09
    s.z_inner_o_delta_10_Sadd = 0.1 / 1.09
    s.delta_10_Sadd = 1.09
    # ------------------------------
    s.k1_eta_50_Sadd, s.energy_50_Sadd = read_json_xy(d50 / "mid-layer-u-spectra(Ue_50).json")
    s.k1_eps_S_50_Sadd, s.E12_50_Sadd = read_json_xy(d50 / "inner-layer-uw-cospectra(Ue_50).json")
    s.k1_eta_compen_inner_50_Sadd, s.compen_energy_inner_50_Sadd = read_json_xy(
        d50 / "inner-layer-compend-u-spectra(Ue_50).json")
    s.k1_eta_E22_compen_inner_50_Sadd, s.compen_E22_inner_50_Sadd = read_json_xy(
        d50 / "inner-layer-compend-w-spectra(Ue_50).json")
    s.k1_eta_compen_mid_50_Sadd, s.compen_energy_mid_50_Sadd = read_json_xy(
        d50 / "mid-layer-compend-u-spectra(Ue_50).json")
    s.z_U_o_Ue_50_Sadd, s.U_o_Ue_50_Sadd = read_json_xy(d50 / "profile_meanvelocity.json")
    s.z_uu_o_Ue2_50_Sadd, s.uu_o_Ue2_50_Sadd = read_json_xy(d50 / "profile_uvar.json")
    s.z_ww_o_Ue2_50_Sadd, s.ww_o_Ue2_50_Sadd = read_json_xy(d50 / "profile_wvar.json")
    s.z_uw_o_Ue2_50_Sadd, s.uw_o_Ue2_50_Sadd = read_json_xy(d50 / "profile_uwcovar.json")
    s.r_o_eta_50_Sadd, s.compen_D11_50_Sadd = read_json_xy(d50 / "mid-layer-D11(Ue_50).json")

    s.Ue_50_Sadd = 50
    s.u_tau_50_Sadd = 2.325
    s.z_0_50_Sadd = 3.8e-4
    s.z_mid_o_delta_50_Sadd = 0.4 / 1.09
    s.z_inner_o_delta_50_Sadd = 0.1 / 1.09
    s.delta_50_Sadd = 1.09
    return s
