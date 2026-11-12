# SPDX-License-Identifier: CC-BY-NC-SA-4.0
# Copyright (c) 2026 ehsan010@umn.edu

"""
Plot Figures 5, 6 and 7 of the paper for a generated field saved by main.py,
compared with the Saddoughi data (the figures of plots.m).

Two ways to use this file:

1. Interactive (like MATLAB sections): in VS Code click "Run Cell" above a
   ``# %%`` line (or press Shift+Enter). Run the first cell once; the loaded
   object stays in memory. Then run any figure cell as often as you like.
   Changes saved in ``sgvf/plots.py`` are picked up automatically, so you
   can edit a figure and rerun its cell without reading the file again.

2. As a script: ``python run_plots.py`` or "plots" in the Run and Debug panel
   draws every figure, saves them to SAVE_DIR and opens the windows.
"""

# %% Setup and load (slow, run once per session)
from pathlib import Path

import matplotlib.pyplot as plt
from IPython import get_ipython

from sgvf import StochasticGeneration, read_saddoughi
from sgvf import plots

# get_ipython() returns the running IPython/Jupyter session (VS Code's
# Interactive Window), or None when the file is run as a normal script.
_ipython = get_ipython()
INTERACTIVE = _ipython is not None
if INTERACTIVE and not any(e.endswith("autoreload") for e in _ipython.extension_manager.loaded):
    # Reload sgvf automatically when its files are saved.
    _ipython.run_line_magic("load_ext", "autoreload")
    _ipython.run_line_magic("autoreload", "2")

try:
    HERE = Path(__file__).resolve().parent
except NameError:  # the interactive window runs in the folder of this file
    HERE = Path.cwd()

# Generated object saved by main.py
PKL_FILE = HERE / "Gen_sample.pkl"
# Folder where the figures are saved as PNG files (None: do not save)
SAVE_DIR = HERE / "figures"

# Reading the file is skipped if it is already in memory, so rerunning this
# cell is cheap. Delete Gen_sample (del Gen_sample) to force a new read.
if "Gen_sample" not in globals():
    print(f"Loading {PKL_FILE} ...")
    Gen_sample = StochasticGeneration.load(PKL_FILE)
sadd = read_saddoughi()

if not INTERACTIVE:
    # Running as a script: draw, save and show every figure once and stop.
    figures = plots.plot_all(Gen_sample, sadd, save_dir=SAVE_DIR, show=True)
    print(f"{len(figures)} figures saved to {SAVE_DIR}")
    raise SystemExit

# %% Mean velocity and second moments (Figure 5 of the paper)
fig = plots.plot_moments(Gen_sample, sadd)
fig.savefig(SAVE_DIR / "plot_moments.png", dpi=600)

# %% Spectra and structure function (Figure 6 of the paper)
fig = plots.plot_spectra(Gen_sample, sadd)
fig.savefig(SAVE_DIR / "plot_spectra.png", dpi=600)

# %% Co-spectrum E12 (Figure 7 of the paper)
fig = plots.plot_E12(Gen_sample, sadd)
fig.savefig(SAVE_DIR / "plot_E12.png", dpi=600)

# %% Compensated co-spectrum: a -7/3 range is a plateau (C = 0.15)
fig = plots.plot_E12_compensated(Gen_sample, sadd, fit_range=(1, 10))
fig.savefig(SAVE_DIR / "plot_E12_compensated.png", dpi=600)

# %% Co-spectrum E12 with a power-law fit of the smoothed line (slope vs -7/3)
fig = plots.plot_E12(Gen_sample, sadd, fit_range=(5.5, 22), sadd_fit_range=(1, float("inf")))
fig.savefig(SAVE_DIR / "plot_E12_fitt.png", dpi=600)

# %% Distribution of u'w' with elevation
fig = plots.plot_uw_distribution(Gen_sample)
fig.savefig(SAVE_DIR / "plot_uw_distribution.png", dpi=600)

# %% In-plane divergence du/dx + dw/dz versus elevation
fig = plots.plot_divergence(Gen_sample)
fig.savefig(SAVE_DIR / "plot_divergence.png", dpi=600)

# %% Save all figures to SAVE_DIR
figures = plots.plot_all(Gen_sample, sadd, save_dir=SAVE_DIR, dpi=600)
plt.close("all")  # saved to files; do not display all of them here
print(f"{len(figures)} figures saved to {SAVE_DIR}")
