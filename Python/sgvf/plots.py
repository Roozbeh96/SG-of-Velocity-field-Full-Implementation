# SPDX-License-Identifier: CC-BY-NC-SA-4.0
# Copyright (c) 2026 ehsan010@umn.edu

"""Figures of the paper (Python version of ``plots.m``), drawn by ``run_plots.py``.

- :func:`plot_moments`: mean velocity and second moments (Figure 5).
- :func:`plot_spectra`: E11 spectra and D11 with the premultiplied E11 inset
  (Figure 6).
- :func:`plot_E12`: co-spectrum E12 (Figure 7), optionally with a power-law fit.
- :func:`plot_E12_compensated`: E12 compensated with the -7/3 power law.
- :func:`plot_uw_distribution`: distribution of u'w' with elevation.
- :func:`plot_divergence`: normalized in-plane divergence du/dx + dw/dz with elevation.

The generated field is compared with the Saddoughi data (U_inf = 50 [m/s]),
read with :func:`sgvf.read_data.read_saddoughi`. :func:`plot_all` draws all
figures. Indices are 0-based.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch
from matplotlib.ticker import LogLocator, NullFormatter
from matplotlib.transforms import blended_transform_factory

from ._utils import first_index, frame_cov, frame_mean, frame_var, log_indices, movmean
from .read_data import read_saddoughi

KAPPA = 0.39
RED, BLUE, BLACK, GREEN = (1.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 0.0, 0.0), (0,1,0)
GOLD = (0.93, 0.69, 0.13)
CYAN = (0.0, 1.0, 1.0)
LIGHT_PURPLE = (0.75, 0.45, 1.0)
LIGHT_GOLD = (1.00, 0.74, 0.17)
#: Light face colour of the Saddoughi markers.
SADD_FACE = (1.00, 0.81, 0.81)

# Computer Modern for all text, as MATLAB's LaTeX interpreter: math in the cm
# math fonts, plain text in cmr10 (shipped with matplotlib). cmr10 has no minus
# sign, so tick labels are written as math.
plt.rcParams.update({
    "mathtext.fontset": "cm",
    "font.family": "serif",
    "font.serif": ["cmr10"],
    "axes.formatter.use_mathtext": True,
    "axes.unicode_minus": False,
})


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------
def _figure(width, height):
    """Figure sized from a MATLAB ``Position`` (pixels)."""
    return plt.figure(figsize=(width / 96, height / 96))


def _grid(ax, fontsize=13):
    ax.grid(True, which="major", alpha=0.35)
    ax.grid(True, which="minor", alpha=0.15)
    ax.minorticks_on()
    ax.tick_params(labelsize=fontsize)


def _legend(fig, ax, position, fontsize=10, frame=False):
    """Legend placed at a MATLAB normalized figure ``Position``."""
    x, y, w, h = position
    leg = ax.legend(loc="upper left", bbox_to_anchor=(x, y + h), bbox_transform=fig.transFigure,
                    fontsize=fontsize, frameon=frame)
    return leg


def _logx(x):
    """x for a log axis: values <= 0 (e.g. k1 = 0) become NaN and are not drawn,
    as MATLAB does on logarithmic axes (matplotlib would clip them instead)."""
    x = np.asarray(x, dtype=float)
    return np.where(x > 0, x, np.nan)


def _ind(obj, value):
    return first_index(np.asarray(obj.z) / obj.delta >= value)


def _textbox(fig, x, y, w, h, text, fontsize=12, color="k"):
    """MATLAB ``annotation('textbox', [x y w h])``: text at the top left of the box."""
    fig.text(x + 0.005, y + h, text, fontsize=fontsize, color=color, va="top", ha="left")


def _arrow(fig, x, y, dx, dy, both=False, color="k", lw=1.0, head=8):
    """MATLAB ``annotation('arrow'/'doublearrow', [x y dx dy])`` in figure coordinates."""
    style = "<|-|>" if both else "-|>"
    fig.add_artist(FancyArrowPatch((x, y), (x + dx, y + dy), transform=fig.transFigure,
                                   arrowstyle=style, mutation_scale=head, color=color, lw=lw))


def _eps_eta_sadd(obj):
    epsilon_Sadd = obj.u_tau ** 3 / (KAPPA * np.asarray(obj.z))
    eta_Sadd = (obj.nu ** 3 / epsilon_Sadd) ** 0.25
    return epsilon_Sadd, eta_Sadd


# ----------------------------------------------------------------------------
# figures
# ----------------------------------------------------------------------------
def plot_moments(obj, sadd):
    """Mean velocity and second moments of the HRVF compared with Saddoughi
    (U_inf = 50 [m/s]): (a) mean velocity, (b) streamwise variance versus
    z/delta, (c) vertical variance, (d) Reynolds shear stress (Figure 5 of the
    paper)."""
    z = np.asarray(obj.z)
    u2 = obj.u_tau ** 2
    li = log_indices(z.size, 30)
    li_b = log_indices(z.size - 50, 30)
    U = frame_mean(obj.Gen_u_HRVF).mean(axis=1)
    uu = frame_var(obj.Gen_u_HRVF).mean(axis=1)
    ww = frame_var(obj.Gen_w_HRVF).mean(axis=1)
    uw = -frame_cov(obj.Gen_u_HRVF, obj.Gen_w_HRVF).mean(axis=1)
    s2 = sadd.Ue_50_Sadd ** 2 / sadd.u_tau_50_Sadd ** 2
    sadd_label = r"Saddoughi($U_{\infty}=50$[m/s])"
    gen_label = r"Gen($U_{\infty}=50$[m/s])"

    fig = _figure(798, 902)
    positions = [[0.0794, 0.7761, 0.8997, 0.2139], [0.0802, 0.4992, 0.8990, 0.2095],
                 [0.0789, 0.2794, 0.9002, 0.1682], [0.0815, 0.0401, 0.8977, 0.1927]]
    axs = [fig.add_axes(p) for p in positions]

    def sadd_points(ax, x, y):
        ax.plot(x, y, linestyle="none", marker="o", markersize=8, markeredgecolor="k",
                markerfacecolor=SADD_FACE, markeredgewidth=1.5, label=sadd_label)

    def gen_line(ax, x, y):
        ax.plot(x, y, color=RED, lw=2, marker="x", markersize=8, markeredgecolor="k",
                markeredgewidth=2, label=gen_label)

    # (a) Mean velocity
    ax = axs[0]
    zz = np.arange(30, 3001)
    ax.semilogy(1 / KAPPA * np.log(zz), zz, color=GOLD, lw=2.5, label="_nolegend_")
    sadd_points(ax, sadd.U_o_Ue_50_Sadd * sadd.Ue_50_Sadd / sadd.u_tau_50_Sadd,
                sadd.z_U_o_Ue_50_Sadd / sadd.z_0_50_Sadd)
    gen_line(ax, U[li] / obj.u_tau, z[li] / obj.z_0)
    ax.set_xlim(8, 22)
    ax.set_ylim(5, 2e4)
    ax.set_yticks([1e1, 1e2, 1e3, 1e4])
    ax.set_xlabel(r"$\overline{U}/u_{\tau}$", fontsize=15)
    ax.set_ylabel(r"$z/z_{0}$", fontsize=15)
    _legend(fig, ax, [0.1621, 0.9268, 0.2818, 0.0524], fontsize=12, frame=True)

    # (b) Streamwise variance versus z/delta
    ax = axs[1]
    zz = np.logspace(-2, 0, 20)
    ax.plot(zz, -1.26 * np.log(zz) + 1.0, color=GOLD, lw=2.5, label=r"$1.00-1.26\ln(z/\delta)$")
    sadd_points(ax, sadd.z_uu_o_Ue2_50_Sadd / sadd.delta_50_Sadd, sadd.uu_o_Ue2_50_Sadd * s2)
    gen_line(ax, z[li_b] / obj.delta, uu[li_b] / u2)
    ax.set_xscale("log")
    ax.set_xlim(1e-3, 1.335)
    ax.set_ylim(0, 8)
    ax.set_xticks([1e-3, 1e-2, 1e-1, 1e0])
    ax.set_xlabel(r"$z/\delta$", fontsize=15)
    ax.set_ylabel(r"$\overline{u^{\prime 2}}/u_{\tau}^{2}$", fontsize=15)

    # (c) Vertical variance
    ax = axs[2]
    sadd_points(ax, sadd.z_ww_o_Ue2_50_Sadd / sadd.z_0_50_Sadd, sadd.ww_o_Ue2_50_Sadd * s2)
    gen_line(ax, z[li] / obj.z_0, ww[li] / u2)
    ax.set_ylim(0, 2)
    ax.set_ylabel(r"$\overline{w^{\prime 2}}/u_{\tau}^{2}$", fontsize=15)

    # (d) Reynolds shear stress
    ax = axs[3]
    sadd_points(ax, sadd.z_uw_o_Ue2_50_Sadd / sadd.z_0_50_Sadd, -sadd.uw_o_Ue2_50_Sadd * s2)
    gen_line(ax, z[li] / obj.z_0, uw[li] / u2)
    ax.set_ylim(-0.5, 1.6)
    ax.set_yticks([-0.5, 0, 0.5, 1, 1.5])
    ax.set_ylabel(r"$-\overline{u^{\prime}w^{\prime}}/u_{\tau}^{2}$", fontsize=15)

    for ax in axs[2:]:
        ax.set_xscale("log")
        ax.set_xlim(10, 5000)
        ax.set_xticks([1e1, 1e2, 1e3])
        ax.set_xlabel(r"$z/z_{0}$", fontsize=15)
    axs[3].set_yticklabels(["-0.5", "0", "0.5", "1", "1.5"])
    for ax in axs:
        _grid(ax, 15)
        # y-labels at the same distance from the figure edge, as in the paper
        ax.yaxis.set_label_coords(0.04, 0.5, transform=blended_transform_factory(
            fig.transFigure, ax.transAxes))
    for ax in axs[1:]:
        # x-labels on the tick-label row, between the ticks (no room below)
        ax.xaxis.labelpad = -20
    for label, (x, y) in zip("abcd", [(0.0847, 0.9481), (0.0759, 0.6665),
                                      (0.0772, 0.4049), (0.0847, 0.1898)]):
        _textbox(fig, x, y - 0.006, 0.0522, 0.0458, f"({label})", fontsize=15)
    return fig


def plot_spectra(obj, sadd, num_points=50, num_points_dots=1000, num_points_dots_inset=500):
    """Spectra and structure function of the HRVF + VorX (Figure 6 of the paper).

    (a) compensated E11 at z/delta = 0.09, (b) E11 in Kolmogorov scaling at
    z/delta = 0.14, (c) D11 at z/delta = 0.05, 0.09 and 0.14 with the
    premultiplied E11 at the same elevations as an inset. The lines are
    sampled at ``num_points`` log-spaced wavenumbers. The noisy spectra are
    drawn as dots (at ``num_points_dots`` wavenumbers in (a) and
    ``num_points_dots_inset`` in the inset) with their moving average as a
    line: ``movmean([2 2])`` in (a) and ``movmean([1 1])`` in the inset; (b)
    is not smoothed.
    """
    z = np.asarray(obj.z)
    k1 = obj.k1
    eps, eta = _eps_eta_sadd(obj)
    i05, i09, i14 = (_ind(obj, v) for v in (0.05, 0.0919, 0.14))
    li = log_indices(k1.size, num_points)
    li_dots = log_indices(k1.size, num_points_dots)
    li_dots_inset = log_indices(k1.size, num_points_dots_inset)

    def gen_label(i):
        return rf"Gen($U_{{\infty}}=50$[m/s]) HRVF + VorX $z/\delta$ = {z[i] / obj.delta:.2f}"

    def sadd_label(zod):
        return rf"Saddoughi($\mathrm{{U}}_{{\infty}}=50$[m/s]) z$/\delta$ = {zod:.2f}"

    def noisy(ax, x_dots, y_dots, x, y, mov, color, ls, label):
        """Raw values as dots and the moving average movmean(y, [mov mov]) as a line."""
        ax.plot(_logx(x_dots), y_dots, linestyle="none", marker=".", markersize=2, color=color,
                label="_nolegend_")
        ax.plot(_logx(x), movmean(y, mov, mov), linestyle=ls, lw=2, color=color, label=label)

    def legend(ax, position):
        x, y, w, h = position
        ax.legend(loc="center", bbox_to_anchor=(x + w / 2, y + h / 2), bbox_transform=fig.transFigure,
                  fontsize=10, frameon=True, handlelength=2.5, borderpad=0.4, labelspacing=0.3,
                  fancybox=False, edgecolor="k")

    def ylabel_at(ax, x_fig, y_axes=0.5):
        """Place the y-label at a fixed figure x (right edge of the rotated text)."""
        ax.yaxis.set_label_coords(x_fig, y_axes, transform=blended_transform_factory(
            fig.transFigure, ax.transAxes))

    fig = _figure(984, 620)

    # (a) Compensated E11 at z/delta = 0.09
    ax = fig.add_axes([0.0681, 0.6226, 0.5504, 0.3538])
    ax.semilogx([1e-5, 1e5], [1.5 * 18 / 55] * 2, "k-.", lw=2, label="_nolegend_")
    ax.plot(sadd.k1_eta_compen_inner_50_Sadd, sadd.compen_energy_inner_50_Sadd, linestyle="none",
            marker="v", markersize=9, markeredgecolor="k", markerfacecolor=SADD_FACE,
            markeredgewidth=1.5, label=sadd_label(sadd.z_inner_o_delta_50_Sadd))
    # dots at num_points_dots wavenumbers, smoothed line at num_points
    ax.plot(_logx(k1[li_dots] * 0.055e-3), obj.E11[i09, li_dots] * k1[li_dots] ** (5 / 3) * 342 ** (-2 / 3),
            linestyle="none", marker=".", markersize=2, color=RED, label="_nolegend_")
    ax.plot(_logx(k1[li] * 0.055e-3), movmean(obj.E11[i09, li] * k1[li] ** (5 / 3) * 342 ** (-2 / 3), 2, 2),
            "-", lw=2, color=RED, marker="x", markersize=8, markeredgecolor="k",
            markeredgewidth=2, label=gen_label(i09))
    ax.set_xlim(1e-5, 1)
    ax.set_ylim(0, 1.5)
    ax.set_yticks([0, 0.5, 1, 1.5])
    ax.set_yticklabels(["0", "0.5", "1", "1.5"])
    ax.set_xlabel(r"$k_{1}\eta$", fontsize=15)
    ax.set_ylabel(r"$\varepsilon^{-2/3}k_{1}^{5/3}E_{11}(k_{1})$", fontsize=15)
    _grid(ax, 15)
    legend(ax, [0.0785, 0.8694, 0.3059, 0.0984])
    ylabel_at(ax, 0.032)

    # (b) E11 in Kolmogorov scaling at z/delta = 0.14
    ax = fig.add_axes([0.6694, 0.0978, 0.3203, 0.8892])
    ax.set_zorder(3)
    kk = np.arange(3e-4, 1 + 1e-12, 1e-5)
    ax.loglog(kk, kk ** (-5 / 3), lw=2, color=LIGHT_GOLD, label="_nolegend_")
    ax.loglog([5e-7, 5], [1, 1], "k--", lw=1.5, label="_nolegend_")
    ax.loglog([1, 1], [1e-8, 10], "k--", lw=1.5, label="_nolegend_")
    ax.plot(sadd.k1_eta_50_Sadd, sadd.energy_50_Sadd, linestyle="none", marker="o", markersize=9,
            markeredgecolor="k", markerfacecolor=SADD_FACE, markeredgewidth=1.5,
            label=sadd_label(sadd.z_mid_o_delta_50_Sadd))
    ax.loglog(_logx(k1[li] * eta[i14]), obj.E11[i14, li] / (eps[i14] * obj.nu ** 5) ** 0.25, "-",
              lw=2, color=BLUE, marker="x", markersize=8, markeredgecolor="k",
              markeredgewidth=2, label=gen_label(i14))
    ax.set_ylim(1e-4, 3.2e6)
    ax.set_yticks([1e0, 1e6])
    ax.set_xticks([1e-4, 1e0])
    ax.set_xlabel(r"$k_{1}\eta$", fontsize=15)
    ax.set_ylabel(r"$E_{11}(\mathrm{k}_{1})/(\varepsilon \nu^5)^{1/4}$", fontsize=15)
    _grid(ax, 15)
    for axis in (ax.xaxis, ax.yaxis):  # after _grid: minorticks_on resets the locators
        axis.set_minor_locator(LogLocator(base=10, subs=np.arange(1, 10), numticks=100))
        axis.set_minor_formatter(NullFormatter())
    _textbox(fig, 0.7993, 0.8618, 0.0847, 0.0802, r"$(k_{1}\eta)^{-5/3}$", fontsize=15)
    legend(ax, [0.6774, 0.1150, 0.3059, 0.0984])
    ylabel_at(ax, 0.66, 0.6)  # above the 10^0 tick, as in the paper

    # (c) D11 at z/delta = 0.05, 0.09 and 0.14
    ax = fig.add_axes([0.0681, 0.0926, 0.5504, 0.4283])
    D = obj.D11.mean(axis=2)
    styles = [(i05, BLACK, ":"), (i09, RED, "-"), (i14, BLUE, "--")]
    for i, c, ls in styles:
        ax.semilogx(obj.r11 / z[i], D[i] / obj.u_tau ** 2, ls, lw=2.5, color=c, label=gen_label(i))
    xx = np.logspace(0.05, 0.5, 25)
    ax.plot(xx, 2 * 1.26 * np.log(xx) + 4.2, "-", lw=2.5, color=GOLD, label="_nolegend_")
    ax.set_xlim(0.1, 20)
    ax.set_ylim(0, 14)
    ax.set_yticks(np.arange(0, 15, 2))
    ax.set_xlabel(r"$r/z$", fontsize=15)
    ax.set_ylabel(r"$D_{11}(r)/u_{\tau}^2$", fontsize=15)
    _grid(ax, 15)
    legend(ax, [0.3074, 0.1083, 0.3059, 0.0888])
    ylabel_at(ax, 0.032)
    _arrow(fig, 0.5697, 0.4565, 0.0003, -0.1565)
    fig.text(0.5697, 0.4565, r"Increasing $z$", fontsize=15, ha="center", va="bottom")

    # inset: premultiplied E11 at the same elevations
    axi = fig.add_axes([0.1011, 0.3199, 0.2905, 0.1630])
    for i, c, ls in [(i05, BLACK, ":"), (i09, RED, "-"), (i14, BLUE, "-.")]:
        # dots at num_points_dots_inset wavenumbers, smoothed line at num_points
        d = li_dots_inset
        noisy(axi, k1[d] * z[i], obj.E11[i, d] * k1[d] / obj.u_tau ** 2,
              k1[li] * z[i], obj.E11[i, li] * k1[li] / obj.u_tau ** 2, 1, c, ls, "_nolegend_")
    axi.set_xscale("log")
    axi.set_xlim(1e-2, 1e2)
    axi.set_ylim(0, 2)
    axi.set_yticks([0, 1, 2])
    axi.xaxis.tick_top()
    axi.yaxis.tick_right()
    axi.yaxis.set_label_position("left")
    axi.set_xlabel(r"$k_{1}z$", fontsize=13, labelpad=2)
    axi.set_ylabel(r"$k_{1}E_{11}(k_{1})/u_{\tau}^{2}$", fontsize=13, labelpad=2)
    _grid(axi, 13)
    _arrow(fig, 0.1363, 0.42, 0.0468, -0.0748)
    fig.text(0.1363, 0.42, r"Increasing $z$", fontsize=13, ha="center", va="bottom")

    for label, (x, y) in zip("acb", [(-0.0032, 0.9403), (0.0080, 0.5225), (0.9341, 0.9131)]):
        _textbox(fig, x, y, 0.0414, 0.0630, f"({label})", fontsize=15)
    return fig


def fit_power_law(x, y, fit_range):
    """Least-squares line log10(y) = n*log10(x) + b through the points with
    fit_range[0] <= x <= fit_range[1] (non-positive and NaN values skipped).

    Returns a dict with the slope ``n``, its standard error ``n_err``, the
    intercept ``b``, the coefficient of determination ``r2``, the number of
    points ``npts`` and the x interval of the fitted points ``x_min``, ``x_max``.
    """
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    m = (x >= fit_range[0]) & (x <= fit_range[1]) & (x > 0) & (y > 0)
    m &= np.isfinite(x) & np.isfinite(y)
    lx, ly = np.log10(x[m]), np.log10(y[m])
    (n, b), cov = np.polyfit(lx, ly, 1, cov=True)
    r2 = 1 - np.sum((ly - (n * lx + b)) ** 2) / np.sum((ly - ly.mean()) ** 2)
    return dict(n=n, n_err=np.sqrt(cov[0, 0]), b=b, r2=r2, npts=int(m.sum()),
                x_min=x[m].min(), x_max=x[m].max())


def plot_E12(obj, sadd, num_points=50, num_points_dots=1000, fit_range=None, sadd_fit_range=None):
    """Co-spectrum E12 of the HRVF + VorX at z/delta = 0.09 normalized with the
    dissipation and the shear rate (Figure 7 of the paper).

    Uses epsilon = 342 and S = 83.5 at z/delta = 0.09 (Saddoughi). The raw
    values are drawn as dots at ``num_points_dots`` log-spaced wavenumbers and
    their moving average ``movmean([2 2])`` at ``num_points`` wavenumbers as a
    line with markers (both start at ``logspace(0.3, ...)``).

    With ``fit_range = (x_min, x_max)`` a power law is fitted in log-log space
    (:func:`fit_power_law`) to the moving-average line and to the Saddoughi
    data with x_min <= k1 eps^1/2 S^-3/2 <= x_max; the fitted lines are drawn
    dashed and the slopes are printed and compared with -7/3.
    ``sadd_fit_range`` sets a different range for the Saddoughi data
    (default: ``fit_range``).
    """
    k1 = obj.k1
    i09 = _ind(obj, 0.0919)
    li = log_indices(k1.size, num_points, start=0.3)
    li_dots = log_indices(k1.size, num_points_dots, start=0.3)
    eps_, S_ = 342.0, 83.5
    xscale = eps_ ** 0.5 * S_ ** (-3 / 2)
    yscale = (eps_ ** 3 * S_ ** (-5)) ** 0.5

    fig = _figure(577, 386)
    ax = fig.add_axes([0.155, 0.1528, 0.8057, 0.8057])
    kk = np.arange(1e-1, 1e3 + 1e-9, 1e1)
    ax.loglog(kk, kk ** (-7 / 3), lw=2, color=LIGHT_GOLD, linestyle="-", label="_nolegend_")
    ax.plot(sadd.k1_eps_S_50_Sadd, sadd.E12_50_Sadd, linestyle="none", marker="v", markersize=9,
            markeredgecolor="k", markerfacecolor=SADD_FACE, markeredgewidth=1.5,
            label=rf"Saddoughi($U_{{\infty}}=50$[m/s]) $z/\delta$ = {sadd.z_inner_o_delta_50_Sadd:.2f}")
    # dots at num_points_dots wavenumbers, smoothed line at num_points
    ax.plot(_logx(k1[li_dots] * xscale), np.abs(obj.E12[i09, li_dots]) / yscale, linestyle="none",
            marker=".", markersize=2, color=RED, label="_nolegend_")
    x_line = _logx(k1[li] * xscale)
    y_line = movmean(np.abs(obj.E12[i09, li]) / yscale, 2, 2)
    ax.plot(x_line, y_line, "-", lw=2,
            color=RED, marker="x", markersize=8, markeredgecolor="k", markeredgewidth=2,
            label=rf"Gen($U_{{\infty}}=50$[m/s]) $z/\delta$ = {obj.z[i09] / obj.delta:.2f}")
    if fit_range is not None:
        sadd_fit_range = fit_range if sadd_fit_range is None else sadd_fit_range
        print(f"Power-law fit of k1 eps^1/2 S^-3/2 (-7/3 = {-7 / 3:.3f}):")
        for name, xs, ys, c, rng in (
                ("Saddoughi", sadd.k1_eps_S_50_Sadd, sadd.E12_50_Sadd, LIGHT_PURPLE, sadd_fit_range),
                ("Gen", x_line, y_line, CYAN, fit_range)):
            fit = fit_power_law(xs, ys, rng)
            print(f"  {name} in [{rng[0]:g}, {rng[1]:g}]: slope = {fit['n']:.3f} +/- {fit['n_err']:.3f}, R^2 = {fit['r2']:.3f}, "
                  f"{fit['npts']} points in [{fit['x_min']:.3g}, {fit['x_max']:.3g}]")
            xf = np.array([fit["x_min"], fit["x_max"]])
            ax.plot(xf, 10 ** fit["b"] * xf ** fit["n"], "--", lw=2.5, color=c, zorder=10,
                    label=rf"{name} fit: slope = ${fit['n']:.2f} \pm {fit['n_err']:.2f}$")
    ax.set_xlim(1e-3, 1e2)
    ax.set_ylim(1e-5, 1e2)
    ax.set_yticks([1e-4, 1e-2, 1e0, 1e2])
    ax.set_xlabel(r"$k_{1}\epsilon^{1/2}S^{-3/2}$", fontsize=15)
    ax.set_ylabel(r"$-E_{12}(k_{1})/(\varepsilon^{3} S^{-5})^{1/2}$", fontsize=15)
    _grid(ax, 15)
    # minor ticks and grid lines at 1-9 x 10^n in every decade (also the unlabelled ones)
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(1, 10), numticks=100))
    ax.yaxis.set_minor_formatter(NullFormatter())  # labels only at the major ticks
    ax.yaxis.set_label_coords(0.045, 0.5, transform=blended_transform_factory(
        fig.transFigure, ax.transAxes))
    _textbox(fig, 0.5775, 0.8616, 0.2405, 0.0824, r"$(k_{1}\epsilon^{1/2}S^{-3/2})^{-7/3}$",
             fontsize=15)
    if fit_range is None:
        x, y, w, h = 0.1613, 0.1885, 0.5107, 0.1225
        ax.legend(loc="center", bbox_to_anchor=(x + w / 2, y + h / 2), bbox_transform=fig.transFigure,
                  fontsize=12, frameon=True, fancybox=False, edgecolor="k", borderpad=0.4,
                  labelspacing=0.3)
    else:
        # four entries: smaller legend in the empty lower-left corner of the axes
        ax.legend(loc="lower left", bbox_to_anchor=(0.01, 0.01), fontsize=10, frameon=True,
                  fancybox=False, edgecolor="k", borderpad=0.3, labelspacing=0.2,
                  handlelength=1.8, framealpha=0.9)
    return fig


def plot_E12_compensated(obj, sadd, num_points=50, fit_range=(1.0, 10.0),
                         C=0.15):
    """Compensated co-spectrum (k1 eps^1/2 S^-3/2)^(7/3) * |E12|/(eps^3 S^-5)^(1/2)
    at z/delta = 0.09, with the same normalization and smoothing as plot_E12.

    A -7/3 range shows as a plateau; its level is compared with C = 0.15
    (Wyngaard & Cote 1972, quoted by Saddoughi & Veeravalli 1994). The range
    ``fit_range`` of k1 eps^1/2 S^-3/2 is shaded and the mean of the smoothed
    curve in it is printed.
    """
    k1 = obj.k1
    i09 = _ind(obj, 0.0919)
    li = log_indices(k1.size, num_points, start=0.3)
    eps_, S_ = 342.0, 83.5
    xscale = eps_ ** 0.5 * S_ ** (-3 / 2)
    yscale = (eps_ ** 3 * S_ ** (-5)) ** 0.5

    def compensated(ind):
        x = k1[ind] * xscale
        return _logx(x), x ** (7 / 3) * np.abs(obj.E12[i09, ind]) / yscale

    x_line, y_line = compensated(li)
    y_line = movmean(y_line, 2, 2)
    x_sadd = np.asarray(sadd.k1_eps_S_50_Sadd)
    y_sadd = x_sadd ** (7 / 3) * np.asarray(sadd.E12_50_Sadd)

    in_range = (x_line >= fit_range[0]) & (x_line <= fit_range[1])
    plateau = np.nanmean(y_line[in_range])
    print(f"Mean compensated co-spectrum for {fit_range[0]:g} <= k1 eps^1/2 S^-3/2 <= "
          f"{fit_range[1]:g}: {plateau:.3f} (C = {C:g})")

    fig = _figure(577, 386)
    ax = fig.add_axes([0.1263, 0.1528, 0.8344, 0.8057])
    ax.axvspan(*fit_range, color=LIGHT_GOLD, alpha=0.35, lw=0, label="_nolegend_")
    ax.axhline(C, color="k", lw=1.5, ls=":", label=rf"$C = {C:g}$")
    ax.semilogx(x_sadd, y_sadd, linestyle="none", marker="v", markersize=9,
                markeredgecolor="k", markerfacecolor=SADD_FACE, markeredgewidth=1.5,
                label=rf"Saddoughi($U_{{\infty}}=50$[m/s]) $z/\delta$ = {sadd.z_inner_o_delta_50_Sadd:.2f}")
    ax.plot(x_line, y_line, "-", lw=2, color=RED, marker="x", markersize=8, markeredgecolor="k",
            markeredgewidth=2, label=rf"Gen($U_{{\infty}}=50$[m/s]) $z/\delta$ = {obj.z[i09] / obj.delta:.2f}")
    ax.set_xlim(1e-2, 1e2)
    ax.set_ylim(0, max(0.5, 1.5 * np.nanmax(y_line[in_range])))
    ax.set_xlabel(r"$k_{1}\epsilon^{1/2}S^{-3/2}$", fontsize=15)
    ax.set_ylabel(r"$(k_{1}\epsilon^{1/2}S^{-3/2})^{7/3}\,[-E_{12}/(\varepsilon^{3} S^{-5})^{1/2}]$",
                  fontsize=13)
    _grid(ax, 15)
    ax.legend(loc="upper left", fontsize=10, frameon=True, fancybox=False, edgecolor="k")
    return fig


def uw_pdfs(obj, bins=None):
    """PDF of the instantaneous u'w'/u_tau^2 at every elevation of the HRVF.

    u' and w' are fluctuations about the mean over the frames at each (z, x),
    as in the second moments. Returns ``(centers, pdf, mean)`` with ``pdf`` of
    shape ``(len(z), len(centers))`` and ``mean`` the mean u'w'/u_tau^2 per
    elevation.
    """
    bins = np.linspace(-8, 8, 321) if bins is None else bins
    u2 = obj.u_tau ** 2
    mean_u = frame_mean(obj.Gen_u_HRVF)
    mean_w = frame_mean(obj.Gen_w_HRVF)
    nz = obj.Gen_u_HRVF.shape[0]
    pdf = np.empty((nz, bins.size - 1))
    mean = np.empty(nz)
    for r in range(nz):
        uw = ((obj.Gen_u_HRVF[r] - mean_u[r][:, None])
              * (obj.Gen_w_HRVF[r] - mean_w[r][:, None])).ravel() / u2
        pdf[r], _ = np.histogram(uw, bins=bins, density=True)
        mean[r] = uw.mean()
    return 0.5 * (bins[1:] + bins[:-1]), pdf, mean


def plot_uw_distribution(obj, sadd=None, z_over_delta=(0.02, 0.03, 0.05, 0.07, 0.09, 0.12, 0.16, 0.20, 0.24),
                         xlim=(-4, 4)):
    """Distribution of u'w'/u_tau^2 with elevation, from the wall (z_min) upward.

    (a) PDFs at the elevations ``z_over_delta`` (the lowest elevation of the
    field is always included), stacked by elevation and coloured by z/delta;
    the circle marks the mean of each PDF. (b) log10 of the PDF at every
    elevation, with the mean profile (white line) and the elevations of (a)
    (dashed lines). ``sadd`` is not used (kept for ``plot_all``).
    """
    z = np.asarray(obj.z)
    zod = z / obj.delta
    centers, pdf, mean = uw_pdfs(obj)
    rows = sorted({0, *(i for i in (_ind(obj, v) for v in z_over_delta) if i is not None)})
    cmap = plt.get_cmap("viridis")
    norm = plt.matplotlib.colors.LogNorm(zod[rows[0]], zod[rows[-1]])

    fig = _figure(1100, 620)

    # (a) ridgeline of the PDFs
    ax = fig.add_axes([0.085, 0.12, 0.37, 0.82])
    spacing = 1.0
    for j, r in enumerate(rows):
        y = pdf[r] / pdf[r].max() * 1.6 * spacing  # scaled so neighbouring PDFs overlap a little
        base = j * spacing
        color = cmap(norm(zod[r]))
        ax.fill_between(centers, base, base + y, color=color, alpha=0.75, lw=0, zorder=len(rows) - j)
        ax.plot(centers, base + y, color="k", lw=1, zorder=len(rows) - j)
        ym = base + np.interp(mean[r], centers, y)
        ax.plot(mean[r], ym, "o", color="w", markeredgecolor="k", markersize=6, zorder=len(rows) + 1)
    ax.axvline(0, color="k", ls="--", lw=1, zorder=0)
    ax.set_yticks([j * spacing for j in range(len(rows))])
    ax.set_yticklabels([f"{zod[r]:.3f}" for r in rows])
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.2, (len(rows) - 1) * spacing + 2.0)
    ax.set_xlabel(r"$u^{\prime}w^{\prime}/u_{\tau}^{2}$", fontsize=15)
    ax.set_ylabel(r"$z/\delta$", fontsize=15)
    ax.tick_params(labelsize=13)
    ax.grid(True, axis="x", alpha=0.3)
    _textbox(fig, 0.075, 0.94, 0.05, 0.05, "(a)", fontsize=15)

    # (b) log10 PDF at every elevation
    ax = fig.add_axes([0.54, 0.12, 0.36, 0.82])
    with np.errstate(divide="ignore"):
        logpdf = np.log10(pdf)
    levels = np.arange(-5, 0.51, 0.25)
    cs = ax.contourf(centers, zod, np.clip(logpdf, levels[0], levels[-1]), levels=levels,
                     cmap="magma", extend="neither")
    ax.plot(mean, zod, color="w", lw=2.5, label=r"$\overline{u^{\prime}w^{\prime}}/u_{\tau}^{2}$")
    for r in rows:
        ax.axhline(zod[r], color="w", ls="--", lw=0.8, alpha=0.7)
    ax.axvline(0, color="w", ls=":", lw=1)
    ax.set_yscale("log")
    ax.set_xlim(*xlim)
    ax.set_ylim(zod[0], zod[-1])
    ticks = [t for t in (0.02, 0.03, 0.05, 0.1, 0.2) if zod[0] <= t <= zod[-1]]
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{t:g}" for t in ticks])
    ax.yaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_xlabel(r"$u^{\prime}w^{\prime}/u_{\tau}^{2}$", fontsize=15)
    ax.set_ylabel(r"$z/\delta$", fontsize=15)
    ax.tick_params(labelsize=13)
    ax.legend(loc="upper left", fontsize=12, framealpha=0.8)
    cax = fig.add_axes([0.915, 0.12, 0.015, 0.82])
    cb = fig.colorbar(cs, cax=cax)
    cb.set_label(r"$\log_{10}$ PDF", fontsize=14)
    cb.ax.tick_params(labelsize=12)
    _textbox(fig, 0.545, 0.94, 0.05, 0.05, "(b)", fontsize=15)
    return fig


def divergence_stats(obj):
    """Statistics of the in-plane divergence of the fluctuating HRVF per elevation.

    With a = du'/dx, d = dw'/dz, div = a + d and the fluctuating vorticity
    omega' = dw'/dx - du'/dz (u' and w' are fluctuations about the mean over
    the frames at each (z, x); gradients with ``np.gradient``), returns a dict
    of arrays over z:

    - ``rho``: correlation coefficient <a d> / sqrt(<a^2> <d^2>)
    - ``R``: <div^2> / (<a^2> + <d^2>)
    - ``div_omega``: rms(div) / rms(omega')
    - ``rms_ratio``: rms(d) / rms(a)

    Averages <.> are over all x and all frames at each elevation.
    """
    dx = obj.Gen_x_HRVF[1] - obj.Gen_x_HRVF[0]
    dz = obj.z[1] - obj.z[0]
    mean_u = frame_mean(obj.Gen_u_HRVF)
    mean_w = frame_mean(obj.Gen_w_HRVF)
    nz = obj.Gen_u_HRVF.shape[0]
    saa, sdd, sad, soo = (np.zeros(nz) for _ in range(4))
    for S in range(obj.Gen_u_HRVF.shape[2]):
        dudz, dudx = np.gradient(obj.Gen_u_HRVF[:, :, S] - mean_u, dz, dx)
        dwdz, dwdx = np.gradient(obj.Gen_w_HRVF[:, :, S] - mean_w, dz, dx)
        om = dwdx - dudz
        saa += (dudx ** 2).sum(axis=1)
        sdd += (dwdz ** 2).sum(axis=1)
        sad += (dudx * dwdz).sum(axis=1)
        soo += (om ** 2).sum(axis=1)
    return dict(rho=sad / np.sqrt(saa * sdd), R=(saa + sdd + 2 * sad) / (saa + sdd),
                div_omega=np.sqrt((saa + sdd + 2 * sad) / soo), rms_ratio=np.sqrt(sdd / saa))


def plot_divergence(obj, sadd=None):
    """Profiles of the normalized in-plane divergence du/dx + dw/dz versus z/delta.

    (a) correlation coefficient of du'/dx and dw'/dz, (b) <div^2>/(<(du'/dx)^2>
    + <(dw'/dz)^2>), (c) rms(div)/rms(omega'). Reference lines: strictly 2D
    divergence-free field (rho = -1, R = 0, div = 0), plane of 3D locally
    isotropic turbulence (rho = -1/2, R = 1/2 and rms(div)/rms(omega') =
    1/sqrt(5), from continuity and isotropy), and gradients with no continuity
    constraint (rho = 0, R = 1). ``sadd`` is not used.
    """
    st = divergence_stats(obj)
    zod = np.asarray(obj.z) / obj.delta
    # the first and last rows use one-sided differences in z; leave them out
    inner = slice(1, -1)
    zod = zod[inner]
    st = {k: v[inner] for k, v in st.items()}
    refs = [(-1.0, 0.0, "2D divergence-free", BLUE, "--"),
            (-0.5, 0.5, "plane of 3D isotropic turbulence", GOLD, "-"),
            (0.0, 1.0, "no continuity", BLACK, ":")]

    fig = _figure(1150, 430)
    axs = [fig.add_axes([0.065 + j * 0.325, 0.15, 0.26, 0.68]) for j in range(3)]

    ax = axs[0]
    for rho_ref, _, label, c, ls in refs:
        ax.axhline(rho_ref, color=c, ls=ls, lw=2, label=label)
    ax.plot(zod, st["rho"], color=RED, lw=2.5, label="Gen HRVF + VorX")
    ax.set_ylim(-1.1, 0.3)
    ax.set_ylabel(r"$\rho\left(\partial u^{\prime}/\partial x,\ \partial w^{\prime}/\partial z\right)$",
                  fontsize=15)

    ax = axs[1]
    for _, r_ref, label, c, ls in refs:
        ax.axhline(r_ref, color=c, ls=ls, lw=2, label=label)
    ax.plot(zod, st["R"], color=RED, lw=2.5, label="Gen HRVF + VorX")
    ax.set_ylim(-0.1, 1.3)
    ax.set_ylabel(r"$R = \langle \mathrm{div}^2\rangle / \left(\langle a^2\rangle+\langle d^2\rangle\right)$",
                  fontsize=14)

    ax = axs[2]
    # isotropy: <div^2> = <(dv/dy)^2> = eps/(15 nu), <omega_y^2> = eps/(3 nu)
    ax.axhline(0.0, color=BLUE, ls="--", lw=2)
    ax.axhline(1 / np.sqrt(5), color=GOLD, ls="-", lw=2)
    ax.plot(zod, st["div_omega"], color=RED, lw=2.5, label="Gen HRVF + VorX")
    ax.set_ylim(-0.05, max(1.0, 1.1 * np.nanmax(st["div_omega"])))
    ax.set_ylabel(r"$\mathrm{div}_{rms}/\omega^{\prime}_{rms}$", fontsize=15)

    for ax, label in zip(axs, "abc"):
        ax.set_xscale("log")
        ax.set_xlim(zod[0], zod[-1])
        ax.set_xlabel(r"$z/\delta$", fontsize=15)
        _grid(ax, 13)
        ax.text(0.03, 0.95, f"({label})", transform=ax.transAxes, fontsize=15, va="top")
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.995), ncol=4,
               fontsize=12, frameon=False)
    fig.text(0.5, 0.875, r"$a=\partial u^{\prime}/\partial x,\ d=\partial w^{\prime}/\partial z,\ "
             r"\mathrm{div}=a+d$", ha="center", fontsize=12)
    return fig


ALL_PLOTS = (plot_moments, plot_spectra, plot_E12, plot_E12_compensated, plot_uw_distribution,
             plot_divergence)


def plot_all(obj, sadd=None, save_dir=None, show=False, fmt="png", dpi=200):
    """Draw Figures 5, 6 and 7, the distribution of u'w' and the divergence profiles.

    Figures are saved to ``save_dir`` (if given) as ``<function name>.<fmt>``.
    Returns a dictionary {function name: figure}.
    """
    sadd = read_saddoughi() if sadd is None else sadd
    if save_dir is not None:
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
    figures = {}
    for func in ALL_PLOTS:
        fig = func(obj, sadd)
        figures[func.__name__] = fig
        if save_dir is not None:
            fig.savefig(save_dir / f"{func.__name__}.{fmt}", dpi=dpi)
    if show:
        plt.show()
    return figures
