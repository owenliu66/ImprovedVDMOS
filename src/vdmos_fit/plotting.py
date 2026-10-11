"""Plot helpers for comparing the improved model against Level-1 and data."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from .data import CurveSet
from .model import vdmos_id, vdmos_id_level1
from .params import DEFAULT, Params

__all__ = ["output_curves", "transfer_curves", "gds_continuity", "overlay_data", "mtriode_shape"]


def _vov(vgs, vth0=3.0, ksub=0.05):
    z = np.clip((vgs - vth0) / ksub, -700.0, 700.0)
    return ksub * np.log1p(np.exp(z))


def mtriode_shape(
    vgs=5.0,
    m_vals=(0.5, 1.0, 2.0, 4.0),
    *,
    kp=2.0,
    vth0=3.0,
    title="Level-1 mtriode is a shape knob only",
    ax=None,
):
    """Level-1 Id(Vds) for several mtriode values at fixed Kp.

    Reproduces the LTspice behaviour: every curve saturates at the SAME
    Id,sat = 0.5*Kp*Vov^2 while Ron = 1/(m*Kp*Vov) and the knee at Vov/m move
    together. You cannot set Ron, knee shape and saturation level independently,
    which is the actual Level-1 limitation this project addresses.
    """
    ax = ax or plt.subplots(figsize=(7, 4.5))[1]
    vds = np.linspace(0.0, 6.0, 600)
    for m in m_vals:
        ax.plot(
            vds,
            vdmos_id_level1(vds, vgs, kp=kp, vth0=vth0, m_triode=m, lam=0.0),
            label=f"mtriode={m:g}",
        )
    ax.set_xlabel("Vds [V]")
    ax.set_ylabel("Id [A]")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    return ax


def output_curves(
    p: Params = DEFAULT,
    vgs_vals=(4, 6, 8, 10, 12),
    vds=None,
    *,
    level1: bool = True,
    title: str = "Output characteristics",
    ax=None,
):
    """Family of Id(Vds) curves, with Level-1 overlaid for comparison.

    Both are C1 continuous; the visible difference is the saturation exponent
    (Level-1 is locked to the square law, the improved model is not).
    """
    vds = np.linspace(0.0, 40.0, 600) if vds is None else np.asarray(vds, float)
    ax = ax or plt.subplots(figsize=(7, 4.5))[1]
    for vg in vgs_vals:
        ax.plot(vds, vdmos_id(vds, vg, p), label=f"improved Vgs={vg:g}")
        if level1:
            ax.plot(vds, vdmos_id_level1(vds, vg), "--", alpha=0.5, label=f"L1 Vgs={vg:g}")
    ax.set_xlabel("Vds [V]")
    ax.set_ylabel("Id [A]")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    return ax


def transfer_curves(
    p: Params = DEFAULT,
    vgs=None,
    vds_vals=(0.05, 1.0, 5.0, 20.0),
    *,
    title: str = "Transfer curves",
    ax=None,
):
    vgs = np.linspace(0.0, 15.0, 400) if vgs is None else np.asarray(vgs, float)
    ax = ax or plt.subplots(figsize=(7, 4.5))[1]
    for vd in vds_vals:
        ax.plot(vgs, vdmos_id(vd, vgs, p), label=f"Vds={vd:g} V")
    ax.set_xlabel("Vgs [V]")
    ax.set_ylabel("Id [A]")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    return ax


def gds_continuity(
    p: Params = DEFAULT,
    vgs_vals=(6.0, 10.0),
    *,
    title: str = "Output conductance across the knee",
    ax=None,
):
    """gds = dId/dVds across the knee, improved model vs Level-1.

    Both models are C1 continuous -- LTspice places the Level-1 knee at
    Vov/mtriode precisely to guarantee that. The point of this plot is that the
    improved model's knee *shape* is a fitted degree of freedom (theta, vsat)
    while Level-1's is locked to the square law.
    """
    ax = ax or plt.subplots(figsize=(7, 4.5))[1]
    h = 1e-4
    for vg in vgs_vals:
        vov = _vov(vg)
        vds = np.linspace(0.5 * vov, 2.0 * vov, 400)
        g_new = (vdmos_id(vds + h, vg, p) - vdmos_id(vds - h, vg, p)) / (2 * h)
        g_l1 = (vdmos_id_level1(vds + h, vg) - vdmos_id_level1(vds - h, vg)) / (2 * h)
        ax.semilogy(vds, np.abs(g_new) + 1e-9, "-", label=f"improved Vgs={vg:g}")
        ax.semilogy(vds, np.abs(g_l1) + 1e-9, "--", alpha=0.7, label=f"Level-1 Vgs={vg:g}")
    ax.set_xlabel("Vds [V]")
    ax.set_ylabel("|gds| [S]")
    ax.set_title(title)
    ax.grid(alpha=0.3, which="both")
    return ax


def overlay_data(data: CurveSet, p: Params = DEFAULT, *, title="Fit vs measured", ax=None):
    """Measured points vs model curves, grouped by Vgs."""
    ax = ax or plt.subplots(figsize=(7, 4.5))[1]
    for vg in data.vgs_grid():
        m = np.isclose(data.vgs, vg)
        vd, meas = data.vds[m], data.ids[m]
        order = np.argsort(vd)
        ax.plot(vd[order], meas[order], "o", ms=4, alpha=0.6, label=f"data Vgs={vg:g}")
        vds = np.linspace(vd.min(), vd.max(), 200)
        ax.plot(vds, vdmos_id(vds, vg, p), "-", alpha=0.8)
    ax.set_xlabel("Vds [V]")
    ax.set_ylabel("Id [A]")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    return ax
