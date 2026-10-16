"""Plot helpers for comparing the improved model against Level-1 and data."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from .data import CurveSet
from .model import vdmos_id, vdmos_id_level1
from .params import DEFAULT, Params

__all__ = ["output_curves", "transfer_curves", "gds_continuity", "overlay_data", "mtriode_shape"]


def _vov(v_gs, v_th0=3.0, v_sub=0.05):
    z = np.clip((v_gs - v_th0) / v_sub, -700.0, 700.0)
    return v_sub * np.log1p(np.exp(z))


def mtriode_shape(
    v_gs=5.0,
    m_vals=(0.5, 1.0, 2.0, 4.0),
    *,
    k_p=2.0,
    v_th0=3.0,
    title="Level-1 mtriode is a shape knob only",
    ax=None,
):
    """Level-1 i_d(v_ds) at several mtriode values: same i_d,sat, Ron and knee
    move together (shape knob only; see README, "Why").
    """
    ax = ax or plt.subplots(figsize=(7, 4.5))[1]
    v_ds = np.linspace(0.0, 6.0, 600)
    for m in m_vals:
        ax.plot(
            v_ds,
            vdmos_id_level1(v_ds, v_gs, k_p=k_p, v_th0=v_th0, m_triode=m, k_lambda=0.0),
            label=f"mtriode={m:g}",
        )
    ax.set_xlabel("v_ds [V]")
    ax.set_ylabel("i_d [A]")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    return ax


def output_curves(
    p: Params = DEFAULT,
    vgs_vals=(4, 6, 8, 10, 12),
    v_ds=None,
    *,
    level1: bool = True,
    title: str = "Output characteristics",
    ax=None,
):
    """Family of i_d(v_ds) curves with Level-1 overlaid; both C1, the visible
    difference is the saturation exponent.
    """
    v_ds = np.linspace(0.0, 40.0, 600) if v_ds is None else np.asarray(v_ds, float)
    ax = ax or plt.subplots(figsize=(7, 4.5))[1]
    for vg in vgs_vals:
        ax.plot(v_ds, vdmos_id(v_ds, vg, p), label=f"improved v_gs={vg:g}")
        if level1:
            ax.plot(v_ds, vdmos_id_level1(v_ds, vg), "--", alpha=0.5, label=f"L1 v_gs={vg:g}")
    ax.set_xlabel("v_ds [V]")
    ax.set_ylabel("i_d [A]")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    return ax


def transfer_curves(
    p: Params = DEFAULT,
    v_gs=None,
    vds_vals=(0.05, 1.0, 5.0, 20.0),
    *,
    title: str = "Transfer curves",
    ax=None,
):
    v_gs = np.linspace(0.0, 15.0, 400) if v_gs is None else np.asarray(v_gs, float)
    ax = ax or plt.subplots(figsize=(7, 4.5))[1]
    for vd in vds_vals:
        ax.plot(v_gs, vdmos_id(vd, v_gs, p), label=f"v_ds={vd:g} V")
    ax.set_xlabel("v_gs [V]")
    ax.set_ylabel("i_d [A]")
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
    """gds = di_d/dv_ds across the knee, improved vs Level-1. Both C1; the
    improved model's knee shape is a fitted degree of freedom (k_theta, v_vsat).
    """
    ax = ax or plt.subplots(figsize=(7, 4.5))[1]
    h = 1e-4
    for vg in vgs_vals:
        v_ov = _vov(vg)
        v_ds = np.linspace(0.5 * v_ov, 2.0 * v_ov, 400)
        g_new = (vdmos_id(v_ds + h, vg, p) - vdmos_id(v_ds - h, vg, p)) / (2 * h)
        g_l1 = (vdmos_id_level1(v_ds + h, vg) - vdmos_id_level1(v_ds - h, vg)) / (2 * h)
        ax.semilogy(v_ds, np.abs(g_new) + 1e-9, "-", label=f"improved v_gs={vg:g}")
        ax.semilogy(v_ds, np.abs(g_l1) + 1e-9, "--", alpha=0.7, label=f"Level-1 v_gs={vg:g}")
    ax.set_xlabel("v_ds [V]")
    ax.set_ylabel("|gds| [S]")
    ax.set_title(title)
    ax.grid(alpha=0.3, which="both")
    return ax


def overlay_data(data: CurveSet, p: Params = DEFAULT, *, title="Fit vs measured", ax=None):
    """Measured points vs model curves, grouped by v_gs."""
    ax = ax or plt.subplots(figsize=(7, 4.5))[1]
    for vg in data.v_gs_grid():
        msk = np.isclose(data.v_gs, vg)
        vd, meas = data.v_ds[msk], data.i_d[msk]
        order = np.argsort(vd)
        ax.plot(vd[order], meas[order], "o", ms=4, alpha=0.6, label=f"data v_gs={vg:g}")
        v_ds = np.linspace(vd.min(), vd.max(), 200)
        ax.plot(v_ds, vdmos_id(v_ds, vg, p), "-", alpha=0.8)
    ax.set_xlabel("v_ds [V]")
    ax.set_ylabel("i_d [A]")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    return ax
