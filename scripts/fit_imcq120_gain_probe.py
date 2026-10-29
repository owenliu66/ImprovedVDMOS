"""Probe (NOT part of the package model): what mechanism the IMCQ120R004M2H data
actually wants, and how much it gains.

HISTORICAL PROBE (session 18) -- the gain it proposed landed in the package in
session 19 (k_beta_trap) and the OFFICIAL IMCQ record is now the session-21
joint fit (scripts/fit_imcq120_joint.py, reference.py::IMCQ120_JOINT). Kept as
the standing demonstration of the super-linear implied-conductance diagnostic
and the gain's effect size. The report it referenced,
out/IMCQ120_TEMPERATURE_REPORT.md, was deleted session 21 (findings carried
into out/IMCQ120_JOINT_FIT_REPORT.md section 8).

Finding that motivated this script (see out/IMCQ120_JOINT_FIT_REPORT.md
section 8 item 6, "fit-quality limit"): the datasheet's own on-curves imply a
channel conductance that grows SUPER-linearly with gate overdrive. Reading the
first (near-origin, ohmic) point of each 25 C curve and subtracting the r_d
drop gives the implied channel overdrive v_ov = g_ch / k_p:

    v_gs      7      8     10     12     15     18
    v_ov    0.24   0.72   2.09   3.83   6.99  10.85  V
    dv/dv   0.48   0.69   0.87   1.05   1.29        V/V  (rises PAST 1)

The improved model's channel conductance is k_p * v_ov_eff with
0 <= dv_ov_eff/dv_ov <= 1 (the charge cap enforces exactly this bound to keep the
package-R solver bracket invariant dg/di_d <= -1). A softplus v_ov(v_gs) is
convex but its slope saturates at 1, so it can NEVER produce slope > 1. That is
why the baseline fit needs v_sub ~ 0.83-0.99 V (near its 1.0 bound) and still
leaves +18 to +38 % residuals on the near-origin points of the v_gs = 8/10/12
curves and -9 to -11 % on the high-v_gs tails.

The mechanism that fixes it is the INVERSE of the charge cap: a trap-filling /
inversion-mobility-recovery gain on the overdrive,

    v_ov_eff = v_ov * (1 + k_beta * v_ov)        (k_beta [1/V], 0 = off)

physically: 4H-SiC channel mobility rises as the surface potential deepens and
interface traps fill, so the conductance grows faster than (v_gs - V_th). This
is the SiC counterpart of the super-junction charge cap the ST8L65 needed --
same knob family, opposite sign.

Result (relative weighting, k_lambda = 0, r_s = 0, cap off):

    25 C : RMSE 8.53 % -> 4.61 %   k_beta = 0.0569 1/V
    175 C: RMSE 6.13 % -> 3.31 %   k_beta = 0.0757 1/V  (+33 %, a NEW TC)

    Ron chord @185.2 A: 3.701 / 4.729 mOhm (25 C, data 3.70 / 4.70)
                        8.808 / 9.685 mOhm (175 C, data 8.90 / 9.43)
    Vth spec point: 57.3 mA (25 C) / 58.2 mA (175 C) vs 58.3 mA
    g_fs @185.2 A/20 V: 128.8 S (25 C) / 126.7 S (175 C) vs 122 S spec

CAVEATS before promoting this into the package:
  1. dv_ov_eff/dv_ov = 1 + 2*k_beta*v_ov reaches 2.4 (25 C) and 3.0 (175 C) at
     v_gs = 18 -- it BREAKS the 0 <= dv_ov_eff/dv_ov <= 1 bound the package-R
     solver bracket proof relies on. With r_s = 0 (this Kelvin-source device) the
     bracket is still valid: g(i) = f(v_ds - r_d*i, v_gs) - i has
     dg/di = -r_d * df/dv_ds - 1 <= -1 regardless of the gain, because the gain
     acts on the GATE derivative, not the drain one. With r_s > 0 (any silicon
     part) the proof must be re-derived. This script uses plain vectorized
     bisection, which is safe because f stays monotone decreasing in i.
  2. The gain makes the high-v_gs transfer MORE convex, which overshoots the
     datasheet transfer figure at V_DS = 20 V (model 3650 A vs figure ~2300 A at
     v_gs = 18). The output-curve data cannot see that region (the v_gs = 18
     curve exits the 1600 A plot frame at ~6.4 V). Bending those curves over is
     exactly what k_lambda does -- so the transfer characteristic (V_DS = 20 V,
     25/175/200 C, datasheet p. 9) should be digitized and k_lambda re-freed in
     a JOINT fit before quoting the high-v_gs / high-v_ds corner.
  3. v_th0 with the gain: 6.06 -> 4.77 V, slope -8.62 mV/K (further from the
     datasheet -6.67) and the offset vs datasheet Vth is no longer constant
     (+1.86 / +1.57 V). k_p's TC also flips sign vs the baseline fit
     (ratio 0.84, beta = +0.44, vs 1.09 / -0.22) -- a direct demonstration that
     the k_p temperature coefficient is model-form dependent, not measurable.

    python scripts/fit_imcq120_gain_probe.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from vdmos_fit import CurveSet, load_csv  # noqa: E402

OUT = ROOT / "out"
NAMES = ("k_p", "v_th0", "v_sub", "k_beta", "k_kneesharp", "k_eta", "v_dbsat")
LOG_IDX = (0, 2, 3, 6)  # packed in log space


def channel_current(v_ds, v_gs, k_p, v_th0, v_sub, k_beta, k_ksh, k_eta, v_dbsat):
    """Improved-model channel with the trap-filling gain; cap/CLM/kink off."""
    v_th = v_th0 - k_eta * v_ds * np.tanh(v_ds / v_dbsat)
    v_ov = v_sub * np.logaddexp(0.0, (v_gs - v_th) / v_sub)
    v_ov_eff = v_ov * (1.0 + k_beta * v_ov)          # <-- the new mechanism
    v_knee = np.maximum(v_ov_eff, 1e-3)
    v_ds_eff = v_ds / (1.0 + (np.maximum(v_ds, 1e-12) / v_knee) ** k_ksh) ** (1 / k_ksh)
    return k_p * v_ov_eff * (v_ds_eff - v_ds_eff**2 / (2 * v_knee))


def drain_current(v_ds, v_gs, p, r_d, iters=70):
    """Package-r_d solve by vectorized bisection (r_s = 0: Kelvin source)."""
    v_ds = np.atleast_1d(np.asarray(v_ds, float))
    if r_d <= 0:
        return channel_current(v_ds, v_gs, *p)
    lo = np.zeros_like(v_ds)
    hi = v_ds / r_d
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        f = channel_current(v_ds - r_d * mid, v_gs, *p) - mid
        pos = f > 0
        lo = np.where(pos, mid, lo)
        hi = np.where(pos, hi, mid)
    return 0.5 * (lo + hi)


def unpack(x):
    p = np.array([np.exp(x[0]), x[1], np.exp(x[2]), np.exp(x[3]), x[4], x[5], np.exp(x[6])])
    return p, x[7]


def residuals(x, data, rmax):
    p, r_d = unpack(x)
    model = drain_current(data.v_ds, data.v_gs, p, r_d)
    return (model - data.i_d) / np.maximum(np.abs(data.i_d), 1e-3)


LO = np.array([np.log(1.0), 3.0, np.log(0.02), np.log(1e-4), 0.1, 1e-5, np.log(0.1), 0.0])
HI = np.array([np.log(1e4), 10.0, np.log(1.0), np.log(2.0), 40.0, 0.05, np.log(10.0), 0.0])
SEEDS = [(50.0, 1e-3, 0.8), (50.0, 0.05, 0.8), (120.0, 0.15, 0.5), (300.0, 0.3, 0.3)]


def fit_gain(data: CurveSet, rmax: float, v_th_start: float):
    hi = HI.copy()
    hi[7] = rmax
    best = None
    for k_p0, beta0, v_sub0 in SEEDS:
        x0 = np.array([np.log(k_p0), v_th_start, np.log(v_sub0), np.log(beta0),
                       1.0, 0.013, np.log(1.4), min(1.9e-3, rmax * 0.6)])
        sol = least_squares(lambda x: residuals(x, data, rmax), x0, bounds=(LO, hi),
                            method="trf", x_scale="jac", max_nfev=2500)
        if best is None or sol.cost < best.cost:
            best = sol
    p, r_d = unpack(best.x)
    rel = residuals(best.x, data, rmax)
    return p, r_d, float(np.sqrt(np.mean(rel**2)))


def main() -> None:
    OUT.mkdir(exist_ok=True)
    d25 = load_csv(ROOT / "data" / "IMCQ120R004M2H_IV_25C.csv")
    d175 = load_csv(ROOT / "data" / "IMCQ120R004M2H_IV_175C.csv")

    print("implied channel overdrive from the near-origin ohmic points (25 C):")
    for vg in (7, 8, 10, 12, 15, 18):
        m = np.isclose(d25.v_gs, vg)
        o = int(np.argmin(d25.v_ds[m]))
        v_ds, i_d = d25.v_ds[m][o], d25.i_d[m][o]
        g_ch = i_d / (v_ds - 1.89e-3 * i_d)
        print(f"   v_gs={vg:2d}: v_ds={v_ds:5.3f} V i_d={i_d:7.1f} A"
              f"  g_ch={g_ch:7.1f} A/V  -> v_ov = g_ch/k_p = {g_ch / 50.9:5.2f} V")

    p25, rd25, rm25 = fit_gain(d25, 3.7e-3, 6.6)
    p175, rd175, rm175 = fit_gain(d175, 8.9e-3, 5.6)
    print(f"\nGAIN fit 25 C : RMSE {rm25 * 100:.2f} %   (baseline cap-free/lam=0: 8.53 %)")
    print(f"GAIN fit 175 C: RMSE {rm175 * 100:.2f} %   (baseline cap-free/lam=0: 6.13 %)")
    print("\nparameter drift 25 -> 175 C:")
    for n, a, b in zip(NAMES, p25, p175, strict=True):
        print(f"   {n:<12s} {a:10.5g} -> {b:10.5g}   ratio {b / a:7.3f}")
    print(f"   {'r_d':<12s} {rd25 * 1e3:10.4f} -> {rd175 * 1e3:10.4f} mOhm"
          f"   ratio {rd175 / rd25:7.3f}")
    print(f"   v_th0 slope = {(p175[1] - p25[1]) / 150 * 1e3:.2f} mV/K"
          "   (datasheet -6.67 mV/K)")
    print(f"   k_beta slope = {(p175[3] - p25[3]) / 150 * 1e3:.3f} mV^-1/K")
    for tag, p in (("25 C", p25), ("175 C", p175)):
        der = 1 + 2 * p[3] * p[2] * np.logaddexp(0.0, (18 - p[1]) / p[2])
        print(f"   {tag}: dv_ov_eff/dv_ov at v_gs=18 = {der:.2f}"
              "   (>1 breaks the package-R bracket proof if r_s > 0)")

    for tag, p, rd, d, vth in (("25C", p25, rd25, d25, 4.2), ("175C", p175, rd175, d175, 3.2)):
        rel = (drain_current(d.v_ds, d.v_gs, p, rd) - d.i_d) / d.i_d
        print(f"\n{tag} per-curve RMSE %:", {vg: round(float(np.sqrt(np.mean(
            rel[np.isclose(d.v_gs, vg)] ** 2)) * 100), 2) for vg in d.v_gs_grid()})
        off = float(drain_current(np.array([vth]), np.array([vth]), p, rd)[0])
        print(f"   Vth spec point: {off * 1e3:.2f} mA vs 58.3 mA")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4))
    for d, p, rd, c, lab in ((d25, p25, rd25, "C0", "25 C"), (d175, p175, rd175, "C1", "175 C")):
        for vg in d.v_gs_grid():
            m = np.isclose(d.v_gs, vg)
            if m.sum() == 1:
                continue
            vd = np.linspace(0.0, d.v_ds[m].max(), 300)
            for ax in axes:
                ax.plot(d.v_ds[m], d.i_d[m], "o", ms=4, alpha=0.55,
                        color=c, label=f"{lab}, v_gs={vg:g}" if ax is axes[0] else None)
                ax.plot(vd, drain_current(vd, vg, p, rd), "--", color=c, alpha=0.9)
    axes[0].set_xlim(0, 21)
    axes[1].set_xlim(0, 3)
    axes[1].set_ylim(0, 120)
    for ax, t in zip(axes, ["Gain-model fit vs data (dots)", "Knee zoom (v_ds <= 3 V)"],
                     strict=True):
        ax.set_xlabel("v_ds [V]")
        ax.set_ylabel("i_d [A]")
        ax.set_title(t)
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=7, ncol=2)
    fig.savefig(OUT / "imcq120_gain_fit.png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    print("wrote", OUT / "imcq120_gain_fit.png")

    txt = OUT / "imcq120_gain_params.json"
    txt.write_text(json.dumps({"names": list(NAMES),
                               "25C": p25.tolist(), "r_d_25C": rd25, "rmse_25C": rm25,
                               "175C": p175.tolist(), "r_d_175C": rd175, "rmse_175C": rm175},
                              indent=2), encoding="utf-8")
    print("wrote", txt)


if __name__ == "__main__":
    main()
