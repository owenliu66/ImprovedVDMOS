"""Re-run the ST8L65N044M9 fit with the drain-knee softness exponent `k_kneesharp` free.

Session-10 follow-up to scripts/fit_st8l65.py (handover "Next tasks" 2-4):
same fit-B recipe and DEVICE_BOUNDS, but `k_kneesharp` is added to the free tuple
(it is deliberately NOT in DEFAULT_FREE) and started near 2.8, which
reproduces the legacy tanh knee (session-9 knee analysis). `k_kneesharp` is fitted
in LINEAR space with lower bound 0.1 (relaxed from 0.5 after the probe
landed at 0.47) -- 0 is the "legacy tanh" sentinel and must stay reachable
only by leaving k_kneesharp out of the free set, never by the optimizer landing
near 0.

Questions this run answers:
  * Does freeing k_kneesharp beat the previous official 9-param fit (9.97 % fit A /
    5.54 % fit B)?  Yes -- and fit A is now the OFFICIAL reference below.
  * Does k_kneesharp pin at 0.1 or 40 (on_bounds)? Pinned = the data does not want
    the new shape.
  * Does the k_theta<->v_vsat ~84 % correlated pair loosen now that knee width
    has its own knob? Checked with two EXPLORATORY runs on the fit-B data
    that additionally free k_theta and v_vsat: 11-param (no k_kneesharp) vs 12-param
    (with k_kneesharp); compare the k_theta<->v_vsat correlation from identifiability().

OFFICIAL STATUS (promoted 2026-10-13): **fit A of this script, 10-param with
k_kneesharp free, 3.79 % relative RMSE on the full 41-point set**, is the quoted
reference for ST8L65N044M9. It supersedes the 9-param fit B (5.54 %,
scripts/fit_st8l65.py) -- see STATUS.md "Current ST8L65 fit". k_lambda and v_dbsat
report sigma = inf (NOT SENSED): the dataset stops at v_ds = 20 V and never
enters deep saturation, so no CLM/drain-gate information exists to sense.
That is expected, not a fit failure; the knee SHAPE (k_kneesharp + v_ovmax + k_ovsharp) and
r_d+r_s are what the data pins.

Outputs (out/):
  st8l65_nd_params.txt   fitted parameters, same format as the capfit file
  st8l65_nd_probe.txt    sensitivity_report + FitResult.summary() +
                         identifiability().summary() per fit
  st8l65_nd_fit.png      full curves + knee zoom, official fit A vs k_kneesharp fit B
  st8l65_nd_percurve.png i_d vs v_ds per measured v_gs, data vs k_kneesharp fit A vs
                         k_kneesharp fit B (knee-shape check vs Figure 3)

    python scripts/fit_st8l65_nd.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from vdmos_fit import (  # noqa: E402
    ST8L65_FIT_A_FREE,
    ST8L65_OFFICIAL,
    ST8L65_OFFICIAL_RMSE_REL,
    CurveSet,
    Params,
    fit,
    load_csv,
    vdmos_id,
)
from vdmos_fit.fit import FitResult, identifiability, sensitivity_report  # noqa: E402
from vdmos_fit.params import BOUNDS  # noqa: E402
from vdmos_fit.plotting import overlay_data  # noqa: E402

CSV = ROOT / "data" / "ST8L65N044M9_IV.csv"
OUT = ROOT / "out"

# OFFICIAL reference -- lives in the package (vdmos_fit.reference) so the
# test suite can pin it too. This script's fit A is what produced it; the
# drift check in main() confirms each run still reproduces the record
# (OFFICIAL_TOL). Supersedes the 9-param fit B of scripts/fit_st8l65.py.
OFFICIAL = {
    "label": "official fit A (10-param, k_kneesharp free)",
    "rmse_rel": ST8L65_OFFICIAL_RMSE_REL,
    "params": ST8L65_OFFICIAL,
}
OFFICIAL_TOL = 1e-3   # max relative drift allowed before the record is stale

# Same re-anchored start as fit_st8l65.py, plus k_kneesharp near the legacy-tanh value.
START = Params(
    k_p=40.0, v_th0=4.7, v_sub=0.06, v_ovmax=2.0, k_ovsharp=10.0,
    k_theta=0.0, v_vsat=1e4, k_kneesharp=2.8, k_eta=0.0, k_lambda=1e-3,
    g_k=0.0, v_thk=1.6, v_kovsoft=0.1, v_dbsat=10.0, r_d=0.030, r_s=0.005,
)

# 9-param free set + k_kneesharp (explicitly freed -- NOT in DEFAULT_FREE). Same
# tuple as the promoted official record (vdmos_fit.ST8L65_FIT_A_FREE).
FREE = ST8L65_FIT_A_FREE

# Per-device bounds to fit_st8l65.py, EXCEPT k_lambda: the global (1e-12, 1.0)
# bound now governs (session-10 bound widening -- the old 0.004 cap pinned
# k_lambda at its LOWER 1e-6 and hid whether the data wants CLM at all).
# k_kneesharp keeps its global (0.1, 40) linear bound; do NOT override it here.
DEVICE_BOUNDS = {
    "r_d": (0.0, 0.036),
    "r_s": (0.0, 0.036),
}

# Same four basins as fit_st8l65.py, k_kneesharp started near 2.8 (legacy tanh) with
# one deliberately sharper (k_kneesharp=6) and one softer (k_kneesharp=1.5) start.
MULTI_START = [
    dict(k_p=90.0, v_th0=4.73, v_sub=0.06, k_lambda=1e-3, v_ovmax=1.9, k_ovsharp=10.0,
         v_dbsat=0.7, r_d=0.026, r_s=0.005, k_kneesharp=2.8),
    dict(k_p=40.0, v_th0=4.70, v_sub=0.06, k_lambda=1e-3, v_ovmax=2.5, k_ovsharp=5.0,
         v_dbsat=2.0, r_d=0.020, r_s=0.010, k_kneesharp=2.8),
    dict(k_p=150.0, v_th0=4.60, v_sub=0.10, k_lambda=3e-3, v_ovmax=1.2, k_ovsharp=20.0,
         v_dbsat=1.0, r_d=0.030, r_s=0.002, k_kneesharp=1.5),
    dict(k_p=60.0, v_th0=4.86, v_sub=0.03, k_lambda=1e-4, v_ovmax=3.0, k_ovsharp=8.0,
         v_dbsat=5.0, r_d=0.015, r_s=0.015, k_kneesharp=6.0),
]


def multi_fit(data: CurveSet, free=FREE, start_base=START) -> FitResult:
    """Best-of-multi-start; all starts reported so basin-hopping is visible."""
    best = None
    for s in MULTI_START:
        st = start_base.replace(**s)
        res = fit(data, start=st, free=free, weight="relative", max_nfev=8000,
                  bounds_overrides=DEVICE_BOUNDS)
        knee_s = s.get("k_kneesharp", 0.0)
        print(f"    start k_p={s['k_p']:8.4g} v_ovmax={s['v_ovmax']:6.4g}"
              f" k_ovsharp={s['k_ovsharp']:6.4g} v_dbsat={s['v_dbsat']:6.4g}"
              f" k_kneesharp={knee_s:5.3g}"
              f" -> RMSE {res.rmse_rel * 100:6.2f} %"
              + ("" if res.success else "   [stopped]"))
        if best is None or res.rmse_rel < best.rmse_rel:
            best = res
    return best


def report(res: FitResult, data: CurveSet, label: str, free=FREE) -> None:
    p = res.params
    state = "converged" if res.success else "STOPPED"
    print(f"\n=== {label}: RMSE_rel = {res.rmse_rel * 100:.2f} % "
          f"({res.n_points} pts, {res.n_free} free, {state}) ===")
    for name in free:
        print(f"   {name:<12s} = {getattr(p, name):.6g}")
    print(f"   r_d+r_s = {(p.r_d + p.r_s) * 1e3:.2f} mOhm"
          f"   (measured Ron chord at v_gs=10: 34.5 mOhm)")
    plateau = float(vdmos_id(19.7, 10.0, p))
    print(f"   plateau v_gs=10 @ 19.7 V: {plateau:.1f} A   (measured 240.2 A)")
    off = float(vdmos_id(3.7, 3.7, p))
    print(f"   v_gs(th) point (3.7, 3.7): {off * 1e6:.1f} uA   (measured 250 uA)")
    pinned = res.on_bounds()
    if pinned:
        print("   PINNED AT BOUND (not a physical result):")
        for name, which in pinned.items():
            print(f"     {name:<12s} {which}")
    else:
        print("   no parameter pinned at a bound")
    ident = identifiability(res, free)
    print("   " + ident.summary().replace("\n", "\n   "))

    model = vdmos_id(data.v_ds, data.v_gs, p)
    rel = (model - data.i_d) / data.i_d
    print("   per-curve RMSE %:", {
        vg: round(float(np.sqrt(np.mean(rel[np.isclose(data.v_gs, vg)] ** 2)) * 100), 2)
        for vg in data.v_gs_grid()
    })


def kneesharp_bound_probe(data_b: CurveSet, p_ks: Params) -> str:
    """Is the k_kneesharp=0.5 sentinel guard costing anything?

    Refit with the lower bound relaxed to 0.1. If the optimum sits just
    below 0.5 with negligible RMSE gain, the guard is harmless and the
    fitted shape is the softest sensible power-law soft-min (k_kneesharp -> 0
    collapses the current, so going lower is not physical anyway).
    Also exposes knee-cluster degeneracy: different v_dbsat/k_p/k_lambda splits
    reaching the same RMSE mean the SHAPE is pinned, not the mechanism.
    """
    lines = ["", "k_kneesharp lower-bound probe (fit-B data, global lower bound 0.1):"]
    bounds = {**DEVICE_BOUNDS, "k_kneesharp": (0.1, 40.0)}
    for ks0 in (0.5, 0.3):
        res = fit(data_b, start=p_ks.replace(k_kneesharp=ks0), free=FREE,
                  weight="relative", max_nfev=8000, bounds_overrides=bounds)
        q = res.params
        lines.append(
            f"  start k_kneesharp={ks0:g}: RMSE {res.rmse_rel * 100:.2f} %"
            f" | k_kneesharp={q.k_kneesharp:.4g} v_dbsat={q.v_dbsat:.4g}"
            f" v_ovmax={q.v_ovmax:.4g} k_ovsharp={q.k_ovsharp:.4g}"
            f" k_p={q.k_p:.4g} k_lambda={q.k_lambda:.4g}"
            + ("  [on_bounds: " + ",".join(res.on_bounds()) + "]"
               if res.on_bounds() else "")
        )
    return "\n".join(lines)


def klambda_cap_probe(data_b: CurveSet) -> str:
    """Does dropping the old per-device k_lambda cap (1e-6, 0.004) change anything?

    The session-9/10 fits pinned k_lambda at that cap's LOWER bound, so the global
    widening to (1e-12, 1.0) is inert unless the override is removed. Refit
    fit B both ways from the same starts and report RMSE + pinned set.
    """
    lines = ["", "k_lambda cap probe (fit-B data):"]
    for tag, ov in (("k_lambda capped (1e-6, 0.004)",
                     {**DEVICE_BOUNDS, "k_lambda": (1e-6, 0.004)}),
                    ("k_lambda free (global 1e-12, 1.0)", DEVICE_BOUNDS)):
        best = None
        for s in MULTI_START:
            st = START.replace(**s)
            res = fit(data_b, start=st, free=FREE, weight="relative",
                      max_nfev=8000, bounds_overrides=ov)
            if best is None or res.rmse_rel < best.rmse_rel:
                best = res
        q = best.params
        lines.append(
            f"  {tag}: RMSE {best.rmse_rel * 100:.2f} %"
            f" | k_lambda={q.k_lambda:.4g} k_p={q.k_p:.4g} v_dbsat={q.v_dbsat:.4g}"
            f" k_kneesharp={q.k_kneesharp:.4g} v_ovmax={q.v_ovmax:.4g} k_ovsharp={q.k_ovsharp:.4g}"
            + ("  [on_bounds: " + ",".join(best.on_bounds()) + "]"
               if best.on_bounds() else "  [no param pinned]")
        )
    return "\n".join(lines)


def ktheta_vsat_probe(data_b: CurveSet) -> str:
    """A/B the k_theta<->v_vsat correlation with and without k_kneesharp, fit-B data.

    Exploratory only (k_theta/v_vsat freed on top of the 9+k_kneesharp set, seeded in the
    fit-B basin). The ~84 % k_theta<->v_vsat pair was measured on the synthetic
    grid; this checks whether k_kneesharp loosens it on the real data.
    """
    base = OFFICIAL["params"].replace(k_kneesharp=2.8, k_theta=0.02, v_vsat=5.0, k_lambda=1e-3)
    lines = ["", "k_theta<->v_vsat correlation probe (fit-B data, exploratory):"]
    for tag, free in (
        ("11-param, no k_kneesharp", FREE[:-1] + ("k_theta", "v_vsat")),
        ("12-param, with k_kneesharp", FREE + ("k_theta", "v_vsat")),
    ):
        res = fit(data_b, start=base, free=free, weight="relative",
                  max_nfev=8000, bounds_overrides=DEVICE_BOUNDS)
        ident = identifiability(res, free)
        pair = next(((a, b, c) for a, b, c in ident.worst_pairs
                     if {a, b} == {"k_theta", "v_vsat"}), None)
        lines.append(
            f"  {tag}: RMSE {res.rmse_rel * 100:.2f} %"
            f" | k_theta={res.params.k_theta:.4g} v_vsat={res.params.v_vsat:.4g}"
            f" k_kneesharp={res.params.k_kneesharp:.4g}"
            + (f" | corr(k_theta,v_vsat) = {pair[2]:+.3f}" if pair
               else " | k_theta<->v_vsat NOT in |corr|>0.8 list")
        )
        for a, b, c in ident.worst_pairs:
            lines.append(f"      pair {a}<->{b}: {c:+.3f}")
    return "\n".join(lines)


def plot(data: CurveSet, fits: dict[str, Params], path: Path) -> None:
    """Left: full output curves. Right: knee zoom, where k_kneesharp changes the shape."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    on = data.i_d > 1e-3
    for vg in data.v_gs_grid():
        msk = np.isclose(data.v_gs, vg) & on
        if not np.any(msk):
            continue
        vd_data, id_data = data.v_ds[msk], data.i_d[msk]
        vd = np.linspace(0.0, vd_data.max(), 300)
        for ax in axes:
            ax.plot(vd_data, id_data, "o", ms=4, alpha=0.65,
                    label=f"data v_gs={vg:g}" if ax is axes[0] else None)
            for name, p in fits.items():
                ax.plot(vd, vdmos_id(vd, vg, p), "--", alpha=0.85,
                        label=f"{name} v_gs={vg:g}" if ax is axes[0] else None)
    axes[0].set_xlim(0, 21)
    axes[1].set_xlim(0, 3)
    axes[1].set_ylim(0, 80)
    for ax, title in zip(axes, ["Output curves: data (dots) vs fits (dashed)",
                                "Knee zoom (v_ds <= 3 V)"], strict=True):
        ax.set_xlabel("v_ds [V]")
        ax.set_ylabel("i_d [A]")
        ax.set_title(title)
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=7, ncol=2)
    fig.savefig(path, dpi=110, bbox_inches="tight")
    print(f"\nwrote {path}")


def plot_per_curve(data: CurveSet, fits: dict[str, Params], path: Path) -> None:
    """One panel per measured v_gs: digitised points vs each fitted curve.

    Uses the plotting.overlay_data machinery per panel. This is the visual
    knee-shape check against data/ST8L65N044M9_Id_Vds_Vgs.png (Figure 3).
    """
    v_gs = data.v_gs_grid()
    ncol = 4
    nrow = int(np.ceil(len(v_gs) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 3.4 * nrow),
                             squeeze=False)
    for ax, vg in zip(axes.ravel()[:len(v_gs)], v_gs, strict=True):
        msk = np.isclose(data.v_gs, vg)
        vd, meas = data.v_ds[msk], data.i_d[msk]
        order = np.argsort(vd)
        ax.plot(vd[order], meas[order], "o", ms=5, alpha=0.7, label="measured")
        v_ds = np.linspace(0.0, vd.max(), 300)
        for name, p in fits.items():
            ax.plot(v_ds, vdmos_id(v_ds, vg, p), "--", alpha=0.9, label=name)
        ax.set_title(f"v_gs = {vg:g} V")
        ax.set_xlabel("v_ds [V]")
        ax.set_ylabel("i_d [A]")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    for ax in axes.ravel()[len(v_gs):]:
        ax.set_visible(False)
    # Generous vertical gaps: the per-panel titles were colliding with the
    # x-label of the row above and with the suptitle.
    fig.subplots_adjust(hspace=0.72, wspace=0.45, top=0.885, bottom=0.09)
    fig.suptitle("i_d vs v_ds at each measured v_gs: data vs k_kneesharp fits A/B",
                 y=0.995)
    fig.savefig(path, dpi=110, bbox_inches="tight")
    print(f"wrote {path}")


def main() -> None:
    data = load_csv(CSV)
    print(f"{data.label}: {len(data)} points, v_gs grid {data.v_gs_grid()}")
    print(f"free = {FREE}")
    print(f"device bounds overrides = {DEVICE_BOUNDS} (k_kneesharp keeps global (0.1, 40))")

    probe: list[str] = []
    sens = sensitivity_report(START, data, FREE)
    probe.append("pre-fit sensitivity_report (mean |d ln i_d / d ln p| at START):\n"
                 + "\n".join(f"  {n:<12s} {v:9.3g}" for n, v in sens.items()))
    print("\npre-fit sensitivity (start): "
          + "  ".join(f"{n}={v:.3g}" for n, v in sens.items()))

    print("\n--- fit A: full dataset, k_kneesharp free ---")
    res_a = multi_fit(data)
    report(res_a, data, "k_kneesharp A, full data")

    # The official record IS fit A: confirm the run still reproduces it.
    drift = max(abs(getattr(res_a.params, n) - getattr(OFFICIAL["params"], n))
                / max(abs(getattr(OFFICIAL["params"], n)), 1e-12)
                for n in FREE)
    print(f"official-record drift vs this run's fit A: {drift:.2e}"
          + ("  [OK]" if drift < OFFICIAL_TOL else
             f"  [STALE: > {OFFICIAL_TOL:g} -- re-promote fit A]"))

    print("\n--- fit B: WITHOUT the v_gs = 5 curve, k_kneesharp free ---")
    data_b = data.subset(~np.isclose(data.v_gs, 5.0))
    res_b = multi_fit(data_b)
    report(res_b, data_b, "k_kneesharp B, no v_gs=5")

    probe_txt = klambda_cap_probe(data_b)
    print(probe_txt)
    probe.append(probe_txt)

    print("\nprevious official (9-param, k_kneesharp=0): fit A 9.97 % | fit B 5.54 %")
    print(f"OFFICIAL now = fit A here:  fit A {res_a.rmse_rel * 100:.2f} %"
          f" | fit B (no v_gs=5) {res_b.rmse_rel * 100:.2f} %")

    probe_txt = ktheta_vsat_probe(data_b)
    print(probe_txt)
    probe.append(probe_txt)
    probe_txt = kneesharp_bound_probe(data_b, res_b.params)
    print(probe_txt)
    probe.append(probe_txt)

    for label, res in (("fit A (full data, k_kneesharp free)", res_a),
                       ("fit B (no v_gs=5, k_kneesharp free)", res_b)):
        probe.append(f"\n--- {label} ---\n{res.summary()}\n"
                     + identifiability(res, FREE).summary())

    OUT.mkdir(exist_ok=True)

    # Fitted parameters, same format as the capfit/st8l65_fitted_params files.
    txt = OUT / "st8l65_nd_params.txt"
    head = (res_a.summary()
            + "\n\n--- without v_gs = 5 curve ---\n" + res_b.summary())

    def pinned(res: FitResult) -> str:
        pb = res.on_bounds()
        return "none" if not pb else ", ".join(f"{k} {v}" for k, v in pb.items())

    tail = (
        f"\n\nRMSE_rel: full {res_a.rmse_rel * 100:.2f} %"
        f" | no-Vgs5 {res_b.rmse_rel * 100:.2f} %"
        f"\nSuperseded reference (9-param, k_kneesharp=0, scripts/fit_st8l65.py):"
        " full 9.97 % | no-Vgs5 5.54 %."
        f"\nk_kneesharp: fit A {res_a.params.k_kneesharp:.4g}"
        f" | fit B {res_b.params.k_kneesharp:.4g}"
        f" (bounds {BOUNDS['k_kneesharp']}). Interior (not on_bounds) = the data wants a"
        " knee far softer than the legacy tanh (k_kneesharp = 2.8) and softer than the"
        " old 0.5 guard allowed."
        f"\nStill riding a bound -- fit A: {pinned(res_a)}"
        f" | fit B: {pinned(res_b)}. Rank, cond and the knee-shape cluster"
        " correlations (k_p, v_ovmax, k_ovsharp, k_kneesharp, v_dbsat, k_lambda)"
        " are in st8l65_nd_probe.txt:"
        " the shape is pinned, the mechanism split is not."
        "\nrd+r_s is expected to sit BELOW the measured Rds(on) chord (34.5 mOhm):"
        " only the SUM is identifiable (r_d <-> r_s ~ -0.999)."
        "\nv_sub is anchored by the single v_gs(th) point (3.7 V, 250 uA) via the"
        " subthreshold channel -- anchor, not slope."
    )
    txt.write_text(head + tail)
    print(f"wrote {txt}")

    probe_file = OUT / "st8l65_nd_probe.txt"
    probe_file.write_text("\n\n".join(probe) + "\n")
    print(f"wrote {probe_file}")

    plot(data, {"official_A": OFFICIAL["params"], "knee_B": res_b.params},
         OUT / "st8l65_nd_fit.png")
    # Per-curve comparison: k_kneesharp fit A vs k_kneesharp fit B (official fit B NOT plotted).
    plot_per_curve(data, {"knee_A": res_a.params, "knee_B": res_b.params},
                   OUT / "st8l65_nd_percurve.png")

    # Plain overlay via plotting.overlay_data (combined family plot).
    fig, ax = plt.subplots(figsize=(8, 5.5))
    overlay_data(data, res_b.params, title="k_kneesharp fit B vs measured (overlay_data)", ax=ax)
    fig.savefig(OUT / "st8l65_nd_overlay.png", dpi=110, bbox_inches="tight")
    print(f"wrote {OUT / 'st8l65_nd_overlay.png'}")


if __name__ == "__main__":
    main()
