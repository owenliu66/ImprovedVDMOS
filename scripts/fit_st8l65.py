"""Fit the ST8L65N044M9 (650 V, 36-44 mOhm super-junction FET) static data.

**SUPERSEDED as the official reference (2026-10-13)**: the quoted ST8L65 fit is
now fit A of `scripts/fit_st8l65_nd.py` (10-param with the drain-knee softness
exponent `k_kneesharp` free, 3.79 % vs this script's 5.54 %). Kept as the 9-param
k_kneesharp = 0 baseline and for the kappa -> cap history; do not quote its numbers as
the device model.

Gate-overdrive charge cap (v_ovmax, k_ovsharp) fit -- the mechanism that replaced the
removed `kappa` (drain-bias-dependent series resistance; see STATUS.md and
the model docstring). The cap keeps r_d a plain resistor and reproduces the
bimodal gate exponent alpha (1.7 -> 0.07 across v_gs 5 -> 10).

Free set and per-device bounds (see the STATUS.md fit table):
  free  = k_p, v_th0, v_sub, k_lambda, r_d, r_s, v_dbsat, v_ovmax, k_ovsharp
  fixed = k_theta = 0, v_vsat = 1e4 (off), k_eta = 0, g_k = 0, v_thk = 1.6,
          v_kovsoft = 0.1, v_kneemin = 1e-3   (all measured dead on this data --
          the cap subsumes k_theta/v_vsat, and 0-20 V data has no DIBL or kink
          evidence)
  bounds: r_d, r_s <= 36 mOhm (datasheet Rds(on)); k_lambda <= 0.004 1/V
          (0.038 is absurd CLM for a 650 V part; the dataset stops at
          v_ds = 20 V and never enters deep saturation)

Two fits are run: the full 41-point set, and one WITHOUT the v_gs = 5 curve.
That curve is physically unreproducible by any C1 model (i_d = 3.937 A at
BOTH v_ds = 1.127 and 5.153 requires beta = 0 over a 4.6x voltage interval,
then beta = 0.118; see data/README.md) and dominates the overall RMSE. The
comparison shows how much of the reported error is that curve rather than
model inadequacy.

Watch v_dbsat: with k_lambda capped, the output slope must come from the channel k_lambda
term alone, and v_dbsat is barely sensed (sigma_rel ~ 110-150) -- the knee shape
is not pinned by this data; check the fit against data/ST8L65N044M9_Id_Vds_Vgs.png
(datasheet Figure 3).

    python scripts/fit_st8l65.py

Dataset note: the (0 V, 650 V, 1 uA) off-state point was removed because the
channel model has no HV leakage element (real leakage waits for the body
diode). The remaining off point (v_gs = v_ds = 3.7 V, 250 uA) is the v_gs(th)
spec condition and anchors v_sub/v_th0 through the subthreshold channel.
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

from vdmos_fit import CurveSet, Params, fit, load_csv, vdmos_id  # noqa: E402
from vdmos_fit.fit import FitResult, identifiability  # noqa: E402

CSV = ROOT / "data" / "ST8L65N044M9_IV.csv"
OUT = ROOT / "out"

# Re-anchored start for a ~36 mOhm, v_th ~ 4.7 V part. DEFAULT is a 4 mOhm
# low-VT device and must never be used as a start for a different one.
START = Params(
    k_p=40.0, v_th0=4.7, v_sub=0.06, v_ovmax=2.0, k_ovsharp=10.0,
    k_theta=0.0, v_vsat=1e4, k_eta=0.0, k_lambda=1e-3,
    g_k=0.0, v_thk=1.6, v_kovsoft=0.1, v_dbsat=10.0, r_d=0.030, r_s=0.005,
)

# 9-param free set (see the STATUS.md fit table). v_ovmax/k_ovsharp are freed
# explicitly here (they are NOT in DEFAULT_FREE: the cap is off for the
# low-VT anchor device).
FREE = ("k_p", "v_th0", "v_sub", "k_lambda", "r_d", "r_s", "v_dbsat", "v_ovmax", "k_ovsharp")

# Per-device bounds: r_d/r_s each capped at the datasheet Rds(on) = 36 mOhm
# (global BOUNDS stay at 0.5 Ohm), k_lambda capped at 0.004 1/V -- the
# unconstrained fit wanted 0.038, which is absurd CLM for a 650 V part.
DEVICE_BOUNDS = {
    "r_d": (0.0, 0.036),
    "r_s": (0.0, 0.036),
    "k_lambda": (1e-6, 0.004),
}

# Multi-start around the known basin (k_p~92, v_ovmax~1.9, k_ovsharp~10, v_dbsat~0.65)
# plus deliberately different basins for v_dbsat/v_ovmax/k_ovsharp, since v_dbsat is barely
# sensed (sigma_rel = 148) and one bad v_dbsat/v_ovmax start stalled at 14 % RMSE.
MULTI_START = [
    dict(k_p=90.0, v_th0=4.73, v_sub=0.06, k_lambda=1e-3, v_ovmax=1.9, k_ovsharp=10.0,
         v_dbsat=0.7, r_d=0.026, r_s=0.005),
    dict(k_p=40.0, v_th0=4.70, v_sub=0.06, k_lambda=1e-3, v_ovmax=2.5, k_ovsharp=5.0,
         v_dbsat=2.0, r_d=0.020, r_s=0.010),
    dict(k_p=150.0, v_th0=4.60, v_sub=0.10, k_lambda=3e-3, v_ovmax=1.2, k_ovsharp=20.0,
         v_dbsat=1.0, r_d=0.030, r_s=0.002),
    dict(k_p=60.0, v_th0=4.86, v_sub=0.03, k_lambda=1e-4, v_ovmax=3.0, k_ovsharp=8.0,
         v_dbsat=5.0, r_d=0.015, r_s=0.015),
]


def multi_fit(data: CurveSet) -> FitResult:
    """Best-of-multi-start; all starts are reported so basin-hopping is visible."""
    best = None
    for s in MULTI_START:
        st = START.replace(**s)
        res = fit(data, start=st, free=FREE, weight="relative", max_nfev=8000,
                  bounds_overrides=DEVICE_BOUNDS)
        print(f"    start k_p={s['k_p']:8.4g} v_ovmax={s['v_ovmax']:6.4g}"
              f" k_ovsharp={s['k_ovsharp']:6.4g} v_dbsat={s['v_dbsat']:6.4g}"
              f" r_d={s['r_d'] * 1e3:5.2f}m"
              f" -> RMSE {res.rmse_rel * 100:6.2f} %"
              + ("" if res.success else "   [stopped]"))
        if best is None or res.rmse_rel < best.rmse_rel:
            best = res
    return best


def report(res: FitResult, data, label: str) -> None:
    p = res.params
    state = "converged" if res.success else "STOPPED"
    print(f"\n=== {label}: RMSE_rel = {res.rmse_rel * 100:.2f} % "
          f"({res.n_points} pts, {res.n_free} free, {state}) ===")
    for name in FREE:
        print(f"   {name:<12s} = {getattr(p, name):.6g}")
    print(f"   r_d+r_s = {(p.r_d + p.r_s) * 1e3:.2f} mOhm"
          f"   (measured Ron chord at v_gs=10: 34.5 mOhm)")
    plateau = float(vdmos_id(19.7, 10.0, p))
    print(f"   plateau v_gs=10 @ 19.7 V: {plateau:.1f} A   (measured 240.2 A)")
    off = float(vdmos_id(3.7, 3.7, p))
    print(f"   v_gs(th) point (3.7, 3.7): {off * 1e6:.1f} uA   (measured 250 uA)"
          "   <- anchors v_sub/v_th0 via the subthreshold channel")
    pinned = res.on_bounds()
    if pinned:
        print("   PINNED AT BOUND (not a physical result):")
        for name, which in pinned.items():
            print(f"     {name:<12s} {which}")
    ident = identifiability(res, FREE)
    print("   " + ident.summary().replace("\n", "\n   "))

    model = vdmos_id(data.v_ds, data.v_gs, p)
    rel = (model - data.i_d) / data.i_d
    print("   per-curve RMSE %:", {
        vg: round(float(np.sqrt(np.mean(rel[np.isclose(data.v_gs, vg)] ** 2)) * 100), 2)
        for vg in data.v_gs_grid()
    })


def plot(data, fits: dict[str, Params], path: Path) -> None:
    """Left: full output curves. Right: knee zoom, where the fits disagree most."""
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


def main() -> None:
    data = load_csv(CSV)
    print(f"{data.label}: {len(data)} points, v_gs grid {data.v_gs_grid()}")
    print(f"free = {FREE}")
    print(f"device bounds overrides = {DEVICE_BOUNDS}")

    print("\n--- fit A: full dataset (k_lambda <= 0.004, r_d/r_s <= 36 mOhm) ---")
    res_a = multi_fit(data)
    report(res_a, data, "cap A, full data")

    print("\n--- fit B: WITHOUT the v_gs = 5 curve ---")
    data_b = data.subset(~np.isclose(data.v_gs, 5.0))
    res_b = multi_fit(data_b)
    report(res_b, data_b, "cap A, no v_gs=5")

    def group_rmse(p, d, mask):
        rel = (vdmos_id(d.v_ds, d.v_gs, p) - d.i_d) / d.i_d
        return float(np.sqrt(np.mean(rel[mask] ** 2)) * 100)

    print(f"\nVgs>=6 only: full-data fit {group_rmse(res_a.params, data, data.v_gs >= 6):.2f} %"
          f" | no-Vgs5 fit {group_rmse(res_b.params, data, data.v_gs >= 6):.2f} %")
    print("watch v_dbsat: with k_lambda capped the output slope must come from the channel"
          " (k_lambda term only); a small v_dbsat"
          " with large sigma_rel means the knee shape is not pinned by this data"
          " -- check against data/ST8L65N044M9_Id_Vds_Vgs.png (Figure 3).")

    OUT.mkdir(exist_ok=True)
    plot(data, {"full": res_a.params, "no_vgs5": res_b.params},
         OUT / "st8l65_capfit.png")
    txt = OUT / "st8l65_capfit_params.txt"
    head = (res_a.summary()
            + "\n\n--- without v_gs = 5 curve ---\n" + res_b.summary())
    tail = (
        f"\n\nRMSE_rel: full {res_a.rmse_rel * 100:.2f} %"
        f" | no-Vgs5 {res_b.rmse_rel * 100:.2f} %"
        "\nrd+r_s is expected to sit BELOW the measured Rds(on) chord (34.5 mOhm):"
        " Rds(on) = r_d + r_s + intrinsic channel resistance. Only the SUM is"
        " identifiable (r_d <-> r_s ~ -0.999); do not quote the split as physical."
        "\nv_sub is anchored by the single v_gs(th) point (3.7 V, 250 uA) via the"
        " subthreshold channel -- the anchor, not the slope, is what this data"
        " pins; a near-threshold transfer curve would pin the swing."
    )
    txt.write_text(head + tail)
    print(f"wrote {txt}")


if __name__ == "__main__":
    main()
