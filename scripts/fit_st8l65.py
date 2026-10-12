"""Fit the ST8L65N044M9 (650 V, 36-44 mOhm super-junction FET) static data.

Gate-overdrive charge cap (kcap, m) fit -- the mechanism that replaced the
removed `kappa` (drain-bias-dependent series resistance; see HANDOVER.md
section 5 and the model docstring). The cap keeps rd a plain resistor and
reproduces the bimodal gate exponent alpha (1.7 -> 0.07 across Vgs 5 -> 10).

Free set and per-device bounds follow HANDOVER section 4:
  free  = kp, vth0, ksub, lam, gds, rd, rs, vk, kcap, m
  fixed = theta = 0, vsat = 1e4 (off), eta = 0, gk = 0, vthk = 1.6,
          kk = 0.1, vov_floor = 1e-3   (all measured dead on this data,
          section 3e -- the cap subsumes theta/vsat, and 0-20 V data has
          no DIBL or kink evidence)
  bounds: rd, rs <= 36 mOhm (datasheet Rds(on)); lam <= 0.004 1/V
          (0.038 is absurd CLM for a 650 V part; the dataset stops at
          Vds = 20 V and never enters deep saturation)

Two fits are run (HANDOVER section 2, decision 5): the full 41-point set,
and one WITHOUT the Vgs = 5 curve. That curve is physically unreproducible
by any C1 model (Id = 3.937 A at BOTH Vds = 1.127 and 5.153 requires
beta = 0 over a 4.6x voltage interval, then beta = 0.118; section 3f) and
dominates the overall RMSE. The comparison shows how much of the reported
error is that curve rather than model inadequacy.

Watch vk: with lam capped, the output slope must come from gds/vk, and the
session-4 fit picked vk = 0.65 V (sigma_rel = 148, barely sensed) producing
a knee narrower than the datasheet figure shows (section 7 item 1).

    python scripts/fit_st8l65.py

Dataset note: the (0 V, 650 V, 1 uA) off-state point was removed because the
channel model has no HV leakage term -- one gds cannot give both 1 uA at 650 V
and 250 uA at 3.7 V (1.5 nS vs 68 uS). The remaining off point (Vgs = Vds =
3.7 V, 250 uA) anchors gds directly.
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

# Re-anchored start for a ~36 mOhm, Vth ~ 4.7 V part. DEFAULT is a 4 mOhm
# low-VT device and must never be used as a start for a different one.
GDS_OFF = 250e-6 / 3.7  # from the single remaining off-state point
START = Params(
    kp=40.0, vth0=4.7, ksub=0.06, kcap=2.0, m=10.0,
    theta=0.0, vsat=1e4, eta=0.0, lam=1e-3, gds=GDS_OFF,
    gk=0.0, vthk=1.6, kk=0.1, vk=10.0, rd=0.030, rs=0.005,
)

# HANDOVER section 4 free set. kcap/m are freed explicitly here (they are NOT
# in DEFAULT_FREE: the cap is off for the low-VT anchor device).
FREE = ("kp", "vth0", "ksub", "lam", "gds", "rd", "rs", "vk", "kcap", "m")

# Per-device bounds: rd/rs each capped at the datasheet Rds(on) = 36 mOhm
# (global BOUNDS stay at 0.5 Ohm), lam capped at 0.004 1/V -- the
# unconstrained fit wanted 0.038, which is absurd CLM for a 650 V part.
DEVICE_BOUNDS = {
    "rd": (0.0, 0.036),
    "rs": (0.0, 0.036),
    "lam": (1e-6, 0.004),
}

# Multi-start around the session-4 basin (kp~92, kcap~1.9, m~10, vk~0.65)
# plus deliberately different basins for vk/kcap/m, since vk is barely
# sensed (sigma_rel = 148) and one bad vk/kcap start stalled at 14 % RMSE.
MULTI_START = [
    dict(kp=90.0, vth0=4.73, ksub=0.06, lam=1e-3, kcap=1.9, m=10.0,
         vk=0.7, rd=0.026, rs=0.005),
    dict(kp=40.0, vth0=4.70, ksub=0.06, lam=1e-3, kcap=2.5, m=5.0,
         vk=2.0, rd=0.020, rs=0.010),
    dict(kp=150.0, vth0=4.60, ksub=0.10, lam=3e-3, kcap=1.2, m=20.0,
         vk=1.0, rd=0.030, rs=0.002),
    dict(kp=60.0, vth0=4.86, ksub=0.03, lam=1e-4, kcap=3.0, m=8.0,
         vk=5.0, rd=0.015, rs=0.015),
]


def multi_fit(data: CurveSet) -> FitResult:
    """Best-of-multi-start; all starts are reported so basin-hopping is visible."""
    best = None
    for s in MULTI_START:
        st = START.replace(**s)
        res = fit(data, start=st, free=FREE, weight="relative", max_nfev=8000,
                  bounds_overrides=DEVICE_BOUNDS)
        print(f"    start kp={s['kp']:8.4g} kcap={s['kcap']:6.4g} m={s['m']:6.4g}"
              f" vk={s['vk']:6.4g} rd={s['rd'] * 1e3:5.2f}m"
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
        print(f"   {name:<8s} = {getattr(p, name):.6g}")
    print(f"   rd+rs = {(p.rd + p.rs) * 1e3:.2f} mOhm"
          f"   (measured Ron chord at Vgs=10: 34.5 mOhm)")
    plateau = float(vdmos_id(19.7, 10.0, p))
    print(f"   plateau Vgs=10 @ 19.7 V: {plateau:.1f} A   (measured 240.2 A)")
    pinned = res.on_bounds()
    if pinned:
        print("   PINNED AT BOUND (not a physical result):")
        for name, which in pinned.items():
            print(f"     {name:<8s} {which}")
    ident = identifiability(res, FREE)
    print("   " + ident.summary().replace("\n", "\n   "))

    model = vdmos_id(data.vds, data.vgs, p)
    rel = (model - data.ids) / data.ids
    print("   per-curve RMSE %:", {
        vg: round(float(np.sqrt(np.mean(rel[np.isclose(data.vgs, vg)] ** 2)) * 100), 2)
        for vg in data.vgs_grid()
    })


def plot(data, fits: dict[str, Params], path: Path) -> None:
    """Left: full output curves. Right: knee zoom, where the fits disagree most."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    on = data.ids > 1e-3
    for vg in data.vgs_grid():
        m = np.isclose(data.vgs, vg) & on
        if not np.any(m):
            continue
        vd_data, id_data = data.vds[m], data.ids[m]
        vd = np.linspace(0.0, vd_data.max(), 300)
        for ax in axes:
            ax.plot(vd_data, id_data, "o", ms=4, alpha=0.65,
                    label=f"data Vgs={vg:g}" if ax is axes[0] else None)
            for name, p in fits.items():
                ax.plot(vd, vdmos_id(vd, vg, p), "--", alpha=0.85,
                        label=f"{name} Vgs={vg:g}" if ax is axes[0] else None)
    axes[0].set_xlim(0, 21)
    axes[1].set_xlim(0, 3)
    axes[1].set_ylim(0, 80)
    for ax, title in zip(axes, ["Output curves: data (dots) vs fits (dashed)",
                                "Knee zoom (Vds <= 3 V)"], strict=True):
        ax.set_xlabel("Vds [V]")
        ax.set_ylabel("Id [A]")
        ax.set_title(title)
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=7, ncol=2)
    fig.savefig(path, dpi=110, bbox_inches="tight")
    print(f"\nwrote {path}")


def main() -> None:
    data = load_csv(CSV)
    print(f"{data.label}: {len(data)} points, Vgs grid {data.vgs_grid()}")
    print(f"gds anchored from off point (3.7 V, 250 uA) = {GDS_OFF * 1e6:.2f} uS")
    print(f"free = {FREE}")
    print(f"device bounds overrides = {DEVICE_BOUNDS}")

    print("\n--- fit A: full dataset (cap variant A, lam <= 0.004, rd/rs <= 36 mOhm) ---")
    res_a = multi_fit(data)
    report(res_a, data, "cap A, full data")

    print("\n--- fit B: WITHOUT the Vgs = 5 curve (section 3f: no C1 model fits it) ---")
    data_b = data.subset(~np.isclose(data.vgs, 5.0))
    res_b = multi_fit(data_b)
    report(res_b, data_b, "cap A, no Vgs=5")

    def group_rmse(p, d, mask):
        rel = (vdmos_id(d.vds, d.vgs, p) - d.ids) / d.ids
        return float(np.sqrt(np.mean(rel[mask] ** 2)) * 100)

    print(f"\nVgs>=6 only: full-data fit {group_rmse(res_a.params, data, data.vgs >= 6):.2f} %"
          f" | no-Vgs5 fit {group_rmse(res_b.params, data, data.vgs >= 6):.2f} %")
    print("watch vk: with lam capped the output slope must come from gds/vk; a small vk"
          " with large sigma_rel means the knee shape is not pinned by this data"
          " (HANDOVER section 7 item 1).")

    OUT.mkdir(exist_ok=True)
    plot(data, {"full": res_a.params, "no_vgs5": res_b.params},
         OUT / "st8l65_capfit.png")
    txt = OUT / "st8l65_capfit_params.txt"
    head = (res_a.summary()
            + "\n\n--- without Vgs = 5 curve ---\n" + res_b.summary())
    tail = (
        f"\n\nRMSE_rel: full {res_a.rmse_rel * 100:.2f} %"
        f" | no-Vgs5 {res_b.rmse_rel * 100:.2f} %"
        "\nrd+rs is expected to sit BELOW the measured Rds(on) chord (34.5 mOhm):"
        " Rds(on) = rd + rs + intrinsic channel resistance. Only the SUM is"
        " identifiable (rd <-> rs ~ -0.999); do not quote the split as physical."
        "\nksub is free but carries no information on this dataset (one off-state"
        " point only) -- expect it to drift to a bound with sigma_rel = inf."
    )
    txt.write_text(head + tail)
    print(f"wrote {txt}")


if __name__ == "__main__":
    main()
