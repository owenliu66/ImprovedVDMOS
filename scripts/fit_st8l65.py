"""Fit the ST8L65N044M9 (650 V, 36-44 mOhm super-junction FET) static data.

Compares the plain rd/rs embedding (kappa = 0) against the drain-bias-dependent
drain resistance (kappa > 0, super-junction JFET/RESURF effect), and reports
identifiability + per-curve residuals.

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
    kp=40.0, vth0=4.7, ksub=0.06, theta=0.10, vsat=100.0, eta=0.01,
    lam=1e-6, gds=GDS_OFF, gk=0.0, vthk=1.6, kk=0.1, vk=10.0,
    rd=0.034, rs=0.005, kappa=0.0,
)
BASE_FREE = ("kp", "vth0", "ksub", "theta", "eta", "gds", "rd", "rs")
DRIFT_FREE = BASE_FREE + ("kappa",)
# Everything the model has, except vov_floor (a numerical floor, never fitted).
ALL_FREE = (
    "kp", "vth0", "ksub", "theta", "vsat", "eta", "lam", "gds", "gk",
    "vthk", "kk", "vk", "rd", "rs", "kappa",
)

# Wider search box for the unconstrained run: the physical defaults in BOUNDS
# are what pinned ksub/theta/vsat/lam in the constrained fits, so a "no
# constraints" fit must widen them or it just reproduces the same corners.
WIDE_BOUNDS = {
    "kp": (1e-3, 1e5),
    "vth0": (-2.0, 12.0),
    "ksub": (1e-4, 1.0),
    "theta": (1e-8, 20.0),
    "vsat": (1e-3, 1e4),
    "eta": (1e-8, 1.0),
    "lam": (1e-10, 20.0),
    "gds": (1e-14, 10.0),
    "gk": (1e-15, 1e3),
    "vthk": (-20.0, 20.0),
    "kk": (1e-4, 20.0),
    "vk": (1e-3, 1e4),
    "rd": (0.0, 0.5),
    "rs": (0.0, 0.5),
    "kappa": (0.0, 20.0),
}

# vsat / lam are pinned at "off" (100 V / 1e-6 1/V): the drift resistance, not
# channel velocity saturation or CLM, sets the output characteristic here.
MULTI_START = [
    dict(kp=40.0, theta=0.10, rd=0.034, rs=0.005, vth0=4.70, ksub=0.06, eta=0.01, kappa=0.07),
    dict(kp=100.0, theta=0.50, rd=0.030, rs=0.010, vth0=4.60, ksub=0.12, eta=0.03, kappa=0.10),
    dict(kp=20.0, theta=0.05, rd=0.036, rs=0.002, vth0=4.86, ksub=0.03, eta=0.005, kappa=0.03),
    dict(kp=300.0, theta=1.50, rd=0.025, rs=0.015, vth0=4.50, ksub=0.20, eta=0.02, kappa=0.15),
    dict(kp=60.0, theta=0.20, rd=0.034, rs=0.005, vth0=4.75, ksub=0.08, eta=0.001, kappa=0.30),
]
# Extra starts that explore the parts of the box the constrained fits could not
# reach (large vsat/lam, active kink, small vk) -- used only by the free fit.
EXTRA_START = [
    dict(kp=500.0, theta=0.02, vsat=3.0, lam=0.5, rd=0.020, rs=0.010, vth0=4.40,
         ksub=0.10, eta=0.05, kappa=0.05, gk=0.5, vthk=4.0, kk=0.3, vk=2.0),
    dict(kp=15.0, theta=5.0, vsat=1e3, lam=1e-4, rd=0.040, rs=0.001, vth0=5.00,
         ksub=0.01, eta=1e-4, kappa=1.0, gk=1e-6, vthk=0.5, kk=0.02, vk=50.0),
    dict(kp=1e3, theta=0.5, vsat=10.0, lam=0.05, rd=0.010, rs=0.030, vth0=4.20,
         ksub=0.30, eta=0.20, kappa=0.01, gk=5.0, vthk=6.0, kk=1.0, vk=0.5),
]


def multi_fit(data: CurveSet, free: tuple[str, ...], *, wide: bool = False) -> FitResult:
    """Best-of-multi-start; all starts are reported so basin-hopping is visible."""
    best = None
    starts = MULTI_START + (EXTRA_START if wide else [])
    for s in starts:
        st = START.replace(**{k: v for k, v in s.items() if k != "kappa"})
        if "kappa" in free:
            st = st.replace(kappa=s["kappa"])
        res = fit(data, start=st, free=free, weight="relative", max_nfev=8000,
                  bounds_overrides=WIDE_BOUNDS if wide else None)
        print(f"    start kp={s['kp']:8.4g} theta={s['theta']:7.4g} kappa={s['kappa']:6.4g}"
              f" vth0={s['vth0']:5.2f} -> RMSE {res.rmse_rel * 100:6.2f} %"
              + ("" if res.success else "   [stopped]"))
        if best is None or res.rmse_rel < best.rmse_rel:
            best = res
    return best


def report(res: FitResult, data, free: tuple[str, ...], label: str) -> None:
    p = res.params
    state = "converged" if res.success else "STOPPED"
    print(f"\n=== {label}: RMSE_rel = {res.rmse_rel * 100:.2f} % "
          f"({res.n_points} pts, {res.n_free} free, {state}) ===")
    for name in free:
        print(f"   {name:<8s} = {getattr(p, name):.6g}")
    print(f"   rd+rs = {(p.rd + p.rs) * 1e3:.2f} mOhm")
    if p.kappa > 0:
        print(f"   drift plateau Id_max = 1/(rd*kappa) = {1.0 / (p.rd * p.kappa):.1f} A"
              f"   (R_d at 19.7 V = {p.rd * (1 + p.kappa * 19.7) * 1e3:.1f} mOhm)")
    pinned = res.on_bounds()
    if pinned:
        print("   PINNED AT BOUND (not a physical result):")
        for name, which in pinned.items():
            print(f"     {name:<8s} {which}")
    ident = identifiability(res, free)
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

    print("\n--- control: kappa = 0 (plain rd/rs embedding) ---")
    res0 = multi_fit(data, BASE_FREE)
    report(res0, data, BASE_FREE, "kappa = 0")

    print("\n--- drift: kappa free (SJ JFET/RESURF drain resistance) ---")
    res1 = multi_fit(data, DRIFT_FREE)
    report(res1, data, DRIFT_FREE, "kappa free")

    print("\n--- FREE: all 15 parameters, widened bounds (no hand-pinning) ---")
    res2 = multi_fit(data, ALL_FREE, wide=True)
    report(res2, data, ALL_FREE, "fully free")

    p0, p1, p2 = res0.params, res1.params, res2.params
    def group_rmse(p, mask):
        rel = (vdmos_id(data.vds, data.vgs, p) - data.ids) / data.ids
        return float(np.sqrt(np.mean(rel[mask] ** 2)) * 100)
    print(f"\nimprovement: kappa=0 {res0.rmse_rel * 100:.2f} %"
          f" -> kappa free {res1.rmse_rel * 100:.2f} %"
          f" -> fully free {res2.rmse_rel * 100:.2f} %")
    print(f"Vgs>=6 only: kappa=0 {group_rmse(p0, data.vgs >= 6):.2f} %"
          f" -> kappa free {group_rmse(p1, data.vgs >= 6):.2f} %"
          f" -> fully free {group_rmse(p2, data.vgs >= 6):.2f} %")

    OUT.mkdir(exist_ok=True)
    plot(data, {"k=0": p0, "k>0": p1, "free": p2}, OUT / "st8l65_driftfit.png")
    txt = OUT / "st8l65_fitted_params.txt"
    head = res2.summary()
    tail = (
        f"\n\nRMSE_rel: kappa=0 {res0.rmse_rel * 100:.2f} %"
        f" | kappa free {res1.rmse_rel * 100:.2f} %"
        f" | fully free {res2.rmse_rel * 100:.2f} %"
        "\nrd+rs is expected to sit BELOW the measured Rds(on) chord (34.5 mOhm):"
        " Rds(on) = rd + rs + intrinsic channel resistance, and rd/rs are the"
        " package (wirebonds, pins) only."
        "\nCAVEAT: kappa multiplies rd, which is defined as package-only; the"
        " element that physically depletes is the drift region. A separate"
        " rdrift(kappa) element would be cleaner -- see STATUS.md."
    )
    txt.write_text(head + tail)
    print(f"wrote {txt}")


if __name__ == "__main__":
    main()
