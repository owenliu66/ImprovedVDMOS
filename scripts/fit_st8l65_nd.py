"""Re-run the ST8L65N044M9 fit with the drain-knee softness exponent `nd` free.

Session-10 follow-up to scripts/fit_st8l65.py (handover "Next tasks" 2-4):
same fit-B recipe and DEVICE_BOUNDS, but `nd` is added to the free tuple
(it is deliberately NOT in DEFAULT_FREE) and started near 2.8, which
reproduces the legacy tanh knee (session-9 knee analysis). `nd` is fitted
in LINEAR space with lower bound 0.5 -- 0 is the "legacy tanh" sentinel and
must stay reachable only by leaving nd out of the free set, never by the
optimizer landing near 0.

Questions this run answers:
  * Does freeing nd beat the official 5.54 % (fit B) / 9.97 % (fit A)?
  * Does nd pin at 0.5 or 40 (on_bounds)? Pinned = the data does not want
    the new shape.
  * Does the theta<->vsat ~84 % correlated pair loosen now that knee width
    has its own knob? Checked with two EXPLORATORY runs on the fit-B data
    that additionally free theta and vsat: 11-param (no nd) vs 12-param
    (with nd); compare the theta<->vsat correlation from identifiability().

Outputs (out/):
  st8l65_nd_params.txt   fitted parameters, same format as the capfit file
  st8l65_nd_probe.txt    sensitivity_report + FitResult.summary() +
                         identifiability().summary() per fit
  st8l65_nd_fit.png      full curves + knee zoom, nd fit vs official fit B
  st8l65_nd_percurve.png Id vs Vds per measured Vgs, data vs nd fit vs
                         official fit B (knee-shape check vs Figure 3)

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

from vdmos_fit import CurveSet, Params, fit, load_csv, vdmos_id  # noqa: E402
from vdmos_fit.fit import FitResult, identifiability, sensitivity_report  # noqa: E402
from vdmos_fit.plotting import overlay_data  # noqa: E402

CSV = ROOT / "data" / "ST8L65N044M9_IV.csv"
OUT = ROOT / "out"

# Official 9-param reference (STATUS.md table, scripts/fit_st8l65.py).
OFFICIAL = {
    "label": "official fit B (9-param, nd=0)",
    "rmse_rel": 0.0554,
    "params": Params(
        kp=50.5805, vth0=4.44589, ksub=0.183585, kcap=2.90428, m=10.2408,
        theta=0.0, vsat=1e4, eta=0.0, lam=0.004, gk=0.0, vthk=1.6, kk=0.1,
        vk=5.5617, rd=0.0298458, rs=0.000657247,
    ),
}

# Same re-anchored start as fit_st8l65.py, plus nd near the legacy-tanh value.
START = Params(
    kp=40.0, vth0=4.7, ksub=0.06, kcap=2.0, m=10.0,
    theta=0.0, vsat=1e4, nd=2.8, eta=0.0, lam=1e-3,
    gk=0.0, vthk=1.6, kk=0.1, vk=10.0, rd=0.030, rs=0.005,
)

# 9-param free set + nd (explicitly freed -- NOT in DEFAULT_FREE).
FREE = ("kp", "vth0", "ksub", "lam", "rd", "rs", "vk", "kcap", "m", "nd")

# Identical per-device bounds to fit_st8l65.py. nd keeps its global
# (0.5, 40) linear bound; do NOT override it here.
DEVICE_BOUNDS = {
    "rd": (0.0, 0.036),
    "rs": (0.0, 0.036),
    "lam": (1e-6, 0.004),
}

# Same four basins as fit_st8l65.py, nd started near 2.8 (legacy tanh) with
# one deliberately sharper (nd=6) and one softer (nd=1.5) start.
MULTI_START = [
    dict(kp=90.0, vth0=4.73, ksub=0.06, lam=1e-3, kcap=1.9, m=10.0,
         vk=0.7, rd=0.026, rs=0.005, nd=2.8),
    dict(kp=40.0, vth0=4.70, ksub=0.06, lam=1e-3, kcap=2.5, m=5.0,
         vk=2.0, rd=0.020, rs=0.010, nd=2.8),
    dict(kp=150.0, vth0=4.60, ksub=0.10, lam=3e-3, kcap=1.2, m=20.0,
         vk=1.0, rd=0.030, rs=0.002, nd=1.5),
    dict(kp=60.0, vth0=4.86, ksub=0.03, lam=1e-4, kcap=3.0, m=8.0,
         vk=5.0, rd=0.015, rs=0.015, nd=6.0),
]


def multi_fit(data: CurveSet, free=FREE, start_base=START) -> FitResult:
    """Best-of-multi-start; all starts reported so basin-hopping is visible."""
    best = None
    for s in MULTI_START:
        st = start_base.replace(**s)
        res = fit(data, start=st, free=free, weight="relative", max_nfev=8000,
                  bounds_overrides=DEVICE_BOUNDS)
        nd_s = s.get("nd", 0.0)
        print(f"    start kp={s['kp']:8.4g} kcap={s['kcap']:6.4g} m={s['m']:6.4g}"
              f" vk={s['vk']:6.4g} nd={nd_s:5.3g}"
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
        print(f"   {name:<8s} = {getattr(p, name):.6g}")
    print(f"   rd+rs = {(p.rd + p.rs) * 1e3:.2f} mOhm"
          f"   (measured Ron chord at Vgs=10: 34.5 mOhm)")
    plateau = float(vdmos_id(19.7, 10.0, p))
    print(f"   plateau Vgs=10 @ 19.7 V: {plateau:.1f} A   (measured 240.2 A)")
    off = float(vdmos_id(3.7, 3.7, p))
    print(f"   Vgs(th) point (3.7, 3.7): {off * 1e6:.1f} uA   (measured 250 uA)")
    pinned = res.on_bounds()
    if pinned:
        print("   PINNED AT BOUND (not a physical result):")
        for name, which in pinned.items():
            print(f"     {name:<8s} {which}")
    else:
        print("   no parameter pinned at a bound")
    ident = identifiability(res, free)
    print("   " + ident.summary().replace("\n", "\n   "))

    model = vdmos_id(data.vds, data.vgs, p)
    rel = (model - data.ids) / data.ids
    print("   per-curve RMSE %:", {
        vg: round(float(np.sqrt(np.mean(rel[np.isclose(data.vgs, vg)] ** 2)) * 100), 2)
        for vg in data.vgs_grid()
    })


def nd_bound_probe(data_b: CurveSet, p_nd: Params) -> str:
    """Is the nd=0.5 sentinel guard costing anything?

    Refit with the lower bound relaxed to 0.1. If the optimum sits just
    below 0.5 with negligible RMSE gain, the guard is harmless and the
    fitted shape is the softest sensible power-law soft-min (nd -> 0
    collapses the current, so going lower is not physical anyway).
    Also exposes knee-cluster degeneracy: different vk/kp/lam splits
    reaching the same RMSE mean the SHAPE is pinned, not the mechanism.
    """
    lines = ["", "nd lower-bound probe (fit-B data, bound relaxed to (0.1, 40)):"]
    bounds = {**DEVICE_BOUNDS, "nd": (0.1, 40.0)}
    for nd0 in (0.5, 0.3):
        res = fit(data_b, start=p_nd.replace(nd=nd0), free=FREE,
                  weight="relative", max_nfev=8000, bounds_overrides=bounds)
        q = res.params
        lines.append(
            f"  start nd={nd0:g}: RMSE {res.rmse_rel * 100:.2f} %"
            f" | nd={q.nd:.4g} vk={q.vk:.4g} kcap={q.kcap:.4g} m={q.m:.4g}"
            f" kp={q.kp:.4g} lam={q.lam:.4g}"
            + ("  [on_bounds: " + ",".join(res.on_bounds()) + "]"
               if res.on_bounds() else "")
        )
    return "\n".join(lines)


def theta_vsat_probe(data_b: CurveSet) -> str:
    """A/B the theta<->vsat correlation with and without nd, fit-B data.

    Exploratory only (theta/vsat freed on top of the 9+nd set, seeded in the
    fit-B basin). The ~84 % theta<->vsat pair was measured on the synthetic
    grid; this checks whether nd loosens it on the real data.
    """
    base = OFFICIAL["params"].replace(nd=2.8, theta=0.02, vsat=5.0, lam=1e-3)
    lines = ["", "theta<->vsat correlation probe (fit-B data, exploratory):"]
    for tag, free in (
        ("11-param, no nd", FREE[:-1] + ("theta", "vsat")),
        ("12-param, with nd", FREE + ("theta", "vsat")),
    ):
        res = fit(data_b, start=base, free=free, weight="relative",
                  max_nfev=8000, bounds_overrides=DEVICE_BOUNDS)
        ident = identifiability(res, free)
        pair = next(((a, b, c) for a, b, c in ident.worst_pairs
                     if {a, b} == {"theta", "vsat"}), None)
        lines.append(
            f"  {tag}: RMSE {res.rmse_rel * 100:.2f} %"
            f" | theta={res.params.theta:.4g} vsat={res.params.vsat:.4g}"
            f" nd={res.params.nd:.4g}"
            + (f" | corr(theta,vsat) = {pair[2]:+.3f}" if pair
               else " | theta<->vsat NOT in |corr|>0.8 list")
        )
        for a, b, c in ident.worst_pairs:
            lines.append(f"      pair {a}<->{b}: {c:+.3f}")
    return "\n".join(lines)


def plot(data: CurveSet, fits: dict[str, Params], path: Path) -> None:
    """Left: full output curves. Right: knee zoom, where nd changes the shape."""
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


def plot_per_curve(data: CurveSet, fits: dict[str, Params], path: Path) -> None:
    """One panel per measured Vgs: digitised points vs each fitted curve.

    Uses the plotting.overlay_data machinery per panel. This is the visual
    knee-shape check against data/ST8L65N044M9_Id_Vds_Vgs.png (Figure 3).
    """
    vgs = data.vgs_grid()
    ncol = 4
    nrow = int(np.ceil(len(vgs) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 3.2 * nrow),
                             squeeze=False)
    for ax, vg in zip(axes.ravel()[:len(vgs)], vgs, strict=True):
        m = np.isclose(data.vgs, vg)
        vd, meas = data.vds[m], data.ids[m]
        order = np.argsort(vd)
        ax.plot(vd[order], meas[order], "o", ms=5, alpha=0.7, label="measured")
        vds = np.linspace(0.0, vd.max(), 300)
        for name, p in fits.items():
            ax.plot(vds, vdmos_id(vds, vg, p), "--", alpha=0.9, label=name)
        ax.set_title(f"Vgs = {vg:g} V")
        ax.set_xlabel("Vds [V]")
        ax.set_ylabel("Id [A]")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    for ax in axes.ravel()[len(vgs):]:
        ax.set_visible(False)
    fig.suptitle("Id vs Vds at each measured Vgs: data vs nd fit vs official fit B")
    fig.savefig(path, dpi=110, bbox_inches="tight")
    print(f"wrote {path}")


def main() -> None:
    data = load_csv(CSV)
    print(f"{data.label}: {len(data)} points, Vgs grid {data.vgs_grid()}")
    print(f"free = {FREE}")
    print(f"device bounds overrides = {DEVICE_BOUNDS} (nd keeps global (0.5, 40))")

    probe: list[str] = []
    sens = sensitivity_report(START, data, FREE)
    probe.append("pre-fit sensitivity_report (mean |d ln Id / d ln p| at START):\n"
                 + "\n".join(f"  {n:<8s} {v:9.3g}" for n, v in sens.items()))
    print("\npre-fit sensitivity (start): "
          + "  ".join(f"{n}={v:.3g}" for n, v in sens.items()))

    print("\n--- fit A: full dataset, nd free ---")
    res_a = multi_fit(data)
    report(res_a, data, "nd A, full data")

    print("\n--- fit B: WITHOUT the Vgs = 5 curve, nd free ---")
    data_b = data.subset(~np.isclose(data.vgs, 5.0))
    res_b = multi_fit(data_b)
    report(res_b, data_b, "nd B, no Vgs=5")

    print(f"\nofficial reference: fit A 9.97 % | fit B {OFFICIAL['rmse_rel'] * 100:.2f} %")
    print(f"nd fit:             fit A {res_a.rmse_rel * 100:.2f} %"
          f" | fit B {res_b.rmse_rel * 100:.2f} %")

    probe_txt = theta_vsat_probe(data_b)
    print(probe_txt)
    probe.append(probe_txt)
    probe_txt = nd_bound_probe(data_b, res_b.params)
    print(probe_txt)
    probe.append(probe_txt)

    for label, res in (("fit A (full data, nd free)", res_a),
                       ("fit B (no Vgs=5, nd free)", res_b)):
        probe.append(f"\n--- {label} ---\n{res.summary()}\n"
                     + identifiability(res, FREE).summary())

    OUT.mkdir(exist_ok=True)

    # Fitted parameters, same format as the capfit/st8l65_fitted_params files.
    txt = OUT / "st8l65_nd_params.txt"
    head = (res_a.summary()
            + "\n\n--- without Vgs = 5 curve ---\n" + res_b.summary())
    tail = (
        f"\n\nRMSE_rel: full {res_a.rmse_rel * 100:.2f} %"
        f" | no-Vgs5 {res_b.rmse_rel * 100:.2f} %"
        f"\nofficial 9-param reference (nd=0): full 9.97 % | no-Vgs5 5.54 %."
        "\nnd pinned at LOWER 0.5 = the data wants a knee SOFTER than the legacy"
        " tanh (nd=2.8); the relaxed-bound probe lands at nd=0.47 for ~0.03 pt"
        " gain, so the 0.5 sentinel guard costs nothing. The knee-shape cluster"
        " (kp, kcap, m, nd, vk, lam) is strongly entangled: the shape is pinned,"
        " the mechanism split is not (see st8l65_nd_probe.txt)."
        "\nrd+rs is expected to sit BELOW the measured Rds(on) chord (34.5 mOhm):"
        " only the SUM is identifiable (rd <-> rs ~ -0.999)."
        "\nksub is anchored by the single Vgs(th) point (3.7 V, 250 uA) via the"
        " subthreshold channel -- anchor, not slope."
    )
    txt.write_text(head + tail)
    print(f"wrote {txt}")

    probe_file = OUT / "st8l65_nd_probe.txt"
    probe_file.write_text("\n\n".join(probe) + "\n")
    print(f"wrote {probe_file}")

    fits = {"official_B": OFFICIAL["params"], "nd_B": res_b.params}
    plot(data, fits, OUT / "st8l65_nd_fit.png")
    plot_per_curve(data, fits, OUT / "st8l65_nd_percurve.png")

    # Plain overlay via plotting.overlay_data (combined family plot).
    fig, ax = plt.subplots(figsize=(8, 5.5))
    overlay_data(data, res_b.params, title="nd fit B vs measured (overlay_data)", ax=ax)
    fig.savefig(OUT / "st8l65_nd_overlay.png", dpi=110, bbox_inches="tight")
    print(f"wrote {OUT / 'st8l65_nd_overlay.png'}")


if __name__ == "__main__":
    main()
