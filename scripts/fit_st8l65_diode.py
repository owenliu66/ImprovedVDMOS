"""Extract the ST8L65N044M9 body-diode parameters (LTspice diode DC branch).

Data: data/ST8L65N044M9_diode.csv -- 30 points, I_sd 5-55 A at -55/25/150 degC
(datasheet body-diode forward characteristics).

PROTOCOL (locked after the probe campaign, README "Body diode"):

  * Fit A (OFFICIAL, session 15): free = i_s, r_diode, v_eg, i_kf, a_ikf,
    a_rdiode, b_rdiode; k_n FIXED at 2.0 (physical diffusion limit), k_xti
    fixed at 3.0 (LTspice default), i_sr = 0. Four multi-starts, r_diode
    capped at 20 mOhm (DEVICE_BOUNDS). RMSE 3.91 %.
  * Why the TCs (a_ikf/a_rdiode/b_rdiode, LTspice Tikf/Trs1/Trs2 laws) are
    in the primary set: TC-off fit A (session-12 record) sits at 4.78 %,
    linear-TC-only (session-14 record) at 4.15 %. The Trs2 promotion was a
    deliberate user override: the gain is modest and the five-knob
    r_diode/i_kf/a_ikf/a_rdiode/b_rdiode block is strongly correlated
    (|corr| up to 0.985, sigma(i_kf) 15) -- quote the CURVE, never the
    split. The stubborn 150 degC residual is a ~12 mV digitization error on
    the 5 A point (exponential region -- r_diode carries only ~10 mV of the
    562 mV drop), not missing series-R physics; excluding that point the
    150 degC rms is ~3.8 % for both the 6- and 7-param fits.
    (identifiability prints rank 6/7: an artifact of b_rdiode's ~1e4 column
    norm inflating the relative dead tolerance; unit-column SVD is 7/7.)
    Baselines live on as probes: tc_off_probe reproduces 4.78 % (session
    12), trs2_off_probe reproduces 4.15 % (session 14).
  * Why k_n fixed: freeing it buys ~1 % RMSE and lands at k_n ~ 1.3-2.1
    (the valley's direction moves with the TCs) with corr(i_s, k_n) ~ +1.0
    and a reverse-blocking floor that drifts by ORDERS of magnitude -- the
    ideality/prefactor split is a degenerate valley. Fix k_n from physics
    (diffusion limit ~2), not from this data.
  * Why i_kf in the primary set: with i_kf pinned off (sentinel 0.0) the fit
    degrades badly. The high-injection knee carries real curvature in the
    5-55 A range.
  * i_sr/k_nr stay off: freeing them changes nothing (unsensed above 5 A).
  * Trs2 (b_rdiode, quadratic r_diode TC) promoted session 15 (user
    override): see the Fit A bullet above for the caveats. See trs2_off_probe
    for the 4.15 % linear-TC baseline.
  * i_s is only meaningful PAIRED with k_n: the data has no low-current tail,
    so sigma_rel(i_s) ~ 1 (order of magnitude). Quote the curve, not i_s.

Temperature law (LTspice Xti/Eg, divisor N per branch) is MANDATORY, not
optional: a 25 degC-only fit predicts -55 degC out-of-sample at ~260 % error.
v_eg lands at ~1.10-1.12 eV -- silicon's bandgap, genuinely sensed
(sigma_rel ~ 0.06).

Outputs (out/):
  st8l65_diode_params.txt   fitted parameters + per-temperature RMSE
  st8l65_diode_probe.txt    sensitivity + identifiability + probes
  st8l65_diode_fit.png      data vs fit at the three temperatures + residuals

    python scripts/fit_st8l65_diode.py
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
    ST8L65_DIODE,
    ST8L65_DIODE_FREE,
    ST8L65_DIODE_RMSE_REL,
    DiodeCurveSet,
    Params,
    diode_sensitivity_report,
    fit_diode,
    load_diode_csv,
    vdmos_diode,
)
from vdmos_fit.fit import FitResult, identifiability  # noqa: E402

CSV = ROOT / "data" / "ST8L65N044M9_diode.csv"
OUT = ROOT / "out"

# OFFICIAL reference -- lives in the package (vdmos_fit.reference) so the test
# suite can pin it too. main() re-derives fit A and checks drift (OFFICIAL_TOL).
OFFICIAL_TOL = 1e-3

FREE = ST8L65_DIODE_FREE  # ("i_s", "r_diode", "v_eg", "i_kf", "a_ikf", "a_rdiode")

# Base start: k_n at the physical diffusion limit, k_xti at the LTspice
# default, i_sr off. Channel parameters are irrelevant to the diode fit.
START = Params(
    i_s=1e-8, k_n=2.0, i_kf=20.0, r_diode=0.002, t_nom=25.0, k_xti=3.0, v_eg=1.11,
)

# Per-device bound: the measured dV/dI slope of the diode data is 3.4-11.7
# mOhm and r_diode carries NO channel current, so 20 mOhm is a generous
# ceiling; the global 0.5 Ohm bound admits a runaway basin (r_diode ~ 0.34
# Ohm with similar RMSE) that is not a physical diode.
DEVICE_BOUNDS = {"r_diode": (0.0, 0.020)}

MULTI_START = [
    dict(i_s=1e-8, r_diode=0.002, v_eg=1.11, i_kf=20.0,
         a_ikf=0.0, a_rdiode=0.0, b_rdiode=0.0),
    dict(i_s=1e-6, r_diode=0.005, v_eg=1.20, i_kf=1.0,
         a_ikf=0.004, a_rdiode=0.004, b_rdiode=2e-5),
    dict(i_s=1e-12, r_diode=0.0005, v_eg=0.30, i_kf=1e4,
         a_ikf=-0.002, a_rdiode=0.010, b_rdiode=-2e-5),
    dict(i_s=1e-4, r_diode=0.001, v_eg=1.10, i_kf=50.0,
         a_ikf=0.010, a_rdiode=0.001, b_rdiode=1e-4),
]


def multi_fit(data: DiodeCurveSet, free=FREE, base=START) -> FitResult:
    """Best-of-multi-start; all starts reported so basin-hopping is visible."""
    best = None
    for s in MULTI_START:
        st = base.replace(**s)
        res = fit_diode(data, start=st, free=free, weight="relative",
                        max_nfev=8000, bounds_overrides=DEVICE_BOUNDS)
        print(f"    start i_s={s['i_s']:8.3g} r_diode={s['r_diode']*1e3:6.3g}m"
              f" v_eg={s['v_eg']:6.3g} i_kf={s['i_kf']:8.3g}"
              f" a=({s['a_ikf']:+.3f},{s['a_rdiode']:+.3f})"
              f" b={s['b_rdiode']:+.1e}"
              f" -> RMSE {res.rmse_rel * 100:6.2f} %"
              + ("" if res.success else "   [stopped]"))
        if best is None or res.rmse_rel < best.rmse_rel:
            best = res
    return best


def reverse_floor(p: Params, v_rev: float = 1.0, t_c: float = 25.0) -> float:
    """Reverse-blocking current magnitude at -v_rev V drain (diode anode-cathode)."""
    return abs(float(vdmos_diode(-v_rev, p, t_c=t_c)))


def report(res: FitResult, data: DiodeCurveSet, label: str, free=FREE) -> None:
    p = res.params
    state = "converged" if res.success else "STOPPED"
    print(f"\n=== {label}: RMSE_rel = {res.rmse_rel * 100:.2f} % "
          f"({res.n_points} pts, {res.n_free} free, {state}) ===")
    for name in free:
        print(f"   {name:<12s} = {getattr(p, name):.6g}")
    print(f"   fixed: k_n={p.k_n:g} k_xti={p.k_xti:g} t_nom={p.t_nom:g} "
          f"i_sr={p.i_sr:g}")
    model = vdmos_diode(data.v_sd, p, t_c=data.temp_c)
    rel = (model - data.i_sd) / data.i_sd
    for t in data.temp_grid():
        msk = np.isclose(data.temp_c, t)
        print(f"   {t:6.1f} C: rms {np.sqrt(np.mean(rel[msk] ** 2)) * 100:5.2f} %"
              f"  range [{rel[msk].min() * 100:+.1f} %, {rel[msk].max() * 100:+.1f} %]")
    print(f"   reverse floor @ -1 V, 25 C: {reverse_floor(p):.3g} A")
    pinned = res.on_bounds()
    if pinned:
        print("   PINNED AT BOUND (not a physical result):")
        for name, which in pinned.items():
            print(f"     {name:<12s} {which}")
    else:
        print("   no parameter pinned at a bound")
    ident = identifiability(res, free)
    print("   " + ident.summary().replace("\n", "\n   "))


def kn_sweep_probe(data: DiodeCurveSet) -> str:
    """Why k_n is FIXED at 2.0: sweep it as a fixed value, i_kf off.

    With i_kf off, low k_n cannot fit the curvature at all (41 % at k_n = 1)
    and the fit only improves as k_n -> 3 -- i.e. ideality absorbs the
    high-injection bending. That is the wrong mechanism: it inflates the
    reverse floor by orders of magnitude (printed per row).
    """
    lines = ["", "k_n sensitivity sweep (free = i_s, r_diode, v_eg; k_xti=3, i_kf off):"]
    for kn in (1.0, 1.5, 2.0, 2.5, 3.0):
        base = START.replace(k_n=kn, i_kf=0.0)
        res = fit_diode(data, start=base, free=("i_s", "r_diode", "v_eg"),
                        weight="relative", max_nfev=8000,
                        bounds_overrides=DEVICE_BOUNDS)
        p = res.params
        lines.append(
            f"  k_n={kn:3.1f}: RMSE {res.rmse_rel * 100:6.2f} %"
            f" | i_s={p.i_s:.3g} r_diode={p.r_diode * 1e3:.3g}m v_eg={p.v_eg:.4g}"
            f" | reverse floor {reverse_floor(p):.2g} A"
        )
    return "\n".join(lines)


def free_kn_probe(data: DiodeCurveSet) -> str:
    """Free k_n on top of fit A: does the data want a different k_n, and at what cost?

    Expect a degenerate valley: with the TCs free the fit lands at k_n ~ 1.3
    (it landed at 2.07 in the TC-off protocol) with i_s <-> k_n correlation
    ~ +1.000 and a reverse floor that drifts by ORDERS of magnitude (nA here,
    0.1-0.4 mA TC-off) -- the direction moves, the unconstrained floor does
    not. The ideality/high-injection/TC split is a valley; only the PAIR
    (k_n, i_s) and the curve are pinned. k_n stays fixed from physics.
    """
    lines = ["", "free-k_n probe (fit A + k_n + k_xti free):"]
    free = FREE + ("k_n", "k_xti")
    res = fit_diode(data, start=START, free=free, weight="relative",
                    max_nfev=8000, bounds_overrides=DEVICE_BOUNDS)
    p = res.params
    ident = identifiability(res, free)
    pair = next(((a, b, c) for a, b, c in ident.worst_pairs if {a, b} == {"i_s", "k_n"}), None)
    lines.append(
        f"  RMSE {res.rmse_rel * 100:.2f} % | k_n={p.k_n:.4g} i_s={p.i_s:.3g}"
        f" i_kf={p.i_kf:.4g} r_diode={p.r_diode * 1e3:.3g}m"
        f" | reverse floor {reverse_floor(p):.2g} A"
        + (f" | corr(i_s,k_n) = {pair[2]:+.3f}" if pair else " | i_s<->k_n not in |corr|>0.8")
    )
    return "\n".join(lines)


def mechanism_probes(data: DiodeCurveSet, p_a: Params) -> str:
    """i_kf off and i_sr/k_nr free -- the two mechanisms kept out of fit A."""
    lines = ["", "mechanism probes (fit A data):"]
    # i_kf = 0.0 is the OFF sentinel; the 1e-3 lower bound is NOT off (an
    # active knee at 1 mA crushes the multi-ampere current).
    res = fit_diode(data, start=p_a.replace(i_kf=0.0),
                    free=("i_s", "r_diode", "v_eg"), weight="relative",
                    max_nfev=8000, bounds_overrides=DEVICE_BOUNDS)
    lines.append(
        f"  i_kf pinned off (k_n=2): RMSE {res.rmse_rel * 100:.2f} %"
        f"  <- why i_kf is in the primary set"
    )
    res = fit_diode(data, start=p_a.replace(i_sr=1e-7),
                    free=FREE + ("i_sr", "k_nr"), weight="relative",
                    max_nfev=8000, bounds_overrides=DEVICE_BOUNDS)
    ident = identifiability(res, FREE + ("i_sr", "k_nr"))
    lines.append(
        f"  i_sr/k_nr freed: RMSE {res.rmse_rel * 100:.2f} %"
        f" | i_sr={res.params.i_sr:.3g} k_nr={res.params.k_nr:.3g}"
        f" | sigma_rel(i_sr)={ident.sigma_rel['i_sr']:.2g}"
        "  <- unsensed, keep off"
    )
    return "\n".join(lines)


def tc_off_probe(data: DiodeCurveSet, p_a: Params) -> str:
    """TCs pinned off (the session-12 record protocol): why they are in the
    primary set now. Linear TCs (LTspice Tikf/Trs1) at 0 = no temperature
    dependence; the fit degrades and the per-temperature residuals tilt.
    """
    lines = ["", "TC-off probe (fit A minus a_ikf/a_rdiode -- the session-12 "
                  "record protocol):"]
    res = fit_diode(data, start=p_a.replace(a_ikf=0.0, a_rdiode=0.0),
                    free=("i_s", "r_diode", "v_eg", "i_kf"), weight="relative",
                    max_nfev=8000, bounds_overrides=DEVICE_BOUNDS)
    model = vdmos_diode(data.v_sd, res.params, t_c=data.temp_c)
    rel = (model - data.i_sd) / data.i_sd
    per_t = "  ".join(
        f"{t:+.0f}C {np.sqrt(np.mean(rel[np.isclose(data.temp_c, t)] ** 2)) * 100:.2f}%"
        for t in data.temp_grid())
    lines.append(
        f"  TC off: RMSE {res.rmse_rel * 100:.2f} % | per-T rms: {per_t}"
        "  <- why the TCs are in the primary set"
    )
    return "\n".join(lines)


def trs2_off_probe(data: DiodeCurveSet, p_a: Params) -> str:
    """b_rdiode pinned at 0 (the session-14 linear-TC record protocol): the
    Trs2 baseline. Reproduces 4.15 %; the promoted 7-param fit A gains
    ~0.24 % absolute on top of it.
    """
    lines = ["", "Trs2-off probe (fit A minus b_rdiode -- the session-14 "
                  "linear-TC record protocol):"]
    free = tuple(n for n in FREE if n != "b_rdiode")
    res = fit_diode(data, start=p_a.replace(b_rdiode=0.0), free=free,
                    weight="relative", max_nfev=8000,
                    bounds_overrides=DEVICE_BOUNDS)
    model = vdmos_diode(data.v_sd, res.params, t_c=data.temp_c)
    rel = (model - data.i_sd) / data.i_sd
    per_t = "  ".join(
        f"{t:+.0f}C {np.sqrt(np.mean(rel[np.isclose(data.temp_c, t)] ** 2)) * 100:.2f}%"
        for t in data.temp_grid())
    lines.append(
        f"  b_rdiode pinned 0: RMSE {res.rmse_rel * 100:.2f} % | per-T rms: {per_t}"
        "  <- session-14 baseline; Trs2 gain is modest (5 A/150 C point is "
        "digitization, not Rs curvature)"
    )
    return "\n".join(lines)


def plot(data: DiodeCurveSet, p: Params, path: Path) -> None:
    """Left: measured points + fit curves at the three temperatures.
    Right: relative residuals vs i_sd, one marker per temperature."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    model = vdmos_diode(data.v_sd, p, t_c=data.temp_c)
    rel = (model - data.i_sd) / data.i_sd
    styles = {-55.0: "s", 25.0: "o", 150.0: "^"}
    for t in data.temp_grid():
        msk = np.isclose(data.temp_c, t)
        i_meas, v_meas = data.i_sd[msk], data.v_sd[msk]
        order = np.argsort(i_meas)
        axes[0].plot(v_meas[order], i_meas[order], styles[t], ms=5, alpha=0.7,
                     label=f"data {t:g} C")
        v = np.linspace(v_meas.min(), v_meas.max(), 200)
        axes[0].plot(v, vdmos_diode(v, p, t_c=t), "--", alpha=0.9,
                     label=f"fit {t:g} C")
        axes[1].plot(i_meas[order], rel[order] * 100, styles[t], ms=5,
                     alpha=0.8, label=f"{t:g} C")
    axes[0].set_xlabel("v_sd [V]")
    axes[0].set_ylabel("i_sd [A]")
    axes[0].set_title("Body diode: data vs fit A")
    axes[0].legend(fontsize=7)
    axes[0].grid(alpha=0.3)
    axes[1].axhline(0.0, color="k", lw=0.8)
    axes[1].set_xlabel("i_sd [A]")
    axes[1].set_ylabel("relative residual [%]")
    axes[1].set_title("Residuals (fit A)")
    axes[1].legend(fontsize=7)
    axes[1].grid(alpha=0.3)
    fig.savefig(path, dpi=110, bbox_inches="tight")
    print(f"wrote {path}")


def main() -> None:
    data = load_diode_csv(CSV)
    print(f"{data.label}: {len(data)} points, temperatures {data.temp_grid()} C, "
          f"i_sd {data.i_sd.min():.1f}-{data.i_sd.max():.1f} A")
    print(f"free = {FREE}  (k_n fixed {START.k_n:g}, k_xti fixed {START.k_xti:g}, i_sr off)")
    print(f"device bounds overrides = {DEVICE_BOUNDS}")

    probe: list[str] = []
    sens = diode_sensitivity_report(START, data, FREE)
    probe.append("pre-fit diode_sensitivity_report (mean |d ln i_sd / d ln p| at START):\n"
                 + "\n".join(f"  {n:<12s} {v:9.3g}" for n, v in sens.items()))
    print("\npre-fit sensitivity (start): "
          + "  ".join(f"{n}={v:.3g}" for n, v in sens.items()))

    print(f"\n--- fit A (OFFICIAL protocol): {', '.join(FREE)} free ---")
    res_a = multi_fit(data)
    report(res_a, data, "diode fit A, all 30 points")

    drift = max(abs(getattr(res_a.params, n) - getattr(ST8L65_DIODE, n))
                / max(abs(getattr(ST8L65_DIODE, n)), 1e-12) for n in FREE)
    print(f"official-record drift vs this run's fit A: {drift:.2e}"
          + ("  [OK]" if drift < OFFICIAL_TOL else
             f"  [STALE: > {OFFICIAL_TOL:g} -- re-promote fit A]"))
    print(f"record RMSE {ST8L65_DIODE_RMSE_REL * 100:.2f} % | this run "
          f"{res_a.rmse_rel * 100:.2f} %")

    for txt in (kn_sweep_probe(data), free_kn_probe(data),
                mechanism_probes(data, res_a.params), tc_off_probe(data, res_a.params),
                trs2_off_probe(data, res_a.params)):
        print(txt)
        probe.append(txt)

    probe.append(f"\n--- fit A ---\n{res_a.summary()}\n"
                 + identifiability(res_a, FREE).summary())

    OUT.mkdir(exist_ok=True)
    txt = OUT / "st8l65_diode_params.txt"
    txt.write_text(res_a.summary()
                   + f"\n\nRMSE_rel = {res_a.rmse_rel * 100:.2f} % on 30 points "
                     f"(-55/25/150 C, 5-55 A).\n"
                     "k_n fixed at 2.0 (physical diffusion limit; freeing it buys ~0.1 %, "
                     "lands at k_n=2.07 with corr(i_s,k_n)=+1.000 -- see probe).\n"
                     "i_s is only meaningful paired with k_n (sigma_rel ~ 1: no low-current "
                     "tail in the data). v_eg ~ 1.09 eV = silicon bandgap, sensed.\n"
                     "r_diode is the diode-path series resistance; the diode bypasses r_d/r_s.\n"
                     "r_diode <-> i_kf corr +0.95: quote the curve, not the split.\n"
                     "a_ikf/a_rdiode/b_rdiode: LTspice Tikf/Trs1/Trs2 TCs (Tnom-"
                     "referenced, LTspice-verified); baselines: TC-off 4.78 %, "
                     "linear-TC-only 4.15 %. The five-knob TC/knee block is "
                     "strongly correlated (|corr| up to 0.985) -- quote the "
                     "curve, never the split.\n")
    print(f"wrote {txt}")

    probe_file = OUT / "st8l65_diode_probe.txt"
    probe_file.write_text("\n\n".join(probe) + "\n")
    print(f"wrote {probe_file}")

    plot(data, res_a.params, OUT / "st8l65_diode_fit.png")


if __name__ == "__main__":
    main()
