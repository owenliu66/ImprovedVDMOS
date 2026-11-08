"""Extract the ISC040N10NM8 body-diode parameters (LTspice diode DC branch).

Data: data/ISC040N10NM8_Diode.csv -- 26 points, I_sd 1.2-490 A at 25/175 degC
(datasheet body-diode forward characteristics). GOTCHA: the temperature
column is `T_j (C)`, NOT `T` -> load_diode_csv() (expects `T`) FAILS here;
this script carries a local loader (same pattern as fit_isc040_channel.py's
spec-CSV loaders).

OFFICIAL RECORD (promoted this session): ``reference.py`` ``ISC040_DIODE``.

PROTOCOL (locked after the session-27 probe campaign):
  * Fit A (OFFICIAL): free = i_s, k_n, r_diode, i_kf, v_eg, a_rdiode;
    k_xti fixed at 3.0 (LTspice default), i_sr = 0 (recombination branch
    unsensed -- freeing it buys NOTHING, session-27 reproduction),
    k_nr = 2.0, t_nom = 25, a_ikf = 0, b_rdiode = 0. Four multi-starts,
    all land in the same basin. RMSE 0.514 % (25 C 0.445 %, 175 C 0.575 %).
    Jacobian rank 5/6: the i_s/k_n pair is a degenerate valley (corr +0.999).
  * KEY FINDING: k_n = 2.0 (the ST8L65 fixed value) is WRONG here. With
    k_n fixed at 2.0 the fit lands ~35 %; freeing it drops the TC-off fit to
    3.31 % and k_n lands ~1.05 (near-ideal diffusion, unlike the SJ part).
    i_s <-> k_n corr +0.999 (degenerate prefactor pair -- quote the PAIR,
    never i_s alone; i_s sigma_rel ~12).
  * KEY FINDING: the LINEAR r_diode TC (a_rdiode, LTspice Trs1) is the
    decisive temperature knob: 3.31 % -> 0.51 %, strongly identified
    (sigma_rel ~0.013). The quadratic b_rdiode and the knee TC a_ikf buy
    NOTHING with only 2 temperatures and go degenerate (a_rdiode <-> b_rdiode
    corr -1.0) -> DO NOT free them (ST8L65 had 3 temperatures; that is what
    made Trs2 identifiable there).
  * i_kf (high-injection knee, ~357 A) and v_eg (~1.11 eV, silicon) are
    genuinely sensed. r_diode <-> i_kf corr +0.994: quote the curve.
  * r_diode ~ 0.96 mOhm is ~10x SMALLER than the recovery bench's old
    R1 = 10 mOhm; the full subckt reconciles them (R1 = r_diode(T) with the
    Trs1 law; recovery is R1-insensitive -- measured, out/recov_r1_sens.cir).

    python scripts/fit_isc040_diode.py
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from vdmos_fit import (  # noqa: E402
    Params,
    diode_sensitivity_report,
    fit_diode,
    identifiability,
    vdmos_diode,
)
from vdmos_fit.data import DiodeCurveSet  # noqa: E402

OUT = ROOT / "out"
CSV = ROOT / "data" / "ISC040N10NM8_Diode.csv"

FREE = ("i_s", "k_n", "r_diode", "i_kf", "v_eg", "a_rdiode")

BASE = Params(i_s=1e-12, k_n=1.1, i_sr=0.0, k_nr=2.0, i_kf=100.0, a_ikf=0.0,
              r_diode=1e-3, a_rdiode=0.0, b_rdiode=0.0, t_nom=25.0,
              k_xti=3.0, v_eg=1.11)

# Per-device bound: the diode-path series R is ~1 mOhm here (the fit's own
# answer); 20 mOhm is the ST8L65 precedent ceiling that keeps the solver out
# of the runaway basin (r_diode ~ 0.3 Ohm with similar RMSE, not a diode).
DEVICE_BOUNDS = {"r_diode": (0.0, 0.020)}

STARTS = [
    dict(i_s=1e-12, k_n=1.2, r_diode=1e-3, i_kf=100, v_eg=1.11, a_rdiode=0.0),
    dict(i_s=1e-10, k_n=1.05, r_diode=5e-4, i_kf=1000, v_eg=1.0, a_rdiode=0.002),
    dict(i_s=1e-14, k_n=1.01, r_diode=2e-3, i_kf=50, v_eg=1.2, a_rdiode=0.001),
    dict(i_s=3.7e-12, k_n=1.05, r_diode=0.96e-3, i_kf=350, v_eg=1.11,
         a_rdiode=0.001),
]


def load_diode_csv_isc040(path: Path) -> DiodeCurveSet:
    """ISC040 diode CSV, header `index,V_sd,I_sd,T_j (C)` (local loader)."""
    rows = list(csv.DictReader(path.open(newline="", encoding="utf-8")))
    return DiodeCurveSet(
        np.array([float(r["I_sd"]) for r in rows]),
        np.array([float(r["V_sd"]) for r in rows]),
        np.array([float(r["T_j (C)"]) for r in rows]),
        label=path.stem,
    )


def multi_fit(data: DiodeCurveSet, free=FREE, base=BASE) -> object:
    best = None
    for s in STARTS:
        # Only override parameters that are actually free: a probe that FIXES
        # a parameter (leaves it out of `free`) must keep the base's value --
        # a start dict must never silently re-introduce it.
        r = fit_diode(data, start=base.replace(**{k: v for k, v in s.items()
                                                  if k in free}),
                      free=free, weight="relative", max_nfev=8000,
                      bounds_overrides=DEVICE_BOUNDS)
        print(f"   start i_s={s['i_s']:8.1e} k_n={s['k_n']:5.3g}"
              f" r_diode={s['r_diode']*1e3:6.3g}m i_kf={s['i_kf']:6.3g}"
              f" a=({s['a_rdiode']:+.4f}) -> RMSE {r.rmse_rel*100:6.3f}%"
              f" {'conv' if r.success else 'stop'}")
        if best is None or r.rmse_rel < best.rmse_rel:
            best = r
    return best


def per_t_rel(data: DiodeCurveSet, p: Params) -> tuple[np.ndarray, dict]:
    model = vdmos_diode(data.v_sd, p, t_c=data.temp_c)
    rel = (model - data.i_sd) / data.i_sd
    per = {int(t): float(np.sqrt(np.mean(rel[np.isclose(data.temp_c, t)] ** 2)))
           for t in data.temp_grid()}
    return rel, per


def probe(data: DiodeCurveSet, label: str, free, base: Params) -> str:
    r = multi_fit(data, free=free, base=base)
    _, per = per_t_rel(data, r.params)
    return (f"  {label}: RMSE {r.rmse_rel*100:6.2f} %"
            f" | per-T rms: {per}"
            f" | k_n={r.params.k_n:.4g} i_s={r.params.i_s:.3g}"
            f" r_diode={r.params.r_diode*1e3:.3g}m"
            f" a_rdiode={r.params.a_rdiode:+.5f}")


def plot(data: DiodeCurveSet, p: Params, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4))
    for t, mk in ((25.0, "o"), (175.0, "s")):
        msk = np.isclose(data.temp_c, t)
        axes[0].plot(data.v_sd[msk], data.i_sd[msk], mk, ms=5, alpha=0.7,
                     label=f"data {t:g} C")
        v = np.linspace(data.v_sd[msk].min(), data.v_sd[msk].max(), 300)
        axes[0].plot(v, vdmos_diode(v, p, t_c=t), "--", alpha=0.85)
    axes[0].set_xlabel("v_sd [V]")
    axes[0].set_ylabel("i_sd [A]")
    axes[0].set_title("ISC040 body diode: data vs fit (ISC040_DIODE)")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3)
    rel, _ = per_t_rel(data, p)
    for t, mk in ((25.0, "o"), (175.0, "s")):
        msk = np.isclose(data.temp_c, t)
        axes[1].plot(data.v_sd[msk], rel[msk] * 100, mk, ms=5,
                     label=f"{t:g} C")
    axes[1].axhline(0.0, color="k", lw=0.8)
    axes[1].set_xlabel("v_sd [V]")
    axes[1].set_ylabel("relative error [%]")
    axes[1].set_title("Residuals (overall RMSE 0.51 %)")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3)
    fig.savefig(path, dpi=110, bbox_inches="tight")
    plt.close(fig)
    print("wrote", path)


def main() -> None:
    OUT.mkdir(exist_ok=True)
    data = load_diode_csv_isc040(CSV)
    print(f"ISC040 diode set: {len(data)} points, T = {data.temp_grid()}")

    print("\npre-fit sensitivity (mean |d ln i_sd / d ln p| at BASE):")
    for n, v in diode_sensitivity_report(BASE, data, FREE).items():
        print(f"  {n:9s} {v:8.3g}")

    print("\n=== FIT A (OFFICIAL): free = i_s, k_n, r_diode, i_kf, v_eg, "
          "a_rdiode ===")
    res = multi_fit(data)
    p = res.params
    print(res.summary())
    print("on_bounds:", res.on_bounds())
    print(identifiability(res, FREE).summary())
    rel, per = per_t_rel(data, p)
    print(f"overall RMSE {res.rmse_rel*100:.3f} % | per-T: {per}")

    print("\n=== PROTOCOL PROBES (why this free set) ===")
    lines = []
    # (1) k_n pinned at the ST8L65 value 2.0: wrong for this part.
    lines.append(probe(data, "k_n FIXED 2.0 (ST8L65 value), TC on        ",
                       ("i_s", "r_diode", "v_eg", "i_kf", "a_rdiode"),
                       BASE.replace(k_n=2.0)))
    # (2) TC off: k_n free but a_rdiode pinned 0 -> the 3.31 % baseline.
    lines.append(probe(data, "TC off (a_rdiode pinned 0), k_n free      ",
                       ("i_s", "k_n", "r_diode", "i_kf", "v_eg"), BASE))
    # (3) free b_rdiode + a_ikf on top of fit A: no gain, degenerate.
    lines.append(probe(data, "fit A + b_rdiode + a_ikf free             ",
                       FREE + ("b_rdiode", "a_ikf"), p))
    # (4) recombination branch free: worse (4.5 %).
    lines.append(probe(data, "fit A + i_sr free                         ",
                       FREE + ("i_sr",), p.replace(i_sr=1e-8)))
    print("\n".join(lines))

    rec = {
        "part": "ISC040N10NM8",
        "params": {n: getattr(p, n) for n in FREE},
        "fixed": {"k_xti": p.k_xti, "i_sr": p.i_sr, "k_nr": p.k_nr,
                  "t_nom": p.t_nom, "a_ikf": p.a_ikf, "b_rdiode": p.b_rdiode},
        "free": list(FREE),
        "rmse_rel": res.rmse_rel,
        "rmse_per_t": per,
        "probes": lines,
    }
    (OUT / "isc040_diode_params.json").write_text(json.dumps(rec, indent=2),
                                                  encoding="utf-8")
    print("\nwrote", OUT / "isc040_diode_params.json")
    plot(data, p, OUT / "isc040_diode_fit.png")


if __name__ == "__main__":
    main()
