"""EXPERIMENT (not a promotion): DIBL soft cap on the ISC040N10NM8 channel.

Problem. The ISC040 channel is fit from LOW-voltage data (output curves stop
near V_ds ~ 5 V, plus Rds(on) chords). The DIBL shift k_eta*V_ds*k_db grows
without bound in V_ds, so the fitted k_eta is an unconstrained linear
extrapolation: the promoted record already drives v_th negative above ~130 V
(v_th(100 V) = 1.16 V, v_th(150 V) = -0.41 V), and the two new off-state
points in ISC040N10NM8_IV_25C.csv (V_gs = 0, V_ds = 80 V, I_d = 1e-7 A @25 C
and 1e-5 A @125 C) are exponentially sensitive to v_th(80 V) -- refitting with
them pins k_eta even harder and makes the 100 V+ extrapolation worse.

Fix under test. Wrap the DIBL shift in the SAME tunable soft cap as the gate
overdrive charge cap:

    v_th = v_th0 - softcap(k_eta*V_ds*k_db ; v_diblmax, k_diblsharp)

v_diblmax ceilings the threshold SHIFT, so v_th >= v_th0 - v_diblmax for every
V_ds; v_diblmax = 0 is the uncapped model, bit for bit. New names: v_diblmax,
k_diblsharp.

This script refits the 25 C channel (the official stage-1 free set) twice --
cap OFF (v_diblmax = 0) and cap ON (v_diblmax, k_diblsharp free) -- on the SAME
data (25 C on-state IV + Rds(on) chords + the new 25 C off-state point), and
compares: fitted k_eta, the threshold-shift cap, the on-state RMSE (must not
degrade -- the cap is free), and the off-state / high-V_ds extrapolation.

    python scripts/experiment_dibl_cap.py
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
from scipy.optimize import brentq, least_squares  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from vdmos_fit import Params, fit, identifiability  # noqa: E402
from vdmos_fit.data import CurveSet  # noqa: E402
from vdmos_fit.model import _soft_cap, vdmos_id  # noqa: E402

OUT = ROOT / "out"
I_TEST = 25.0

# ---------------------------------------------------------------------------
# data loading (IV CSV carries a T_j column -> filter to 25 C here)
# ---------------------------------------------------------------------------


def load_iv_with_T(path: Path) -> tuple[CurveSet, CurveSet, tuple]:
    """Split the IV CSV into (25 C on-state, 25 C off-state V_gs=0, 125 C off)."""
    rows = list(csv.DictReader(path.open(newline="", encoding="utf-8")))
    on_vg, on_vd, on_id = [], [], []
    off_vg, off_vd, off_id = [], [], []
    off125 = None
    for r in rows:
        T = float(r["T_j"])
        vg, vd, idv = float(r["V_gs"]), float(r["V_ds"]), float(r["I_d"])
        if abs(T - 125.0) < 1e-6:
            off125 = (vg, vd, idv)
            continue
        if vg == 0.0:
            off_vg.append(vg)
            off_vd.append(vd)
            off_id.append(idv)
        else:
            on_vg.append(vg)
            on_vd.append(vd)
            on_id.append(idv)
    return (
        CurveSet(np.array(on_vg), np.array(on_vd), np.array(on_id), "iv25_on"),
        CurveSet(np.array(off_vg), np.array(off_vd), np.array(off_id), "off25"),
        off125,
    )


def load_rdson_id(path: Path) -> CurveSet:
    rows = list(csv.DictReader(path.open(newline="", encoding="utf-8")))
    vg, vd, idv = [], [], []
    for r in rows:
        i = float(r["I_d"])
        R = float(r["R_dson (mOhm)"]) * 1e-3
        vg.append(float(r["V_gs"]))
        vd.append(R * i)
        idv.append(i)
    return CurveSet(np.array(vg), np.array(vd), np.array(idv), "rdson_id")


def load_rdson_vgs25(path: Path) -> CurveSet:
    rows = list(csv.DictReader(path.open(newline="", encoding="utf-8")))
    vg, vd, idv = [], [], []
    for r in rows:
        if abs(float(r["T_j (C)"]) - 25.0) > 1e-6:
            continue
        R = float(r["R_dson (mOhm)"]) * 1e-3
        vg.append(float(r["V_gs"]))
        vd.append(R * I_TEST)
        idv.append(I_TEST)
    return CurveSet(np.array(vg), np.array(vd), np.array(idv), "rv25")


def combine(*sets: CurveSet) -> CurveSet:
    return CurveSet(
        np.concatenate([s.v_gs for s in sets]),
        np.concatenate([s.v_ds for s in sets]),
        np.concatenate([s.i_d for s in sets]),
        "combined25",
    )


def load_rdson_vgs_all(path: Path) -> list[tuple]:
    """Rds(on)-vs-Vgs CSV -> (T_c, V_gs, chord_ohm) rows at ALL temperatures."""
    rows = list(csv.DictReader(path.open(newline="", encoding="utf-8")))
    return [(float(r["T_j (C)"]), float(r["V_gs"]),
             float(r["R_dson (mOhm)"]) * 1e-3) for r in rows]


def load_vth(path: Path) -> list[tuple]:
    """Vth-vs-Tj CSV -> (T_c, V_th, I_d) rows."""
    rows = list(csv.DictReader(path.open(newline="", encoding="utf-8")))
    return [(float(r["T_j (C)"]), float(r["V_th"]), float(r["I_d"])) for r in rows]


# ---------------------------------------------------------------------------
# fits
# ---------------------------------------------------------------------------

FREE25 = ("k_p", "v_th0", "v_sub", "v_ovmax", "k_ovsharp", "k_kneesharp",
          "k_eta", "v_dbsat", "r_d")
FREE_CAP = FREE25 + ("v_diblmax", "k_diblsharp")

BASE25 = Params(k_p=222.0, v_th0=4.29, v_sub=0.244, v_ovmax=5.8, k_ovsharp=2.2,
                k_theta=0.0, v_vsat=1e4, k_kneesharp=1.8, k_eta=0.031,
                k_lambda=0.0, g_k=0.0, v_dbsat=0.53, r_d=0.00248, r_s=0.0,
                r_sk=0.0)

MULTI25 = [
    dict(k_p=222, v_ovmax=5.8, k_ovsharp=2.2, k_kneesharp=1.8, v_dbsat=0.53),
    dict(k_p=400, v_ovmax=3.0, k_ovsharp=6.0, k_kneesharp=1.0, v_dbsat=1.0),
    dict(k_p=150, v_ovmax=8.0, k_ovsharp=1.5, k_kneesharp=2.5, v_dbsat=0.3),
    dict(k_p=600, v_ovmax=2.0, k_ovsharp=10.0, k_kneesharp=0.6, v_dbsat=2.0),
    dict(k_p=300, v_ovmax=4.5, k_ovsharp=3.5, k_kneesharp=1.5, v_dbsat=0.8),
]
# cap-start multipliers: spread v_diblmax / k_diblsharp so the cap basin is sampled
MULTI_CAP = [
    dict(v_diblmax=2.0, k_diblsharp=10.0),
    dict(v_diblmax=1.0, k_diblsharp=3.0),
    dict(v_diblmax=3.5, k_diblsharp=20.0),
    dict(v_diblmax=0.5, k_diblsharp=40.0),
    dict(v_diblmax=5.0, k_diblsharp=6.0),
]


def multi_fit(D: CurveSet, free, extra_starts):
    best = None
    for s in MULTI25:
        start = BASE25.replace(**s)
        starts = [start.replace(**e) for e in extra_starts] if extra_starts else [start]
        for st in starts:
            r = fit(D, start=st, free=free, weight="relative",
                    max_nfev=9000, bounds_overrides={"r_d": (0.0, 3.4e-3)})
            if best is None or r.rmse_rel < best.rmse_rel:
                best = r
    return best


# ---------------------------------------------------------------------------
# diagnostics
# ---------------------------------------------------------------------------


def vth_of(v_ds, p: Params) -> float:
    k_db = np.tanh(v_ds / p.v_dbsat)
    raw = p.k_eta * v_ds * k_db
    return p.v_th0 - _soft_cap(raw, p.v_diblmax, p.k_diblsharp)


def vth_crossing(p: Params) -> float | None:
    """V_ds where the (model) threshold first reaches 0, or None if never <= 200 V."""
    vs = np.linspace(0.0, 200.0, 40001)
    t = np.array([vth_of(v, p) for v in vs])
    neg = np.where(t <= 0.0)[0]
    return float(vs[neg[0]]) if neg.size else None


def on_rmse(p: Params, on: CurveSet) -> float:
    m = vdmos_id(on.v_ds, on.v_gs, p)
    rel = (m - on.i_d) / np.maximum(np.abs(on.i_d), 1e-3)
    return float(np.sqrt(np.mean(rel ** 2)))


def off_current(p: Params, v_ds: float) -> float:
    return float(vdmos_id(v_ds, 0.0, p))


# ---------------------------------------------------------------------------
# TC stage (official protocol, 25 C base LOCKED): free = s_vth, beta_p, a_rd,
# b_rd on Rds(on)-vs-Vgs at all 4 T (chord residual as current at R*25 A) +
# Vth-vs-Tj at both currents (Vth-VOLTAGE residual). v_sub/v_ovmax/k_ovsharp/
# k_kneesharp/k_eta/v_dbsat AND the cap knobs (v_diblmax/k_diblsharp) are
# T-independent, carried through p25.as_dict() unchanged.
# ---------------------------------------------------------------------------

T0K = 298.15
TC_START = np.array([-6.5e-3, 1.5, 6.0e-3, 2.5e-5])
TC_LO = np.array([-0.03, 0.0, -5e-3, -1e-4])
TC_HI = np.array([0.0, 8.0, 5e-2, 1e-3])


def params_at(p25: Params, tc: float, s_vth: float, beta_p: float,
              a_rd: float, b_rd: float) -> Params:
    dt = tc - 25.0
    tk = tc + 273.15
    fac = 1.0 + a_rd * dt + b_rd * dt * dt
    return Params(**{**p25.as_dict(),
                     "k_p": p25.k_p * (tk / T0K) ** (-beta_p),
                     "v_th0": p25.v_th0 + s_vth * dt,
                     "r_d": p25.r_d * max(fac, 1e-6)})


def model_ron(p25, tc, vg, s, b, a, bb) -> float:
    p = params_at(p25, tc, s, b, a, bb)
    v = brentq(lambda x: float(vdmos_id(x, vg, p)) - I_TEST, 1e-5, 20.0)
    return v / I_TEST


def model_vth(p25, tc, itest, s, b, a, bb) -> float:
    p = params_at(p25, tc, s, b, a, bb)
    return brentq(lambda vg: float(vdmos_id(vg, vg, p)) - itest, 0.5, 8.0)


def fit_tc(p25: Params, rv: list, vt: list) -> dict:
    def resid(x):
        s, b, a, bb = x
        r1 = [(model_ron(p25, t, vg, s, b, a, bb) - R) / R for t, vg, R in rv]
        r2 = [(model_vth(p25, t, it, s, b, a, bb) - vth) / vth
              for t, vth, it in vt]
        return np.array(r1 + r2)

    sol = least_squares(resid, TC_START, bounds=(TC_LO, TC_HI), method="trf",
                        x_scale="jac", max_nfev=6000, ftol=1e-13, xtol=1e-13,
                        gtol=1e-13)
    s, b, a, bb = sol.x
    rel = resid(sol.x)
    n1 = len(rv)
    return {"x": sol.x, "s_vth": s, "beta_p": b, "a_rd": a, "b_rd": bb,
            "rmse_rdson": float(np.sqrt(np.mean(rel[:n1] ** 2))),
            "rmse_vth": float(np.sqrt(np.mean(rel[n1:] ** 2))),
            "rmse_overall": float(np.sqrt(np.mean(rel ** 2)))}


def tc_report(tag: str, p25: Params, tcres: dict, off125) -> dict:
    s, b, a, bb = tcres["x"]
    p125 = params_at(p25, 125.0, s, b, a, bb)
    vg, vd, i125 = off125
    pred = float(vdmos_id(vd, vg, p125))
    ratio = pred / i125
    decades = np.log10(max(pred, 1e-30) / i125)
    print(f"\n=== TC + 125 C off-state: {tag} ===")
    print(f"  TC coeffs: s_vth {s:.5g}  beta_p {b:.4g}  a_rd {a:.5g}  b_rd {bb:.5g}")
    print(f"  TC RMSE: rdson {tcres['rmse_rdson']*100:.2f}%  "
          f"vth {tcres['rmse_vth']*100:.2f}%  overall {tcres['rmse_overall']*100:.2f}%")
    print(f"  v_th0(125C) {p125.v_th0:.4g}  k_p(125C) {p125.k_p:.4g}  "
          f"r_d(125C) {p125.r_d*1e3:.4g} mO")
    print(f"  v_th(V_ds=80, 125C) = {vth_of(80.0, p125):.4g}")
    print(f"  125 C off-state i_d(V_gs=0, V_ds=80): model {pred:.3g} A  "
          f"vs data {i125:.3g} A")
    print(f"    -> model/data = {ratio:.3g}  ({decades:+.2f} decades; "
          f"{'UNDER' if ratio < 1 else 'OVER'}-predicts by {abs(decades):.2f} decades)")
    return {"tag": tag, "s_vth": s, "beta_p": b, "a_rd": a, "b_rd": bb,
            "rmse_rdson": tcres["rmse_rdson"], "rmse_vth": tcres["rmse_vth"],
            "rmse_tc": tcres["rmse_overall"], "vth0_125": p125.v_th0,
            "vth_80_125": vth_of(80.0, p125),
            "off125_model": pred, "off125_data": i125,
            "ratio": ratio, "decades": float(decades)}


# ---------------------------------------------------------------------------
# plot: off-state leakage i_d(V_ds) at V_gs = 0, both fits
# ---------------------------------------------------------------------------


def plot_leak(p_off: Params, p_cap: Params, off_pts, fname: str) -> None:
    v_ds = np.linspace(0.0, 150.0, 601)
    i_off = np.maximum(np.asarray(vdmos_id(v_ds, 0.0, p_off)), 1e-12)
    i_cap = np.maximum(np.asarray(vdmos_id(v_ds, 0.0, p_cap)), 1e-12)
    fig, ax = plt.subplots(figsize=(8.5, 5.6))
    ax.semilogy(v_ds, i_off, label="DIBL cap OFF (uncapped)")
    ax.semilogy(v_ds, i_cap, label="DIBL cap ON (v_diblmax, k_diblsharp)")
    for _vg, vd, idv in off_pts:
        ax.plot(vd, idv, "k*", ms=14, label=f"data V_gs=0, {vd:g} V, {idv:g} A")
    ax.axhline(1e-7, ls=":", c="gray", lw=0.8)
    ax.set_xlabel("V_ds [V]")
    ax.set_ylabel("off-state i_d at V_gs = 0 [A]")
    ax.set_title("ISC040N10NM8 off-state leakage: uncapped vs soft-capped DIBL")
    ax.set_xlim(0, 150)
    ax.set_ylim(1e-12, 1e3)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)
    fig.savefig(OUT / fname, dpi=110, bbox_inches="tight")
    plt.close(fig)
    print("wrote", OUT / fname)


# ---------------------------------------------------------------------------


def report(tag: str, res, D: CurveSet, on25: CurveSet) -> dict:
    p = res.params
    print(f"\n=== {tag} ===")
    print(f"  converged {res.success}   overall RMSE {res.rmse_rel*100:.2f}%   "
          f"on-state-only RMSE {on_rmse(p, on25)*100:.2f}%")
    print(f"  k_p {p.k_p:.4g}  v_th0 {p.v_th0:.4g}  v_sub {p.v_sub:.4g}  "
          f"v_ovmax {p.v_ovmax:.4g}  k_ovsharp {p.k_ovsharp:.4g}")
    print(f"  k_kneesharp {p.k_kneesharp:.4g}  v_dbsat {p.v_dbsat:.4g}  "
          f"r_d {p.r_d*1e3:.4g} mO")
    print(f"  k_eta {p.k_eta:.5g}  v_diblmax {p.v_diblmax:.5g}  "
          f"k_diblsharp {p.k_diblsharp:.5g}")
    print("  on_bounds:", res.on_bounds())
    print("  v_th(V_ds):", {v: round(vth_of(v, p), 3)
                            for v in (5, 20, 40, 80, 100, 150)})
    cross = vth_crossing(p)
    print(f"  v_th reaches 0 at V_ds = "
          f"{'never (<=200 V)' if cross is None else f'{cross:.1f} V'}")
    print("  off-state i_d(V_gs=0):", {f"{v}V": f"{off_current(p, v):.3g}"
                                       for v in (80, 100, 150)},
          " (data @80V = 1e-7)")
    return {"tag": tag, "params": p.as_dict(), "rmse": res.rmse_rel,
            "on_rmse": on_rmse(p, on25), "k_eta": p.k_eta,
            "v_diblmax": p.v_diblmax, "k_diblsharp": p.k_diblsharp,
            "vth_100": vth_of(100.0, p), "vth_150": vth_of(150.0, p),
            "vth_zero_at": cross,
            "off80": off_current(p, 80.0), "off100": off_current(p, 100.0),
            "off150": off_current(p, 150.0)}


def main() -> None:
    OUT.mkdir(exist_ok=True)
    on25, off25, off125 = load_iv_with_T(ROOT / "data" / "ISC040N10NM8_IV_25C.csv")
    rid = load_rdson_id(ROOT / "data" / "ISC040N10NM8_Rdson_Id_25C.csv")
    rv25 = load_rdson_vgs25(ROOT / "data" / "ISC040N10NM8_Rdson_Vgs.csv")
    D25 = combine(on25, rid, rv25, off25)
    print(f"25 C fit set: IV-on {len(on25)} + RdsonId {len(rid)} + "
          f"RdsonVgs25 {len(rv25)} + off25 {len(off25)} = {len(D25)} points")
    print(f"125 C off-state point (excluded from the 25 C fit): {off125}")

    res_off = multi_fit(D25, FREE25, None)
    res_cap = multi_fit(D25, FREE_CAP, MULTI_CAP)

    r_off = report("CAP OFF (v_diblmax = 0)", res_off, D25, on25)
    r_cap = report("CAP ON (v_diblmax, k_diblsharp free)", res_cap, D25, on25)

    print("\n=== COMPARISON ===")
    print(f"  on-state RMSE:   cap OFF {r_off['on_rmse']*100:.2f}%   "
          f"cap ON {r_cap['on_rmse']*100:.2f}%   (cap must not hurt the fit)")
    print(f"  k_eta:           cap OFF {r_off['k_eta']:.5g}   "
          f"cap ON {r_cap['k_eta']:.5g}")
    print(f"  v_th @100 V:     cap OFF {r_off['vth_100']:+.3f}   "
          f"cap ON {r_cap['vth_100']:+.3f}")
    print(f"  v_th @150 V:     cap OFF {r_off['vth_150']:+.3f}   "
          f"cap ON {r_cap['vth_150']:+.3f}")
    print(f"  v_th=0 at V_ds:  cap OFF {r_off['vth_zero_at']}   "
          f"cap ON {r_cap['vth_zero_at']}")
    print(f"  off i_d @100 V:  cap OFF {r_off['off100']:.3g}   "
          f"cap ON {r_cap['off100']:.3g}")
    print(f"  off i_d @150 V:  cap OFF {r_off['off150']:.3g}   "
          f"cap ON {r_cap['off150']:.3g}")

    print("\ncap-ON identifiability:")
    print(identifiability(res_cap, FREE_CAP).summary())

    # ---- TC stage on each locked 25 C base, then the 125 C off-state check ----
    rv_all = load_rdson_vgs_all(ROOT / "data" / "ISC040N10NM8_Rdson_Vgs.csv")
    vt = load_vth(ROOT / "data" / "ISC040N10NM8_Vth_Tj.csv")
    tc_off = fit_tc(res_off.params, rv_all, vt)
    tc_cap = fit_tc(res_cap.params, rv_all, vt)
    r_tc_off = tc_report("CAP OFF base", res_off.params, tc_off, off125)
    r_tc_cap = tc_report("CAP ON base", res_cap.params, tc_cap, off125)

    off_pts = [(0.0, 80.0, 1e-7)]
    plot_leak(res_off.params, res_cap.params, off_pts, "isc040_dibl_cap_leak.png")

    rec = {"part": "ISC040N10NM8", "experiment": "dibl_soft_cap",
           "n_points": len(D25), "cap_off": r_off, "cap_on": r_cap,
           "tc_cap_off": r_tc_off, "tc_cap_on": r_tc_cap}
    (OUT / "isc040_dibl_cap_experiment.json").write_text(
        json.dumps(rec, indent=2), encoding="utf-8")
    print("\nwrote", OUT / "isc040_dibl_cap_experiment.json")


if __name__ == "__main__":
    main()
