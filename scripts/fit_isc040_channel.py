"""Fit the ISC040N10NM8 (100 V trench super-junction, non-Kelvin) channel.

OFFICIAL RECORD (promoted this session): ``reference.py`` ``ISC040_CHANNEL``.
This closes the HANDOVER queue item that was BLOCKED ON DATA (the ISC040
channel record -- data/ now carries ISC040N10NM8_IV_25C.csv + the Rdson/Vth
spec CSVs).

User-directed protocol (three notes from the handover):
  1. NO transfer curve (datasheet V_ds conditions for the I_d-V_gs figure are
     vague) -> the channel is fit from OUTPUT curves + Rds(on) chords only.
     Consequence: the gate-loop degeneration (r_s) is NOT identifiable and the
     data REJECTS it (r_s -> 0 from every start, monotonically worse as r_s
     grows -- the same finding as the IMCQ Kelvin study, session 22). The DC
     record carries the whole current-loop series R in r_d; r_s = r_sk = 0.
  2. NOT a Kelvin package -> r_sk = 0 (no kelvin tap; the gate loop returns at
     the power source). The 4-terminal subckt keeps the K pin for the driver
     return; the kelvin split is inert at r_s = r_sk = 0.
  3. Fit 25 C FIRST and LOCK, then fit temperature coefficients. There is NO
     high-temperature I-V data -- the TC comes from Rds(on)-vs-Vgs at
     -55/25/100/175 C (I_test = 25 A) and Vth-vs-Tj (520 uA / 52 uA).

Model family: the SJ charge-cap form (v_ovmax / k_ovsharp free, k_theta = 0,
v_vsat = 1e4 off, k_lambda = 0, g_k = 0) -- the same family as the ST8L65N044M9
super-junction record. It beats the velocity-saturation form (1.33 % vs 1.73 %
25 C RMSE) and is better conditioned. The cap is genuinely sensed: it ceilings
the V_gs = 15 channel conductance (without it the Ron chord is ~13 % low).

Stage 1 (25 C, locked): free = k_p, v_th0, v_sub, v_ovmax, k_ovsharp,
  k_kneesharp, k_eta, v_dbsat, r_d.  Data = IV_25C (69) + Rds(on)-vs-Id_25C
  (55, converted to V_ds = R*I chords) + Rds(on)-vs-Vgs 25 C subset (8, at
  I_test = 25 A).  Relative weighting.  r_d <= datasheet 25 C chord (3.4 mO).
  Fixed: k_theta = 0, v_vsat = 1e4, k_lambda = 0, g_k = 0, r_s = r_sk = 0.

Stage 2 (TC, 25 C locked): free = s_vth, beta_p, a_rd, b_rd.  TC laws:
  k_p(T)   = k_p25 * (T_K/298.15)^(-beta_p)
  v_th0(T) = v_th0_25 + s_vth*(T - 25)
  r_d(T)   = r_d25*(1 + a_rd*(T-25) + b_rd*(T-25)^2)   [quadratic, Trs1/Trs2
             law form -- NEVER linear (goes negative below ~-150 C); vertex
             ~-130 C, valid over the -55..175 C data span]
  v_sub, v_ovmax, k_ovsharp, k_kneesharp, k_eta, v_dbsat: T-independent.
  Data = Rds(on)-vs-Vgs at all 4 T (chord residual as CURRENT at R_data*25 A)
  + Vth-vs-Tj at both currents (residual on the Vth VOLTAGE: solve the model
  gate voltage at V_ds = V_gs that gives I_d, vs the datasheet Vth -- this
  fits the Vth(T) curve and its slope, avoiding the exponential-current
  sensitivity that the single-softplus subthreshold cannot reproduce).

Identifiability doctrine (quote the CURVES + TC laws, never the split):
  * r_d <-> r_s corr -1.0 (only the sum sensed); gate degeneration rejected.
  * k_p <-> v_th0 <-> v_ovmax <-> r_d cluster (|corr| up to 0.86): the composite
    channel conductance k_p*v_ov_eff2(v_gs) is the record.
  * beta_p weakly determined (sigma ~ 4): the k_p(T)/r_d(T)/v_th0(T) TC split is
    model-form dependent; beta_p is interpretable ONLY because Rds(on)-Vgs(T)
    pins r_d(T) at high V_gs (the r_d-dominated floor) leaving the low-V_gs
    channel shape for beta_p -- the IMCQ doctrine.

    python scripts/fit_isc040_channel.py
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

from vdmos_fit import Params, fit, identifiability, load_csv, vdmos_id  # noqa: E402
from vdmos_fit.data import CurveSet  # noqa: E402

OUT = ROOT / "out"
T0K = 298.15
I_TEST = 25.0  # Rds(on)-vs-Vgs test current [A] (verified: Vgs=15.7 -> 3.33 mO
               # matches IV_25C row Vgs=15,Vds=0.08325,I=25)

# ---------------------------------------------------------------------------
# data loading
# ---------------------------------------------------------------------------


def load_rdson_id(path: Path) -> CurveSet:
    """Rds(on)-vs-Id CSV `index,I_d,R_dson (mOhm),V_gs` -> ohmic IV points."""
    rows = list(csv.DictReader(path.open(newline="", encoding="utf-8")))
    vg, vd, idv = [], [], []
    for r in rows:
        i = float(r["I_d"])
        R = float(r["R_dson (mOhm)"]) * 1e-3
        vg.append(float(r["V_gs"]))
        vd.append(R * i)
        idv.append(i)
    return CurveSet(np.array(vg), np.array(vd), np.array(idv), "rdson_id")


def load_rdson_vgs(path: Path, Tsel: float | None = None) -> list[tuple]:
    """Rds(on)-vs-Vgs CSV -> (T_c, V_gs, chord_ohm) rows (I_test = 25 A)."""
    rows = list(csv.DictReader(path.open(newline="", encoding="utf-8")))
    out = []
    for r in rows:
        T = float(r["T_j (C)"])
        if Tsel is not None and abs(T - Tsel) > 1e-6:
            continue
        out.append((T, float(r["V_gs"]), float(r["R_dson (mOhm)"]) * 1e-3))
    return out


def load_vth(path: Path) -> list[tuple]:
    """Vth-vs-Tj CSV -> (T_c, V_th, I_d) rows."""
    rows = list(csv.DictReader(path.open(newline="", encoding="utf-8")))
    return [(float(r["T_j (C)"]), float(r["V_th"]), float(r["I_d"])) for r in rows]


def combine(*sets: CurveSet) -> CurveSet:
    return CurveSet(
        np.concatenate([s.v_gs for s in sets]),
        np.concatenate([s.v_ds for s in sets]),
        np.concatenate([s.i_d for s in sets]),
        "combined25",
    )


# ---------------------------------------------------------------------------
# stage 1: 25 C channel fit
# ---------------------------------------------------------------------------

FREE25 = ("k_p", "v_th0", "v_sub", "v_ovmax", "k_ovsharp", "k_kneesharp",
          "k_eta", "v_dbsat", "r_d")

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


def fit_25c(D: CurveSet) -> object:
    best = None
    for s in MULTI25:
        r = fit(D, start=BASE25.replace(**s), free=FREE25, weight="relative",
                max_nfev=8000, bounds_overrides={"r_d": (0.0, 3.4e-3)})
        print(f"   start k_p={s['k_p']:5g} v_ovmax={s['v_ovmax']:5g}"
              f" k_ovsh={s['k_ovsharp']:5g} -> RMSE {r.rmse_rel*100:5.2f}%"
              f" {'conv' if r.success else 'stop'}")
        if best is None or r.rmse_rel < best.rmse_rel:
            best = r
    return best


# ---------------------------------------------------------------------------
# stage 2: temperature coefficients (25 C locked)
# ---------------------------------------------------------------------------

TCLAW_FREE = ("s_vth", "beta_p", "a_rd", "b_rd")
TC_START = np.array([-6.5e-3, 1.5, 6.0e-3, 2.5e-5])
TC_LO = np.array([-0.03, 0.0, -5e-3, -1e-4])
TC_HI = np.array([0.0, 8.0, 5e-2, 1e-3])


def params_at(p25: Params, tc: float, s_vth: float, beta_p: float,
              a_rd: float, b_rd: float) -> Params:
    """Full Params at temperature tc [degC] through the record TC laws."""
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
    jac = sol.jac
    col = np.linalg.norm(jac, axis=0)
    tol = col.max() * 1e-5
    cov = np.linalg.pinv(jac.T @ jac)
    sd = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    sigma = {n: (float("inf") if c < tol else float(v))
             for n, c, v in zip(TCLAW_FREE, col, sd, strict=True)}
    return {"x": sol.x, "s_vth": s, "beta_p": b, "a_rd": a, "b_rd": bb,
            "rmse_rdson": float(np.sqrt(np.mean(rel[:n1] ** 2))),
            "rmse_vth": float(np.sqrt(np.mean(rel[n1:] ** 2))),
            "rmse_overall": float(np.sqrt(np.mean(rel ** 2))),
            "sigma_rel": sigma}


# ---------------------------------------------------------------------------
# acceptance checks (independent specs, never fitted)
# ---------------------------------------------------------------------------


def acceptance(p25, tcres, D25, iv25) -> dict:
    s, b, a, bb = tcres["x"]
    out = {}
    # Ron chords @25 A (Vgs 15/10/8) at 25 C
    out["ron25"] = {vg: model_ron(p25, 25.0, vg, s, b, a, bb) * 1e3
                    for vg in (15.0, 10.0, 8.0)}
    # Ron @Vgs=15 across T
    out["ron_vgs15_T"] = {int(t): model_ron(p25, t, 15.7, s, b, a, bb) * 1e3
                          for t in (-55.0, 25.0, 100.0, 175.0)}
    # Vth(T) at 520 uA
    out["vth520_T"] = {int(t): model_vth(p25, t, 5.2e-4, s, b, a, bb)
                       for t in (-52.4, 22.5, 107.5, 172.9)}
    # g_fs-like: dI/dVgs at Vgs=10, Vds=2 (a transfer slope we DO have via
    # the Rdson_Vgs shape) -- reported for reference only
    p = p25
    h = 0.01
    out["gm_vgs10_vds2"] = (float(vdmos_id(2.0, 10.0 + h, p))
                            - float(vdmos_id(2.0, 10.0 - h, p))) / (2 * h)
    return out


# ---------------------------------------------------------------------------
# plots
# ---------------------------------------------------------------------------


def plot_iv(d: CurveSet, p: Params, title: str, fname: str) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4))
    for vg in d.v_gs_grid():
        msk = np.isclose(d.v_gs, vg)
        if msk.sum() == 1:
            axes[0].plot(d.v_ds[msk], d.i_d[msk], "k*", ms=10)
            continue
        vd_data, id_data = d.v_ds[msk], d.i_d[msk]
        vd = np.linspace(0.0, vd_data.max(), 300)
        for ax in axes:
            ax.plot(vd_data, id_data, "o", ms=4, alpha=0.6,
                    label=f"v_gs={vg:g}" if ax is axes[0] else None)
            ax.plot(vd, vdmos_id(vd, vg, p), "--", alpha=0.85)
    axes[0].set_xlim(0, 5.2)
    axes[1].set_xlim(0, 1.0)
    axes[1].set_ylim(0, 130)
    for ax, t in zip(axes, [title, title + " (knee/ohmic zoom)"], strict=True):
        ax.set_xlabel("v_ds [V]")
        ax.set_ylabel("i_d [A]")
        ax.set_title(t)
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=7, ncol=2, title="dots = data, dashed = fit")
    fig.savefig(OUT / fname, dpi=110, bbox_inches="tight")
    plt.close(fig)
    print("wrote", OUT / fname)


def plot_tc(p25, tcres, rv, vt) -> None:
    s, b, a, bb = tcres["x"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4))
    # Rds(on) vs Vgs at each T
    temps = sorted({t for t, _, _ in rv})
    for t in temps:
        pts = [(vg, R) for tt, vg, R in rv if abs(tt - t) < 1e-6]
        pts.sort()
        axes[0].plot([p[0] for p in pts], [p[1] * 1e3 for p in pts], "o", ms=4,
                     label=f"data {t:g} C")
        vgs = np.linspace(min(p[0] for p in pts), max(p[0] for p in pts), 60)
        mod = [model_ron(p25, t, vg, s, b, a, bb) * 1e3 for vg in vgs]
        axes[0].plot(vgs, mod, "--", alpha=0.85)
    axes[0].set_xlabel("V_gs [V]")
    axes[0].set_ylabel("Rds(on) chord @25 A [m$\\Omega$]")
    axes[0].set_title("Rds(on) vs V_gs: data vs TC model")
    axes[0].legend(fontsize=7)
    axes[0].grid(alpha=0.3)
    # Vth vs T at both currents
    for cur, mk in ((5.2e-4, "o"), (5.2e-5, "s")):
        pts = [(t, vth) for t, vth, it in vt if abs(it - cur) < 1e-9]
        pts.sort()
        axes[1].plot([p[0] for p in pts], [p[1] for p in pts], mk, ms=5,
                     label=f"data {cur*1e6:g} uA")
        Ts = np.linspace(-55, 175, 60)
        mod = [model_vth(p25, t, cur, s, b, a, bb) for t in Ts]
        axes[1].plot(Ts, mod, "--", alpha=0.85)
    axes[1].set_xlabel("T_j [degC]")
    axes[1].set_ylabel("V_th [V]  (V_ds = V_gs)")
    axes[1].set_title("Threshold vs T_j: data vs TC model")
    axes[1].legend(fontsize=7)
    axes[1].grid(alpha=0.3)
    fig.savefig(OUT / "isc040_channel_tc.png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    print("wrote", OUT / "isc040_channel_tc.png")


# ---------------------------------------------------------------------------


def main() -> None:
    OUT.mkdir(exist_ok=True)
    iv = load_csv(ROOT / "data" / "ISC040N10NM8_IV_25C.csv")
    rid = load_rdson_id(ROOT / "data" / "ISC040N10NM8_Rdson_Id_25C.csv")
    rv25 = load_rdson_vgs(ROOT / "data" / "ISC040N10NM8_Rdson_Vgs.csv", Tsel=25.0)
    rv_all = load_rdson_vgs(ROOT / "data" / "ISC040N10NM8_Rdson_Vgs.csv")
    vt = load_vth(ROOT / "data" / "ISC040N10NM8_Vth_Tj.csv")
    rv25_set = CurveSet(np.array([r[1] for r in rv25]),
                        np.array([r[2] * I_TEST for r in rv25]),
                        np.array([I_TEST] * len(rv25)), "rv25")
    D25 = combine(iv, rid, rv25_set)
    print(f"25 C set: IV {len(iv)} + RdsonId {len(rid)} + RdsonVgs25 {len(rv25)}"
          f" = {len(D25)} points")

    print("\n=== STAGE 1: 25 C channel (locked) ===")
    r25 = fit_25c(D25)
    p25 = r25.params
    print(r25.summary())
    print("on_bounds:", r25.on_bounds())
    print(identifiability(r25, FREE25).summary())

    print("\n=== STAGE 2: temperature coefficients (25 C locked) ===")
    tcres = fit_tc(p25, rv_all, vt)
    for n, v in zip(TCLAW_FREE, tcres["x"], strict=True):
        print(f"  {n:8s} = {v:.6g}   sigma_rel {tcres['sigma_rel'][n]:.3g}")
    print(f"  RMSE rdson {tcres['rmse_rdson']*100:.2f}%  vth {tcres['rmse_vth']*100:.2f}%"
          f"  overall {tcres['rmse_overall']*100:.2f}%")

    print("\n=== ACCEPTANCE (independent specs) ===")
    acc = acceptance(p25, tcres, D25, iv)
    print("  Ron@25A 25C [mO]:", {k: round(v, 3) for k, v in acc["ron25"].items()},
          " (data 15:3.33 10:3.52 8:3.87)")
    print("  Ron@Vgs15 vs T [mO]:", {k: round(v, 3) for k, v in acc["ron_vgs15_T"].items()},
          " (data -55:2.05 25:3.33 100:5.23 175:8.00)")
    print("  Vth@520uA vs T [V]:", {k: round(v, 3) for k, v in acc["vth520_T"].items()},
          " (data -52:3.37 22:3.05 107:2.58 173:2.07)")

    # record
    rec = {
        "part": "ISC040N10NM8",
        "p25": p25.as_dict(),
        "tc": {"s_vth": tcres["s_vth"], "beta_p": tcres["beta_p"],
               "a_rd": tcres["a_rd"], "b_rd": tcres["b_rd"]},
        "free25": list(FREE25), "tc_free": list(TCLAW_FREE),
        "rmse25": r25.rmse_rel, "rmse_tc": tcres["rmse_overall"],
        "rmse_rdson": tcres["rmse_rdson"], "rmse_vth": tcres["rmse_vth"],
    }
    (OUT / "isc040_channel_params.json").write_text(json.dumps(rec, indent=2),
                                                    encoding="utf-8")
    print("\nwrote", OUT / "isc040_channel_params.json")

    plot_iv(iv, p25, "ISC040N10NM8 output curves, 25 C", "isc040_iv_25C_fit.png")
    plot_tc(p25, tcres, rv_all, vt)


if __name__ == "__main__":
    main()
