"""Head-to-head: ISC040 25 C channel fit, mobility knee form (k_kneeform=0)
vs the new linear knee-scale form (k_kneeform=1, v_knee = k_theta*v_ov_eff2).

Both use the same free set and data; k_theta is bounded per form
(-1/v_ovmax for mobility, 0 for scale). Reports RMSE, off-current at 80 V,
and the saturation gate exponent alpha = d ln i_d,sat / d ln v_ov_eff2.
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402

from vdmos_fit import Params, fit, vdmos_id  # noqa: E402
from vdmos_fit.data import CurveSet  # noqa: E402
from vdmos_fit.model import _soft_cap  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
I_TEST = 25.0

FREE = ("k_p", "v_th0", "v_sub", "v_ovmax", "k_ovsharp", "k_theta", "v_vsat",
        "k_kneesharp", "k_eta", "v_diblmax", "k_diblsharp", "v_dbsat", "r_d")

BASE = dict(k_p=280.0, v_th0=4.45, v_sub=0.242, v_ovmax=9.0, k_ovsharp=3.3,
            k_theta=1e-6, v_vsat=5.0, k_kneesharp=1.7, k_eta=40.0,
            v_diblmax=7.0, k_diblsharp=0.15, k_lambda=0.0, g_k=0.0,
            v_thk=1.6, v_kovsoft=0.1, v_dbsat=0.33, v_kneemin=1e-3,
            r_d=2.9e-3, r_s=0.0, r_sk=0.0)


def load_iv_T(path):
    rows = list(csv.DictReader(path.open(newline="", encoding="utf-8")))
    on = [[], [], []]
    off = [[], [], []]
    for r in rows:
        if abs(float(r["T_j"]) - 25.0) > 1e-6:
            continue
        vg, vd, idv = float(r["V_gs"]), float(r["V_ds"]), float(r["I_d"])
        tgt = off if vg == 0.0 else on
        tgt[0].append(vg)
        tgt[1].append(vd)
        tgt[2].append(idv)
    return (CurveSet(*(np.array(a) for a in on), label="on"),
            CurveSet(*(np.array(a) for a in off), label="off"))


def load_rid(path):
    rows = list(csv.DictReader(path.open(newline="", encoding="utf-8")))
    vg, vd, idv = [], [], []
    for r in rows:
        i = float(r["I_d"])
        R = float(r["R_dson (mOhm)"]) * 1e-3
        vg.append(float(r["V_gs"]))
        vd.append(R * i)
        idv.append(i)
    return CurveSet(np.array(vg), np.array(vd), np.array(idv), "rid")


def load_rv25(path):
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


def combine(*s):
    return CurveSet(np.concatenate([x.v_gs for x in s]),
                    np.concatenate([x.v_ds for x in s]),
                    np.concatenate([x.i_d for x in s]), "D25")


on25, off25 = load_iv_T(ROOT / "data" / "ISC040N10NM8_IV_25C.csv")
D = combine(on25, load_rid(ROOT / "data" / "ISC040N10NM8_Rdson_Id_25C.csv"),
            load_rv25(ROOT / "data" / "ISC040N10NM8_Rdson_Vgs.csv"), off25)
print(f"fit set: {len(D)} points\n")


def v_ov_eff2(p, vg, vds):
    k_db = np.tanh(vds / p.v_dbsat)
    v_th = p.v_th0 - _soft_cap(p.k_eta * vds * k_db, p.v_diblmax, p.k_diblsharp)
    z = np.clip((vg - v_th) / p.v_sub, -700, 700)
    v_ov = p.v_sub * np.log1p(np.exp(z))
    return _soft_cap(v_ov, p.v_ovmax, p.k_ovsharp)


def alpha(p, vg=10.0, vds=5.0):
    h = 0.02
    a = np.log(float(vdmos_id(vds, vg - h, p)))
    b = np.log(float(vdmos_id(vds, vg + h, p)))
    va = np.log(v_ov_eff2(p, vg - h, vds))
    vb = np.log(v_ov_eff2(p, vg + h, vds))
    return (b - a) / (vb - va)


def run(tag, kneeform, kt_start, kt_lo, kt_hi):
    best = None
    for kts in kt_start:
        st = Params().replace(**{**BASE, "k_kneeform": kneeform, "k_theta": kts})
        r = fit(D, start=st, free=FREE, weight="relative", max_nfev=15000,
                bounds_overrides={"r_d": (0.0, 3.4e-3), "k_theta": (kt_lo, kt_hi)})
        if best is None or r.rmse_rel < best.rmse_rel:
            best = r
    p = best.params
    print(f"=== {tag} ===  RMSE {best.rmse_rel*100:.3f}%  conv {best.success}")
    print(f"  k_theta {p.k_theta:.5g}  v_vsat {p.v_vsat:.4g}  v_ovmax {p.v_ovmax:.4g}"
          f"  k_kneesharp {p.k_kneesharp:.4g}  v_dbsat {p.v_dbsat:.4g}")
    print(f"  off i_d(V_gs=0,80V) = {float(vdmos_id(80.0, 0.0, p)):.3g} A   "
          f"(data 1e-7)")
    print(f"  saturation gate exponent alpha(v_gs=10,V_ds=5) = {alpha(p):.3f}")
    print(f"  on_bounds: {best.on_bounds()}")
    return best


# mobility form: k_theta free down to -1/v_ovmax (dynamic from start v_ovmax)
lo_mob = -1.0 / BASE["v_ovmax"]
run("MOBILITY form (k_kneeform=0, k_theta free to -1/v_ovmax)",
    0.0, [1e-6], lo_mob, 5.0)
# scale form: k_theta a non-negative multiplier
run("SCALE form (k_kneeform=1, v_knee = k_theta*v_ov_eff2)",
    1.0, [1.0], 0.0, 5.0)
