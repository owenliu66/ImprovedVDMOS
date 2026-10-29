"""JOINT multi-temperature fit for the IMCQ120R004M2H (CoolSiC 1200 V, 4 mOhm).

OFFICIAL RECORD (promoted session 21, user go-ahead): ``reference.py``
``IMCQ120_JOINT``, pinned by ``tests/test_imcq120_official.py``. Every run of
this script re-checks the record: sensed-parameter drift + record-RMSE
re-evaluation (OK/STALE -- session-19 lesson: drift in UNSENSED directions is
solver-trajectory noise, not staleness). The session-18
``out/IMCQ120_TEMPERATURE_REPORT.md`` was deleted (session 21); its durable
findings are carried in ``out/IMCQ120_JOINT_FIT_REPORT.md`` section 8.

HANDOVER PART 1: one residual vector over BOTH IV sets + the digitized transfer
characteristic (V_DS = 20 V, 25/175/200 C) + Rds(on)-Tj (185.2 A, Vgs = 15/18)
+ Vth-Tj (58.3 mA), with the TC LAWS as the only extra unknowns instead of
independent per-temperature fits.

User-directed protocol changes vs the HANDOVER spec:
  * r_d(T) is a QUADRATIC in T(degC) referenced to 25 C (LTspice Trs1/Trs2 law
    form), NOT the power law / exponential of the (deleted) temperature
    report. The digitized Rds(on)-Tj chords fit that quadratic to 0.07 % RMSE
    over the full -46..196 C span (power law: 9.5 %, exponential: 5.9 %) -- the
    quadratic curvature IS the cold-end flattening that killed the two-
    temperature power-law extrapolation (-15 % at -46 C).
  * k_theta is FREED (the transfer data senses the knee position: v_knee ->
    1/k_theta bends the high-v_gs transfer over from quadratic toward linear).
  * k_lambda re-freed (HANDOVER): the transfer points at V_DS = 20 V are
    exactly what senses it through k_db = tanh(v_ds/v_dbsat).

TC laws (extra unknowns beyond the 25 C Params):
  k_p(T)          = k_p25 * (T_K/298.15)^(-beta_p)      [power law, exponent free]
  v_th0(T)        = v_th0_25 + s_vth*(T - 25)           [linear, slope free]
  v_sub           = const                               [measured invariant]
  k_beta_trap(T)  = max(k_bt25 + s_bkt*(T - 25), 0)     [linear; Arrhenius variant]
  r_d(T)          = r_d25*(1 + a_rd*(T-25) + b_rd*(T-25)^2)   [QUADRATIC, user]
  r_d(T) variant  = r_d25*(T_K/298.15)^gamma_rd         [power law, comparison]
  k_theta, k_kneesharp, k_eta, v_dbsat, k_lambda: T-independent.
  Fixed: r_s = 0 (Kelvin source), g_k = 0 (no trench kink), v_ovmax = 0 (cap
  off -- not sensed on this SiC part, wrong sign), v_vsat = 1e4 (off).

Identifiability doctrine (HANDOVER 1.3): r_d(T) is pinned EXTERNALLY by the
Rds(on)-Tj data inside the joint residual vector (22 dedicated chord points),
NOT by the IV curves -- only then is the k_p(T) exponent beta_p interpretable
(the k_p/r_d TC split is otherwise model-form dependent, corr +0.82..+0.88).

Transfer-curve data note: the 175 C and 200 C curves are indistinguishable
below V_gs = 12 V (digitization limit); they cross at 12.025 V / 999.1 A and
that crossing point is carried in the CSV for BOTH curves. Residuals exist
only where data exists.

Acceptance checks (independent specs, never fitted -- doctrine: quote curves
and spec checks, never a single RMSE):
  g_fs = 122 S @ 185.2 A / 20 V; ZTC crossing ~14.2 V / 1500 A (25 vs 175 C);
  Ron chords 3.70 / 4.70 mOhm @25 C (Vgs 18/15); Vth spec 58.3 mA;
  Rds(on)-Tj over the FULL -46..195 C span; the high-v_gs corner at
  v_gs = 18 / V_DS = 20 must come down from the 3650 A gain-fit overshoot
  toward the digitized ~2070 A (175 C).

Variants (ablations, one basin each, started from the primary optimum):
  arrhenius  : k_beta_trap(T) Arrhenius instead of linear
  rd_powerlaw: r_d(T) power law instead of quadratic
  no_ktheta  : k_theta pinned at its 1e-4 floor (what freeing it buys)
  no_klambda : k_lambda pinned at 1e-12 (what re-freeing it buys)
  no_transfer: transfer + Rdson + Vth groups dropped (what the spec data buys)

    python scripts/fit_imcq120_joint.py

Long-running: launch detached and poll (120 s tool wall), e.g.
    subprocess.Popen([sys.executable, script], stdout=log, stderr=log,
                     creationflags=DETACHED_PROCESS)
"""

from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.optimize import brentq, least_squares  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from vdmos_fit import Params, load_csv, vdmos_id  # noqa: E402
from vdmos_fit.reference import (  # noqa: E402
    IMCQ120_JOINT,
    IMCQ120_JOINT_RMSE_REL,
)

OUT = ROOT / "out"
T0K = 298.15  # 25 C reference for all power laws
I_TEST = 185.2  # Rds(on) / g_fs test current [A]
I_VTH = 0.0583  # Vth spec current [A]

# --------------------------------------------------------------------------
# data loading (spec CSVs use the digitizer export style; script-local loaders
# per the established pattern -- data/README.md, IMCQ spec CSV entries)
# --------------------------------------------------------------------------


def load_transfer(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transfer CSV `index,V_gs,I_d,condition` -> (v_gs, T_c, i_d).

    condition = "V_ds=20V,T_j=25C". The 175/200 C curves share the crossing
    point (12.025 V, 999.1 A) -- it appears once per curve, keep both rows.
    """
    with path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    vg, tc, id_ = [], [], []
    for r in rows:
        cond = r["condition"]
        m_tj = re.search(r"T_j\s*=\s*(-?[\d.]+)\s*C", cond)
        m_vd = re.search(r"V_ds\s*=\s*(-?[\d.]+)\s*V", cond)
        if not (m_tj and m_vd):
            raise ValueError(f"unparsed condition {cond!r}")
        if float(m_vd.group(1)) != 20.0:
            raise ValueError(f"transfer fit assumes V_ds = 20 V, got {cond!r}")
        vg.append(float(r["V_gs"]))
        tc.append(float(m_tj.group(1)))
        id_.append(float(r["I_d"]))
    return np.array(vg), np.array(tc), np.array(id_)


def load_rdson(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Rds(on)-Tj CSV `index,Tj (C),R_dson (mOhm),V_gs` -> (T_c, V_gs, chord Ohm)."""
    with path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    t = np.array([float(r["Tj (C)"]) for r in rows])
    r_ohm = np.array([float(r["R_dson (mOhm)"]) for r in rows]) * 1e-3
    vg = np.array([float(r["V_gs"]) for r in rows])
    return t, vg, r_ohm


def load_vth(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Vth-Tj CSV `index,T_j (C),V_th,I_d` -> (T_c, V_th)."""
    with path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    t = np.array([float(r["T_j (C)"]) for r in rows])
    v = np.array([float(r["V_th"]) for r in rows])
    ids = {float(r["I_d"]) for r in rows}
    if ids != {I_VTH}:
        raise ValueError(f"Vth CSV current {ids} != spec {I_VTH}")
    return t, v


# --------------------------------------------------------------------------
# TC laws: 25 C base dict + temperature -> full Params
# --------------------------------------------------------------------------

# name -> (log?, lo, hi) in the optimisation vector
FREE_SPEC: dict[str, tuple[bool, float, float]] = {
    "k_p25": (True, 1e-2, 1e4),
    "beta_p": (False, -2.0, 8.0),          # k_p ~ T^-beta_p
    "v_th0_25": (False, 0.0, 10.0),
    "s_vth": (False, -0.03, 0.0),          # dVth0/dT [V/K]
    "v_sub": (True, 5e-3, 1.0),
    "k_bt25": (False, 0.0, 2.0),           # k_beta_trap(25 C)
    "s_bkt": (False, -1e-3, 1e-2),         # dk_beta_trap/dT [1/V/K] (linear law)
    "ea_bkt": (False, -500.0, 3000.0),     # Arrhenius variant [K]
    "k_theta": (True, 1e-4, 5.0),
    "k_ksh": (False, 0.1, 40.0),           # k_kneesharp (linear space)
    "k_eta": (True, 1e-5, 0.05),
    "v_dbsat": (True, 0.1, 1300.0),
    "k_lambda": (True, 1e-12, 1.0),
    "r_d25": (False, 0.0, 3.7e-3),         # <= datasheet 25 C chord
    "a_rd": (False, -5e-3, 5e-2),          # r_d TC linear  [1/K]
    "b_rd": (False, -1e-4, 1e-3),          # r_d TC quadratic [1/K^2]
    "gamma_rd": (False, 0.5, 6.0),         # r_d power-law variant
}

PRIMARY_FREE = (
    "k_p25", "beta_p", "v_th0_25", "s_vth", "v_sub", "k_bt25", "s_bkt",
    "k_theta", "k_ksh", "k_eta", "v_dbsat", "k_lambda",
    "r_d25", "a_rd", "b_rd",
)

VARIANTS: dict[str, dict] = {
    "primary": {"free": PRIMARY_FREE, "rd_law": "quad", "bkt_law": "lin",
                "groups": ("iv", "transfer", "rdson", "vth")},
    "arrhenius": {"free": tuple(n for n in PRIMARY_FREE if n != "s_bkt") + ("ea_bkt",),
                  "rd_law": "quad", "bkt_law": "arr",
                  "groups": ("iv", "transfer", "rdson", "vth")},
    "rd_powerlaw": {"free": tuple(n for n in PRIMARY_FREE if n not in ("a_rd", "b_rd"))
                    + ("gamma_rd",),
                    "rd_law": "pow", "bkt_law": "lin",
                    "groups": ("iv", "transfer", "rdson", "vth")},
    "no_ktheta": {"free": tuple(n for n in PRIMARY_FREE if n != "k_theta"),
                  "rd_law": "quad", "bkt_law": "lin",
                  "groups": ("iv", "transfer", "rdson", "vth")},
    "no_klambda": {"free": tuple(n for n in PRIMARY_FREE if n != "k_lambda"),
                   "rd_law": "quad", "bkt_law": "lin",
                   "groups": ("iv", "transfer", "rdson", "vth")},
    "no_transfer": {"free": PRIMARY_FREE, "rd_law": "quad", "bkt_law": "lin",
                    "groups": ("iv",)},
}

START = {
    "k_p25": 25.7, "beta_p": 0.4, "v_th0_25": 6.06, "s_vth": -7.0e-3,
    "v_sub": 0.76, "k_bt25": 0.057, "s_bkt": 1.25e-4, "ea_bkt": 250.0,
    "k_theta": 0.02, "k_ksh": 1.0, "k_eta": 0.013, "v_dbsat": 1.4,
    "k_lambda": 1e-3, "r_d25": 1.73e-3, "a_rd": 4.3e-3, "b_rd": 3.2e-5,
    "gamma_rd": 3.3,
}

# multi-start seeds for the primary fit (channel-shape knobs only; the rest
# from START). Mirrors the gain-probe seed families.
SEEDS = [
    {"k_p25": 25.7, "k_bt25": 0.057, "v_sub": 0.76, "k_theta": 0.01,
     "k_lambda": 1e-4, "k_ksh": 1.0},
    {"k_p25": 50.0, "k_bt25": 0.03, "v_sub": 0.50, "k_theta": 0.05,
     "k_lambda": 1e-2, "k_ksh": 1.5},
    {"k_p25": 120.0, "k_bt25": 0.10, "v_sub": 0.30, "k_theta": 0.02,
     "k_lambda": 1e-3, "k_ksh": 0.7},
    {"k_p25": 300.0, "k_bt25": 0.06, "v_sub": 0.85, "k_theta": 1e-3,
     "k_lambda": 1e-6, "k_ksh": 2.0},
    {"k_p25": 15.0, "k_bt25": 0.08, "v_sub": 0.60, "k_theta": 0.10,
     "k_lambda": 3e-2, "k_ksh": 1.0},
    {"k_p25": 10.0, "k_bt25": 0.45, "v_sub": 0.79, "k_theta": 0.11,
     "k_lambda": 1e-6, "k_ksh": 1.4},
]


def params_at(base: dict, t_c: float, rd_law: str, bkt_law: str) -> Params:
    """Full Params at temperature t_c [degC] from the 25 C base + TC laws."""
    tk = t_c + 273.15
    dt = t_c - 25.0
    if bkt_law == "lin":
        kbt = max(base["k_bt25"] + base["s_bkt"] * dt, 0.0)
    else:  # Arrhenius, rises with T for ea > 0
        kbt = base["k_bt25"] * np.exp(base["ea_bkt"] * (1.0 / T0K - 1.0 / tk))
    if rd_law == "quad":
        fac = 1.0 + base["a_rd"] * dt + base["b_rd"] * dt * dt
        r_d = base["r_d25"] * max(fac, 1e-6)  # Trs1/Trs2 clamp precedent
    else:
        r_d = base["r_d25"] * (tk / T0K) ** base["gamma_rd"]
    return Params(
        k_p=base["k_p25"] * (tk / T0K) ** (-base["beta_p"]),
        v_th0=base["v_th0_25"] + base["s_vth"] * dt,
        v_sub=base["v_sub"],
        k_beta_trap=kbt,
        k_theta=base["k_theta"],
        v_vsat=1e4,
        k_kneesharp=base["k_ksh"],
        k_eta=base["k_eta"],
        k_lambda=base["k_lambda"],
        g_k=0.0,
        v_dbsat=base["v_dbsat"],
        r_d=r_d,
        r_s=0.0,
    )


# --------------------------------------------------------------------------
# joint residual vector
# --------------------------------------------------------------------------


class JointData:
    """All four data groups, pre-grouped by temperature for vectorized eval."""

    def __init__(self, include: tuple[str, ...]) -> None:
        iv25 = load_csv(ROOT / "data" / "IMCQ120R004M2H_IV_25C.csv")
        iv175 = load_csv(ROOT / "data" / "IMCQ120R004M2H_IV_175C.csv")
        vg_tr, t_tr, i_tr = load_transfer(ROOT / "data" / "IMCQ120R004M2H_Transfer.csv")
        t_rd, vg_rd, r_rd = load_rdson(ROOT / "data" / "IMCQ120R004M2H_Rdson_Tj_185.2A.csv")
        t_vt, v_vt = load_vth(ROOT / "data" / "IMCQ120R004M2H_Vth_Tj.csv")

        # flat (v_ds, v_gs, i_ref, T, group) rows
        rows = []
        if "iv" in include:
            for d, t_c in ((iv25, 25.0), (iv175, 175.0)):
                rows += [(vd, vg, i, t_c, "iv")
                         for vg, vd, i in zip(d.v_gs, d.v_ds, d.i_d, strict=True)]
        if "transfer" in include:
            rows += [(20.0, vg, i, t, "transfer")
                     for vg, t, i in zip(vg_tr, t_tr, i_tr, strict=True)]
        if "rdson" in include:
            # chord residual as CURRENT at the measured test voltage:
            # i_model(R_data * 185.2, V_gs, T) vs 185.2 A (relative)
            rows += [(r_ * I_TEST, vg, I_TEST, t, "rdson")
                     for t, vg, r_ in zip(t_rd, vg_rd, r_rd, strict=True)]
        if "vth" in include:
            rows += [(v, v, I_VTH, t, "vth") for t, v in zip(t_vt, v_vt, strict=True)]

        self.v_ds = np.array([r[0] for r in rows])
        self.v_gs = np.array([r[1] for r in rows])
        self.i_ref = np.array([r[2] for r in rows])
        self.t_c = np.array([r[3] for r in rows])
        self.group = np.array([r[4] for r in rows])
        self.temps = sorted({round(float(t), 9) for t in self.t_c})
        self.masks = {t: (np.isclose(self.t_c, t)) for t in self.temps}
        self.iv25 = iv25
        self.iv175 = iv175
        self.rd_T, self.rd_vg, self.rd_chord = t_rd, vg_rd, r_rd
        self.vt_T, self.vt_vth = t_vt, v_vt

    def currents(self, base: dict, rd_law: str, bkt_law: str) -> np.ndarray:
        out = np.zeros_like(self.v_ds)
        for t, m in self.masks.items():
            p = params_at(base, t, rd_law, bkt_law)
            out[m] = vdmos_id(self.v_ds[m], self.v_gs[m], p)
        return out

    def residuals(self, base: dict, rd_law: str, bkt_law: str) -> np.ndarray:
        i = self.currents(base, rd_law, bkt_law)
        return (i - self.i_ref) / np.maximum(np.abs(self.i_ref), 1e-3)


def pack_x(base: dict, free: tuple[str, ...]) -> np.ndarray:
    x = []
    for n in free:
        v = base[n]
        x.append(np.log(v) if FREE_SPEC[n][0] else v)
    return np.asarray(x, dtype=float)


def unpack_x(x: np.ndarray, free: tuple[str, ...], base: dict) -> dict:
    out = dict(base)
    for n, xi in zip(free, x, strict=True):
        out[n] = float(np.exp(xi)) if FREE_SPEC[n][0] else float(xi)
    return out


def fit_variant(name: str, start: dict,
                seeds: list[dict] | None = None) -> dict:
    """Fit one variant; the variant's `groups` config selects the data used
    (the no_transfer ablation must actually drop the spec groups)."""
    cfg = VARIANTS[name]
    data = JointData(include=cfg["groups"])
    free = cfg["free"]
    # raw (physical-unit) bounds; the optimisation vector maps log params
    # through log, exactly like vdmos_fit.params.vector_bounds
    raw_lo = np.array([FREE_SPEC[n][1] for n in free])
    raw_hi = np.array([FREE_SPEC[n][2] for n in free])
    # keep a free r_d25 strictly inside the chord bound (0 is the "no series R"
    # sentinel neighbourhood; same doctrine as the per-T r_d cap)
    if "r_d25" in free:
        raw_lo[free.index("r_d25")] = 1e-6
    lo = np.array([np.log(a) if FREE_SPEC[n][0] else a
                   for n, a in zip(free, raw_lo, strict=True)])
    hi = np.array([np.log(b) if FREE_SPEC[n][0] else b
                   for n, b in zip(free, raw_hi, strict=True)])

    def res(x):
        return data.residuals(unpack_x(x, free, start), cfg["rd_law"], cfg["bkt_law"])

    best = None
    starts = [{}] + (seeds or [])
    for s in starts:
        st = {**start, **s}
        x0 = np.clip(pack_x(st, free), lo + 1e-9, hi - 1e-9)
        sol = least_squares(res, x0, bounds=(lo, hi), method="trf",
                            x_scale="jac", max_nfev=8000,
                            ftol=1e-13, xtol=1e-13, gtol=1e-13)
        tag = "converged" if sol.success else "stopped"
        if starts is not None and len(starts) > 1:
            print(f"   [{name}] start k_p25={st['k_p25']:7.4g} k_bt={st['k_bt25']:6.3g}"
                  f" v_sub={st['v_sub']:5.3g} k_theta={st['k_theta']:7.3g}"
                  f" k_lam={st['k_lambda']:8.3g} -> cost {sol.cost:.4e} [{tag}]")
        if best is None or sol.cost < best.cost:
            best = sol
    base = unpack_x(best.x, free, start)
    rel = data.residuals(base, cfg["rd_law"], cfg["bkt_law"])
    per_group = {g: float(np.sqrt(np.mean(rel[data.group == g] ** 2)) * 100)
                 for g in ("iv", "transfer", "rdson", "vth")
                 if np.any(data.group == g)}
    overall = float(np.sqrt(np.mean(rel ** 2)) * 100)
    # identifiability (same semantics as vdmos_fit.fit.identifiability)
    jac = best.jac
    col = np.linalg.norm(jac, axis=0)
    tol = col.max() * 1e-5
    dead = col < tol
    cov = np.linalg.pinv(jac.T @ jac)
    sd = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    sigma = {n: (float("inf") if d else float(s))
             for n, d, s in zip(free, dead, sd, strict=True)}
    ok = ~dead
    pairs = []
    if ok.sum() > 1:
        sdc = sd.copy()
        sdc[~ok] = 1.0
        corr = cov / np.outer(sdc, sdc)
        for i, a in enumerate(free):
            for j, b in enumerate(free):
                if j > i and ok[i] and ok[j] and abs(corr[i, j]) > 0.8:
                    pairs.append((a, b, float(corr[i, j])))
    pinned = {}
    for i, n in enumerate(free):
        v, (a, b) = base[n], (raw_lo[i], raw_hi[i])
        m = max(abs(a), abs(v), 1e-5)
        if abs(v - a) <= 1e-3 * m:
            pinned[n] = f"LOWER {a:g}"
        elif abs(v - b) <= 1e-3 * max(abs(b), abs(v), 1e-5):
            pinned[n] = f"UPPER {b:g}"
    return {"name": name, "base": base, "free": free, "cfg": cfg, "sol": best,
            "data": data,
            "rmse_overall": overall, "rmse_group": per_group,
            "sigma_rel": sigma, "corr_pairs": pairs, "pinned": pinned,
            "residuals": rel}


# --------------------------------------------------------------------------
# acceptance checks (independent specs, never fitted)
# --------------------------------------------------------------------------


def model_ron(base: dict, vg: float, t_c: float, rd_law: str, bkt_law: str) -> float:
    p = params_at(base, t_c, rd_law, bkt_law)
    v = brentq(lambda x: float(vdmos_id(x, vg, p)) - I_TEST, 1e-4, 10.0)
    return v / I_TEST


def acceptance(base: dict, rd_law: str, bkt_law: str) -> dict:
    out = {}
    p25 = params_at(base, 25.0, rd_law, bkt_law)
    p175 = params_at(base, 175.0, rd_law, bkt_law)
    p200 = params_at(base, 200.0, rd_law, bkt_law)

    # g_fs spec: 122 S at 185.2 A / 20 V, 25 C
    vg = brentq(lambda g: float(vdmos_id(20.0, g, p25)) - I_TEST, 5.0, 18.0)
    h = 0.01
    gfs = (float(vdmos_id(20.0, vg + h, p25)) - float(vdmos_id(20.0, vg - h, p25))) / (2 * h)
    out["g_fs_S"] = gfs
    out["g_fs_bias_vg"] = vg

    # ZTC crossings on the transfer (V_DS = 20 V)
    f25_175 = lambda g: (float(vdmos_id(20.0, g, p25))  # noqa: E731
                         - float(vdmos_id(20.0, g, p175)))
    f175_200 = lambda g: (float(vdmos_id(20.0, g, p175))  # noqa: E731
                          - float(vdmos_id(20.0, g, p200)))
    try:
        vz = brentq(f25_175, 8.0, 19.5)
        out["ztc_25_175"] = (vz, float(vdmos_id(20.0, vz, p25)))
    except ValueError:
        out["ztc_25_175"] = None
    try:
        vz2 = brentq(f175_200, 8.0, 19.5)
        out["ztc_175_200"] = (vz2, float(vdmos_id(20.0, vz2, p175)))
    except ValueError:
        out["ztc_175_200"] = None

    # Ron chords (brentq solve at 185.2 A)
    out["ron_25_vg18"] = model_ron(base, 18.0, 25.0, rd_law, bkt_law)
    out["ron_25_vg15"] = model_ron(base, 15.0, 25.0, rd_law, bkt_law)
    out["ron_175_vg18"] = model_ron(base, 18.0, 175.0, rd_law, bkt_law)
    out["ron_175_vg15"] = model_ron(base, 15.0, 175.0, rd_law, bkt_law)

    # Vth spec points (25/175 C IV anchors)
    out["vth_25_mA"] = float(vdmos_id(4.2, 4.2, p25)) * 1e3
    out["vth_175_mA"] = float(vdmos_id(3.2, 3.2, p175)) * 1e3

    # high-v_gs corner at V_DS = 20, v_gs = 18 (was 3650 A in the gain fit)
    out["corner_18V"] = {
        "25C": float(vdmos_id(20.0, 18.0, p25)),
        "175C": float(vdmos_id(20.0, 18.0, p175)),
        "200C": float(vdmos_id(20.0, 18.0, p200)),
    }
    out["corner_data_interp"] = {"25C": 3213.0, "175C": 2070.0, "200C": 1836.0}

    # Rds(on)-Tj full span: chord solve per digitized point
    errs = []
    for t, vg_, r_ in zip(base["_rd_T"], base["_rd_vg"], base["_rd_chord"], strict=True):
        m = model_ron(base, vg_, t, rd_law, bkt_law)
        errs.append((t, vg_, (m - r_) / r_ * 100.0))
    e = np.array([x[2] for x in errs])
    out["rdson_errs"] = errs
    out["rdson_rmse"] = float(np.sqrt(np.mean(e ** 2)))
    out["rdson_max"] = float(np.abs(e).max())
    out["rdson_cold_rmse"] = float(np.sqrt(np.mean(e[[x[0] < 0 for x in errs]] ** 2)))
    out["rdson_warm_rmse"] = float(np.sqrt(np.mean(e[[x[0] >= 0 for x in errs]] ** 2)))
    return out


# --------------------------------------------------------------------------
# plots
# --------------------------------------------------------------------------


def plot_transfer(data: JointData, base: dict, rd_law: str, bkt_law: str) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4))
    colors = {25.0: "C0", 175.0: "C1", 200.0: "C2"}
    for t in (25.0, 175.0, 200.0):
        p = params_at(base, t, rd_law, bkt_law)
        vgs = np.linspace(6.5, 20.0, 200)
        for ax in axes:
            ax.plot(vgs, vdmos_id(20.0, vgs, p), "-", color=colors[t],
                    label=f"model {t:g} C")
    m = data.group == "transfer"
    for t in (25.0, 175.0, 200.0):
        sel = m & np.isclose(data.t_c, t)
        axes[0].plot(data.v_gs[sel], data.i_ref[sel], "o", color=colors[t], ms=6,
                     label=f"data {t:g} C")
        axes[1].plot(data.v_gs[sel], data.i_ref[sel], "o", color=colors[t], ms=6)
    for ax, t in zip(axes, ["Transfer @ V_DS = 20 V: joint fit vs digitized data",
                            "same, log scale (ZTC crossing)"], strict=True):
        ax.set_xlabel("v_gs [V]")
        ax.set_ylabel("i_d [A]")
        ax.set_title(t)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    axes[1].set_yscale("log")
    fig.savefig(OUT / "imcq120_joint_transfer.png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    print("wrote", OUT / "imcq120_joint_transfer.png")


def plot_iv(data: JointData, base: dict, rd_law: str, bkt_law: str) -> None:
    for d, t_c, fname in ((data.iv25, 25.0, "imcq120_joint_iv_25C.png"),
                          (data.iv175, 175.0, "imcq120_joint_iv_175C.png")):
        p = params_at(base, t_c, rd_law, bkt_law)
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
        axes[0].set_xlim(0, 21)
        axes[1].set_xlim(0, 3)
        axes[1].set_ylim(0, 120)
        for ax, ttl in zip(axes, [f"Output curves {t_c:g} C (joint fit)",
                                  f"{t_c:g} C knee zoom"], strict=True):
            ax.set_xlabel("v_ds [V]")
            ax.set_ylabel("i_d [A]")
            ax.set_title(ttl)
            ax.grid(alpha=0.3)
        axes[0].legend(fontsize=7, ncol=2, title="dots = data, dashed = fit")
        fig.savefig(OUT / fname, dpi=110, bbox_inches="tight")
        plt.close(fig)
        print("wrote", OUT / fname)


def plot_rdson(base: dict, rd_law: str, bkt_law: str) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 5.6))
    T = base["_rd_T"]
    for vg, c, mk in ((18.0, "C1", "o"), (15.0, "C2", "s")):
        m = base["_rd_vg"] == vg
        ax.plot(T[m], base["_rd_chord"][m] * 1e3, mk, color=c, ms=5,
                label=f"data Vgs={vg:g}")
        Ts = np.linspace(-50, 200, 80)
        pred = [model_ron(base, vg, t, rd_law, bkt_law) * 1e3 for t in Ts]
        ax.plot(Ts, pred, "-", color=c, label=f"model Vgs={vg:g}")
    ax.axvline(25, color="gray", ls=":")
    ax.axvline(175, color="gray", ls=":")
    ax.set_xlabel("T_j [degC]")
    ax.set_ylabel("Rds(on) chord @ 185.2 A [m$\\Omega$]")
    ax.set_title("Rds(on)-Tj: joint fit with quadratic r_d(T) law, full -46..196 C span")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.savefig(OUT / "imcq120_joint_ron_Tj.png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    print("wrote", OUT / "imcq120_joint_ron_Tj.png")


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------


def report_variant(r: dict) -> None:
    print(f"\n=== variant {r['name']}: overall RMSE {r['rmse_overall']:.2f} %"
          f"  per-group: " + "  ".join(f"{g} {v:.2f}%" for g, v in r["rmse_group"].items()))
    for n in r["free"]:
        s = r["sigma_rel"][n]
        pin = r["pinned"].get(n, "")
        print(f"   {n:<10s} = {r['base'][n]:.6g}   sigma_rel {s:9.3g}"
              + (f"   <-- pinned {pin}" if pin else ""))
    if r["corr_pairs"]:
        print("   |corr|>0.8:", ", ".join(f"{a}~{b} {c:+.2f}"
                                          for a, b, c in r["corr_pairs"]))
    # per-curve IV RMSE
    rel = r["residuals"]
    data = r["data"]
    for tag, d, t_c in (("25C", data.iv25, 25.0), ("175C", data.iv175, 175.0)):
        m = (data.group == "iv") & np.isclose(data.t_c, t_c)
        if not np.any(m):
            continue
        per = {}
        for vg in d.v_gs_grid():
            mm = m & np.isclose(data.v_gs, vg)
            if np.any(mm):
                per[vg] = round(float(np.sqrt(np.mean(rel[mm] ** 2)) * 100), 2)
        print(f"   {tag} per-curve RMSE %: {per}")


def main() -> None:
    OUT.mkdir(exist_ok=True)
    data = JointData(include=("iv", "transfer", "rdson", "vth"))
    counts = {g: int((data.group == g).sum())
              for g in ("iv", "transfer", "rdson", "vth")}
    print(f"joint data: {len(data.v_ds)} pts {counts}, "
          f"{len(data.temps)} distinct temperatures")

    start = dict(START)
    start["_rd_T"], start["_rd_vg"], start["_rd_chord"] = data.rd_T, data.rd_vg, data.rd_chord

    results = {}
    primary = fit_variant("primary", start, seeds=SEEDS)
    results["primary"] = primary
    report_variant(primary)

    acc = acceptance(primary["base"], "quad", "lin")
    print("\n=== acceptance checks (primary) ===")
    print(f"   g_fs @185.2 A/20 V: {acc['g_fs_S']:.1f} S (spec 122, "
          f"{(acc['g_fs_S'] - 122) / 122 * 100:+.1f} %) at v_gs={acc['g_fs_bias_vg']:.2f}")
    if acc["ztc_25_175"]:
        print(f"   ZTC 25/175: v_gs={acc['ztc_25_175'][0]:.2f} V, "
              f"{acc['ztc_25_175'][1]:.0f} A (datasheet ~14.2 V / 1500 A)")
    else:
        print("   ZTC 25/175: NO CROSSING in 8..19.5 V")
    if acc["ztc_175_200"]:
        print(f"   ZTC 175/200: v_gs={acc['ztc_175_200'][0]:.2f} V, "
              f"{acc['ztc_175_200'][1]:.0f} A (digitized crossing 12.03 V / 999 A)")
    print(f"   Ron chords 25 C: {acc['ron_25_vg18'] * 1e3:.3f} / "
          f"{acc['ron_25_vg15'] * 1e3:.3f} mOhm (data 3.70 / 4.70)")
    print(f"   Ron chords 175 C: {acc['ron_175_vg18'] * 1e3:.3f} / "
          f"{acc['ron_175_vg15'] * 1e3:.3f} mOhm (data 8.90 / 9.43)")
    print(f"   Vth spec: {acc['vth_25_mA']:.1f} / {acc['vth_175_mA']:.1f} mA (58.3)")
    print("   corner v_gs=18 @20 V: " + "  ".join(
        f"{k} {v:.0f} A (data ~{acc['corner_data_interp'][k]:.0f})"
        for k, v in acc["corner_18V"].items()))
    print(f"   Rds(on)-Tj: RMSE {acc['rdson_rmse']:.2f}%  max {acc['rdson_max']:.2f}%"
          f"  cold(<0C) {acc['rdson_cold_rmse']:.2f}%  warm {acc['rdson_warm_rmse']:.2f}%")
    print("   per-point Rds(on) err %:", ", ".join(
        f"{t:g}C/{vg:g}V {e:+.2f}" for t, vg, e in acc["rdson_errs"]))

    # --- official-record drift check -------------------------------------
    # SENSED parameters only: unsensed directions (k_eta cap-pin, v_dbsat,
    # k_lambda -- all sigma inf) are solver-trajectory noise and scipy-version
    # sensitive (session-19 A/B lesson); drift there is NOT staleness.
    rec_fields = ("k_p25", "beta_p", "v_th0_25", "s_vth", "v_sub", "k_bt25",
                  "s_bkt", "k_theta", "k_ksh", "k_eta", "v_dbsat", "k_lambda",
                  "r_d25", "a_rd", "b_rd")
    drift = {n: abs(primary["base"][n] - getattr(IMCQ120_JOINT, n))
             / max(abs(getattr(IMCQ120_JOINT, n)), 1e-12)
             for n in primary["free"]
             if n in rec_fields and not np.isinf(primary["sigma_rel"][n])}
    dmax = max(drift.values())
    worst = max(drift, key=lambda k: drift[k])
    print(f"\nrecord drift (sensed only): max rel {dmax:.2e} on {worst}"
          f"  [{'OK' if dmax < 1e-3 else 'STALE -- re-quote the record'}]")
    # record reproducibility through the script's own machinery (cross-checks
    # the reference.py params_at law implementation against this script's)
    rec_base = {n: float(getattr(IMCQ120_JOINT, n)) for n in rec_fields}
    rec_base.update({"ea_bkt": 250.0, "gamma_rd": 3.3})
    rec_rel = data.residuals(rec_base, "quad", "lin")
    rec_rmse = float(np.sqrt(np.mean(rec_rel**2)))
    print(f"record RMSE re-eval: {rec_rmse * 100:.3f} % vs published "
          f"{IMCQ120_JOINT_RMSE_REL * 100:.2f} %"
          f"  [{'OK' if abs(rec_rmse - IMCQ120_JOINT_RMSE_REL) < 2e-3 else 'STALE'}]")

    # variants: local basin from the primary optimum; knobs pinned by the
    # ablation are set to their OFF value first (what FREEING buys, not what
    # the primary value buys)
    base_p = primary["base"]
    for name in ("arrhenius", "rd_powerlaw", "no_ktheta", "no_klambda", "no_transfer"):
        st = {**base_p}
        if name == "rd_powerlaw":
            st["gamma_rd"] = 3.3
        if name == "arrhenius":
            st["ea_bkt"] = 250.0
        if name == "no_ktheta":
            st["k_theta"] = 1e-4      # off value (identity knee)
        if name == "no_klambda":
            st["k_lambda"] = 1e-12    # off value (no CLM)
        if name == "no_transfer":
            r = fit_variant(name, start, seeds=[base_p])
        else:
            r = fit_variant(name, st, seeds=None)
        results[name] = r
        report_variant(r)
        a2 = acceptance(r["base"], r["cfg"]["rd_law"], r["cfg"]["bkt_law"])
        print(f"   [{name}] corner18 175C {a2['corner_18V']['175C']:.0f} A"
              f"  rdson RMSE {a2['rdson_rmse']:.2f}% cold {a2['rdson_cold_rmse']:.2f}%"
              f"  g_fs {a2['g_fs_S']:.1f} S")

    # artifacts
    plot_transfer(data, base_p, "quad", "lin")
    plot_iv(data, base_p, "quad", "lin")
    plot_rdson(base_p, "quad", "lin")

    def ser(r: dict) -> dict:
        b = {k: v for k, v in r["base"].items() if not k.startswith("_")}
        return {"base_25C_and_TC": b, "free": list(r["free"]),
                "rmse_overall_pct": r["rmse_overall"],
                "rmse_group_pct": r["rmse_group"],
                "sigma_rel": r["sigma_rel"], "pinned": r["pinned"],
                "corr_pairs": r["corr_pairs"]}

    payload = {k: ser(r) for k, r in results.items()}
    payload["acceptance_primary"] = {k: (v if not isinstance(v, list) else v)
                                     for k, v in acc.items()}
    payload["record_drift_sensed_max_rel"] = dmax
    payload["record_rmse_reeval_pct"] = rec_rmse * 100
    payload["official_record"] = "reference.py::IMCQ120_JOINT"
    (OUT / "imcq120_joint_params.json").write_text(
        json.dumps(payload, indent=2, default=float), encoding="utf-8")
    print("wrote", OUT / "imcq120_joint_params.json")


if __name__ == "__main__":
    main()
