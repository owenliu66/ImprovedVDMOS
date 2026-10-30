"""Pins the OFFICIAL IMCQ120R004M2H record (vdmos_fit.reference.IMCQ120_JOINT).

Promoted session 21 (user go-ahead): the joint multi-temperature fit of
scripts/fit_imcq120_joint.py -- ONE fit over IV 25/175 C + transfer @20 V
(25/175/200 C) + Rds(on)-Tj chords + Vth-Tj (200 pts, 29 temperatures) with
the TC laws as the only extra unknowns; overall 1.93 % relative RMSE.

Like test_st8l65_official.py these tests do NOT re-run the fit (too slow for
the unit suite): they evaluate the RECORDED TC-law parameter set against all
four measured CSVs and against the datasheet acceptance specs, so a model /
bounds / TC-law change that silently moves the official fit fails here.

Deliberately NOT pinned as physical results (see reference.py record notes):
k_p25/k_bt25 individually (corr -0.95 -- the composite k_p*v_ov_eff2 is the
record), beta_p/s_bkt individually (+0.97), and the inert DIBL/CLM block
(k_eta pinned at the 0.05 cap, v_dbsat/k_lambda sigma = inf).

The two KNOWN OPEN GAPS are pinned as tolerances so they cannot silently
worsen while the Kelvin-aware subckt (queue) is being built: g_fs +6.2 % vs
the 122 S spec, and the 25 C transfer corner -5..-7 % at v_gs 15-19.
Session-22 experiment (out/IMCQ120_KELVIN_SPLIT_REPORT.md): neither gap is
a series-R topology problem -- gate-loop degeneration (r_s inside the gate
loop) is REJECTED by the data, the current-loop split is DC-degenerate, and
the terminal g_fs is pinned by the transfer curve itself (~130 S). The
tolerances below stay exactly as they were; only the explanation changed.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import brentq

from vdmos_fit import (
    IMCQ120_JOINT,
    IMCQ120_JOINT_GROUP_RMSE_REL,
    IMCQ120_JOINT_RMSE_REL,
    load_csv,
    vdmos_id,
)

DATA = Path(__file__).resolve().parents[1] / "data"
I_TEST = 185.2   # Rds(on) / g_fs test current [A]
I_VTH = 0.0583   # Vth spec current [A]


def _rows(name: str) -> list[dict]:
    with (DATA / name).open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _load_transfer() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(v_gs, T_c, i_d) from `index,V_gs,I_d,condition` (V_ds = 20 V fixed)."""
    vg, tc, id_ = [], [], []
    for r in _rows("IMCQ120R004M2H_Transfer.csv"):
        cond = r["condition"].lower()
        assert "v_ds=20v" in cond, f"transfer fit assumes V_ds = 20 V: {cond!r}"
        vg.append(float(r["V_gs"]))
        tc.append(float(cond.split("t_j=")[1].rstrip("c")))
        id_.append(float(r["I_d"]))
    return np.array(vg), np.array(tc), np.array(id_)


def _load_rdson() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(T_c, V_gs, chord Ohm) from `index,Tj (C),R_dson (mOhm),V_gs`."""
    rows = _rows("IMCQ120R004M2H_Rdson_Tj_185.2A.csv")
    return (np.array([float(r["Tj (C)"]) for r in rows]),
            np.array([float(r["V_gs"]) for r in rows]),
            np.array([float(r["R_dson (mOhm)"]) for r in rows]) * 1e-3)


def _load_vth() -> tuple[np.ndarray, np.ndarray]:
    rows = _rows("IMCQ120R004M2H_Vth_Tj.csv")
    assert all(np.isclose(float(r["I_d"]), I_VTH) for r in rows)
    return (np.array([float(r["T_j (C)"]) for r in rows]),
            np.array([float(r["V_th"]) for r in rows]))


@pytest.fixture(scope="module")
def joint():
    """The 200-point joint set, grouped exactly like scripts/fit_imcq120_joint.py."""
    d25 = load_csv(DATA / "IMCQ120R004M2H_IV_25C.csv")
    d175 = load_csv(DATA / "IMCQ120R004M2H_IV_175C.csv")
    vg_tr, t_tr, i_tr = _load_transfer()
    t_rd, vg_rd, r_rd = _load_rdson()
    t_vt, v_vt = _load_vth()

    v_ds, v_gs, i_ref, t_c, group = [], [], [], [], []
    for d, t in ((d25, 25.0), (d175, 175.0)):
        for vg, vd, i in zip(d.v_gs, d.v_ds, d.i_d, strict=True):
            v_ds.append(vd)
            v_gs.append(vg)
            i_ref.append(i)
            t_c.append(t)
            group.append("iv")
    for vg, t, i in zip(vg_tr, t_tr, i_tr, strict=True):
        v_ds.append(20.0)
        v_gs.append(vg)
        i_ref.append(i)
        t_c.append(t)
        group.append("transfer")
    for t, vg, r_ in zip(t_rd, vg_rd, r_rd, strict=True):
        v_ds.append(r_ * I_TEST)
        v_gs.append(vg)
        i_ref.append(I_TEST)
        t_c.append(t)
        group.append("rdson")
    for t, v in zip(t_vt, v_vt, strict=True):
        v_ds.append(v)
        v_gs.append(v)
        i_ref.append(I_VTH)
        t_c.append(t)
        group.append("vth")

    return {
        "v_ds": np.array(v_ds), "v_gs": np.array(v_gs),
        "i_ref": np.array(i_ref), "t_c": np.array(t_c),
        "group": np.array(group),
        "iv25": d25, "iv175": d175,
        "rd_T": t_rd, "rd_vg": vg_rd, "rd_chord": r_rd,
    }


def _model(joint) -> np.ndarray:
    out = np.zeros_like(joint["v_ds"])
    for t in {round(float(x), 9) for x in joint["t_c"]}:
        m = np.isclose(joint["t_c"], t)
        p = IMCQ120_JOINT.params_at(t)
        out[m] = vdmos_id(joint["v_ds"][m], joint["v_gs"][m], p)
    return out


def _chord(vg: float, t_c: float) -> float:
    """Model Rds(on) chord at the datasheet test current (brentq, as in the script)."""
    p = IMCQ120_JOINT.params_at(t_c)
    v = brentq(lambda x: float(vdmos_id(x, vg, p)) - I_TEST, 1e-4, 10.0)
    return v / I_TEST


# --------------------------------------------------------------------------
# RMSE reproduction (record evaluation is deterministic -- no refit)
# --------------------------------------------------------------------------


def test_record_reproduces_joint_rmse(joint):
    rel = (_model(joint) - joint["i_ref"]) / np.maximum(np.abs(joint["i_ref"]), 1e-3)
    rmse = float(np.sqrt(np.mean(rel**2)))
    assert rmse == pytest.approx(IMCQ120_JOINT_RMSE_REL, abs=2e-3)
    assert rmse < 0.021, "official joint fit regressed"


def test_record_per_group_rmse(joint):
    rel = (_model(joint) - joint["i_ref"]) / np.maximum(np.abs(joint["i_ref"]), 1e-3)
    for g, ref in IMCQ120_JOINT_GROUP_RMSE_REL.items():
        r = float(np.sqrt(np.mean(rel[joint["group"] == g] ** 2)))
        assert r == pytest.approx(ref, abs=2e-3), f"{g} group RMSE moved"
    # absolute ceilings (the transfer group carries the digitization limit:
    # 175/200 C curves indistinguishable below 12 V, data/README.md)
    assert float(np.sqrt(np.mean(rel[joint["group"] == "iv"] ** 2))) < 0.020
    assert float(np.sqrt(np.mean(rel[joint["group"] == "transfer"] ** 2))) < 0.040
    assert float(np.sqrt(np.mean(rel[joint["group"] == "rdson"] ** 2))) < 0.010
    assert float(np.sqrt(np.mean(rel[joint["group"] == "vth"] ** 2))) < 0.030


# --------------------------------------------------------------------------
# Acceptance specs (never fitted) -- the two FIXED failures pinned hard
# --------------------------------------------------------------------------


def test_rdson_full_span_cold_end_fixed(joint):
    """Rds(on)-Tj over the FULL -46..196 C span: max |err| < 1.5 %.

    The two-temperature power-law protocol failed at -15.3 % at -46 C; the
    quadratic r_d(T) law (user call) fixed it. Pinning the cold end here is
    what stops a future TC-law change from silently reintroducing that.
    """
    errs = []
    for t, vg, r_ in zip(joint["rd_T"], joint["rd_vg"], joint["rd_chord"],
                         strict=True):
        errs.append(((_chord(vg, t) - r_) / r_) * 100.0)
    errs = np.array(errs)
    assert np.abs(errs).max() < 1.5
    cold = np.abs(errs[joint["rd_T"] < 0.0])
    assert np.sqrt(np.mean(cold**2)) < 1.2, "cold-end Rds(on) TC regressed"


def test_ron_chord_specs(joint):
    assert _chord(18.0, 25.0) * 1e3 == pytest.approx(3.70, rel=0.015)
    assert _chord(15.0, 25.0) * 1e3 == pytest.approx(4.70, rel=0.015)
    assert _chord(18.0, 175.0) * 1e3 == pytest.approx(8.90, rel=0.03)
    assert _chord(15.0, 175.0) * 1e3 == pytest.approx(9.43, rel=0.03)


def test_high_vgs_corner_no_overshoot(joint):
    """The PART 1 fix: v_gs = 18 @ V_DS = 20 must stay near the digitized
    transfer figure (175 C ~2070 A). The pre-joint gain fit predicted 3650 A
    -- the 2400 A ceiling is the regression guard against that overshoot."""
    p175 = IMCQ120_JOINT.params_at(175.0)
    p200 = IMCQ120_JOINT.params_at(200.0)
    p25 = IMCQ120_JOINT.params_at(25.0)
    i175 = float(vdmos_id(20.0, 18.0, p175))
    assert i175 < 2400.0, "high-v_gs corner overshot again (gain-fit was 3650 A)"
    assert i175 == pytest.approx(2070.0, rel=0.04)
    assert float(vdmos_id(20.0, 18.0, p200)) == pytest.approx(1836.0, rel=0.04)
    # 25 C: chord-interpolated digitized value ~3213 A; the -7 % corner
    # deficit is a KNOWN open gap -- a channel-shape compromise, NOT a
    # series-R/Kelvin-path gap (session-22 experiment: gate degeneration
    # makes it strictly worse) -- pin it so it cannot silently worsen.
    assert float(vdmos_id(20.0, 18.0, p25)) == pytest.approx(3213.0, rel=0.10)


def test_ztc_crossings(joint):
    """ZTC 25/175 on the transfer: the model matches the DIGITIZED output-curve
    crossing (~14.9 V); the datasheet-figure ~14.2 V read is itself soft."""
    p25 = IMCQ120_JOINT.params_at(25.0)
    p175 = IMCQ120_JOINT.params_at(175.0)
    vz = brentq(lambda g: float(vdmos_id(20.0, g, p25))
                - float(vdmos_id(20.0, g, p175)), 8.0, 19.5)
    assert 13.5 < vz < 16.5
    assert 1200.0 < float(vdmos_id(20.0, vz, p25)) < 2200.0

    # 175/200 C digitized crossing point (12.025 V, 999.1 A) -- both curves
    # must land within 2 % there (the model separates them slightly more than
    # the data: documented, see out/IMCQ120_JOINT_FIT_REPORT.md).
    p200 = IMCQ120_JOINT.params_at(200.0)
    for p in (p175, p200):
        i = float(vdmos_id(20.0, 12.0253405, p))
        assert i == pytest.approx(999.1, rel=0.02)


def test_vth_spec_points(joint):
    assert float(vdmos_id(4.2, 4.2, IMCQ120_JOINT.params_at(25.0))) * 1e3 \
        == pytest.approx(58.3, rel=0.03)
    assert float(vdmos_id(3.2, 3.2, IMCQ120_JOINT.params_at(175.0))) * 1e3 \
        == pytest.approx(58.3, rel=0.03)


def test_gfs_known_gap_pinned(joint):
    """g_fs @185.2 A / 20 V: 129.5 S vs the 122 S spec (+6.2 %). KNOWN open
    gap -- pinned so it cannot silently worsen. Session 22 showed the
    terminal slope at this bias is pinned by the digitized transfer curve
    itself (~130 S): no series-R topology (r_d/r_s/r_sk, degenerated or
    not) can bring it to 122 without breaking the curve fit -- this is a
    spec-vs-digitized-slope softness, not missing physics."""
    p = IMCQ120_JOINT.params_at(25.0)
    vg = brentq(lambda g: float(vdmos_id(20.0, g, p)) - I_TEST, 5.0, 18.0)
    h = 0.01
    gfs = (float(vdmos_id(20.0, vg + h, p)) - float(vdmos_id(20.0, vg - h, p))) / (2 * h)
    assert gfs == pytest.approx(122.0, rel=0.10)
    assert gfs > 115.0


# --------------------------------------------------------------------------
# TC-law physical invariants over the data span
# --------------------------------------------------------------------------


def test_tc_laws_physical_over_span():
    for t in np.linspace(-46.0, 200.0, 25):
        p = IMCQ120_JOINT.params_at(float(t))
        assert p.k_beta_trap >= 0.0          # hard constraint (model.py guard)
        assert p.r_d > 0.0                   # quadratic never negative in span
        assert p.k_p > 0.0
    # r_d carries the whole Ron TC: x3.2 at 175 C, flattening at the cold end
    p25 = IMCQ120_JOINT.params_at(25.0)
    p175 = IMCQ120_JOINT.params_at(175.0)
    assert p175.r_d / p25.r_d == pytest.approx(3.22, abs=0.1)
    assert IMCQ120_JOINT.params_at(-46.0).r_d < p25.r_d * 0.85  # flattening
    # v_th0 slope -7.16 mV/K (datasheet V_th slope -6.67; fitted offset
    # ~ +1.5 V constant -- NEVER substitute the datasheet V_th for v_th0)
    slope = (p175.v_th0 - p25.v_th0) / 150.0
    assert slope == pytest.approx(-7.16e-3, abs=1e-4)
    # v_sub T-invariant (shape knob, not a thermal voltage)
    assert IMCQ120_JOINT.params_at(-46.0).v_sub == IMCQ120_JOINT.v_sub


def test_model_finite_and_monotone_on_transfer():
    """No NaN/inf; i_d strictly increasing in v_gs at V_DS = 20 V (the gain +
    knee-saturation model must stay monotone -- k_beta_trap >= 0 keeps it)."""
    vgs = np.linspace(6.0, 20.0, 60)
    for t in (25.0, 175.0, 200.0):
        p = IMCQ120_JOINT.params_at(t)
        i_d = np.asarray(vdmos_id(20.0, vgs, p))
        assert np.all(np.isfinite(i_d))
        assert np.all(np.diff(i_d) > 0.0), f"transfer TC went non-monotone at {t:g} C"
        assert float(vdmos_id(0.0, 18.0, p)) == 0.0


def test_record_sanity_and_documented_bound_pins():
    """Physical sanity of the record; the documented bound exceptions are the
    inert DIBL/CLM block (k_eta at the 0.05 cap, k_lambda ~ 0, v_dbsat
    sigma = inf) -- everything that MATTERS is interior."""
    r = IMCQ120_JOINT
    assert 0.0 < r.r_d25 < 3.7e-3, "r_d25 must stay below the 25 C datasheet chord"
    assert r.k_bt25 > 0.0 and r.s_bkt > 0.0     # trap-filling gain rises with T
    assert 0.0 < r.beta_p < 2.0
    assert -0.03 < r.s_vth < 0.0
    assert 0.5 < r.v_sub < 1.0                  # shape knob (near the 1.0 cap)
    assert r.k_theta > 0.0                      # freed session 21: knee position
    assert 0.1 < r.k_ksh < 40.0
    assert r.k_eta == pytest.approx(0.05, rel=1e-3)   # pinned at cap -- INERT
    assert r.k_lambda < 1e-8                    # ~ 0: CLM unnecessary
    assert r.r_s == 0.0                         # Kelvin source
    assert r.v_ovmax == 0.0 and r.g_k == 0.0    # cap off, no trench kink
