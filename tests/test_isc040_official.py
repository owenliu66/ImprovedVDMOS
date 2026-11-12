"""Pins the OFFICIAL ISC040N10NM8 records: channel (reference.ISC040_CHANNEL)
and body diode (reference.ISC040_DIODE).

Promoted session 27. Like test_imcq120_official.py / test_st8l65_official.py
these tests do NOT re-run the fits (too slow for the unit suite): they
evaluate the RECORDED parameter sets against the measured CSVs and against
the datasheet acceptance specs, so a model / bounds / TC-law change that
silently moves an official fit fails here.

Channel protocol (scripts/fit_isc040_channel.py -- three user notes):
NO transfer curve (datasheet V_ds conditions vague) -> output curves +
Rds(on) chords only; NOT Kelvin -> r_s = r_sk = 0 and the data REJECT gate-
loop degeneration (the IMCQ session-22 finding), whole current-loop series
R in r_d; fit 25 C FIRST and LOCK, then TC (no high-T IV data -- TC comes
from Rds(on)-vs-Vgs at -55/25/100/175 C + Vth-vs-Tj at 520 uA / 52 uA).
Model family: SJ charge-cap (v_ovmax/k_ovsharp free -- the cap ceilings the
V_gs=15 conductance; without it the Ron chord is ~13 % low).

Identifiability caveats (quote curves, never split): r_d<->r_s corr -1.0
(only the sum sensed); k_p<->v_th0<->v_ovmax<->r_d cluster |corr| <= 0.86
(the composite channel conductance is the record); beta_p weak (sigma ~ 4,
interpretable only because Rds(on)-Vgs(T) pins r_d(T) at high V_gs);
v_dbsat sigma ~ 77 (inert, data stops at V_ds ~ 5 V).

Diode caveats (scripts/fit_isc040_diode.py): k_n = 2.0 (the ST8L65 fixed
value) is WRONG here -- with k_n free RMSE drops 35 % -> 3.31 %, k_n lands
~1.05 (near-ideal diffusion). i_s<->k_n corr +0.999 (degenerate prefactor
pair -- quote the pair/curve, never i_s alone). The LINEAR r_diode TC
(a_rdiode) is the decisive TC (3.31 % -> 0.51 %); quadratic b_rdiode and
a_ikf buy nothing with only 2 temperatures and go degenerate (a_rdiode
<-> b_rdiode corr -1.0) -> they stay pinned at 0.

NOTE the CSV loaders here are LOCAL: the ISC040 CSVs carry a `T_j (C)`
temperature column, which the package loaders (expecting `T` / `Tj (C)`)
choke on -- see the HANDOVER gotchas.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import brentq

from vdmos_fit import (
    ISC040_CHANNEL,
    ISC040_CHANNEL_RMSE25_REL,
    ISC040_CHANNEL_RMSE_TC_REL,
    ISC040_DIODE,
    ISC040_DIODE_GROUP_RMSE_REL,
    ISC040_DIODE_RMSE_REL,
    ISC040_RECOV,
    load_csv,
    vdmos_diode,
    vdmos_id,
)

DATA = Path(__file__).resolve().parents[1] / "data"
OUT = Path(__file__).resolve().parents[1] / "out"
I_TEST = 25.0   # Rds(on)-vs-Vgs test current [A] (verified: Vgs=15.7 ->
                # 3.33 mOhm matches the IV_25C row Vgs=15, Vds=0.08325, I=25)


def _rows(name: str) -> list[dict]:
    with (DATA / name).open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _load_rdson_id(name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Rds(on)-vs-Id CSV `index,I_d,R_dson (mOhm),V_gs` -> ohmic IV points."""
    rows = _rows(name)
    vg = np.array([float(r["V_gs"]) for r in rows])
    vd = np.array([float(r["R_dson (mOhm)"]) * 1e-3 * float(r["I_d"])
                   for r in rows])
    id_ = np.array([float(r["I_d"]) for r in rows])
    return vg, vd, id_


def _load_rdson_vgs(name: str) -> list[tuple[float, float, float]]:
    """Rds(on)-vs-Vgs CSV -> (T_c, V_gs, chord_ohm) rows (I_test = 25 A)."""
    rows = _rows(name)
    return [(float(r["T_j (C)"]), float(r["V_gs"]),
             float(r["R_dson (mOhm)"]) * 1e-3) for r in rows]


def _load_vth(name: str) -> list[tuple[float, float, float]]:
    """Vth-vs-Tj CSV -> (T_c, V_th, I_d) rows."""
    rows = _rows(name)
    return [(float(r["T_j (C)"]), float(r["V_th"]), float(r["I_d"]))
            for r in rows]


def _load_diode() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Diode CSV `index,V_sd,I_sd,T_j (C)` -> (V_sd, I_sd, T_c)."""
    rows = _rows("ISC040N10NM8_Diode.csv")
    v = np.array([float(r["V_sd"]) for r in rows])
    i = np.array([float(r["I_sd"]) for r in rows])
    t = np.array([float(r["T_j (C)"]) for r in rows])
    return v, i, t


@pytest.fixture(scope="module")
def ch() -> dict:
    """The 25 C set (132 pts) + TC rows, grouped like the fit script."""
    iv = load_csv(DATA / "ISC040N10NM8_IV_25C.csv")
    # The promoted record was fit on the 69 on-state / transfer IV rows
    # (V_gs > 0). The two V_gs = 0 off-state rows (80 V, 25/125 C) added for
    # the DIBL-cap experiment are NOT part of this record's fit data, so exclude
    # them here -- the pin stays against the data the record was actually fit
    # on. Promoting a cap-based record that uses them is a separate decision
    # (see scripts/experiment_dibl_cap.py).
    iv = iv.subset(iv.v_gs != 0.0)
    rid_vg, rid_vd, rid_id = _load_rdson_id("ISC040N10NM8_Rdson_Id_25C.csv")
    rv25 = [r for r in _load_rdson_vgs("ISC040N10NM8_Rdson_Vgs.csv")
            if abs(r[0] - 25.0) < 1e-6]
    v_ds = np.concatenate([iv.v_ds, rid_vd,
                           np.array([r[2] * I_TEST for r in rv25])])
    v_gs = np.concatenate([iv.v_gs, rid_vg, np.array([r[1] for r in rv25])])
    i_ref = np.concatenate([iv.i_d, rid_id, np.array([I_TEST] * len(rv25))])
    return {"v_ds": v_ds, "v_gs": v_gs, "i_ref": i_ref,
            "rd_all": _load_rdson_vgs("ISC040N10NM8_Rdson_Vgs.csv"),
            "vth": _load_vth("ISC040N10NM8_Vth_Tj.csv")}


def _chord(vg: float, t_c: float) -> float:
    """Model Rds(on) chord at I_TEST (brentq, as in the fit script)."""
    p = ISC040_CHANNEL.params_at(t_c)
    v = brentq(lambda x: float(vdmos_id(x, vg, p)) - I_TEST, 1e-5, 20.0)
    return v / I_TEST


def _model_vth(t_c: float, i_test: float) -> float:
    """Model gate voltage at V_ds = V_gs giving i_test (Vth spec reading)."""
    p = ISC040_CHANNEL.params_at(t_c)
    return brentq(lambda vg: float(vdmos_id(vg, vg, p)) - i_test, 0.5, 8.0)


# ---------------------------------------------------------------------------
# Channel: RMSE reproduction (record evaluation is deterministic -- no refit)
# ---------------------------------------------------------------------------


def test_record_reproduces_25c_rmse(ch):
    p25 = ISC040_CHANNEL.params_at(25.0)
    model = vdmos_id(ch["v_ds"], ch["v_gs"], p25)
    rel = (model - ch["i_ref"]) / np.maximum(np.abs(ch["i_ref"]), 1e-3)
    assert len(ch["v_ds"]) == 132, "25 C set must stay IV(69)+RdsonId(55)" \
                                   "+RdsonVgs25(8) like the fit script"
    rmse = float(np.sqrt(np.mean(rel**2)))
    assert rmse == pytest.approx(ISC040_CHANNEL_RMSE25_REL, abs=2e-3)
    assert rmse < 0.016, "official 25 C channel fit regressed"


def test_record_reproduces_tc_rmse(ch):
    """TC residuals exactly as fit_tc(): chord residual + Vth-VOLTAGE
    residual (solve the model gate V at V_ds=V_gs giving I_d, vs datasheet
    Vth -- fits the Vth(T) slope without the exponential-current sensitivity
    the single-softplus subthreshold cannot reproduce)."""
    r1 = [(_chord(vg, t) - R) / R for t, vg, R in ch["rd_all"]]
    r2 = [(_model_vth(t, it) - vth) / vth for t, vth, it in ch["vth"]]
    allr = np.array(r1 + r2)
    rmse = float(np.sqrt(np.mean(allr**2)))
    assert rmse == pytest.approx(ISC040_CHANNEL_RMSE_TC_REL, abs=2e-3)
    assert rmse < 0.032, "official TC fit regressed"
    # absolute ceilings per group (rdson carries the chord data, vth the
    # subthreshold-shape softness of the model form)
    assert float(np.sqrt(np.mean(np.array(r1) ** 2))) < 0.028
    assert float(np.sqrt(np.mean(np.array(r2) ** 2))) < 0.045


# ---------------------------------------------------------------------------
# Acceptance specs (never fitted) -- pinned so they cannot silently worsen
# ---------------------------------------------------------------------------


def test_ron_chord_specs_25c(ch):
    """Ron @ 25 A, 25 C: model 3.35/3.56/3.89 mOhm at Vgs 15/10/8 against
    the digitized data curve read at the same Vgs (3.35/3.59/3.93 mOhm by
    interpolation; the handover quotes the nearest digitized points
    3.33/3.52/3.87 at Vgs 15.74/10.55/8.76)."""
    assert _chord(15.0, 25.0) * 1e3 == pytest.approx(3.35, rel=0.015)
    assert _chord(10.0, 25.0) * 1e3 == pytest.approx(3.56, rel=0.015)
    assert _chord(8.0, 25.0) * 1e3 == pytest.approx(3.89, rel=0.015)
    pts = sorted((vg, R) for t, vg, R in ch["rd_all"] if abs(t - 25.0) < 1e-6)
    vg_d = np.array([p[0] for p in pts])
    r_d = np.array([p[1] for p in pts])
    for vg in (15.0, 10.0, 8.0):
        assert _chord(vg, 25.0) / float(np.interp(vg, vg_d, r_d)) - 1.0 \
            == pytest.approx(0.0, abs=0.02)


def test_ron_vs_T_specs(ch):
    """Ron @ Vgs=15.7 vs T: model 2.12/3.34/5.20/7.76 mOhm at -55/25/100/175
    (data 2.05/3.33/5.23/8.00 at the top digitized Vgs ~15.7 of each T).
    The -55 C +3.2 % and 175 C -3.0 % deviations against the data are the
    model-form softness of the quadratic r_d(T) law over the full span --
    pinned so they cannot silently worsen."""
    model = {t: _chord(15.7, t) * 1e3 for t in (-55.0, 25.0, 100.0, 175.0)}
    for t, ref in ((-55.0, 2.12), (25.0, 3.34), (100.0, 5.20), (175.0, 7.76)):
        assert model[t] == pytest.approx(ref, rel=0.02)
    data = {}
    for t in (-55.0, 25.0, 100.0, 175.0):
        top = max((row for row in ch["rd_all"] if abs(row[0] - t) < 1e-6),
                  key=lambda row: row[1])
        data[t] = top[2]
    for t in (-55.0, 25.0, 100.0, 175.0):
        assert model[t] / (data[t] * 1e3) - 1.0 == pytest.approx(0.0, abs=0.05)
    assert model[-55.0] < model[25.0] < model[100.0] < model[175.0]


def test_vth_spec_points(ch):
    """Vth @ 520 uA vs T: model 3.52/3.06/2.50/2.07 V at -52.4/22.5/107.5/
    172.9 C (data 3.37/3.05/2.58/2.07). s_vth -7.49 mV/K vs the datasheet
    Vth_Tj slope ~ -6.3 mV/K: within model form -- pin the resulting curve."""
    for t, ref in ((-52.4, 3.52), (22.5, 3.06), (107.5, 2.50), (172.9, 2.07)):
        assert _model_vth(t, 5.2e-4) == pytest.approx(ref, rel=0.015)
    slope = (_model_vth(172.9, 5.2e-4) - _model_vth(-52.4, 5.2e-4)) / 225.3
    assert slope == pytest.approx(-6.44e-3, abs=2e-4)


# ---------------------------------------------------------------------------
# TC-law physical invariants over the data span
# ---------------------------------------------------------------------------


def test_tc_laws_physical_over_span():
    for t in np.linspace(-55.0, 175.0, 25):
        p = ISC040_CHANNEL.params_at(float(t))
        assert p.r_d > 0.0        # quadratic never negative in span
        assert p.k_p > 0.0
        assert np.isfinite(p.v_th0)
    p25 = ISC040_CHANNEL.params_at(25.0)
    p175 = ISC040_CHANNEL.params_at(175.0)
    # r_d carries the whole Ron TC: x2.29 at 175 C, flattening at the cold
    # end (quadratic Trs1/Trs2 law, vertex ~-130 C -- NEVER linear)
    assert p175.r_d / p25.r_d == pytest.approx(2.293, abs=0.02)
    assert ISC040_CHANNEL.params_at(-55.0).r_d < p25.r_d * 0.75
    # mobility-like k_p falls with T (beta_p = 2.20 essential: without it the
    # Rdson RMSE goes to 13 %, the channel-dominated low-Vgs end breaks)
    assert p175.k_p < p25.k_p
    assert ISC040_CHANNEL.beta_p == pytest.approx(2.20, abs=0.05)
    # v_th0 slope -7.486 mV/K (datasheet Vth_Tj slope ~ -6.3, within form)
    slope = (p175.v_th0 - p25.v_th0) / 150.0
    assert slope == pytest.approx(-7.486e-3, abs=1e-4)
    # shape knobs T-invariant
    assert ISC040_CHANNEL.params_at(-55.0).v_sub == ISC040_CHANNEL.v_sub
    assert ISC040_CHANNEL.params_at(175.0).v_ovmax == ISC040_CHANNEL.v_ovmax


def test_model_finite_and_monotone_on_transfer():
    """No NaN/inf; i_d strictly increasing in v_gs at V_DS = 2 V; the channel
    B-source is forward-polarity only, so v_ds = 0 gives exactly zero."""
    vgs = np.linspace(6.0, 15.0, 40)
    for t in (-55.0, 25.0, 100.0, 175.0):
        p = ISC040_CHANNEL.params_at(t)
        i_d = np.asarray(vdmos_id(2.0, vgs, p))
        assert np.all(np.isfinite(i_d))
        assert np.all(np.diff(i_d) > 0.0), f"transfer went non-monotone at {t:g} C"
        assert float(vdmos_id(0.0, 12.0, p)) == 0.0


def test_record_sanity_and_documented_pins():
    """Physical sanity of the record. The documented structure: SJ cap form
    ACTIVE (v_ovmax > 0 -- it ceilings the Vgs=15 conductance; without it the
    Ron chord is ~13 % low), k_theta = 0, velocity-sat off (v_vsat = 1e4),
    no CLM, no kink, non-Kelvin (r_s = r_sk = 0, the data reject gate-loop
    degeneration). k_eta interior (NOT at the 0.05 cap like IMCQ120)."""
    r = ISC040_CHANNEL
    assert 0.0 < r.r_d25 < 3.4e-3, "r_d25 must stay below the 25 C datasheet chord"
    assert r.r_s == 0.0 and r.r_sk == 0.0
    assert r.k_theta == 0.0 and r.k_lambda == 0.0 and r.g_k == 0.0
    assert r.v_vsat == pytest.approx(1e4)
    assert r.v_ovmax > 0.0 and r.k_ovsharp > 0.0   # SJ cap sensed (active)
    assert 0.1 < r.v_sub < 1.0
    assert -0.03 < r.s_vth < 0.0
    assert 0.0 < r.beta_p < 8.0                    # interior (TC bounds 0..8)
    assert 0.0 < r.k_eta < 0.05                    # interior, unlike IMCQ pin
    assert 0.1 < r.v_dbsat < 5.0                   # inert (data stops ~5 V)
    assert 0.1 < r.k_kneesharp < 40.0
    assert 0.0 < r.a_rd < 5e-2 and abs(r.b_rd) < 1e-3


# ---------------------------------------------------------------------------
# Body diode record
# ---------------------------------------------------------------------------


def test_diode_record_reproduces():
    """The promoted DC record still reproduces the measured diode CSV.
    26 pts at 25/175 C, 1.18 A .. 490 A; overall RMSE 0.514 %."""
    v, i, t = _load_diode()
    model = vdmos_diode(v, ISC040_DIODE, t_c=t)
    rel = (model - i) / i
    rmse = float(np.sqrt(np.mean(rel**2)))
    assert rmse == pytest.approx(ISC040_DIODE_RMSE_REL, abs=2e-3)
    assert rmse < 0.006, "official diode DC fit regressed"
    for temp, ceil in ((25.0, 0.006), (175.0, 0.007)):
        msk = np.isclose(t, temp)
        r = float(np.sqrt(np.mean(rel[msk] ** 2)))
        assert r == pytest.approx(ISC040_DIODE_GROUP_RMSE_REL[int(temp)],
                                  abs=2e-3)
        assert r < ceil


def test_diode_record_sanity_and_documented_pins():
    """Interior values + the documented fixed structure. KEY FINDING pinned:
    k_n is ~1.05 (near-ideal diffusion), NOT the 2.0 that ST8L65 pins --
    k_n = 2.0 here costs 35 % RMSE. Recombination branch not needed (i_sr =
    0); quadratic b_rdiode / a_ikf stay 0 (2 temperatures only, degenerate:
    a_rdiode<->b_rdiode corr -1.0). The linear a_rdiode is the decisive TC
    (3.31 % -> 0.51 %). r_diode ~ 0.96 mOhm is ~10x SMALLER than the old
    recovery-bench R1 = 10 mOhm -- the subckt reconciles this by putting
    r_diode(T) in R1 with the Trs1 law and Rs = 0 on D1 (single count)."""
    d = ISC040_DIODE
    assert 1.001 < d.k_n < 1.2, "k_n must stay interior and near-ideal"
    assert 1.0 < d.v_eg < 1.2                      # silicon-like, well sensed
    assert d.i_sr == 0.0 and d.k_nr == 2.0         # recombination branch off
    assert d.a_ikf == 0.0 and d.b_rdiode == 0.0    # stay pinned (2 T only)
    assert d.t_nom == 25.0 and d.k_xti == 3.0
    assert 0.5e-3 < d.r_diode < 2e-3
    assert 0.0 < d.a_rdiode < 0.05
    assert d.i_kf > 100.0                          # high-injection knee sensed
    assert abs(float(vdmos_diode(-1.0, d))) < 1e-5  # reverse floor
    # the i_s<->k_n prefactor pair is degenerate (corr +0.999): pin the PAIR
    # through the forward curve, not i_s alone -- 1 V @25 C as a pair anchor
    assert float(vdmos_diode(1.0, d)) == pytest.approx(152.6, rel=0.03)


# ---------------------------------------------------------------------------
# Generated full subckt: the D1 branch must carry the DC-fit set (session-27
# reconciliation). out/ is gitignored -> skip on a fresh clone.
# ---------------------------------------------------------------------------


def test_generated_subckt_d1_branch_matches_record():
    sub = OUT / "ISC040N10NM8_full.subckt"
    if not sub.exists():
        pytest.skip("out/ not present (gitignored deliverable)")
    text = sub.read_text(encoding="utf-8")
    d = ISC040_DIODE
    xrec = [ln for ln in text.splitlines() if ln.startswith("Xrec")]
    assert len(xrec) == 1
    for token in (f"Is={d.i_s:.6g}", f"Rs={d.r_diode:.6g}", f"N={d.k_n:.6g}",
                  f"Ikf={d.i_kf:.6g}", f"Eg={d.v_eg:.6g}",
                  f"Xti={d.k_xti:g}", f"Tnom={d.t_nom:g}",
                  f"Trs1={d.a_rdiode:.6g}"):
        assert token in xrec[0], f"Xrec instance params missing {token!r}"
    assert "Is=1e-06" not in xrec[0] and "N=1.4" not in xrec[0], \
        "old bench diode set must not come back"
    assert ".model MyDiode D(Is=Is,Rs=0,N=N,Ikf=Ikf,Eg=Eg,Xti=Xti,Tnom=Tnom)" \
        in text, "D1 must carry the full DC set with Rs=0 (single series-R \
count: r_diode(T) lives in R1 with the Trs1 law)"
    assert "R1 K Ki {Rs*(1+Trs1*(temp-Tnom))}" in text, \
        "session 30: at r_max = 0 the branch is R1 tied straight to K"
    net = "\n".join(ln for ln in text.splitlines()
                    if not ln.startswith("*"))
    assert "B3" not in net, "B3 only exists with a non-zero r_max (session 30)"
    assert "N001" not in net, "orphan node N001 must go with B3"


def test_generated_subckt_mirror_and_fast_variant():
    """Session-29 pins. (a) The channel expression is MIRROR-SYMMETRIC
    third quadrant: abs(V(dx,sx)) inside, nested-if sign factor outside
    (LTspice has no sign() function) -- one flat expression, no buffered
    helper nodes (measured worse). (b) The FAST variant file exists with
    the SAME subckt name, carries NO recovery network, maps the pole-gauge
    storage onto the plain diode's tt = t_minor + t_pole/0.69, single-
    counts the series R in the model card (no R1), and puts the Cds
    charge-master at the D-S terminals."""
    full = OUT / "ISC040N10NM8_full.subckt"
    fast = OUT / "ISC040N10NM8_fast.subckt"
    if not (full.exists() and fast.exists()):
        pytest.skip("out/ not present (gitignored deliverable)")
    for path in (full, fast):
        bch = [ln for ln in path.read_text(encoding="utf-8").splitlines()
               if ln.startswith("Bch ")]
        assert len(bch) == 1, f"{path.name}: exactly one channel B-source"
        assert "if(V(dx,sx)<0,-1,if(V(dx,sx)>0,1,0))*" in bch[0], \
            f"{path.name}: mirror sign factor missing"
        assert "abs(V(dx,sx))" in bch[0], \
            f"{path.name}: mirror must evaluate the forward expr at |v_ds|"
        assert "nc_vda" not in bch[0], \
            f"{path.name}: buffered helper nodes are a measured regression"
    tx = fast.read_text(encoding="utf-8")
    net = "\n".join(ln for ln in tx.splitlines() if not ln.startswith("*"))
    assert "Xrec" not in net and "DiodeRecov" not in net, \
        "fast variant must carry NO recovery network (header comments aside)"
    assert "D1 s d MyDiode" in tx, "fast: body diode anode=source, cathode=drain"
    tt = ISC040_RECOV.t_minor + ISC040_RECOV.t_pole / 0.69
    card = [ln for ln in tx.splitlines() if ln.startswith(".model MyDiode")]
    assert len(card) == 1
    assert f"tt={tt * 1e9:.6g}n" in card[0], "fast: tt mapping"
    assert f"Rs={ISC040_DIODE.r_diode:.6g}" in card[0] and \
        f"Trs1={ISC040_DIODE.a_rdiode:.6g}" in card[0], \
        "fast: series R single-counted in the model card (Trs1 law)"
    assert "Rs=0," not in card[0], "fast: no Rs=0 split (that is the FULL gauge)"
    assert "Cds D S Q=" in tx, "fast: Cds charge-master at the D-S terminals"
    assert ".subckt ISC040N10NM8 G D S K" in tx, "fast: same drop-in subckt name"
