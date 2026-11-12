"""Phase B (session 30): fitting stages + presets on real ISC040 data.

These are QUALITY pins, not record pins: the official records are pinned
by test_isc040_official.py; refits from generic preset starts may land in
different flat-valley basins, so we assert convergence, RMSE bands, and
the TC stage reproducing the pinned record's acceptance quality.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from vdmos_fit.ingest import IVTable, load_table
from vdmos_fit.presets import PRESETS
from vdmos_fit.reference import ISC040_CHANNEL
from vdmos_fit.stages import (
    channel_rows,
    diode_rows,
    fit_channel_ref,
    fit_channel_tc,
    fit_diode_stage,
    resolve_test_current,
    third_quadrant_rows,
)

DATA = Path(__file__).resolve().parents[1] / "data"


@pytest.fixture(scope="module")
def isc040_tables():
    iv = load_table(DATA / "ISC040N10NM8_IV_25C.csv")
    rd_id = load_table(DATA / "ISC040N10NM8_Rdson_Id_25C.csv")
    rd_vgs = load_table(DATA / "ISC040N10NM8_Rdson_Vgs.csv")
    vth = load_table(DATA / "ISC040N10NM8_Vth_Tj.csv")
    return iv, rd_id, rd_vgs, vth


# ---- row-splitting doctrine ---------------------------------------------

def test_row_split_gate_off_vs_third_quadrant():
    iv = IVTable(vgs=np.array([0.0, 5.0, 5.0, 5.0]),
                 vds=np.array([-0.8, -0.5, 0.1, 0.1]),
                 ids=np.array([-3.0, -20.0, 10.0, 0.0]),
                 tj=np.array([25.0] * 4), gate_known=True, source="t")
    d = diode_rows(iv)
    assert d is not None and len(d.i_sd) == 1      # gate-off reverse only
    assert d.i_sd[0] == 3.0 and d.v_sd[0] == 0.8   # diode convention
    tq = third_quadrant_rows(iv)
    assert tq is not None and len(tq.i_d) == 1     # gate-on reverse
    ch = channel_rows(iv, 25.0)
    assert len(ch.i_d) == 1                        # I > 0 only


def test_resolve_test_current_sources(tmp_path):
    iv = load_table(DATA / "ISC040N10NM8_Rdson_Id_25C.csv")
    it = resolve_test_current(iv, None)
    assert (it > 0).all()
    # session 30 (user data edit): Rdson_Vgs now carries I_d = 25 A in a
    # column -- resolve succeeds straight from the column:
    iv2 = load_table(DATA / "ISC040N10NM8_Rdson_Vgs.csv")
    it2 = resolve_test_current(iv2, None)
    assert np.all(it2 == 25.0)
    assert np.all(resolve_test_current(iv2, 25.0) == 25.0)
    # the unknown-current raise path, pinned on a synthetic file (no I_d
    # column, no condition text):
    p = tmp_path / "nochord.csv"
    p.write_text("index,V_gs,R_dson (mOhm),T_j (C)\n1,10,12.5,25\n",
                 encoding="utf-8")
    spec = load_table(p)
    with pytest.raises(Exception, match="test current unknown"):
        resolve_test_current(spec, None)
    assert np.all(resolve_test_current(spec, 25.0) == 25.0)


# ---- stage 1: channel @ 25 C --------------------------------------------

def test_fit_channel_ref_isc040(isc040_tables):
    iv, rd_id, rd_vgs, _ = isc040_tables
    # chords for stage 1 come from the Rdson-Id file (has I column)
    res = fit_channel_ref(iv, rd_id, PRESETS["trench"], ref_temp=25.0,
                          multistart_limit=2, max_nfev=3000)
    assert res.success
    assert res.rmse_rel < 0.02
    assert 0.0 < res.params.r_d < 3.4e-3
    # the record's fitted k_p magnitude band (composite conductance)
    assert 100.0 < res.params.k_p < 800.0


# ---- stage 2: TC on the PINNED record (deterministic) -------------------

def test_fit_channel_tc_on_pinned_record(isc040_tables):
    _, rd_id, rd_vgs, vth = isc040_tables
    p25 = ISC040_CHANNEL.params_at(25.0)
    tcres = fit_channel_tc(p25, rd_vgs, vth, PRESETS["trench"],
                           rdson_i_test=25.0)
    # the TC stage on the pinned p25 must REPRODUCE the pinned TC record
    # (same data, same law, same bounds): bit-level agreement
    assert tcres["success"]
    for name in ("s_vth", "beta_p", "a_rd", "b_rd"):
        assert tcres[name] == pytest.approx(getattr(ISC040_CHANNEL, name),
                                            rel=1e-6), name
    assert tcres["rmse_rdson"] < 0.03
    assert tcres["rmse_vth"] < 0.05
    assert tcres["s_vth"] < 0.0          # Vth falls with temperature
    assert tcres["a_rd"] > 0.0           # r_d rises with temperature


# ---- stage 3: diode ------------------------------------------------------

def test_fit_diode_stage_isc040():
    iv = load_table(DATA / "ISC040N10NM8_Diode.csv")
    res = fit_diode_stage([iv], PRESETS["trench"], multistart_limit=2,
                          max_nfev=4000)
    assert res is not None and res.success
    assert res.rmse_rel < 0.01          # record: 0.51 %
    assert res.params.a_rdiode > 0.0    # linear TC promoted by the data
    assert res.params.k_n < 1.5         # near-ideal diffusion, not 2.0


def test_sic_preset_skips_diode_without_data():
    assert PRESETS["sic"].diode_free        # SiC diode support ENABLED
    iv = load_table(DATA / "IMCQ120R004M2H_IV_25C.csv")
    assert diode_rows(iv) is None       # all-positive IV: no diode rows
    assert fit_diode_stage([iv], PRESETS["sic"]) is None


def test_v_ovmax_zero_when_not_free():
    """User rule (session 30): the charge-cap term is OFF (v_ovmax = 0)
    whenever v_ovmax is NOT in the free set -- a locked preset start
    value must never leak into the export."""
    iv = load_table(DATA / "ISC040N10NM8_IV_25C.csv")
    rd = load_table(DATA / "ISC040N10NM8_Rdson_Vgs.csv")
    free = tuple(n for n in PRESETS["trench"].channel_free25
                 if n != "v_ovmax")
    res = fit_channel_ref(iv, rd, PRESETS["trench"], free_override=free,
                          multistart_limit=1, max_nfev=2000)
    assert res.params.v_ovmax == 0.0


def test_fit_channel_bounds_override_pins():
    """GUI min/max columns reach the fitter: a tight v_dbsat ceiling
    pins the fit at its upper bound (the optimum sits well above it).
    A k_p ceiling was tried first and starved the package-R solve at
    the 25 A chords -- extreme boxes can make the model unsolvable."""
    iv = load_table(DATA / "ISC040N10NM8_IV_25C.csv")
    rd = load_table(DATA / "ISC040N10NM8_Rdson_Vgs.csv")
    res = fit_channel_ref(iv, rd, PRESETS["trench"],
                          bounds_override={"v_dbsat": (0.1, 0.3)},
                          multistart_limit=1, max_nfev=1500)
    assert res.params.v_dbsat <= 0.3 * 1.001
    assert res.bounds_used["v_dbsat"] == (0.1, 0.3)


def test_fit_channel_ref_no_chords_at_ref_temp():
    """Regression (GUI field report): a digitized Rds(on)-vs-T curve with
    NO row exactly at ref temp (IMCQ: closest 16.0/32.1 C) must not crash
    the r_d upper bound -- the bound is simply skipped."""
    iv = load_table(DATA / "IMCQ120R004M2H_IV_25C.csv")
    rd = load_table(DATA / "IMCQ120R004M2H_Rdson_Tj_185.2A.csv")
    assert not np.any(np.isclose(rd.tj, 25.0))
    res = fit_channel_ref(iv, rd, PRESETS["sic"], multistart_limit=2,
                          max_nfev=4000)
    assert np.isfinite(res.params.r_d) and res.params.r_d > 0.0
    assert res.rmse_rel < 0.20          # regression target: no crash, sane fit


def test_fit_channel_kelvin_ratio_constraint():
    """Regression (GUI field report): freeing r_s AND r_sk must not crash
    on the model guard r_sk <= r_s -- fit() optimises the ratio
    f_sk = r_sk/r_s in [0, 1] (session-22 kelvin doctrine)."""
    iv = load_table(DATA / "IMCQ120R004M2H_IV_25C.csv")
    rd = load_table(DATA / "IMCQ120R004M2H_Rdson_Tj_185.2A.csv")
    free = tuple(PRESETS["sic"].channel_free25) + ("r_s", "r_sk")
    res = fit_channel_ref(iv, rd, PRESETS["sic"], free_override=free,
                          multistart_limit=1, max_nfev=3000)
    assert 0.0 <= res.params.r_sk <= res.params.r_s


def test_start_override_carries_nonpreset_knob():
    """Regression (GUI field report): the channel fit must honour the live
    values as the first start, so a knob the preset never knew about --
    k_kneeform -- survives the fit instead of reverting to the Params()
    default (which silently dropped the GUI's scale-form selection and pinned
    k_theta at 0). Without start_override the official preset default is used."""
    iv = load_table(DATA / "ISC040N10NM8_IV_25C.csv")
    rd = load_table(DATA / "ISC040N10NM8_Rdson_Vgs.csv")
    p = PRESETS["trench"]
    so = {**p.channel_fixed, **p.channel_start, "k_kneeform": 1.0,
          "k_theta": 1.0, "v_vsat": 5.0}
    free = tuple(p.channel_free25) + ("k_theta", "v_vsat")
    res = fit_channel_ref(iv, rd, p, free_override=free, start_override=so,
                          multistart_limit=1, max_nfev=1500)
    assert res.params.k_kneeform == 1.0        # survived, not reset to 0
    res0 = fit_channel_ref(iv, rd, p, free_override=free,
                           multistart_limit=1, max_nfev=1500)
    assert res0.params.k_kneeform == 0.0       # preset default (mobility form)
