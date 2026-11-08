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


def test_resolve_test_current_sources():
    iv = load_table(DATA / "ISC040N10NM8_Rdson_Id_25C.csv")
    it = resolve_test_current(iv, None)
    assert (it > 0).all()
    iv2 = load_table(DATA / "ISC040N10NM8_Rdson_Vgs.csv")
    with pytest.raises(Exception, match="test current unknown"):
        resolve_test_current(iv2, None)
    it2 = resolve_test_current(iv2, 25.0)
    assert np.all(it2 == 25.0)


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
    iv = load_table(DATA / "IMCQ120R004M2H_IV_25C.csv")
    assert diode_rows(iv) is None       # all-positive IV: no diode rows
    assert fit_diode_stage([iv], PRESETS["sic"]) is None
