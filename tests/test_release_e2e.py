"""Phase D (session 30): end-to-end release golden test.

The ISC040 CSVs -> orchestrator -> FAST subckt (no recovery CSV exists for
ISC040 in the repo, so per the export rule the fast variant with a plain
tt-less diode must come out). Slow (real refits, ~2 min): quality bands,
not record equality -- the records themselves are pinned elsewhere.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from vdmos_fit.release import build_part, recov_spec_from_table

DATA = Path(__file__).resolve().parents[1] / "data"
OUT = Path(__file__).resolve().parents[1] / "out"

CSVS = [DATA / n for n in (
    "ISC040N10NM8_IV_25C.csv",
    "ISC040N10NM8_Diode.csv",
    "ISC040N10NM8_Rdson_Id_25C.csv",
    "ISC040N10NM8_Rdson_Vgs.csv",
    "ISC040N10NM8_Vth_Tj.csv",
    "ISC040N10NM8_Cap_Vds.csv",
    "ISC040N10NM8_Vgs_Qg.csv",
)]


@pytest.fixture(scope="module")
def release(tmp_path_factory):
    out = tmp_path_factory.mktemp("release_out")
    res = build_part(CSVS, "ISC040E2E", preset="trench", ref_temp=25.0,
                     rdson_i_test=25.0, out_dir=out, multistart_limit=2,
                     max_nfev=3000)
    yield res
    stale = Path("data") / "ISC040E2E_cds_expr.txt"
    if stale.exists():
        stale.unlink()


def test_record_quality(release):
    assert release.rmse["channel_ref"] < 0.02
    assert release.rmse["diode"] < 0.01
    assert release.rmse["tc"] < 0.05
    assert release.tc["s_vth"] < 0 and release.tc["a_rd"] > 0
    assert release.diode["k_n"] < 1.5
    assert any("FAST export" in w for w in release.warnings)


def test_fast_subckt_tokens(release, tmp_path_factory):
    # the subckt path is recorded on the result
    sub = Path(next(f for f in release.files
                    if f.endswith("_fast.subckt")))
    text = sub.read_text(encoding="utf-8")
    net = "\n".join(ln for ln in text.splitlines()
                    if not ln.startswith("*"))
    assert ".subckt ISC040E2E G D S K" in net
    assert "Bch dx sx I=if(V(dx,sx)<0,-1,if(V(dx,sx)>0,1,0))*" in net
    assert "D1 s d MyDiode" in net
    assert ",tt=" not in net, "no recovery data -> NO tt (user rule)"
    assert "Xrec" not in net and "DiodeRecov" not in net
    assert "Cds D S Q=" in net
    assert "Xcaps G D S isc040e2e_caps" in net
    assert "B3" not in net


def test_record_json_roundtrip(release):
    rec = Path(next(f for f in release.files
                    if f.endswith("_record.json")))
    data = json.loads(rec.read_text(encoding="utf-8"))
    assert data["part"] == "ISC040E2E"
    assert data["channel"]["k_p"] == pytest.approx(release.channel["k_p"])
    assert "scipy" in data["meta"]


def test_recov_spec_uses_csv_currents():
    """Session 30: the recovery fit takes its forward current from the
    CSV's I_sd column (per row), never the RecovSpec 12.5 A hardcode."""
    import numpy as np

    from vdmos_fit.ingest import RecTable
    from vdmos_fit.release import PLACEHOLDER_DIODE
    rt = RecTable(qrr=np.array([247e-9, 36e-9]),
                  trr=np.array([25e-9, 29e-9]),
                  irm=np.array([np.nan, np.nan]),
                  i_f=np.array([12.5, -25.0]),   # I_sd sign-normalized
                  didt=np.array([1000e6, 100e6]),
                  source="x")
    spec = recov_spec_from_table(rt, "P", PLACEHOLDER_DIODE)
    assert spec.i_f == 12.5
    assert spec.i_f_points == (12.5, 25.0)
    for got, want in zip(spec.points, ((1000.0, 247.0, 25.0),
                                       (100.0, 36.0, 29.0)), strict=True):
        for g, w in zip(got, want, strict=True):
            assert g == pytest.approx(w)   # A/us, nC, ns rescaling
    rt2 = RecTable(qrr=np.array([247e-9, 36e-9]),
                   trr=np.array([25e-9, 29e-9]),
                   irm=np.array([np.nan, np.nan]),
                   i_f=np.array([12.5, 12.5]),
                   didt=np.array([1000e6, 100e6]),
                   source="x")
    assert recov_spec_from_table(rt2, "P", PLACEHOLDER_DIODE).i_f_points == ()
