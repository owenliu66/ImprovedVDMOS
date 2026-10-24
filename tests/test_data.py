"""Tests for the CSV I/O contract in data.py (HANDOVER queue item 6, CLOSED).

load_csv() accepts BOTH on-disk header spellings, matched per column,
case-insensitively, order-free, extra columns (the digitizer `index`) ignored:

    legacy SPICE-style   vgs,vds,ids            (what save_csv() writes)
    digitizer/model      index,V_ds,I_d,V_gs    (what curve_digitizer exports)

Before this, the IMCQ digitizer exports (index,V_ds,I_d,V_gs) broke load_csv()
and both IMCQ scripts carried a local load_imcq() workaround -- now gone.
save_csv() still WRITES the legacy spelling (external contract with the
already-digitised datasheet files); the round trip must be bit-identical.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from vdmos_fit import CurveSet, load_csv, save_csv

DATA = Path(__file__).resolve().parents[1] / "data"

_ROWS = [(4.0, 0.05, 0.812), (4.0, 1.0, 15.3), (10.0, 20.0, 240.2)]


def _write(tmp_path: Path, name: str, header: str, rows=_ROWS) -> Path:
    """Write `rows` (v_gs, v_ds, i_d) under `header`, per-column value mapping."""
    path = tmp_path / name
    cols = [h.strip().lower() for h in header.split(",")]
    lines = []
    for n, (g, d, i) in enumerate(rows, start=1):
        vals = {"vgs": g, "v_gs": g, "vds": d, "v_ds": d, "ids": i, "i_d": i, "index": n}
        lines.append(",".join(str(vals[c]) for c in cols))
    path.write_text(f"{header}\n" + "\n".join(lines) + "\n", encoding="utf-8")
    return path


# ---------------------------------------------------------------- spellings
def test_load_csv_legacy_spelling(tmp_path):
    data = load_csv(_write(tmp_path, "legacy.csv", "vgs,vds,ids"))
    assert len(data) == 3
    assert np.array_equal(data.v_gs, [4.0, 4.0, 10.0])
    assert np.array_equal(data.v_ds, [0.05, 1.0, 20.0])
    assert np.array_equal(data.i_d, [0.812, 15.3, 240.2])


def test_load_csv_digitizer_spelling(tmp_path):
    """The curve_digitizer export spelling that used to break load_csv()."""
    data = load_csv(_write(tmp_path, "digitizer.csv", "index,V_gs,V_ds,I_d"))
    assert np.array_equal(data.v_gs, [4.0, 4.0, 10.0])
    assert np.array_equal(data.v_ds, [0.05, 1.0, 20.0])
    assert np.array_equal(data.i_d, [0.812, 15.3, 240.2])


def test_load_csv_both_spellings_identical(tmp_path):
    """Same numbers under either header -> bit-identical CurveSet arrays."""
    rows = [(3.7, 3.7, 2.5e-4), (10.0, 19.7, 240.2)]
    a = load_csv(_write(tmp_path, "a.csv", "vgs,vds,ids", rows))
    b = load_csv(_write(tmp_path, "b.csv", "index,V_ds,I_d,V_gs", rows))
    assert np.array_equal(a.v_gs, b.v_gs)
    assert np.array_equal(a.v_ds, b.v_ds)
    assert np.array_equal(a.i_d, b.i_d)


def test_load_csv_mixed_spelling_and_case(tmp_path):
    """Per-column alias matching: spellings may be mixed, case and spacing free."""
    data = load_csv(_write(tmp_path, "mixed.csv", " Index , Vds , i_d , VGS "))
    assert np.array_equal(data.v_gs, [4.0, 4.0, 10.0])
    assert np.array_equal(data.v_ds, [0.05, 1.0, 20.0])
    assert np.array_equal(data.i_d, [0.812, 15.3, 240.2])


def test_load_csv_error_lists_both_spellings(tmp_path):
    path = tmp_path / "bad.csv"
    path.write_text("index,x,y,condition\n1,0.5,10.0,foo\n", encoding="utf-8")
    with pytest.raises(ValueError, match=r"missing column\(s\)"):
        load_csv(path)
    with pytest.raises(ValueError, match=r"vgs/vds/ids or v_gs/v_ds/i_d"):
        load_csv(path)


# ---------------------------------------------------------------- round trip
def test_save_load_roundtrip_bit_identical(tmp_path):
    """save_csv writes the legacy spelling; the round trip is bit-identical
    (Python float str/repr round-trips exactly, so no tolerance needed)."""
    rng = np.random.default_rng(7)
    data = CurveSet(
        rng.uniform(2.0, 14.0, 50), rng.uniform(0.01, 40.0, 50),
        rng.uniform(1e-3, 300.0, 50), label="roundtrip",
    )
    path = tmp_path / "rt.csv"
    save_csv(data, path)
    assert path.read_text(encoding="utf-8").splitlines()[0] == "vgs,vds,ids"
    back = load_csv(path)
    assert np.array_equal(back.v_gs, data.v_gs)
    assert np.array_equal(back.v_ds, data.v_ds)
    assert np.array_equal(back.i_d, data.i_d)
    assert back.label == "rt"


# ---------------------------------------------------------------- real files
def test_real_imcq_digitizer_csvs_load_natively():
    """The IMCQ IV CSVs (header index,V_ds,I_d,V_gs) load through load_csv()
    with the same shape the old local load_imcq() produced: 78 / 73 points,
    gate grids up to v_gs = 20."""
    d25 = load_csv(DATA / "IMCQ120R004M2H_IV_25C.csv")
    d175 = load_csv(DATA / "IMCQ120R004M2H_IV_175C.csv")
    assert len(d25) == 78
    assert len(d175) == 73
    assert d25.v_gs_grid()[-1] == 20.0
    assert d175.v_gs_grid()[-1] == 20.0
    assert np.all(d25.i_d > 0.0)
    assert d25.label == "IMCQ120R004M2H_IV_25C"


def test_real_st8l65_csv_still_loads():
    """The legacy-spelling official dataset is untouched by the alias change."""
    data = load_csv(DATA / "ST8L65N044M9_IV.csv")
    assert len(data) == 41
    assert data.v_gs_grid() == [3.7, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
