"""Phase A (session 30): the input-contract parser (ingest.py).

Grammar cases from the agreed contract, unit prefixes, reverse-polarity
negation, per-file defaults, typed classification, and a regression pass
over every real data CSV in the repo.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from vdmos_fit.ingest import (
    CapTable,
    IngestError,
    IVTable,
    QgTable,
    RdsSpec,
    RecTable,
    VthSpec,
    load_table,
    load_tables,
    parse_condition,
)

DATA = Path(__file__).resolve().parents[1] / "data"


def _csv(tmp_path: Path, name: str, header: list[str],
         rows: list[list[str]]) -> Path:
    body = "\n".join([",".join(header),
                      *(",".join(r) for r in rows)]) + "\n"
    p = tmp_path / name
    p.write_text(body, encoding="utf-8")
    return p


def _iv(tmp_path: Path, header, rows, ref_temp=25.0) -> IVTable:
    t = load_table(_csv(tmp_path, "x.csv", header, rows), ref_temp=ref_temp)
    assert isinstance(t, IVTable)
    return t


# ---- header grammar ----------------------------------------------------

@pytest.mark.parametrize("hdr", ["V_ds", "V_ds[V]", "V_ds [V]", "Vds"])
def test_voltage_header_spellings_ok(tmp_path, hdr):
    t = _iv(tmp_path, [hdr, "I_d"], [["1.5", "2"]])
    assert t.vds == pytest.approx([1.5])
    assert t.ids == pytest.approx([2.0])


@pytest.mark.parametrize("hdr", ["V_ds[kV]", "V_ds[KV]", "V_ds [kV]"])
def test_kilo_both_cases(tmp_path, hdr):
    t = _iv(tmp_path, [hdr, "I_d"], [["0.0015", "2"]])
    assert t.vds == pytest.approx([1.5])


@pytest.mark.parametrize("hdr", ["I_d[uA]", "I_d[µA]", "I_d[μA]"])
def test_micro_three_spellings(tmp_path, hdr):
    t = _iv(tmp_path, ["V_ds", hdr], [["1", "2"]])
    assert t.ids == pytest.approx([2e-6])


def test_mega_vs_milli(tmp_path):
    t = _iv(tmp_path, ["I_d[mA]"], [["2"]])
    assert t.ids == pytest.approx([2e-3])
    t = _iv(tmp_path, ["I_d[MA]"], [["2"]])
    assert t.ids == pytest.approx([2e6])


def test_bad_prefix_rejected(tmp_path):
    with pytest.raises(IngestError):
        _iv(tmp_path, ["I_d[Ama]"], [["1"]])


# ---- reverse polarity --------------------------------------------------

def test_isd_negated(tmp_path):
    t = _iv(tmp_path, ["I_sd", "V_sd"], [["5.0", "0.8"]])
    assert t.ids == pytest.approx([-5.0])
    assert t.vds == pytest.approx([-0.8])
    assert any("negated" in n for n in t.notes)


def test_vsg_negated(tmp_path):
    t = _iv(tmp_path, ["V_sg", "I_d", "V_ds"], [["-5", "1", "0.1"]])
    assert t.vgs == pytest.approx([5.0])


# ---- temperature -------------------------------------------------------

def test_temp_bracket_styles_and_kelvin(tmp_path):
    for hdr in ("T_j (C)", "Tj (C)", "T_j[C]", "T_j [C]"):
        t = _iv(tmp_path, ["V_ds", "I_d", hdr], [["1", "1", "100"]])
        assert t.tj == pytest.approx([100.0]), hdr
    t = _iv(tmp_path, ["V_ds", "I_d", "T_j[K]"], [["1", "1", "298.15"]])
    assert t.tj == pytest.approx([25.0])
    t = _iv(tmp_path, ["V_ds", "I_d", "T"], [["1", "1", "25"]])
    assert t.tj == pytest.approx([25.0])  # bare T = Celsius, never Kelvin
    assert any("Celsius" in n for n in t.notes)


def test_temp_prefix_rejected(tmp_path):
    with pytest.raises(IngestError):
        _iv(tmp_path, ["V_ds", "I_d", "T_j[mC]"], [["1", "1", "1"]])


def test_temp_default_per_file(tmp_path):
    t = _iv(tmp_path, ["V_ds", "I_d"], [["1", "1"], ["2", "2"]],
            ref_temp=27.0)
    assert np.all(t.tj == 27.0)
    assert any("defaulted" in n for n in t.notes)


# ---- defaults ----------------------------------------------------------

def test_vds_defaults_to_vgs(tmp_path):
    t = _iv(tmp_path, ["V_gs", "I_d"], [["3.3", "0.1"]])
    assert t.vds == pytest.approx([3.3])
    assert any("V_ds = V_gs" in n for n in t.notes)


def test_no_gate_column_means_gate_off(tmp_path):
    t = _iv(tmp_path, ["V_ds", "I_d"], [["0.8", "-3"]])
    assert not t.gate_known
    assert np.all(t.vgs == 0.0)


# ---- classification ----------------------------------------------------

def test_cap_file_split_by_type(tmp_path):
    t = load_table(_csv(tmp_path, "c.csv",
                        ["index", "V_ds", "Cap (pF)", "Type"],
                        [["1", "0.5", "3500", "iss"],
                         ["2", "0.5", "900", "rss"]]))
    assert isinstance(t, list) or hasattr(t, "tables")
    tabs = t.tables if hasattr(t, "tables") else [t]
    kinds = sorted(x.kind for x in tabs)
    assert kinds == ["iss", "rss"]
    assert tabs[0].c[0] == pytest.approx(3500e-12)


def test_qg_rds_vth_classification(tmp_path):
    qg = load_table(_csv(tmp_path, "q.csv", ["Q_g (nC)", "V_gs", "condition"],
                         [["0", "0", "I_d=12.5A, V_ds=20V"]]))
    assert isinstance(qg, QgTable) and qg.qg[0] == pytest.approx(0.0)
    rds = load_table(_csv(tmp_path, "r.csv",
                          ["V_gs", "R_dson (mOhm)", "T_j (C)"],
                          [["5.5", "5.3", "-55"]]))
    assert isinstance(rds, RdsSpec)
    assert rds.rds[0] == pytest.approx(5.3e-3)
    assert rds.tj[0] == pytest.approx(-55.0)
    vth = load_table(_csv(tmp_path, "v.csv", ["T_j (C)", "V_th", "I_d"],
                          [["-52.4", "3.37", "5.2E-04"]]))
    assert isinstance(vth, VthSpec)
    assert vth.i_test[0] == pytest.approx(5.2e-4)


def test_rec_table_and_didt_units(tmp_path):
    hdr = ["Q_rr[nC]", "t_rr[ns]", "I_RM[A]", "I_sd[A]", "di/dt[MA/s]"]
    t = load_table(_csv(tmp_path, "rec.csv", hdr,
                        [["100", "50", "8", "10", "100"]]))
    assert isinstance(t, RecTable)
    assert t.qrr[0] == pytest.approx(100e-9)
    assert t.trr[0] == pytest.approx(50e-9)
    assert t.didt[0] == pytest.approx(100e6)
    # [A/us] is numerically the same as [MA/s]
    t2 = load_table(_csv(tmp_path, "rec2.csv",
                         ["Q_rr[nC]", "t_rr[ns]", "I_RM[A]", "I_sd[A]",
                          "di/dt[A/us]"],
                         [["100", "50", "8", "10", "100"]]))
    assert t2.didt[0] == pytest.approx(t.didt[0])


def test_rec_missing_column_clear_error(tmp_path):
    with pytest.raises(IngestError, match="missing column"):
        load_table(_csv(tmp_path, "rec.csv", ["Q_rr[nC]", "t_rr[ns]"],
                        [["1", "2"]]))


def test_unknown_column_clear_error(tmp_path):
    with pytest.raises(IngestError, match="unknown column"):
        _iv(tmp_path, ["V_ds", "I_d", "banana"], [["1", "1", "1"]])


def test_bad_number_reports_row(tmp_path):
    with pytest.raises(IngestError, match="row 2"):
        _iv(tmp_path, ["V_ds", "I_d"], [["1", "N.A."]])


# ---- condition strings -------------------------------------------------

def test_parse_condition():
    assert parse_condition("V_ds=20V,T_j=25C") == {"vds": 20.0, "tj": 25.0}
    c = parse_condition("I_d=12.5A, V_ds=20V")
    assert c["ids"] == pytest.approx(12.5) and c["vds"] == 20.0


# ---- real repo CSVs (regression) ---------------------------------------

def test_all_repo_data_csvs_parse():
    paths = sorted(DATA.glob("*.csv"))
    assert len(paths) >= 18
    tables, report = load_tables(paths)
    kinds = {type(t).__name__ for t in tables}
    assert {"IVTable", "CapTable", "QgTable", "RdsSpec", "VthSpec"} <= kinds
    # the ST8L65 diode file: reverse names negated -> diode-forward rows land
    # as I_ds < 0, V_ds < 0 with gate off:
    di = load_table(DATA / "ST8L65N044M9_diode.csv")
    assert isinstance(di, IVTable) and not di.gate_known
    assert np.all(di.ids < 0) and np.all(di.vds < 0)
    # ISC040 diode file: T_j (C) header parses (the old loader gotcha)
    di2 = load_table(DATA / "ISC040N10NM8_Diode.csv")
    assert set(di2.tj) == {25.0, 175.0} or set(di2.tj) == {25.0}
    # caps: pF scaled to farads, Type column split
    caps = [t for t in tables if isinstance(t, CapTable)]
    assert caps and all(1e-12 < t.c.mean() < 1e-5 for t in caps)
