"""Session 30 addendum (user rule 2026-11-09): the GUI "Save subckt"
button asks for a destination FOLDER and writes ONLY the device
subckts -- FAST always, FULL when recovery knobs exist. The caps
subckt and the cds expression text file are NOT standalone artifacts:
the caps delta and the Cds charge-master are embedded inline in every
device subckt anyway. Pinned through export_part(caps_files=False,
both_variants=True); the CLI default flags keep the old artifact set.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from vdmos_fit import build_cap_model, load_cap_csv
from vdmos_fit.export_subckt import RecovParams
from vdmos_fit.params import Params
from vdmos_fit.release import PLACEHOLDER_DIODE, export_part

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def caps_model():
    return build_cap_model(
        load_cap_csv(ROOT / "data" / "ISC040N10NM8_Cap_Vds.csv"))


def test_gui_save_fast_only_no_artifacts(tmp_path, caps_model):
    # no recovery -> FAST only, and NOTHING but the device subckt
    res, files = export_part(tmp_path, "PFAST", Params(), PLACEHOLDER_DIODE,
                             None, caps_model, n_sub=4, caps_files=False,
                             both_variants=True)
    assert sorted(Path(f).name for f in files) == ["PFAST_fast.subckt"]
    assert not (tmp_path / "PFAST_caps.subckt").exists()
    assert not any(f.endswith("cds_expr.txt") for f in files)
    assert any("FAST export" in w for w in res.warnings)
    text = Path(files[0]).read_text(encoding="utf-8")
    assert ".subckt PFAST G D S K" in text
    # caps delta + Cds charge-master embedded inline -> self-contained
    assert ".subckt pfast_caps G D S" in text
    assert "Cds D S Q=" in text


def test_gui_save_writes_both_variants(tmp_path, caps_model):
    rec = RecovParams(tau=1e-9, t_minor=5e-9, t_pole=12e-9)
    res, files = export_part(tmp_path, "PBOTH", Params(), PLACEHOLDER_DIODE,
                             rec, caps_model, n_sub=4, caps_files=False,
                             both_variants=True)
    assert sorted(Path(f).name for f in files) == ["PBOTH_fast.subckt",
                                                   "PBOTH_full.subckt"]
    assert not (tmp_path / "PBOTH_caps.subckt").exists()
    assert not any("FAST export" in w for w in res.warnings)
    full = (tmp_path / "PBOTH_full.subckt").read_text(encoding="utf-8")
    fast = (tmp_path / "PBOTH_fast.subckt").read_text(encoding="utf-8")
    assert "Xrec" in full and ".subckt DiodeRecov" in full
    net = "\n".join(ln for ln in fast.splitlines()
                    if not ln.startswith("*"))
    assert "Xrec" not in net and "DiodeRecov" not in net, \
        "the FAST twin carries NO recovery network (header comments aside)"
    # FAST twin maps the pole-gauge storage onto tt (gen_full_subckt rule)
    tt = rec.t_minor + rec.t_pole / 0.69
    assert f"tt={tt * 1e9:.6g}n" in fast
    # same subckt NAME in both files: drop-in swap, never both at once
    assert ".subckt PBOTH G D S K" in full and \
        ".subckt PBOTH G D S K" in fast


def test_cli_default_artifacts_unchanged(tmp_path, caps_model, monkeypatch):
    # default flags keep the old export set (caps file + rule variant);
    # the cds expression txt lands under CWD/data, so run the test there
    monkeypatch.chdir(tmp_path)
    (tmp_path / "data").mkdir()  # export_part writes cds txt to CWD/data
    res, files = export_part(tmp_path, "PCLI", Params(), PLACEHOLDER_DIODE,
                             None, caps_model, n_sub=4)
    assert sorted(Path(f).name for f in files) == ["PCLI_caps.subckt",
                                                   "PCLI_fast.subckt"]
    assert (tmp_path / "PCLI_caps.subckt").exists()
    assert (tmp_path / "data" / "PCLI_cds_expr.txt").exists()
    assert ",tt=" not in (tmp_path / "PCLI_fast.subckt").read_text(
        encoding="utf-8"), "CLI rule: no recovery -> FAST, no tt"
