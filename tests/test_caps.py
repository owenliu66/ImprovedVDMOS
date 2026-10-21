"""Capacitance/charge model tests (CAP_SPEC.md).

Covers: CSV contracts and unit conversion, the delta decomposition, the
charge master (dq/dv == C identity), the c_gd splice (C0/C1 at v_gd = 0),
PCHIP monotonicity through the super-junction cliff, the clamp policy, the
trench graded-junction analytic anchor, the pinned c_gdmax extraction
records (drift checks in the diode-record style), and the LTspice export
(charge-table ramp-sum identity, subckt structure, and -- when LTspice is
installed -- a live parity check of the exported nonlinear caps against
the Python model under PWL ramp drive).
"""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import curve_fit

from vdmos_fit import (
    CAP_ANCHOR_IMCQ120,
    CAP_ANCHOR_ISC040,
    CAP_ANCHOR_ST8L65,
    CAP_IMCQ120,
    CAP_ISC040,
    CAP_PLATEAU_DEFICIT_REL,
    CAP_ST8L65,
    CapModel,
    CapParams,
    QgCurveSet,
    build_cap_model,
    cap_charge_tables,
    export_cap_subckt,
    extract_cmax,
    load_cap_csv,
    load_qg_csv,
    save_cap_csv,
    save_qg_csv,
)
from vdmos_fit.caps import ANCHOR_N, _BranchTable
from vdmos_fit.ltspice import _rampsum_expr

ROOT = Path(__file__).resolve().parents[1]
PARTS = {
    "ST8L65N044M9": {"cap_pts": 46, "qg_pts": 4, "record": CAP_ST8L65,
                     "anchor": CAP_ANCHOR_ST8L65, "r_ds_on": 0.036},
    "IMCQ120R004M2H": {"cap_pts": 36, "qg_pts": 4, "record": CAP_IMCQ120,
                       "anchor": CAP_ANCHOR_IMCQ120, "r_ds_on": 0.004},
    "ISC040N10NM8": {"cap_pts": 28, "qg_pts": 10, "record": CAP_ISC040,
                     "anchor": CAP_ANCHOR_ISC040, "r_ds_on": 0.040},
}


def _model(name: str) -> CapModel:
    return build_cap_model(load_cap_csv(ROOT / "data" / f"{name}_Cap_Vds.csv"))


# ---------------------------------------------------------------------------
# CSV contracts
# ---------------------------------------------------------------------------

def test_loaders_units_and_counts():
    for name, cfg in PARTS.items():
        caps = load_cap_csv(ROOT / "data" / f"{name}_Cap_Vds.csv")
        assert len(caps) == cfg["cap_pts"]
        assert caps.branches() == ["iss", "oss", "rss"]
        # pF -> F conversion landed in a sane range
        assert 1e-14 < caps.c.min() and caps.c.max() < 1e-6
        qg = load_qg_csv(ROOT / "data" / f"{name}_Vgs_Qg.csv")
        assert len(qg) == cfg["qg_pts"]
        assert 1e-12 < qg.q_g.max() < 1e-6  # nC -> C
    # condition parsing: ISC040 carries three test biases
    qg = load_qg_csv(ROOT / "data" / "ISC040N10NM8_Vgs_Qg.csv")
    assert qg.conditions() == [(12.5, 20.0), (12.5, 50.0), (12.5, 80.0)]


def test_cap_qg_save_roundtrip(tmp_path):
    caps = load_cap_csv(ROOT / "data" / "ISC040N10NM8_Cap_Vds.csv")
    save_cap_csv(caps, tmp_path / "c.csv")
    back = load_cap_csv(tmp_path / "c.csv")
    assert np.allclose(back.v_ds, caps.v_ds)
    assert np.allclose(back.c, caps.c, rtol=1e-9)
    assert list(back.branch) == list(caps.branch)
    qg = load_qg_csv(ROOT / "data" / "ISC040N10NM8_Vgs_Qg.csv")
    save_qg_csv(qg, tmp_path / "q.csv")
    back = load_qg_csv(tmp_path / "q.csv")
    assert np.allclose(back.q_g, qg.q_g, rtol=1e-9)
    assert np.allclose(back.v_gs, qg.v_gs)
    assert back.conditions() == qg.conditions()


def test_cap_branch_validation():
    from vdmos_fit.data import CapCurveSet
    with pytest.raises(ValueError, match="unknown capacitance branch"):
        CapCurveSet(np.array([1.0]), np.array([1e-12]), np.array(["mig"]))
    caps = load_cap_csv(ROOT / "data" / "ISC040N10NM8_Cap_Vds.csv")
    only_iss = CapCurveSet(caps.v_ds[caps.branch == "iss"],
                           caps.c[caps.branch == "iss"],
                           caps.branch[caps.branch == "iss"])
    with pytest.raises(ValueError, match="missing branch"):
        from vdmos_fit.caps import CapTables
        CapTables.from_curve_set(only_iss)


# ---------------------------------------------------------------------------
# Decomposition + tables
# ---------------------------------------------------------------------------

def test_decomposition_positive_everywhere():
    for name in PARTS:
        m = _model(name)
        v = np.geomspace(0.05, m.tables.rss.v[-1], 300)
        assert (m.c_gs(v) > 0).all()
        assert (m.c_ds(v) > 0).all()


def test_dq_dv_equals_c_identity():
    """Charge is master: C = dq/dv by construction (central FD check)."""
    for name in PARTS:
        m = _model(name)
        for f_q, f_c, pts in (
            (m.q_gd, m.c_gd, np.array([-300.0, -20.0, -1.0, -0.05, 0.05, 1.0, 5.0, 15.0])),
            (m.q_gs, m.c_gs, np.array([0.5, 3.0, 20.0, 60.0])),
            (m.q_ds, m.c_ds, np.array([0.5, 3.0, 20.0, 60.0])),
        ):
            h = 1e-4 * np.maximum(np.abs(pts), 1.0)
            d = (f_q(pts + h) - f_q(pts - h)) / (2 * h)
            c = f_c(pts)
            assert np.allclose(d, c, rtol=1e-4, atol=1e-12), name


def test_pchip_monotone_through_sj_cliff():
    """No new extrema between nodes (PCHIP on ln C is shape-preserving)."""
    m = _model("ST8L65N044M9")
    t = m.tables.rss
    for i in range(len(t.v) - 1):
        vg = np.linspace(t.v[i], t.v[i + 1], 200)
        ci = t.c_at(vg)
        lo, hi = min(t.c[i], t.c[i + 1]), max(t.c[i], t.c[i + 1])
        assert ci.min() >= lo * (1 - 1e-9) and ci.max() <= hi * (1 + 1e-9)


def test_clamp_policy_beyond_digitised_range():
    m = _model("ST8L65N044M9")
    t = m.tables.rss
    v_end, c_end = t.v[-1], t.c[-1]
    # C clamps flat, q grows linearly with that flat C
    assert np.allclose(t.c_at(v_end * 2.0), c_end)
    q_extra = float(t.q_at(v_end * 2.0) - t.q_at(v_end))
    assert np.isclose(q_extra, c_end * v_end, rtol=1e-9)


# ---------------------------------------------------------------------------
# The c_gd splice
# ---------------------------------------------------------------------------

def test_splice_c0_c1_continuity():
    for name in PARTS:
        m = _model(name)
        a = m.anchor
        eps = 1e-6
        c_lo = float(m.c_gd(np.array([-eps]))[0])
        c_hi = float(m.c_gd(np.array([eps]))[0])
        assert np.isclose(c_lo, a.c_gd0, rtol=1e-6)
        assert np.isclose(c_hi, a.c_gd0, rtol=1e-6)
        # directed slopes at the splice: table side = s_gd0 (quadratic head,
        # exact), exp side = k_agd*(c_gdmax - c_gd0) == s_gd0 by definition
        d_lo = (float(m.c_gd(np.array([0.0]))[0]) - c_lo) / eps
        assert d_lo == pytest.approx(a.s_gd0, rel=1e-3)
        p = m.params
        if np.isfinite(p.k_agd):
            assert np.isclose(p.k_agd * (p.c_gdmax - p.c_gd0), a.s_gd0, rtol=1e-12)


def test_capparams_guards_and_flat_branch():
    with pytest.raises(ValueError):
        CapParams(c_gd0=1e-9, s_gd0=1e-10, c_gdmax=0.5e-9)  # below c_gd0
    flat = CapParams(c_gd0=1e-9, s_gd0=1e-10, c_gdmax=1e-9)
    assert np.isinf(flat.k_agd)
    v = np.linspace(0, 20, 50)
    assert np.allclose(flat.c_gdmax, 1e-9)
    m = _model("ISC040N10NM8")
    flat_model = CapModel(m.tables, m.anchor, flat)
    assert np.allclose(flat_model.c_gd(v), 1e-9)
    assert np.allclose(flat_model.q_gd(v), 1e-9 * v)


# ---------------------------------------------------------------------------
# Analytic anchor: trench graded-junction law
# ---------------------------------------------------------------------------

def test_trench_junction_law_anchor():
    caps = load_cap_csv(ROOT / "data" / "ISC040N10NM8_Cap_Vds.csv")
    v, c = caps.branch_xy("rss")

    def junction(x, c0, vbi, mm):
        return c0 / (1.0 + x / vbi) ** mm

    p, _ = curve_fit(junction, v, c, p0=[c[0] * 1.1, 1.0, 0.4], maxfev=20000)
    c0, vbi, mm = p
    assert 0.6 < mm < 0.85 and 4.0 < vbi < 9.0
    assert np.max(np.abs(junction(v, *p) - c) / c) < 0.12
    # the splice slope implied by the junction law matches the first-4-point
    # anchor convention (independent cross-check, CAP_SPEC)
    m = build_cap_model(caps)
    assert m.anchor.s_gd0 == pytest.approx(mm * c0 / vbi, rel=0.35)


def test_anchor_window_convention_pinned():
    """ANCHOR_N = 4 is the convention; the record must reproduce it exactly."""
    for name, cfg in PARTS.items():
        m = _model(name)
        assert m.anchor.n_points == ANCHOR_N == 4
        assert m.anchor.c_gd0 == pytest.approx(cfg["anchor"].c_gd0, rel=1e-5)
        assert m.anchor.s_gd0 == pytest.approx(cfg["anchor"].s_gd0, rel=1e-5)
        assert m.anchor.b2 == pytest.approx(cfg["anchor"].b2, rel=1e-4)


# ---------------------------------------------------------------------------
# c_gdmax extraction records (drift checks)
# ---------------------------------------------------------------------------

def test_cmax_extraction_matches_records():
    for name, cfg in PARTS.items():
        caps = load_cap_csv(ROOT / "data" / f"{name}_Cap_Vds.csv")
        qg = load_qg_csv(ROOT / "data" / f"{name}_Vgs_Qg.csv")
        m = build_cap_model(caps)
        fits = extract_cmax(m, qg, r_ds_on=cfg["r_ds_on"])
        good = [f for f in fits if f.c_gdmax_raw is not None]
        chosen = float(np.median([f.c_gdmax for f in good])) if good else fits[0].c_gdmax
        assert chosen == pytest.approx(cfg["record"].c_gdmax, rel=1e-6), name
        # bracket [c_gd0, 2*c_gs(0)] always respected
        lo, hi = m.anchor.c_gd0, 2.0 * m.c_gs0()
        assert lo - 1e-15 <= chosen <= hi + 1e-15 * hi
        # plateau deficit POSITIVE for every condition (channel charge,
        # limitation 1) and the pinned worst-condition magnitudes hold
        for f in fits:
            assert f.plateau_deficit > 0.0, name
            assert abs(f.plateau_deficit) / f.qgd_meas <= 0.75
        worst = max(f.plateau_deficit / f.qgd_meas for f in fits)
        assert worst == pytest.approx(CAP_PLATEAU_DEFICIT_REL[name], rel=0.02)
        # Qgs table check within 20 % (limitation 4)
        for f in fits:
            assert abs(f.qgs_model - f.qgs_meas) / f.qgs_meas < 0.20


def test_cmax_clamp_diagnostics():
    # SiC: measured seg-3 charge below the flat-branch prediction -> LOWER
    caps = load_cap_csv(ROOT / "data" / "IMCQ120R004M2H_Cap_Vds.csv")
    qg = load_qg_csv(ROOT / "data" / "IMCQ120R004M2H_Vgs_Qg.csv")
    m = build_cap_model(caps)
    (fit,) = extract_cmax(m, qg, r_ds_on=0.004)
    assert fit.clamp == "LOWER" and fit.c_gdmax_raw is None
    assert fit.c_gdmax == m.anchor.c_gd0

    # synthetic UPPER clamp: inflate the segment-3 charge (x3 keeps the real
    # plateau the max-dQ/dV segment, but pushes the remainder past the
    # upper bracket -> clamp UPPER)
    q, v = qg.condition_xy(185.2, 800.0)
    q_big = np.array([q[0], q[1], q[2], q[2] + 3.0 * (q[3] - q[2])])
    v_big = np.array([v[0], v[1], v[2], v[3]])
    qg_big = QgCurveSet(q_big, v_big, np.full(4, 185.2), np.full(4, 800.0))
    (fit2,) = extract_cmax(m, qg_big, r_ds_on=0.004)
    assert fit2.clamp == "UPPER" and fit2.c_gdmax_raw is None
    assert fit2.c_gdmax == pytest.approx(2.0 * m.c_gs0())


def test_branch_table_rejects_bad_nodes():
    with pytest.raises(ValueError):
        _BranchTable(np.array([1.0, 2.0]), np.array([1e-12, 0.0]))
    with pytest.raises(ValueError):
        _BranchTable(np.array([1.0]), np.array([1e-12]))


# ---------------------------------------------------------------------------
# LTspice export (charge-table ramp sum + live parity)
# ---------------------------------------------------------------------------

def _model_record(name: str) -> CapModel:
    """Model re-spliced at the extracted c_gdmax record (as export_caps does)."""
    m = _model(name)
    rec = PARTS[name]["record"]
    return CapModel(m.tables, m.anchor, rec)


def _eval_rampsum(v: np.ndarray, q: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Evaluate the exported LTspice expression in Python (max -> np.maximum)."""
    expr = _rampsum_expr(v, q).replace("max(0,", "np.maximum(0,")
    return eval(expr, {"np": np}, {"x": np.asarray(x, float)})


def test_charge_tables_monotone_and_rampsum_identity():
    for name in PARTS:
        m = _model_record(name)
        qs = cap_charge_tables(m)
        for _br, (v, q) in qs.items():
            assert np.all(np.diff(v) > 0) and np.all(np.diff(q) > 0), name
            # the exported expression reproduces the piecewise-linear q...
            xt = np.linspace(v[0], v[-1], 60)
            qe = _eval_rampsum(v, q, xt)
            assert np.allclose(qe, np.interp(xt, v, q), rtol=1e-6), name
            # ... and its clamp-C policy beyond the span (end-segment slope,
            # q keeps growing linearly -- not table()'s flat-q clamp).
            m_lo = (q[1] - q[0]) / (v[1] - v[0])
            m_hi = (q[-1] - q[-2]) / (v[-1] - v[-2])
            d = 0.2 * (v[-1] - v[0])
            qlo = _eval_rampsum(v, q, np.array([v[0] - d]))[0]
            qhi = _eval_rampsum(v, q, np.array([v[-1] + d]))[0]
            assert (q[0] - qlo) / d == pytest.approx(m_lo, rel=1e-6), name
            assert (qhi - q[-1]) / d == pytest.approx(m_hi, rel=1e-6), name


def test_export_subckt_structure():
    m = _model_record("ISC040N10NM8")
    txt = export_cap_subckt(m, "isc_caps")
    assert txt.startswith("* ImprovedVDMOS")
    assert ".subckt isc_caps G D S" in txt and txt.rstrip().endswith(".ends isc_caps")
    for dev in ("Cgs G S", "Cgd G D", "Cds D S"):
        assert f"{dev} Q=" in txt
    # continuation lines: content (never a second '+') follows the marker
    for line in txt.splitlines():
        if line.startswith("+"):
            assert not line[1:].lstrip().startswith("+")
    qs = cap_charge_tables(m)
    # one max() per interior knot (slope increments), 10-g significants
    assert txt.count("max(0,") == sum(len(vv) - 2 for vv, _ in qs.values())
    assert re.search(r"\d\.\d{9}[eE]?", txt)


def _ltspice_exe() -> str | None:
    for name in ("LTspice", "ltspice"):
        found = shutil.which(name)
        if found:
            return found
    candidates = [
        Path(r"C:\Program Files\ADI\LTspice\LTspice.exe"),
        Path.home() / r"AppData\Local\Programs\ADI\LTspice\LTspice.exe",
    ]
    for c in candidates:
        if c.exists():
            return str(c)
    return None


LTSPICE = _ltspice_exe()


def _run_ltspice_meas(cir_text: str, tmp_path: Path) -> dict[str, float]:
    tmp_path.mkdir(parents=True, exist_ok=True)
    cir = tmp_path / "parity.cir"
    cir.write_text(cir_text, encoding="utf-8")
    subprocess.run([LTSPICE, "-b", "-w", str(cir)], cwd=tmp_path,
                   capture_output=True, timeout=300)
    log = (tmp_path / "parity.log").read_text()
    out: dict[str, float] = {}
    for m in re.finditer(r"^([a-z]\d+):\s*I\(V\w+\)\s*=\s*([-\d.e+]+)", log, re.MULTILINE):
        out[m.group(1)] = float(m.group(2))
    assert out, f"no .meas results parsed from {tmp_path / 'parity.log'}"
    return out


def _stair_c(qs: dict, branch: str, v: float) -> float:
    vv, qq = qs[branch]
    m = (qq[1:] - qq[:-1]) / (vv[1:] - vv[:-1])
    i = min(max(int(np.searchsorted(vv, v)) - 1, 0), len(m) - 1)
    return float(m[i])


@pytest.mark.skipif(LTSPICE is None, reason="LTspice not installed")
def test_ltspice_caps_parity(tmp_path):
    """Exported nonlinear caps vs the Python model under PWL ramp drive.

    Three drive configurations, each exciting two delta branches (the
    delta cannot isolate one branch by DC wiring alone): A drives G with
    D=S=0 (c_gs + c_gd), B drives G=D with S=0 (c_gs + c_ds), C drives D
    with G=S=0 (c_ds + c_gd(-v)). Expected I(Vdrive) = -dq/dt =
    -C_staircase*v_slope, exact against the EXPORTED staircase (the
    staircase-vs-PCHIP agreement is checked separately at 5 %). The
    maxstep (5 ns -> ~0.05 V Gear-memory window at the 1e7 V/s ramps)
    keeps the order-2 companion's history term from smoothing across the
    staircase C jumps at the segment-midpoint measure points.
    """
    name = "ISC040N10NM8"
    m = _model_record(name)
    qs = cap_charge_tables(m)

    def _seg_mid(branches: tuple[str, ...], v: float) -> float:
        """Midpoint of the common segment of all branch grids around v.

        FIND AT must land strictly INSIDE staircase-C segments: a target
        exactly on a knot (e.g. v = 80, a native rss node) measures the
        solver's blend across the C jump, not either segment's value.
        """
        lo, hi = -np.inf, np.inf
        for br in branches:
            vv = qs[br][0]
            if v < vv[0] or v > vv[-1]:
                continue  # flat-C extension: no knot constraint outside
            i = min(max(int(np.searchsorted(vv, v)) - 1, 0), len(vv) - 2)
            lo, hi = max(lo, vv[i]), min(hi, vv[i + 1])
        return 0.5 * (lo + hi)

    sub = export_cap_subckt(m, "isc_caps")
    S = 1e7
    raw = {"a1": (-30.0, -S), "a2": (-60.0, -S), "a3": (-90.0, -S),
           "a4": (-75.0, 5 * S), "a5": (-55.0, 5 * S), "a6": (-35.0, 5 * S),
           "a7": (15.0, 5 * S),
           "b1": (30.0, S), "b2": (60.0, S), "b3": (90.0, S), "b4": (65.0, -S),
           "c1": (20.0, S), "c2": (50.0, S), "c3": (80.0, S)}
    segs = {"a": ("gs", "gd"), "b": ("gs", "ds"), "c": ("ds", "gd")}
    cases = {k: (_seg_mid(segs[k[0]], v), s) for k, (v, s) in raw.items()}
    # times: config A ramp1 0->-100 in 10u (slope -S), ramp2 -100->+25 in
    # 2.5u (slope +5S); config B 0->+100 in 10u, back at -S; config C 0->+100.
    def t_of(k: str, v: float, s: float) -> float:
        if k[0] == "a":
            return abs(v) / S if s < 0 else 11e-6 + (v + 100.0) / (5 * S)
        if k[0] == "b":
            return v / S if s > 0 else 11e-6 + (100.0 - v) / S
        return v / S
    meas = []
    for k, (v, s) in cases.items():
        src = {"a": "Vgd", "b": "Vg", "c": "Vd"}[k[0]]
        meas.append(f".meas tran {k} FIND I({src}) AT={t_of(k, v, s) * 1e6:.6f}u")
    net = sub + """
XA gd 0 0 isc_caps
Vgd gd 0 pwl(0 0, 10u -100, 11u -100, 13.5u 25, 14u 25)
XB g g 0 isc_caps
Vg g 0 pwl(0 0, 10u 100, 11u 100, 21u 0, 22u 0)
XC 0 d 0 isc_caps
Vd d 0 pwl(0 0, 10u 100, 11u 100)
""" + "\n".join(meas) + "\n.tran 0 22u 0 5n\n.end\n"
    got = _run_ltspice_meas(net, tmp_path)
    f = lambda x: float(np.ravel(x)[0])  # noqa: E731
    for k, (v, s) in cases.items():
        if k[0] == "a":
            c_st = _stair_c(qs, "gs", v) + _stair_c(qs, "gd", v)
            c_sm = f(m.c_gs(v)) + f(m.c_gd(v))
        elif k[0] == "b":
            c_st = _stair_c(qs, "gs", v) + _stair_c(qs, "ds", v)
            c_sm = f(m.c_gs(v)) + f(m.c_ds(v))
        else:
            c_st = _stair_c(qs, "ds", v) + _stair_c(qs, "gd", -v)
            c_sm = f(m.c_ds(v)) + f(m.c_gd(-v))
        assert got[k] == pytest.approx(-s * c_st, rel=2e-3), k  # exported model
        assert c_st == pytest.approx(c_sm, rel=0.05), k  # staircase vs PCHIP
