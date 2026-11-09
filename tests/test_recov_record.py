"""Pinned tests for the recovery record (ISC040_RECOV, pole gauge).

Recovery lives on the LTspice side (DiodeRecov netlist), so the pin is:
(1) the record values themselves, (2) the bench netlist STRUCTURE the
extraction driver renders from them (pole gauge, no dTmax, no .ic --
session 25b/25c doctrine), and (3) the deterministic index-step
decomposition. Live-simulation parity is recorded in
out/recov_extract_isc040_run/extract_report.txt; this test itself needs
no LTspice run (only the module import, which locates LTspice.exe).
"""

import re

import pytest

try:
    from vdmos_fit.extract_recov import (
        ISC040,
        grid_points,
        index_step_block,
        make_netlist,
    )
except RuntimeError:  # LTspice.exe not installed
    pytest.skip("extract_recov requires an LTspice installation",
                allow_module_level=True)

from vdmos_fit import ISC040_RECOV


def test_record_values_pinned():
    """The official pole-gauge fit. tau is a FLAT valley coupled to both
    other knobs (interior optimum ~3e-4): quote the coupled set + RMS,
    never tau alone. These literals are the extraction-run optimum
    (out/recov_extract_isc040_run/extract_results.json)."""
    assert ISC040_RECOV.tau == pytest.approx(3.7892757739719363e-4, rel=1e-12)
    assert ISC040_RECOV.t_minor == pytest.approx(2.2363459877611823e-8, rel=1e-12)
    assert ISC040_RECOV.t_pole == pytest.approx(8.32388859768649e-10, rel=1e-12)
    assert ISC040_RECOV.rms_rel == pytest.approx(0.011132, rel=1e-4)
    assert ISC040_RECOV.targets == ((1000.0, 247.0, 25.0), (100.0, 36.0, 29.0))
    # both datasheet rows hit within 2 % on the intrinsic path
    for (didt, q_t, t_t), (_d, q, t) in zip(ISC040_RECOV.targets,
                                            ISC040_RECOV.achieved,
                                            strict=True):
        assert didt == _d
        assert abs(q / q_t - 1) < 0.02
        assert abs(t / t_t - 1) < 0.02


def test_record_within_search_bounds():
    for nm, v in (("tau", ISC040_RECOV.tau),
                  ("t_minor", ISC040_RECOV.t_minor),
                  ("t_pole", ISC040_RECOV.t_pole)):
        lo, hi = ISC040.bounds[nm]
        assert lo < v < hi, nm


def test_weights_broadcast():
    """Session-30 crash fix: the old fixed 4-tuple weights default made
    strict zip() die for any point count other than 4 (GUI recovery fit
    on a 1- or 2-row datasheet table). Weights now broadcast."""
    from vdmos_fit.extract_recov import _weights
    assert _weights(None, 2) == (1.0, 1.0)
    assert _weights(None, 1) == (1.0,)
    assert _weights((2.0,), 3) == (2.0, 2.0, 2.0)
    assert _weights((1.0, 2.0), 2) == (1.0, 2.0)
    with pytest.raises(ValueError):
        _weights((1.0, 1.0, 1.0, 1.0), 2)


def test_make_netlist_i_f_override():
    """Session 30: the bench forward current comes from the datasheet
    table (I_sd column), not a hardcode -- make_netlist takes a per-bench
    i_f; the spec default (ISC040 record: 12.5 A) stays the fallback."""
    net = make_netlist(ISC040, 1000.0, tau=ISC040_RECOV.tau,
                       t_minor=ISC040_RECOV.t_minor,
                       t_pole=ISC040_RECOV.t_pole)
    assert ".param i_test=12.5 " in net
    net2 = make_netlist(ISC040, 1000.0, tau=ISC040_RECOV.tau,
                        t_minor=ISC040_RECOV.t_minor,
                        t_pole=ISC040_RECOV.t_pole, i_f=3.5)
    assert ".param i_test=3.5 " in net2
    assert ".param i_test=12.5 " not in net2


def test_bench_netlist_pole_gauge_structure():
    """The rendered bench must stay in the pole gauge: q_0 (the deleted
    gauge) and the pre-25b form must never return; the 25c doctrine
    (no dTmax, no .ic) is part of the contract."""
    net = make_netlist(ISC040, 1000.0, tau=ISC040_RECOV.tau,
                       t_minor=ISC040_RECOV.t_minor,
                       t_pole=ISC040_RECOV.t_pole)
    assert "C1 q_j0 0 {tm_p} Rser=0" in net
    assert "E1 N002 0 q_j0 0 1" in net and "R3 q_j N002 1" in net
    assert "C3 q_j 0 {tp_p}" in net
    assert "Roff=1" in net
    assert "q_0" not in net
    assert ".tran 0 {tr*5+50n} 0 convreport" in net   # no dTmax arg
    assert re.search(r"(?m)^\s*\.ic\b", net) is None  # no .ic statement (comment OK)
    assert "method=gear" in net and "solver=alt" in net
    assert f".param tau_p={ISC040_RECOV.tau:.10e}" in net


def test_index_step_decomposition_deterministic():
    """Mixed-radix index sweep: tau fastest, t_pole slowest, corners
    pinned (this mapping is what parse_stepped row order assumes)."""
    pts = grid_points(ISC040)
    na, nb, nc = (ISC040.grid_tau[0], ISC040.grid_tminor[0],
                  ISC040.grid_tpole[0])
    assert len(pts) == na * nb * nc
    assert pts[0][0] == pytest.approx(ISC040.grid_tau[1])       # tau lo
    assert pts[na - 1][0] == pytest.approx(ISC040.grid_tau[2])  # tau hi fastest
    assert pts[-1][2] == pytest.approx(ISC040.grid_tpole[2])    # t_pole hi last
    blk = index_step_block(ISC040)
    first = blk.splitlines()[0]
    assert first.startswith(".step param k list 0 ")
    assert first.endswith(str(na * nb * nc - 1))
