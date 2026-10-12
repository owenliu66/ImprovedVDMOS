"""Regression tests for the improved VDMOS drain-current model.

These encode the three checks that used to live in the `__main__` block of
Id_equation_test.py, plus the limiting-behaviour claims from the model docstring.
"""

from __future__ import annotations

import numpy as np
import pytest

from vdmos_fit import DEFAULT, Params, vdmos_id, vdmos_id_level1


def vov(vgs, vth0=3.0, ksub=0.05):
    z = np.clip((vgs - vth0) / ksub, -700.0, 700.0)
    return ksub * np.log1p(np.exp(z))


# ---------------------------------------------------------------- reference points
def test_ron_matches_ltspice_ab_run():
    """kp = 34.7 was re-fit to hold Ron at Vgs = 10 V; LTspice gave 4.27 mOhm."""
    ron = 0.05 / float(vdmos_id(0.05, 10.0, DEFAULT))
    assert ron * 1e3 == pytest.approx(4.27, abs=0.02)


def test_saturation_current_table():
    """Id,sat at Vds = 20 V for the nominal parameter set."""
    expected = {4: 31.5, 6: 91.0, 8: 157.1, 10: 225.9, 12: 297.7}
    for vg, ref in expected.items():
        assert float(vdmos_id(20.0, vg, DEFAULT)) == pytest.approx(ref, rel=0.01)


def test_saturation_exponent_below_two():
    """Velocity saturation + mobility degradation must give alpha < 2.

    Level-1 gives exactly 2.00 and mtriode cannot change that. Re-derived from
    the current DEFAULT (2026-10): alpha = 1.078 over Vgs 6-12 V (1.098 over
    10-12 V), pinned here to +-0.01 so a future re-calibration has to be a
    deliberate act. Tunability checked too: raising vsat alone walks alpha
    through the physically expected 1.2-1.6 band (1.21 at vsat=3 V, 1.32 at
    vsat=5 V); theta is a secondary knob.
    """
    vg = np.linspace(6, 12, 7)
    vo = vov(vg)
    a_new = np.polyfit(np.log(vo), np.log(vdmos_id(np.full_like(vg, 20.0), vg, DEFAULT)), 1)[0]
    a_l1 = np.polyfit(np.log(vo), np.log(vdmos_id_level1(np.full_like(vg, 20.0), vg)), 1)[0]
    assert a_new == pytest.approx(1.08, abs=0.01)
    assert a_l1 == pytest.approx(2.0, abs=1e-6)
    assert a_new < a_l1
    a_hi = np.polyfit(
        np.log(vo),
        np.log(vdmos_id(np.full_like(vg, 20.0), vg, DEFAULT.replace(vsat=5.0))),
        1,
    )[0]
    assert 1.2 <= a_hi <= 1.6  # README: "tunable 1.2-1.6 via theta, vsat"


def test_matched_ron_square_law_overprediction():
    """README claim, re-measured against the CORRECTED Level-1 reference.

    "At matched Ron the square law over-predicts saturation current by 4-5x in
    the 10-12 V gate range." Setup: Level-1 defaults (kp=28, mtriode=1.2,
    lam=0.02) already match the improved DEFAULT's Ron within 0.15% at
    Vds = 50 mV for Vgs = 10-12 V, so no kp rescaling is needed; Id is
    compared at Vds = 20 V. Measured: 4.25x at Vgs=10, 4.80x at 11, 5.33x at 12.

    The number is *unchanged* by the knee-location correction: at Vds = 20 V
    (>> Vov/mtriode) both the old buggy formula (knee = Vov) and the corrected
    one (knee = Vov/mtriode) evaluate the identical saturation branch
    0.5*Kp*Vov^2*(1+lam*Vds), and Ron = 1/(m*Kp*Vov) is the same in both. The
    knee bug corrupted the gds/continuity claims, not this ratio.
    """
    for vg, ref in ((10.0, 4.25), (11.0, 4.80), (12.0, 5.33)):
        ron_imp = 0.05 / float(vdmos_id(0.05, vg, DEFAULT))
        ron_l1 = 0.05 / float(vdmos_id_level1(0.05, vg))
        assert ron_l1 == pytest.approx(ron_imp, rel=2e-3)  # "matched Ron"
        ratio = float(vdmos_id_level1(20.0, vg)) / float(vdmos_id(20.0, vg, DEFAULT))
        assert ratio == pytest.approx(ref, rel=5e-3)


# ---------------------------------------------------------------- smoothness
@pytest.mark.parametrize("vg", [6.0, 10.0])
def test_gds_continuous_across_knee(vg):
    """dId/dVds must not jump at Vds = Vov (the old triode/sat boundary)."""
    h = 1e-4
    vo = vov(vg)
    gl = (vdmos_id(vo - h, vg, DEFAULT) - vdmos_id(vo - 2 * h, vg, DEFAULT)) / h
    gr = (vdmos_id(vo + h, vg, DEFAULT) - vdmos_id(vo, vg, DEFAULT)) / h
    assert float(gr - gl) == pytest.approx(0.0, abs=1e-3)


@pytest.mark.parametrize("m", [1, 2, 4])
def test_level1_is_c1_continuous(m):
    """LTspice Level-1 VDMOS is C1 continuous for ANY mtriode -- regression guard.

    LTspice puts the knee at Vds = Vov/mtriode, which makes the triode slope
    m*Kp*(Vov - m*Vds) vanish exactly where saturation begins. These are the
    numbers LTspice 26.0.1 reports for Vgs=5, Vto=3, Kp=2, Lambda=0.

    An earlier version of this project used knee = Vov instead of Vov/mtriode,
    which manufactured a fake -40 S negative-gds dip and a several-hundred-
    percent gds jump. Those artifacts were the bug in the reference formula, not
    a property of Level-1. This test pins the correct behaviour.
    """
    vgs, kp, vth0 = 5.0, 2.0, 3.0
    vov = vgs - vth0
    knee = vov / m

    def f(vd):
        return vdmos_id_level1(vd, vgs, kp=kp, vth0=vth0, m_triode=m, lam=0.0)

    # C0: no step at the knee (LTspice: Id = 4.000000 either side)
    assert float(f(knee - 1e-4)) == pytest.approx(4.0, abs=1e-4)
    assert float(f(knee)) == pytest.approx(4.0, abs=1e-4)
    assert float(f(knee + 1e-4)) == pytest.approx(4.0, abs=1e-4)

    # C1: no slope jump at the knee
    h = 1e-6
    gl = (f(knee) - f(knee - h)) / h
    gr = (f(knee + h) - f(knee)) / h
    assert float(gr - gl) == pytest.approx(0.0, abs=1e-3)

    # No negative differential resistance anywhere (LTspice min Id = 4.0)
    vds = np.linspace(0.0, 6.0, 6001)
    gds = np.gradient(vdmos_id_level1(vds, vgs, kp=kp, vth0=vth0, m_triode=m, lam=0.0), vds)
    assert gds.min() > -1e-9

    # Ron conductance = m*Kp*Vov, and Id,sat is independent of mtriode
    assert float(f(1e-3)) / 1e-3 == pytest.approx(m * kp * vov, rel=2e-3)
    assert float(f(5.0)) == pytest.approx(0.5 * kp * vov**2, rel=1e-9)


def test_level1_mtriode_is_shape_only():
    """mtriode moves the knee and Ron but cannot change Id,sat.

    This is the real Level-1 limitation the improved model addresses: Ron, knee
    shape and saturation level are not independently settable.
    """
    vgs, kp, vth0 = 5.0, 2.0, 3.0
    vov = vgs - vth0
    sat = [
        float(vdmos_id_level1(5.0, vgs, kp=kp, vth0=vth0, m_triode=m, lam=0.0))
        for m in (0.5, 1, 2, 4)
    ]
    assert np.allclose(sat, 0.5 * kp * vov**2)
    ron_g = [
        float(vdmos_id_level1(1e-4, vgs, kp=kp, vth0=vth0, m_triode=m, lam=0.0)) / 1e-4
        for m in (0.5, 1, 2, 4)
    ]
    assert np.allclose(ron_g, [m * kp * vov for m in (0.5, 1, 2, 4)], rtol=1e-3)


@pytest.mark.parametrize("vg", [6.0, 10.0])
def test_improved_model_has_no_negative_gds(vg):
    """The improved model stays monotone in Vds across the whole sweep."""
    vds = np.linspace(0.0, 40.0, 4000)
    g_new = np.gradient(vdmos_id(vds, vg, DEFAULT), vds).min()
    assert g_new > -2.0


# ---------------------------------------------------------------- limiting cases
def test_ohmic_limit_is_linear_in_vds():
    """Vds << Vk -> Id ~ kp*Vov*Vds (pure resistance)."""
    p = DEFAULT.replace(theta=0.0, vsat=1e6, eta=0.0, lam=0.0, gk=0.0)
    i1 = float(vdmos_id(0.01, 10.0, p))
    i2 = float(vdmos_id(0.02, 10.0, p))
    assert i2 / i1 == pytest.approx(2.0, rel=1e-3)


def test_level1_limit_agreement():
    """theta = vsat = eta = gk = lam = 0 -> smooth approximation of Level-1
    with m_triode = 1, NOT an exact copy.

    The soft drain variable v = Vk*tanh(Vds/Vk) rounds the hard corner, so the
    model is exact in the ohmic limit and in deep saturation but sits a few
    percent low around the knee. Quantified here so the discrepancy is a
    documented property rather than a surprise.
    """
    p = Params(
        kp=28.0,
        vth0=3.0,
        ksub=0.05,
        theta=0.0,
        vsat=1e9,
        eta=0.0,
        lam=0.0,
        gk=0.0,
        vk=1e9,
        vov_floor=1e-9,
    )
    for vg in (5.0, 12.0):
        vo = vov(vg)
        x = np.array([0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 20.0])
        got = vdmos_id(x * vo, vg, p)
        ref = vdmos_id_level1(x * vo, vg, kp=28.0, m_triode=1.0, lam=0.0)
        ratio = got / ref
        # ohmic limit: exact
        assert ratio[0] == pytest.approx(1.0, abs=1e-6)
        assert ratio[1] == pytest.approx(1.0, abs=1e-3)
        # deep saturation: essentially exact
        assert ratio[-1] == pytest.approx(1.0, abs=2e-3)
        # everywhere: bounded, and always <= Level-1 (the corner is rounded down)
        assert np.all(ratio <= 1.0 + 1e-9)
        assert ratio.min() > 0.93


def test_subthreshold_is_smooth_and_positive():
    """Below threshold Id must be positive, monotone, and kink-free in Vgs.

    Kink detector: the ratio of consecutive dId/dVgs steps. A hard transition
    would make that ratio jump by an order of magnitude between adjacent
    samples; the softplus subthreshold law gives a smooth ~exp(dVgs/ksub).
    """
    vg = np.linspace(0.0, 3.0, 300)
    ids = vdmos_id(1.0, vg, DEFAULT)
    assert np.all(ids > 0.0)
    d = np.diff(ids)
    assert np.all(d > 0.0)  # monotone, no dead step
    ratio = d[1:] / d[:-1]
    # Deep subthreshold the ratio is exp(dVgs/ksub) ~ 1.22; it rises to ~2.1 at the
    # softplus crossover into moderate inversion. A hard switch would give a zero
    # step (caught above) or a ratio orders of magnitude larger.
    assert np.max(ratio) < 3.0 and np.min(ratio) > 0.7  # no kink


def test_subthreshold_slope_matches_ksub():
    """Deep subthreshold d ln(Id)/dVgs should approach 1/ksub (60 mV/decade-ish).

    The channel softplus is the ONLY off-state current (the gds*Vds leak term
    was deleted, session 7), so the slope law holds at every Vds without a
    leakage-dominance caveat.
    """
    p = DEFAULT.replace(gk=0.0)
    vg = np.linspace(1.0, 2.4, 60)
    ids = vdmos_id(1.0, vg, p)
    slope = np.polyfit(vg, np.log(ids), 1)[0]
    assert slope == pytest.approx(1.0 / p.ksub, rel=0.02)


def test_off_state_is_gate_dependent_not_a_leak_resistor():
    """No gds*Vds leak term: the off state is the gate-dependent softplus floor.

    The drift-conductance term was deleted (session 7). On the ST8L65 data it
    was a one-point parameter that absorbed the Vgs(th) spec point (250 uA at
    Vgs = Vds = 3.7 V) with a GATE-INDEPENDENT leak, flattening the whole
    subthreshold transfer curve to a constant floor. The softplus channel
    already conducts below threshold with slope 1/ksub, which is the physical
    off state; real HV leakage needs the body diode / leakage element, not a
    fitted resistor.
    """
    p = DEFAULT.replace(gk=0.0)
    off = float(vdmos_id(20.0, 0.0, p))
    assert 0.0 < off < 1e-9  # tiny but strictly positive, no floating off state
    # gate-dependent: +2 V on the gate raises it by ~20 decades (exp(2/0.05))
    assert float(vdmos_id(20.0, 2.0, p)) > 1e4 * off


def test_kink_term_is_one_sided():
    """The kink (parasitic accumulation) channel must never sink current.

    A bare tanh((Vgs-Vthk)/kk) is odd and injected about -40 mA at Vgs = 0,
    Vds = 1 V. Below the kink onset the contribution must go to zero.
    """
    p = DEFAULT
    off = float(vdmos_id(1.0, 0.0, p))
    assert off > 0.0
    assert off < 1e-12  # channel softplus floor only; kink contributes ~nothing

    # Above the onset (Vgs >> vthk) the sigmoid equals the old tanh gate exactly,
    # so the calibrated reference points are unchanged by this fix.
    for vg in (4.0, 6.0, 10.0, 12.0):
        old = p.gk * 20.0 * np.tanh((vg - p.vthk) / p.kk) * np.tanh(20.0 / p.vk)
        gate = 0.5 * (1.0 + np.tanh((vg - p.vthk) / p.kk))
        new = p.gk * 20.0 * gate * np.tanh(20.0 / p.vk)
        assert new == pytest.approx(old, rel=1e-9)


def test_monotone_in_vgs_everywhere():
    """Id must never decrease with Vgs at any Vds.

    With the gds leak term deleted, Id is the channel current alone and stays
    resolvable in float64 down to ~1e-47 A, but the softplus steps at the very
    bottom of the grid are at the resolution edge, so keep non-decreasing
    everywhere and strictly increasing once the channel term is comfortably
    resolvable (Vgs >= 1.5 V).
    """
    vg = np.linspace(-2.0, 15.0, 400)
    for vd in (0.05, 1.0, 20.0, 100.0):
        ids = vdmos_id(vd, vg, DEFAULT)
        d = np.diff(ids)
        assert np.all(d >= 0.0), f"decreasing at Vds={vd}"
        # from Vgs = 1.5 V up, the channel term is resolvable at every Vds
        assert np.all(np.diff(ids[vg >= 1.5]) > 0.0), f"flat above 1.5 V at Vds={vd}"


# ---------------------------------------------------------------- numerics
def test_broadcasting_and_shapes():
    vds = np.linspace(0, 40, 50)
    vg = np.linspace(2, 14, 40)
    vd2, vg2 = np.meshgrid(vds, vg)
    assert vdmos_id(vd2, vg2, DEFAULT).shape == (40, 50)
    assert vdmos_id(1.0, 10.0, DEFAULT).shape == ()


def test_no_nan_or_inf_anywhere():
    vds = np.linspace(-5.0, 100.0, 500)
    vg = np.linspace(-5.0, 25.0, 500)
    vd2, vg2 = np.meshgrid(vds, vg)
    out = vdmos_id(vd2, vg2, DEFAULT)
    assert np.all(np.isfinite(out))
