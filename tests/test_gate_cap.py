"""Tests for the gate-overdrive charge cap (kcap, m).

    Vov_eff = Vov / (1 + (Vov/kcap)**m)**(1/m)    (kcap = 0 -> identity)

The cap models finite mobile charge in the super-junction pillars: the measured
gate exponent alpha = d ln Id / d ln Vov collapses from ~1.6 to ~0.07 above
Vgs ~ 8 V on the ST8L65N044M9, and no theta/vsat/lam setting can do that.

These tests pin the three properties the cap is *built* to have --
identity at kcap = 0 (legacy anchors bit-identical), 0 <= dVov_eff/dVov <= 1
(the package-R solver bracket invariant dg/dId <= -1 survives), and no
overflow in the log-space evaluation at the fitted regime m ~ 10-40 -- plus
the physical consequences: subthreshold untouched, no NDR with the fitted
parameters, and a synthetic round trip that recovers kcap and m.
"""

from __future__ import annotations

import numpy as np
import pytest

from vdmos_fit import DEFAULT, Params, fit, make_synthetic, vdmos_id
from vdmos_fit.model import _soft_cap

# Best fit on the ST8L65N044M9 datasheet data (HANDOVER section 3b):
# Vov cap + single rd, 9 free, RMSE 6.88 %. (gds kwarg stripped in
# session 7 when the leak term was deleted; the fitted gds was 6.8e-5 S.)
FITTED = Params(
    kp=91.9481, vth0=4.7327, ksub=0.06, kcap=1.90909, m=9.86552,
    theta=0.0, vsat=1e4, eta=0.0, lam=0.0383371,
    gk=0.0, vthk=1.6, kk=0.1, vk=0.653951,
    rd=26.03e-3, rs=4.644e-3,
)


# ---------------------------------------------------------------- identity
def test_kcap_zero_is_bit_identical():
    """Params(kcap=0.0) must equal Params() and reproduce the legacy model
    bit for bit -- both on the intrinsic fast path and through the package-R
    solver. (test_model.py's anchors -- Ron 4.27 mOhm, sat table
    31.5/91.0/157.1/225.9/297.7 A, alpha = 1.078 -- ride DEFAULT.kcap = 0;
    this test closes the explicit Params-level gap.)"""
    assert Params(kcap=0.0) == Params()

    vds = np.linspace(0.0, 60.0, 121)
    vgs = np.linspace(-2.0, 16.0, 91)
    vd2, vg2 = np.meshgrid(vds, vgs)

    ref = vdmos_id(vd2, vg2, DEFAULT)
    assert np.array_equal(vdmos_id(vd2, vg2, Params(kcap=0.0, m=37.0)), ref)

    p = DEFAULT.replace(rd=2e-3, rs=1e-3)
    assert np.array_equal(vdmos_id(vd2, vg2, p), vdmos_id(vd2, vg2, p.replace(kcap=0.0, m=37.0)))


# ---------------------------------------------------------------- cap math
def test_cap_derivative_bounded():
    """0 <= dVov_eff/dVov <= 1 everywhere.

    This is the property that keeps the package-R Newton bracket invariant
    dg/dId <= -1 intact: the cap can only shrink the channel's dI/dVds, never
    blow it up (the failure mode that made the old kappa model need bound (b)
    and a pole-avoidance dance). Checked numerically across the full fitted
    regime: vov in [1e-6, 50], kcap in [0.1, 20], m in [0.5, 40].
    """
    vov = np.logspace(-6.0, np.log10(50.0), 300)
    for kcap in (0.1, 0.5, 2.0, 20.0):
        for m in (0.5, 1.0, 10.0, 40.0):
            h = 1e-6 * np.maximum(vov, 1e-6)
            d = (_soft_cap(vov + h, kcap, m) - _soft_cap(vov - h, kcap, m)) / (2.0 * h)
            assert np.all(d >= -1e-9), f"negative slope at kcap={kcap}, m={m}"
            assert np.all(d <= 1.0 + 1e-9), f"slope > 1 at kcap={kcap}, m={m}"


def test_cap_monotone_with_kcap_asymptote():
    """Non-decreasing in vov, never exceeds kcap, and approaches kcap at the
    analytic rate. The leading term of kcap - Vov_eff is (kcap/m)*(kcap/vov)**m
    (Bernoulli: 1-(1+u**-m)**(-1/m) <= 1/(m*u**m) for u >= 1), so the bound
    carries the 1/m factor -- the m = 0.5 case violates the m >= 1 constant
    kcap*(kcap/vov)**m by up to 1.34x at u = 10. Near the asymptote the true
    error drops below float64 resolution of kcap, so the comparison floors at
    8 eps * kcap (the m = 40 rows sit at u**-m ~ 1e-30 << eps)."""
    eps = np.finfo(float).eps
    for kcap in (0.1, 2.0, 20.0):
        for m in (0.5, 8.0, 40.0):
            vov = np.linspace(0.0, 10.0 * kcap, 501)
            capped = _soft_cap(vov, kcap, m)
            assert np.all(np.diff(capped) >= -8.0 * eps * kcap)
            assert np.all(capped <= kcap + 8.0 * eps * kcap)
            big = vov[vov >= 2.0 * kcap]
            err = kcap - _soft_cap(big, kcap, m)
            # err is computed by cancellation (kcap - capped), so it carries an
            # absolute rounding floor of ~eps*kcap on top of the analytic bound.
            bound = (kcap / m) * (kcap / big) ** m + 8.0 * eps * kcap
            assert np.all(err <= bound)


def test_cap_no_overflow_at_extremes():
    """The log-space form must survive the regime the naive power would not:
    kcap = 0.05, m = 40, vov = 50 gives (vov/kcap)**m = 1e120, and vov = 1e6
    gives ~1e292. No inf/NaN, and the result stays pinned at kcap."""
    capped = _soft_cap(np.array([0.05, 1.0, 50.0]), 0.05, 40.0)
    assert np.all(np.isfinite(capped))
    assert capped[-1] == pytest.approx(0.05, rel=1e-9)

    huge = _soft_cap(np.array([1e6]), 0.05, 40.0)
    assert np.all(np.isfinite(huge))
    assert huge[0] <= 0.05 * (1.0 + 1e-12)  # log-space exp rounds ~1 ulp over

    # End to end through the model, including the package-R solver.
    p = DEFAULT.replace(kcap=0.05, m=40.0, rd=2e-3, rs=1e-3)
    out = vdmos_id(np.linspace(0.0, 100.0, 201), np.linspace(-5.0, 25.0, 101)[:, None], p)
    assert np.all(np.isfinite(out))


# ---------------------------------------------------------------- physics
def test_subthreshold_slope_untouched_by_cap():
    """Deep below threshold vov ~ ksub*exp((Vgs-Vth)/ksub) << kcap, so the cap
    is the identity there for ANY kcap > 0 and d ln Id/dVgs stays 1/ksub."""
    for kcap in (0.05, 1.0, 20.0):
        p = DEFAULT.replace(gk=0.0, kcap=kcap, m=10.0)
        vg = np.linspace(1.0, 2.4, 60)
        ids = vdmos_id(1.0, vg, p)
        slope = np.polyfit(vg, np.log(ids), 1)[0]
        assert slope == pytest.approx(1.0 / p.ksub, rel=0.02)


def test_no_ndr_with_fitted_cap():
    """dId/dVds > 0 for Vds in [0, 200] at Vgs = 10 with the fitted cap and
    rd = 26 mOhm. This is the check the kappa model could not pass by
    construction (its u(Id) pole); the cap is safe because dVov_eff/dVov <= 1
    and every other channel term is monotone in Vds."""
    vds = np.linspace(0.0, 200.0, 4001)
    ids = vdmos_id(vds, 10.0, FITTED)
    assert np.all(np.isfinite(ids))
    assert np.gradient(ids, vds).min() > 0.0


# ---------------------------------------------------------------- fitting
CAPTRUTH = Params(
    kp=30.0, vth0=3.1, ksub=0.05, kcap=2.0, m=8.0,
    theta=0.02, vsat=1.2, eta=0.03, lam=0.015,
    gk=0.0, vthk=1.6, kk=0.1, vk=10.0,
)


def test_round_trip_recovers_cap():
    """Synthetic round trip: data generated with kcap = 2.0, m = 8 must
    recover both within ~10 %. The Vgs grid deliberately reaches Vov ~ 7 V,
    well past the cap, so the collapse in alpha is actually sampled."""
    vgs = np.linspace(2.0, 10.0, 33)
    vds = np.array([0.05, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0])
    data = make_synthetic(CAPTRUTH, vgs, vds, noise=0.005, seed=11)

    free = ("kp", "vth0", "ksub", "kcap", "m", "theta", "vsat", "eta", "lam")
    # Non-free params come from `start`: it must match the truth where params
    # are fixed (gk = 0!), or the kink sigmoid injects amps of phantom current.
    start = DEFAULT.replace(kp=30.0, gk=0.0, kcap=1.0, m=5.0)
    res = fit(data, start=start, free=free, weight="relative", max_nfev=8000)

    assert res.success
    assert res.rmse_rel < 0.02
    assert res.params.kcap == pytest.approx(2.0, rel=0.10)
    assert res.params.m == pytest.approx(8.0, rel=0.10)
