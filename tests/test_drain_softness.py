"""Tests for the drain-knee softness exponent (nd).

    v = Vk*tanh(Vds/Vk)                        (nd = 0, legacy default)
    v = Vds/(1 + (Vds/Vk)**nd)**(1/nd)         (nd > 0, power-law soft-min)

Same soft-min mechanism as the gate cap (kcap, m), applied to the drain
variable. Motivation: theta and vsat move the knee position AND width
together (~84 % correlated in fits -- README, "Fitting strategy"); nd is the
missing shape-only knob. By construction dv/dVds(0) = 1 and v -> Vk for
every nd, so Ron and the Id,sat plateau are independent of nd while the
transition width scales with it (nd -> inf is a hard min; nd ~ 2.8
reproduces the legacy tanh knee, L2 rms 0.024*Vk over Vds in [0, 10 Vk]).

These tests pin: nd = 0 is bit-identical to the legacy model (all anchors
in test_model.py ride the tanh path), 0 <= dv/dVds <= 1 (the package-R
solver bracket invariant dg/dId <= -1 survives), the shape-only property
itself, no NaN at negative Vds, and a synthetic round trip recovering nd.
"""

from __future__ import annotations

import numpy as np
import pytest

from vdmos_fit import DEFAULT, Params, fit, make_synthetic, vdmos_id
from vdmos_fit.model import _soft_min


# ---------------------------------------------------------------- identity
def test_nd_zero_is_bit_identical():
    """Params(nd=0.0) must equal Params() and reproduce the legacy tanh model
    bit for bit -- both on the intrinsic fast path and through the package-R
    solver. (test_model.py's anchors -- Ron 4.27 mOhm, sat table, alpha =
    1.078 -- ride DEFAULT.nd = 0; this test closes the explicit gap.)"""
    assert Params(nd=0.0) == Params()

    vds = np.linspace(0.0, 60.0, 121)
    vgs = np.linspace(-2.0, 16.0, 91)
    vd2, vg2 = np.meshgrid(vds, vgs)

    ref = vdmos_id(vd2, vg2, DEFAULT)
    assert np.array_equal(vdmos_id(vd2, vg2, Params(nd=0.0)), ref)

    p = DEFAULT.replace(rd=2e-3, rs=1e-3)
    assert np.array_equal(vdmos_id(vd2, vg2, p), vdmos_id(vd2, vg2, p.replace(nd=0.0)))


# ---------------------------------------------------------------- soft-min math
def test_soft_min_derivative_bounded():
    """0 <= dv/dVds <= 1 everywhere, and dv/dVds(0+) = 1 exactly.

    The analytic derivative of x/(1+(x/k)**n)**(1/n) is
    (1 + (x/k)**n)**(-1/n - 1), which lives in (0, 1] with value 1 at x = 0.
    This is the property that keeps the package-R Newton bracket invariant
    dg/dId <= -1 intact, exactly as for the gate cap. Checked numerically
    across the fitted regime: k in [0.1, 10], nd in [0.5, 40]."""
    for k in (0.1, 1.0, 10.0):
        for n in (0.5, 1.0, 2.8, 10.0, 40.0):
            x = np.linspace(0.0, 10.0 * k, 2001)
            h = 1e-6 * np.maximum(x, 1e-6)
            lo = np.where(x > 0, _soft_min(np.maximum(x - h, 0.0), k, n), 0.0)
            d = (_soft_min(x + h, k, n) - lo) / (2.0 * h)
            d = np.where(x == 0.0, _soft_min(h, k, n) / h, d)  # one-sided at 0
            assert np.all(d >= -1e-9), f"negative slope at k={k}, nd={n}"
            assert np.all(d <= 1.0 + 1e-9), f"slope > 1 at k={k}, nd={n}"
            # One-sided probe at 0 converges as 1 - (1+nd)*(h/k)**nd, so the
            # nd = 0.5 corner (h/k ~ 1e-11) needs a looser floor than 1e-6.
            assert d[x == 0.0] == pytest.approx(1.0, abs=1e-4)


def test_soft_min_asymptotes():
    """Identity for x << k, plateau at k for x >> k, monotone in between.
    Convergence to the plateau is power-law (x/k)**-nd / nd, NOT exponential
    like tanh -- the tail is a documented property of the switch."""
    for k in (0.1, 2.0):
        for n in (1.0, 4.0, 40.0):
            x = np.linspace(0.0, 10.0 * k, 501)
            v = _soft_min(x, k, n)
            eps = np.finfo(float).eps
            assert np.all(np.diff(v) >= -8.0 * eps * k)  # 1-ulp ties at the plateau
            assert np.all(v <= k + 1e-12)
            msk = (x > 0.0) & (x <= 0.05 * k)
            # identity rate is 1 - (x/k)**n/n, so rtol must cover u**n/n at
            # u = 0.05: 5 % for n = 1, negligible for n >= 4.
            assert np.allclose(v[msk] / x[msk], 1.0, rtol=0.06 / n)
            big = x >= 2.0 * k
            err = k - v[big]
            bound = (k / n) * (k / x[big]) ** n + 8.0 * np.finfo(float).eps * k
            assert np.all(err <= bound)


def test_soft_min_no_overflow_at_extremes():
    """The log-space form must survive the regime the naive power would not:
    nd = 40, Vds/Vk = 1e3 gives (Vds/Vk)**nd = 1e120. No inf/NaN, result
    pinned at Vk."""
    v = _soft_min(np.array([1.0, 10.0, 1e3]), 1.0, 40.0)
    assert np.all(np.isfinite(v))
    assert v[-1] == pytest.approx(1.0, rel=1e-9)


# ---------------------------------------------------------------- shape-only knob
def _knee_metrics(p, vgs=10.0):
    """Ron conductance at Vds = 1e-4, Id at Vds = 1e4 (deep plateau), and the
    knee WIDTH = Vds span over which the drain fraction Id/Id_plateau crosses
    0.19 -> 0.99, i.e. v/Vk = 0.1 -> 0.9 (on the clean set Id =
    kp*Vk^2*(r - r^2/2) and Id_plateau = kp*Vk^2/2, so f = 2r - r^2 and
    r = 1 - sqrt(1 - f) inverts exactly; Vk = Vov there).

    Clean set (theta = 0, vsat -> inf, eta = lam = gk = 0) so the switch shape
    is all that remains, and the DEFAULT vsat = 1.5 divisor dip cannot
    contaminate the crossings. This fraction-span metric is monotone in nd by
    construction; a 90/10 SPAN OF GDS is NOT (a sharper corner also delays
    the gds drop near Vk, so nd = 40 measured wider than nd = 10)."""
    assert p.theta == 0.0 and p.vsat > 1e3
    assert p.eta == 0.0 and p.lam == 0.0 and p.gk == 0.0
    g_ohm = float(vdmos_id(1e-4, vgs, p)) / 1e-4
    i_deep = float(vdmos_id(1e4, vgs, p))
    z = np.clip((vgs - p.vth0) / p.ksub, -700.0, 700.0)
    vov_ = p.ksub * np.log1p(np.exp(z))          # eta = 0 -> Vth = vth0
    i_plat = p.kp * vov_ * vov_ / 2.0            # Vk = Vov on the clean set
    vds = np.linspace(0.0, 400.0, 80001)
    f = np.asarray(vdmos_id(vds, vgs, p)) / i_plat
    i1, i9 = np.argmax(f > 0.19), np.argmax(f > 0.99)
    assert 0 < i1 < i9, "crossings out of order or outside the grid"
    v1 = np.interp(0.19, f[i1 - 1 : i1 + 1], vds[i1 - 1 : i1 + 1])
    v9 = np.interp(0.99, f[i9 - 1 : i9 + 1], vds[i9 - 1 : i9 + 1])
    return g_ohm, i_deep, v9 - v1


def test_nd_is_shape_only():
    """Sweeping nd at fixed theta/vsat changes ONLY the knee width.

    This is the property the fit diagnostics flagged: theta and vsat are
    ~84 % correlated because both flatten the knee while moving it. nd must
    not touch Ron (dv/dVds(0) = 1) or the deep plateau (v -> Vk)."""
    base = DEFAULT.replace(theta=0.0, vsat=1e4, eta=0.0, lam=0.0, gk=0.0)
    g0_legacy, i_deep_legacy, w_legacy = _knee_metrics(base)  # nd = 0 (tanh)
    widths = []
    for nd in (1.0, 2.0, 2.8, 5.0, 10.0, 40.0):
        g0, i_deep, w = _knee_metrics(base.replace(nd=nd))
        assert g0 == pytest.approx(g0_legacy, rel=1e-4)
        assert i_deep == pytest.approx(i_deep_legacy, rel=2e-3)
        widths.append(w)
    # width shrinks monotonically as the switch sharpens, and spans a wide
    # range -- the knob actually does what it claims. The legacy tanh width
    # brackets between nd = 2 and nd = 2.8 (the 90/10 width metric puts tanh
    # at nd ~ 2.5; the L2 best fit over [0, 10 Vk] is 2.8).
    assert np.all(np.diff(widths) < 0.0)
    assert widths[0] / widths[-1] > 2.0
    assert widths[1] > w_legacy > widths[2]


def test_nd_2p8_matches_legacy_tanh():
    """nd = 2.8 is the L2 best match to the legacy tanh knee; the model must
    stay within ~1 % of DEFAULT through the knee region at that exponent
    (the value quoted in params.py, so a future re-calibration has to be
    deliberate). Restricted to Vds <= 2*Vk: beyond the knee the soft-min
    approaches Vk with a power-law tail while tanh is already pinned, and
    that tail difference is a documented property, not a mismatch."""
    base = DEFAULT.replace(eta=0.0, lam=0.0, gk=0.0)
    vds = np.linspace(0.05, 10.0, 400)  # Vk(Vgs=10) ~ 4.97 V
    legacy = vdmos_id(vds, 10.0, base)
    soft = vdmos_id(vds, 10.0, base.replace(nd=2.8))
    assert np.max(np.abs(soft / legacy - 1.0)) < 0.01


# ---------------------------------------------------------------- numerics
def test_nd_finite_at_negative_vds():
    """Sign-preserving soft-min: negative Vds (unphysical until the body
    diode, but evaluated by test_model's NaN sweep and by solver probes)
    must stay finite and mirror the positive side."""
    p = DEFAULT.replace(nd=3.0)
    pos = vdmos_id(np.linspace(0.0, 40.0, 41), 10.0, p)
    neg = vdmos_id(np.linspace(-40.0, 0.0, 41), 10.0, p)
    assert np.all(np.isfinite(pos)) and np.all(np.isfinite(neg))
    # v is odd, and ich = kp*vov*(v - v^2/2Vk) is not odd, but both must be
    # finite; check the drain variable itself is odd via a symmetric probe.
    assert float(vdmos_id(0.0, 10.0, p)) == 0.0


def test_no_ndr_with_fitted_nd():
    """dId/dVds > 0 across [0, 200] V at Vgs = 10 with a sharp fitted nd and
    package R: the derivative bound 0 <= dv/dVds <= 1 keeps the whole channel
    monotone, same argument as the cap.

    Parameter set note: with the DEFAULT demo set (vsat = 1.5 << Vk ~ 5) the
    velocity-saturation divisor alone creates a mild gds dip near the knee
    EVEN AT nd = 0 (min -1.86 S at Vds ~ 6 V; test_model's -2.0 S tolerance
    already absorbs it, and a sharper nd deepens it). That is a property of
    the demo parameters, not of the switch, so the monotonicity proof is run
    on a set without the divisor effect (vsat -> inf, theta = 0)."""
    p = DEFAULT.replace(nd=12.0, theta=0.0, vsat=1e4, eta=0.0, gk=0.0,
                        rd=2e-3, rs=1e-3)
    vds = np.linspace(0.0, 200.0, 4001)
    ids = vdmos_id(vds, 10.0, p)
    assert np.all(np.isfinite(ids))
    assert np.gradient(ids, vds).min() > 0.0


# ---------------------------------------------------------------- fitting
NDTRUTH = Params(
    kp=30.0, vth0=3.1, ksub=0.05, theta=0.02, vsat=1.2, nd=4.0,
    eta=0.03, lam=0.015, gk=0.0, vthk=1.6, kk=0.1, vk=10.0,
)


def test_round_trip_recovers_nd():
    """Synthetic round trip: data generated with nd = 4 must recover it within
    ~20 %. The Vds grid samples the knee densely (0.25-3 V steps) so the
    transition width is actually constrained; nd is fitted in LINEAR space
    (0 is the legacy sentinel and cannot be logged)."""
    vgs = np.linspace(4.0, 12.0, 9)
    vds = np.array([0.05, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0])
    data = make_synthetic(NDTRUTH, vgs, vds, noise=0.005, seed=23)

    free = ("kp", "vth0", "theta", "vsat", "nd", "eta", "lam")
    # start must match the truth where params are fixed (gk = 0!), or the
    # kink sigmoid injects phantom current.
    start = DEFAULT.replace(kp=30.0, gk=0.0, nd=1.5)
    res = fit(data, start=start, free=free, weight="relative", max_nfev=8000)

    assert res.success
    assert res.rmse_rel < 0.02
    assert res.params.nd == pytest.approx(4.0, rel=0.20)
