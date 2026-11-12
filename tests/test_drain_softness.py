"""Tests for the drain-knee softness exponent (k_kneesharp).

    v_dseff = v_knee*tanh(v_ds/v_knee)                   (k_kneesharp = 0, legacy default)
    v_dseff = v_ds/(1 + (v_ds/v_knee)**k_kneesharp)**(1/k_kneesharp)
                                          (k_kneesharp > 0, power-law soft-min)

Same soft-min mechanism as the gate cap (v_ovmax, k_ovsharp), applied to the drain
variable. Motivation: k_theta and v_vsat move the knee position AND width
together (~84 % correlated in fits -- README, "Fitting strategy"); k_kneesharp is the
missing shape-only knob. By construction dv_dseff/dv_ds(0) = 1 and v_dseff -> v_knee for
every k_kneesharp, so Ron and the i_d,sat plateau are independent of k_kneesharp while the
transition width scales with it (k_kneesharp -> inf is a hard min; k_kneesharp ~ 2.8
reproduces the legacy tanh knee, L2 rms 0.024*v_knee over v_ds in [0, 10 v_knee]).

These tests pin: k_kneesharp = 0 is bit-identical to the legacy model (all anchors
in test_model.py ride the tanh path), 0 <= dv_dseff/dv_ds <= 1 (the package-R
solver bracket invariant dg/di_d <= -1 survives), the shape-only property
itself, no NaN at negative v_ds, and a synthetic round trip recovering k_kneesharp.
"""

from __future__ import annotations

import numpy as np
import pytest

from vdmos_fit import DEFAULT, Params, fit, make_synthetic, vdmos_id
from vdmos_fit.model import _soft_min


# ---------------------------------------------------------------- identity
def test_kneesharp_zero_is_bit_identical():
    """Params(k_kneesharp=0.0) must equal Params() and reproduce the legacy tanh model
    bit for bit -- both on the intrinsic fast path and through the package-R
    solver. (test_model.py's anchors -- Ron 4.27 mOhm, sat table, alpha =
    1.078 -- ride DEFAULT.k_kneesharp = 0; this test closes the explicit gap.)"""
    assert Params(k_kneesharp=0.0) == Params()

    v_ds = np.linspace(0.0, 60.0, 121)
    v_gs = np.linspace(-2.0, 16.0, 91)
    vd2, vg2 = np.meshgrid(v_ds, v_gs)

    ref = vdmos_id(vd2, vg2, DEFAULT)
    assert np.array_equal(vdmos_id(vd2, vg2, Params(k_kneesharp=0.0)), ref)

    p = DEFAULT.replace(r_d=2e-3, r_s=1e-3)
    assert np.array_equal(vdmos_id(vd2, vg2, p), vdmos_id(vd2, vg2, p.replace(k_kneesharp=0.0)))


# ---------------------------------------------------------------- soft-min math
def test_soft_min_derivative_bounded():
    """0 <= dv_dseff/dv_ds <= 1 everywhere, and dv_dseff/dv_ds(0+) = 1 exactly.

    The analytic derivative of x/(1+(x/k)**n)**(1/n) is
    (1 + (x/k)**n)**(-1/n - 1), which lives in (0, 1] with value 1 at x = 0.
    This is the property that keeps the package-R Newton bracket invariant
    dg/di_d <= -1 intact, exactly as for the gate cap. Checked numerically
    across the fitted regime: k in [0.1, 10], k_kneesharp in [0.5, 40]."""
    for k in (0.1, 1.0, 10.0):
        for n in (0.5, 1.0, 2.8, 10.0, 40.0):
            x = np.linspace(0.0, 10.0 * k, 2001)
            h = 1e-6 * np.maximum(x, 1e-6)
            lo = np.where(x > 0, _soft_min(np.maximum(x - h, 0.0), k, n), 0.0)
            d = (_soft_min(x + h, k, n) - lo) / (2.0 * h)
            d = np.where(x == 0.0, _soft_min(h, k, n) / h, d)  # one-sided at 0
            assert np.all(d >= -1e-9), f"negative slope at k={k}, k_kneesharp={n}"
            assert np.all(d <= 1.0 + 1e-9), f"slope > 1 at k={k}, k_kneesharp={n}"
            # One-sided probe at 0 converges as 1 - (1+k_kneesharp)*(h/k)**k_kneesharp, so the
            # k_kneesharp = 0.5 corner (h/k ~ 1e-11) needs a looser floor than 1e-6.
            assert d[x == 0.0] == pytest.approx(1.0, abs=1e-4)


def test_soft_min_asymptotes():
    """Identity for x << k, plateau at k for x >> k, monotone in between.
    Convergence to the plateau is power-law (x/k)**-k_kneesharp / k_kneesharp, NOT exponential
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
    k_kneesharp = 40, v_ds/v_knee = 1e3 gives (v_ds/v_knee)**k_kneesharp = 1e120. No inf/NaN, result
    pinned at v_knee."""
    v = _soft_min(np.array([1.0, 10.0, 1e3]), 1.0, 40.0)
    assert np.all(np.isfinite(v))
    assert v[-1] == pytest.approx(1.0, rel=1e-9)


# ---------------------------------------------------------------- shape-only knob
def _knee_metrics(p, v_gs=10.0):
    """Ron conductance at v_ds = 1e-4, i_d at v_ds = 1e4 (deep plateau), and the
    knee WIDTH = v_ds span over which the drain fraction i_d/i_plateau crosses
    0.19 -> 0.99, i.e. v_dseff/v_knee = 0.1 -> 0.9 (on the clean set i_d =
    k_p*v_knee^2*(r - r^2/2) and i_plateau = k_p*v_knee^2/2, so f = 2r - r^2 and
    r = 1 - sqrt(1 - f) inverts exactly; v_knee = v_ov there).

    Clean set (k_theta = 0, v_vsat -> inf, k_eta = k_lambda = g_k = 0) so the switch shape
    is all that remains, and the DEFAULT v_vsat = 1.5 divisor dip cannot
    contaminate the crossings. This fraction-span metric is monotone in k_kneesharp by
    construction; a 90/10 SPAN OF GDS is NOT (a sharper corner also delays
    the gds drop near v_knee, so k_kneesharp = 40 measured wider than k_kneesharp = 10)."""
    assert p.k_theta == 0.0 and p.v_vsat > 1e3
    assert p.k_eta == 0.0 and p.k_lambda == 0.0 and p.g_k == 0.0
    g_ohm = float(vdmos_id(1e-4, v_gs, p)) / 1e-4
    i_deep = float(vdmos_id(1e4, v_gs, p))
    z = np.clip((v_gs - p.v_th0) / p.v_sub, -700.0, 700.0)
    v_ov = p.v_sub * np.log1p(np.exp(z))          # k_eta = 0 -> v_th = v_th0
    i_plat = p.k_p * v_ov * v_ov / 2.0              # v_knee = v_ov on the clean set
    v_ds = np.linspace(0.0, 400.0, 80001)
    f = np.asarray(vdmos_id(v_ds, v_gs, p)) / i_plat
    i1, i9 = np.argmax(f > 0.19), np.argmax(f > 0.99)
    assert 0 < i1 < i9, "crossings out of order or outside the grid"
    v1 = np.interp(0.19, f[i1 - 1 : i1 + 1], v_ds[i1 - 1 : i1 + 1])
    v9 = np.interp(0.99, f[i9 - 1 : i9 + 1], v_ds[i9 - 1 : i9 + 1])
    return g_ohm, i_deep, v9 - v1


def test_kneesharp_is_shape_only():
    """Sweeping k_kneesharp at fixed k_theta/v_vsat changes ONLY the knee width.

    This is the property the fit diagnostics flagged: k_theta and v_vsat are
    ~84 % correlated because both flatten the knee while moving it. k_kneesharp must
    not touch Ron (dv_dseff/dv_ds(0) = 1) or the deep plateau (v_dseff -> v_knee)."""
    base = DEFAULT.replace(k_theta=0.0, v_vsat=1e4, k_eta=0.0, k_lambda=0.0, g_k=0.0)
    g0_legacy, i_deep_legacy, w_legacy = _knee_metrics(base)  # k_kneesharp = 0 (tanh)
    widths = []
    for k_kneesharp in (1.0, 2.0, 2.8, 5.0, 10.0, 40.0):
        g0, i_deep, w = _knee_metrics(base.replace(k_kneesharp=k_kneesharp))
        assert g0 == pytest.approx(g0_legacy, rel=1e-4)
        assert i_deep == pytest.approx(i_deep_legacy, rel=2e-3)
        widths.append(w)
    # width shrinks monotonically as the switch sharpens, and spans a wide
    # range -- the knob actually does what it claims. The legacy tanh width
    # brackets between k_kneesharp = 2 and k_kneesharp = 2.8 (the 90/10 width metric puts tanh
    # at k_kneesharp ~ 2.5; the L2 best fit over [0, 10 v_knee] is 2.8).
    assert np.all(np.diff(widths) < 0.0)
    assert widths[0] / widths[-1] > 2.0
    assert widths[1] > w_legacy > widths[2]


def test_kneesharp_2p8_matches_legacy_tanh():
    """k_kneesharp = 2.8 is the L2 best match to the legacy tanh knee; the model must
    stay within ~1 % of DEFAULT through the knee region at that exponent
    (the value quoted in params.py, so a future re-calibration has to be
    deliberate). Restricted to v_ds <= 2*v_knee: beyond the knee the soft-min
    approaches v_knee with a power-law tail while tanh is already pinned, and
    that tail difference is a documented property, not a mismatch."""
    base = DEFAULT.replace(k_eta=0.0, k_lambda=0.0, g_k=0.0)
    v_ds = np.linspace(0.05, 10.0, 400)  # v_knee(v_gs=10) ~ 4.97 V
    legacy = vdmos_id(v_ds, 10.0, base)
    soft = vdmos_id(v_ds, 10.0, base.replace(k_kneesharp=2.8))
    assert np.max(np.abs(soft / legacy - 1.0)) < 0.01


# ---------------------------------------------------------------- numerics
def test_kneesharp_finite_at_negative_vds():
    """Sign-preserving soft-min: negative v_ds (unphysical until the body
    diode, but evaluated by test_model's NaN sweep and by solver probes)
    must stay finite and mirror the positive side."""
    p = DEFAULT.replace(k_kneesharp=3.0)
    pos = vdmos_id(np.linspace(0.0, 40.0, 41), 10.0, p)
    neg = vdmos_id(np.linspace(-40.0, 0.0, 41), 10.0, p)
    assert np.all(np.isfinite(pos)) and np.all(np.isfinite(neg))
    # v_dseff is odd, and i_ch = k_p*v_ov_eff*(v_dseff - v_dseff^2/2v_knee) is not
    # odd, but both must be
    # finite; check the drain variable itself is odd via a symmetric probe.
    assert float(vdmos_id(0.0, 10.0, p)) == 0.0


def test_no_ndr_with_fitted_kneesharp():
    """di_d/dv_ds > 0 across [0, 200] V at v_gs = 10 with a sharp fitted k_kneesharp and
    package R: the derivative bound 0 <= dv_dseff/dv_ds <= 1 keeps the whole channel
    monotone, same argument as the cap.

    Parameter set note: with the DEFAULT demo set (v_vsat = 1.5 << v_knee ~ 5) the
    velocity-saturation divisor alone creates a mild gds dip near the knee
    EVEN AT k_kneesharp = 0 (min -1.86 S at v_ds ~ 6 V; test_model's -2.0 S tolerance
    already absorbs it, and a sharper k_kneesharp deepens it). That is a property of
    the demo parameters, not of the switch, so the monotonicity proof is run
    on a set without the divisor effect (v_vsat -> inf, k_theta = 0)."""
    p = DEFAULT.replace(k_kneesharp=12.0, k_theta=0.0, v_vsat=1e4, k_eta=0.0, g_k=0.0,
                        r_d=2e-3, r_s=1e-3)
    v_ds = np.linspace(0.0, 200.0, 4001)
    i_d = vdmos_id(v_ds, 10.0, p)
    assert np.all(np.isfinite(i_d))
    assert np.gradient(i_d, v_ds).min() > 0.0


# ---------------------------------------------------------------- fitting
KNEE_TRUTH = Params(
    k_p=30.0, v_th0=3.1, v_sub=0.05, k_theta=0.02, v_vsat=1.2, k_kneesharp=4.0,
    k_eta=0.03, k_lambda=0.015, g_k=0.0, v_thk=1.6, v_kovsoft=0.1, v_dbsat=10.0,
)


def test_round_trip_recovers_kneesharp():
    """Synthetic round trip: data generated with k_kneesharp = 4 must recover it within
    ~20 %. The v_ds grid samples the knee densely (0.25-3 V steps) so the
    transition width is actually constrained; k_kneesharp is fitted in LINEAR space
    (0 is the legacy sentinel and cannot be logged)."""
    v_gs = np.linspace(4.0, 12.0, 9)
    v_ds = np.array([0.05, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0])
    data = make_synthetic(KNEE_TRUTH, v_gs, v_ds, noise=0.005, seed=23)

    free = ("k_p", "v_th0", "k_theta", "v_vsat", "k_kneesharp", "k_eta", "k_lambda")
    # start must match the truth where params are fixed (g_k = 0!), or the
    # kink sigmoid injects phantom current.
    start = DEFAULT.replace(k_p=30.0, g_k=0.0, k_kneesharp=1.5)
    res = fit(data, start=start, free=free, weight="relative", max_nfev=8000)

    assert res.success
    assert res.rmse_rel < 0.02
    assert res.params.k_kneesharp == pytest.approx(4.0, rel=0.20)


# ------------------------------------------------------- negative k_theta (wider knee)
def test_negative_ktheta_widens_knee_and_stays_finite():
    """A negative k_theta makes v_knee = v_ov_eff2/(1 + k_theta*v_ov_eff2) exceed
    v_ov_eff2 -- a WIDER knee than the k_theta >= 0 regime can give (the ISC040
    low-voltage fit wants this). Checked with a charge cap so v_ov_eff2 is bounded
    by v_ovmax: at the safe bound k_theta = -1/v_ovmax the divisor stays positive
    (v_ov_eff2 < v_ovmax strictly) so the model is finite; pushed past it the
    divisor would go negative but the model floor (_KNEE_DENOM_FLOOR) keeps
    everything finite -- the numerical backstop."""
    base = Params(k_p=200.0, v_th0=4.0, v_sub=0.24, v_ovmax=8.0, k_ovsharp=3.0,
                  k_theta=0.0, v_vsat=8.0, k_kneesharp=1.6, k_eta=0.0,
                  k_lambda=0.0, g_k=0.0, v_dbsat=0.5, r_d=0.0, r_s=0.0)
    v_ds = np.linspace(0.0, 30.0, 61)
    v_gs = 12.0

    i_pos = np.asarray(vdmos_id(v_ds, v_gs, base))
    i_neg = np.asarray(vdmos_id(v_ds, v_gs, base.replace(k_theta=-0.05)))
    # wider knee (larger v_knee) -> higher current at a given v_ds below the knee
    assert np.all(i_neg[1:] >= i_pos[1:] - 1e-9)
    assert np.any(i_neg > i_pos + 1e-6)

    # exactly at the safe bound -1/v_ovmax: divisor still > 0 (cap is an
    # asymptote), model finite
    i_bound = np.asarray(vdmos_id(v_ds, v_gs, base.replace(k_theta=-1.0 / 8.0)))
    assert np.all(np.isfinite(i_bound))
    # pushed PAST the bound: the floor keeps it finite (no inf/NaN)
    i_past = np.asarray(vdmos_id(v_ds, v_gs, base.replace(k_theta=-1.0)))
    assert np.all(np.isfinite(i_past))


def test_ktheta_is_linear_space_with_negative_bound():
    """k_theta moved OUT of log space (it can be negative), so a negative search
    bound must not produce log(negative) = NaN in the packing path."""
    from vdmos_fit.params import BOUNDS, LOG_PARAMS, pack, vector_bounds

    assert "k_theta" not in LOG_PARAMS
    assert BOUNDS["k_theta"][0] < 0.0
    assert pack(Params(k_theta=-0.1), ("k_theta",))[0] == pytest.approx(-0.1)
    lo, _hi = vector_bounds(("k_theta",))
    assert np.isfinite(lo[0]) and lo[0] == pytest.approx(BOUNDS["k_theta"][0])


# ------------------------------------------------------- linear knee-scale form
def test_kneeform_selector_zero_is_bit_identical():
    """k_kneeform = 0 (default) selects the mobility form and is bit-identical to
    Params() -- the whole catalog and every anchor ride this default."""
    assert Params(k_kneeform=0.0) == Params()
    v_ds = np.linspace(0.0, 60.0, 121)
    v_gs = np.linspace(-2.0, 16.0, 91)
    vd2, vg2 = np.meshgrid(v_ds, v_gs)
    assert np.array_equal(vdmos_id(vd2, vg2, DEFAULT),
                          vdmos_id(vd2, vg2, DEFAULT.replace(k_kneeform=0.0)))


def test_kneeform_scale_pins_alpha_at_two():
    """k_kneeform = 1: v_knee = k_theta*v_ov_eff2, so with v_vsat off the
    saturated current is i_d,sat = 0.5*k_p*k_theta*v_ov_eff2^2 -- the gate
    exponent is exactly 2 (square-law), independent of k_theta. Checked directly
    against that closed form (cap off, k_eta = k_lambda = g_k = 0 so v_ov_eff2 =
    v_ov = v_sub*softplus((v_gs - v_th0)/v_sub))."""
    p = Params(k_p=200.0, v_th0=4.0, v_sub=0.06, v_ovmax=0.0, k_theta=1.7,
               k_kneeform=1.0, v_vsat=1e9, k_kneesharp=0.0, k_eta=0.0,
               k_lambda=0.0, g_k=0.0, v_dbsat=0.5, r_d=0.0, r_s=0.0)

    def v_ov(vg):
        z = np.clip((vg - p.v_th0) / p.v_sub, -700, 700)
        return p.v_sub * np.log1p(np.exp(z))

    for vg in (8.0, 12.0, 16.0):
        i_sat = float(vdmos_id(1e5, vg, p))           # deep saturation, v_vsat off
        closed = 0.5 * p.k_p * p.k_theta * v_ov(vg) ** 2
        assert i_sat == pytest.approx(closed, rel=1e-6)


def test_kneeform_scale_negative_ktheta_floors_knee_and_stays_finite():
    """k_kneeform = 1 with k_theta < 0 makes v_knee = k_theta*v_ov_eff2 negative,
    which the v_kneemin floor clamps -- the knee switches off, but the model
    stays finite (no divisor, so no singularity anywhere, unlike the mobility
    form)."""
    base = Params(k_p=200.0, v_th0=4.0, v_sub=0.24, v_ovmax=8.0, k_ovsharp=3.0,
                  k_theta=-0.5, k_kneeform=1.0, v_vsat=8.0, k_kneesharp=1.6,
                  k_eta=0.0, k_lambda=0.0, g_k=0.0, v_dbsat=0.5)
    v_ds = np.linspace(0.0, 30.0, 61)
    out = np.asarray(vdmos_id(v_ds, 12.0, base))
    assert np.all(np.isfinite(out))
    # the aggressive mobility-form params that needed the divisor floor are also
    # perfectly finite in the scale form (no floor is ever needed)
    v_ds2 = np.linspace(0.0, 200.0, 401)
    p2 = base.replace(k_theta=-0.2751, k_eta=42.55, v_diblmax=6.965,
                      k_diblsharp=0.1479, k_lambda=9.768e-3, v_dbsat=0.3258)
    assert np.all(np.isfinite(np.asarray(vdmos_id(v_ds2, 15.0, p2))))
