"""Tests for the gate-overdrive charge cap (v_ovmax, k_ovsharp).

    v_ov_eff = v_ov / (1 + (v_ov/v_ovmax)**k_ovsharp)**(1/k_ovsharp)    (v_ovmax = 0 -> identity)

The cap models finite mobile charge in the super-junction pillars: the measured
gate exponent alpha = d ln i_d / d ln v_ov collapses from ~1.6 to ~0.07 above
v_gs ~ 8 V on the ST8L65N044M9, and no k_theta/v_vsat/k_lambda setting can do that.

These tests pin the three properties the cap is *built* to have --
identity at v_ovmax = 0 (legacy anchors bit-identical), 0 <= dv_ov_eff/dv_ov <= 1
(the package-R solver bracket invariant dg/di_d <= -1 survives), and no
overflow in the log-space evaluation at the fitted regime k_ovsharp ~ 10-40 -- plus
the physical consequences: subthreshold untouched, no NDR with the fitted
parameters, and a synthetic round trip that recovers v_ovmax and k_ovsharp.
"""

from __future__ import annotations

import numpy as np
import pytest

from vdmos_fit import DEFAULT, Params, fit, make_synthetic, vdmos_id
from vdmos_fit.model import _soft_cap

# Best fit on the ST8L65N044M9 datasheet data (HANDOVER section 3b):
# v_ov cap + single r_d, 9 free, RMSE 6.88 %. (gds kwarg stripped in
# session 7 when the leak term was deleted; the fitted gds was 6.8e-5 S.)
FITTED = Params(
    k_p=91.9481, v_th0=4.7327, v_sub=0.06, v_ovmax=1.90909, k_ovsharp=9.86552,
    k_theta=0.0, v_vsat=1e4, k_eta=0.0, k_lambda=0.0383371,
    g_k=0.0, v_thk=1.6, v_kovsoft=0.1, v_dbsat=0.653951,
    r_d=26.03e-3, r_s=4.644e-3,
)


# ---------------------------------------------------------------- identity
def test_cap_zero_is_bit_identical():
    """Params(v_ovmax=0.0) must equal Params() and reproduce the legacy model
    bit for bit -- both on the intrinsic fast path and through the package-R
    solver. (test_model.py's anchors -- Ron 4.27 mOhm, sat table
    31.5/91.0/157.1/225.9/297.7 A, alpha = 1.078 -- ride DEFAULT.v_ovmax = 0;
    this test closes the explicit Params-level gap.)"""
    assert Params(v_ovmax=0.0) == Params()

    v_ds = np.linspace(0.0, 60.0, 121)
    v_gs = np.linspace(-2.0, 16.0, 91)
    vd2, vg2 = np.meshgrid(v_ds, v_gs)

    ref = vdmos_id(vd2, vg2, DEFAULT)
    assert np.array_equal(vdmos_id(vd2, vg2, Params(v_ovmax=0.0, k_ovsharp=37.0)), ref)

    p = DEFAULT.replace(r_d=2e-3, r_s=1e-3)
    p_off = p.replace(v_ovmax=0.0, k_ovsharp=37.0)
    assert np.array_equal(vdmos_id(vd2, vg2, p), vdmos_id(vd2, vg2, p_off))


# ---------------------------------------------------------------- cap math
def test_cap_derivative_bounded():
    """0 <= dv_ov_eff/dv_ov <= 1 everywhere.

    This is the property that keeps the package-R Newton bracket invariant
    dg/di_d <= -1 intact: the cap can only shrink the channel's di_d/dv_ds, never
    blow it up (the failure mode that made the old kappa model need bound (b)
    and a pole-avoidance dance). Checked numerically across the full fitted
    regime: v_ov in [1e-6, 50], v_ovmax in [0.1, 20], k_ovsharp in [0.5, 40].
    """
    v_ov = np.logspace(-6.0, np.log10(50.0), 300)
    for v_ovmax in (0.1, 0.5, 2.0, 20.0):
        for k_ovsharp in (0.5, 1.0, 10.0, 40.0):
            h = 1e-6 * np.maximum(v_ov, 1e-6)
            d = (_soft_cap(v_ov + h, v_ovmax, k_ovsharp)
                 - _soft_cap(v_ov - h, v_ovmax, k_ovsharp)) / (2.0 * h)
            assert np.all(d >= -1e-9), f"negative slope at v_ovmax={v_ovmax}, k_ovsharp={k_ovsharp}"
            assert np.all(d <= 1.0 + 1e-9), f"slope > 1 at v_ovmax={v_ovmax}, k_ovsharp={k_ovsharp}"


def test_cap_monotone_with_v_ovmax_asymptote():
    """Non-decreasing in v_ov, never exceeds v_ovmax, and approaches v_ovmax at the
    analytic rate. The leading term of v_ovmax - v_ov_eff is
    (v_ovmax/k_ovsharp)*(v_ovmax/v_ov)**k_ovsharp (Bernoulli:
    1-(1+u**-k_ovsharp)**(-1/k_ovsharp) <= 1/(k_ovsharp*u**k_ovsharp) for u >= 1), so the bound
    carries the 1/k_ovsharp factor -- the k_ovsharp = 0.5 case violates the k_ovsharp >= 1 constant
    v_ovmax*(v_ovmax/v_ov)**k_ovsharp by up to 1.34x at u = 10. Near the asymptote the true
    error drops below float64 resolution of v_ovmax, so the comparison floors at
    8 eps * v_ovmax (the k_ovsharp = 40 rows sit at u**-k_ovsharp ~ 1e-30 << eps)."""
    eps = np.finfo(float).eps
    for v_ovmax in (0.1, 2.0, 20.0):
        for k_ovsharp in (0.5, 8.0, 40.0):
            v_ov = np.linspace(0.0, 10.0 * v_ovmax, 501)
            capped = _soft_cap(v_ov, v_ovmax, k_ovsharp)
            assert np.all(np.diff(capped) >= -8.0 * eps * v_ovmax)
            assert np.all(capped <= v_ovmax + 8.0 * eps * v_ovmax)
            big = v_ov[v_ov >= 2.0 * v_ovmax]
            err = v_ovmax - _soft_cap(big, v_ovmax, k_ovsharp)
            # err is computed by cancellation (v_ovmax - capped), so it carries an
            # absolute rounding floor of ~eps*v_ovmax on top of the analytic bound.
            bound = (v_ovmax / k_ovsharp) * (v_ovmax / big) ** k_ovsharp + 8.0 * eps * v_ovmax
            assert np.all(err <= bound)


def test_cap_no_overflow_at_extremes():
    """The log-space form must survive the regime the naive power would not:
    v_ovmax = 0.05, k_ovsharp = 40, v_ov = 50 gives (v_ov/v_ovmax)**k_ovsharp
    = 1e120, and v_ov = 1e6
    gives ~1e292. No inf/NaN, and the result stays pinned at v_ovmax."""
    capped = _soft_cap(np.array([0.05, 1.0, 50.0]), 0.05, 40.0)
    assert np.all(np.isfinite(capped))
    assert capped[-1] == pytest.approx(0.05, rel=1e-9)

    huge = _soft_cap(np.array([1e6]), 0.05, 40.0)
    assert np.all(np.isfinite(huge))
    assert huge[0] <= 0.05 * (1.0 + 1e-12)  # log-space exp rounds ~1 ulp over

    # End to end through the model, including the package-R solver.
    p = DEFAULT.replace(v_ovmax=0.05, k_ovsharp=40.0, r_d=2e-3, r_s=1e-3)
    out = vdmos_id(np.linspace(0.0, 100.0, 201), np.linspace(-5.0, 25.0, 101)[:, None], p)
    assert np.all(np.isfinite(out))


# ---------------------------------------------------------------- physics
def test_subthreshold_slope_untouched_by_cap():
    """Deep below threshold v_ov ~ v_sub*exp((v_gs-v_th)/v_sub) << v_ovmax, so the cap
    is the identity there for ANY v_ovmax > 0 and d ln i_d/dv_gs stays 1/v_sub."""
    for v_ovmax in (0.05, 1.0, 20.0):
        p = DEFAULT.replace(g_k=0.0, v_ovmax=v_ovmax, k_ovsharp=10.0)
        vg = np.linspace(1.0, 2.4, 60)
        i_d = vdmos_id(1.0, vg, p)
        slope = np.polyfit(vg, np.log(i_d), 1)[0]
        assert slope == pytest.approx(1.0 / p.v_sub, rel=0.02)


def test_no_ndr_with_fitted_cap():
    """di_d/dv_ds > 0 for v_ds in [0, 200] at v_gs = 10 with the fitted cap and
    r_d = 26 mOhm. This is the check the kappa model could not pass by
    construction (its u(i_d) pole); the cap is safe because dv_ov_eff/dv_ov <= 1
    and every other channel term is monotone in v_ds."""
    v_ds = np.linspace(0.0, 200.0, 4001)
    i_d = vdmos_id(v_ds, 10.0, FITTED)
    assert np.all(np.isfinite(i_d))
    assert np.gradient(i_d, v_ds).min() > 0.0


# ---------------------------------------------------------------- fitting
CAPTRUTH = Params(
    k_p=30.0, v_th0=3.1, v_sub=0.05, v_ovmax=2.0, k_ovsharp=8.0,
    k_theta=0.02, v_vsat=1.2, k_eta=0.03, k_lambda=0.015,
    g_k=0.0, v_thk=1.6, v_kovsoft=0.1, v_dbsat=10.0,
)


def test_round_trip_recovers_cap():
    """Synthetic round trip: data generated with v_ovmax = 2.0, k_ovsharp = 8 must
    recover both within ~10 %. The v_gs grid deliberately reaches v_ov ~ 7 V,
    well past the cap, so the collapse in alpha is actually sampled."""
    v_gs = np.linspace(2.0, 10.0, 33)
    v_ds = np.array([0.05, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0])
    data = make_synthetic(CAPTRUTH, v_gs, v_ds, noise=0.005, seed=11)

    free = ("k_p", "v_th0", "v_sub", "v_ovmax", "k_ovsharp",
            "k_theta", "v_vsat", "k_eta", "k_lambda")
    # Non-free params come from `start`: it must match the truth where params
    # are fixed (g_k = 0!), or the kink sigmoid injects amps of phantom current.
    start = DEFAULT.replace(k_p=30.0, g_k=0.0, v_ovmax=1.0, k_ovsharp=5.0)
    res = fit(data, start=start, free=free, weight="relative", max_nfev=8000)

    assert res.success
    assert res.rmse_rel < 0.02
    assert res.params.v_ovmax == pytest.approx(2.0, rel=0.10)
    assert res.params.k_ovsharp == pytest.approx(8.0, rel=0.10)
