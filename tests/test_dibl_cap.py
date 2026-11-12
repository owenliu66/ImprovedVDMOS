"""Tests for the tunable DIBL-shift soft cap (v_diblmax, k_diblsharp).

    v_th = v_th0 - softcap(k_eta*V_ds*k_db ; v_diblmax, k_diblsharp)
    softcap(x; c, n) = x / (1 + (x/c)**n)**(1/n)      (c = 0 -> identity)

The raw DIBL shift k_eta*V_ds*k_db grows without bound in V_ds, so a fit
anchored on low-voltage data extrapolates an aggressive linear DIBL that drives
v_th negative at high V_ds (the ISC040N10NM8 100 V failure mode: the uncapped
record hits v_th = 0 near 140 V and off-state leakage runs away to amps). The
cap -- the SAME power-law soft-min as the gate-charge cap -- ceilings the shift
at v_diblmax, so v_th >= v_th0 - v_diblmax for every V_ds.

These tests pin the three properties the cap is built to have --
identity at v_diblmax = 0 (legacy anchors bit-identical), monotonicity of the
shift (0 <= dshift/draw <= 1, so v_th stays non-increasing and the package-R
bracket proof survives), and no overflow in the log-space evaluation at the
fitted regime k_diblsharp ~ 10-40 -- plus the physical consequence: a capped
fit keeps v_th positive and off-state leakage bounded at high V_ds where the
uncapped fit runs away.
"""

from __future__ import annotations

import numpy as np
import pytest

from vdmos_fit import DEFAULT, Params, vdmos_id
from vdmos_fit.model import _soft_cap


# ---------------------------------------------------------------- identity
def test_dibl_cap_zero_is_bit_identical():
    """Params(v_diblmax=0.0) must equal Params() and reproduce the pre-cap
    model bit for bit -- intrinsic fast path and through the package-R solver."""
    assert Params(v_diblmax=0.0) == Params()

    v_ds = np.linspace(0.0, 200.0, 201)
    v_gs = np.linspace(-2.0, 16.0, 91)
    vd2, vg2 = np.meshgrid(v_ds, v_gs)

    ref = vdmos_id(vd2, vg2, DEFAULT)
    assert np.array_equal(
        vdmos_id(vd2, vg2, Params(v_diblmax=0.0, k_diblsharp=37.0)), ref)

    p = DEFAULT.replace(r_d=2e-3, r_s=1e-3)
    p_off = p.replace(v_diblmax=0.0, k_diblsharp=37.0)
    assert np.array_equal(vdmos_id(vd2, vg2, p), vdmos_id(vd2, vg2, p_off))


# ---------------------------------------------------------------- cap math
def test_dibl_shift_monotone_in_v_ds():
    """The capped shift is non-decreasing in V_ds (soft-min derivative
    0 <= dshift/draw <= 1), so v_th stays non-increasing and every channel
    partial stays >= 0 (the package-R bracket invariant). Checked across the
    fitted regime v_diblmax in [0.05, 20], k_diblsharp in [0.5, 40]."""
    v_ds = np.linspace(0.0, 200.0, 2001)
    for v_diblmax in (0.05, 1.0, 5.0, 20.0):
        for k_diblsharp in (0.5, 1.0, 10.0, 40.0):
            p = DEFAULT.replace(k_eta=0.05, v_dbsat=1.0,
                                v_diblmax=v_diblmax, k_diblsharp=k_diblsharp)
            k_db = np.tanh(v_ds / p.v_dbsat)
            raw = p.k_eta * v_ds * k_db
            shift = _soft_cap(raw, v_diblmax, k_diblsharp)
            assert np.all(np.diff(shift) >= -1e-12)
            # shift never exceeds the asymptote v_diblmax
            assert np.all(shift <= v_diblmax + 1e-9)
            # and never exceeds the raw (uncapped) shift
            assert np.all(shift <= raw + 1e-12)


def test_dibl_cap_derivative_bounded():
    """0 <= d(softcap)/d(raw) <= 1 everywhere -- the property the package-R
    bracket proof relies on (the cap can only shrink the DIBL sensitivity to
    V_ds, never amplify it)."""
    raw = np.logspace(-6.0, np.log10(50.0), 300)
    for v_diblmax in (0.05, 1.0, 5.0, 20.0):
        for k_diblsharp in (0.5, 1.0, 10.0, 40.0):
            h = 1e-6 * np.maximum(raw, 1e-6)
            d = (_soft_cap(raw + h, v_diblmax, k_diblsharp)
                 - _soft_cap(raw - h, v_diblmax, k_diblsharp)) / (2.0 * h)
            assert np.all(d >= -1e-9)
            assert np.all(d <= 1.0 + 1e-9)


def test_dibl_cap_no_overflow_at_extremes():
    """Log-space form survives the regime the naive power would not:
    v_diblmax = 0.05, k_diblsharp = 40, raw = 50 -> (raw/v_diblmax)**40 ~ 1e120.
    End to end through the model including the package-R solver."""
    capped = _soft_cap(np.array([0.05, 1.0, 50.0]), 0.05, 40.0)
    assert np.all(np.isfinite(capped))
    assert capped[-1] == pytest.approx(0.05, rel=1e-9)

    p = DEFAULT.replace(k_eta=0.5, v_dbsat=1.0, v_diblmax=0.05,
                        k_diblsharp=40.0, r_d=2e-3, r_s=1e-3)
    out = vdmos_id(np.linspace(0.0, 400.0, 401),
                   np.linspace(-5.0, 25.0, 101)[:, None], p)
    assert np.all(np.isfinite(out))


# ---------------------------------------------------------------- physics
def test_dibl_cap_keeps_vth_positive_and_leakage_bounded():
    """The cap's whole point: with the shift ceiling v_diblmax < v_th0, v_th
    stays positive for every V_ds and off-state leakage (V_gs = 0) stays
    bounded, whereas the uncapped fit runs away. Uses a deliberately aggressive
    k_eta (the low-voltage-fit failure mode) at v_th0 = 4 V."""
    v_ds = np.linspace(0.0, 400.0, 801)

    p_uncapped = Params(k_p=200.0, v_th0=4.0, v_sub=0.24, v_ovmax=5.8,
                        k_ovsharp=2.2, k_theta=0.0, v_vsat=1e4, k_kneesharp=1.8,
                        k_eta=0.05, k_lambda=0.0, g_k=0.0, v_dbsat=0.5,
                        v_diblmax=0.0)
    p_capped = p_uncapped.replace(v_diblmax=1.5, k_diblsharp=10.0)

    def vth(pp, v):
        k_db = np.tanh(v / pp.v_dbsat)
        return pp.v_th0 - _soft_cap(pp.k_eta * v * k_db,
                                    pp.v_diblmax, pp.k_diblsharp)

    # uncapped: v_th goes negative at high V_ds (the failure mode)
    assert vth(p_uncapped, 400.0) < 0.0
    # capped: v_th >= v_th0 - v_diblmax > 0 everywhere
    assert vth(p_capped, 400.0) == pytest.approx(4.0 - 1.5, rel=1e-3)
    assert all(vth(p_capped, v) > 0.0 for v in (100.0, 200.0, 400.0))

    # off-state leakage at V_gs = 0: bounded for the capped model, runs away
    # for the uncapped one (v_th < 0 -> the channel spuriously turns on)
    i_unc = np.asarray(vdmos_id(v_ds, 0.0, p_uncapped))
    i_cap = np.asarray(vdmos_id(v_ds, 0.0, p_capped))
    assert np.all(np.isfinite(i_cap))
    assert i_cap.max() < 1e-3
    assert i_unc.max() > 1.0   # runaway: amps of phantom off-current
