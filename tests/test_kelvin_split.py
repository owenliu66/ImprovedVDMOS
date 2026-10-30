"""Tests for the Kelvin-aware source split (r_sk, names.xlsx "sk").

Terminal equations (module docstring of model.py):

    v_ds' = v_ds - (r_d + r_s)*i_d        current loop: ALL of i_d
    v_gs' = v_gs - (r_s - r_sk)*i_d       gate loop: returns at the KELVIN tap

r_sk is the resistance between the kelvin tap and the power source pin (the
downstream segment of r_s, carrying only gate current -- zero in the DC
model). Pins: r_sk = 0 bit-identical to the legacy r_s semantics; r_sk = r_s
(gate loop clean) equals the r_s = 0 model with r_d -> r_d + r_s; the
intermediate split is self-consistent with the (r_d + r_s, r_s - r_sk)
fixed point and matches independent brentq roots; the 0 <= r_sk <= r_s guard;
the bracket survives r_g > 0 with the trap-filling gain active (the IMCQ
record's dv_ov_eff2/dv_ov_eff = 12.5 regime); and the fitting-space
decisions (linear space, bounds, not in DEFAULT_FREE).
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.optimize import brentq

from vdmos_fit import (
    BOUNDS,
    DEFAULT,
    DEFAULT_FREE,
    LOG_PARAMS,
    Params,
    vdmos_id,
)
from vdmos_fit.model import _channel_current

# IMCQ-like set: trap-filling gain ACTIVE (gate gain 1 + 2*k_bt*v_ov >> 1)
# -- the regime the HANDOVER says must be re-checked before any r_s > 0 use.
GAINED = Params(
    k_p=9.3, v_th0=5.75, v_sub=0.79, k_beta_trap=0.47, k_theta=0.109,
    v_vsat=1e4, k_kneesharp=1.46, k_eta=0.05, k_lambda=0.0, g_k=0.0,
    v_dbsat=320.0, r_d=2.0e-3, r_s=1.0e-3, r_sk=0.4e-3,
)


# ---------------------------------------------------------------- identity
def test_r_sk_zero_bit_identical_legacy():
    """r_sk = 0 must reproduce the legacy r_s semantics bit for bit --
    Params equality and the full grid through the package-R solver."""
    assert Params(r_sk=0.0) == Params()
    assert DEFAULT.r_sk == 0.0

    legacy = Params(k_p=91.9, v_th0=4.73, v_sub=0.06, v_ovmax=1.909,
                    k_ovsharp=9.87, k_theta=0.0, v_vsat=1e4, k_eta=0.0,
                    k_lambda=0.0383, g_k=0.0, v_dbsat=0.654,
                    r_d=26.03e-3, r_s=4.644e-3)
    split = legacy.replace(r_sk=0.0)
    vd2, vg2 = np.meshgrid(np.linspace(0.0, 200.0, 201),
                           np.linspace(0.0, 14.0, 71))
    assert np.array_equal(vdmos_id(vd2, vg2, legacy),
                          vdmos_id(vd2, vg2, split))


def test_ideal_kelvin_equals_clean_gate_loop():
    """r_sk = r_s (ideal kelvin tap): the gate loop is NOT degenerated, so
    the model must equal the r_s = 0 model with the whole source path moved
    into the drain leg -- bit-identical (identical solver path: same r_tot,
    r_g = 0)."""
    p_k = Params(k_p=30.0, v_th0=3.1, v_sub=0.05, k_theta=0.02, v_vsat=1.2,
                 k_eta=0.03, k_lambda=0.015, g_k=0.0, v_dbsat=10.0,
                 r_d=2e-3, r_s=1e-3, r_sk=1e-3)
    p_d = p_k.replace(r_d=p_k.r_d + p_k.r_s, r_s=0.0, r_sk=0.0)
    vd2, vg2 = np.meshgrid(np.linspace(0.0, 40.0, 101),
                           np.linspace(0.0, 14.0, 71))
    assert np.array_equal(vdmos_id(vd2, vg2, p_k), vdmos_id(vd2, vg2, p_d))


# ---------------------------------------------------------------- split physics
def test_split_self_consistency_and_brentq():
    """Intermediate split: the returned i_d satisfies the (r_d + r_s,
    r_s - r_sk) fixed point, and matches independent brentq roots."""
    p = GAINED
    r_tot = p.r_d + p.r_s
    r_g = p.r_s - p.r_sk

    vd2, vg2 = np.meshgrid(np.linspace(0.0, 40.0, 101),
                           np.linspace(0.0, 20.0, 81))
    i_d = vdmos_id(vd2, vg2, p)
    assert np.all(np.isfinite(i_d)) and np.all(i_d >= 0.0)
    resid = np.abs(_channel_current(vd2 - r_tot * i_d, vg2 - r_g * i_d, p) - i_d)
    mask = i_d != 0.0
    assert np.all(resid[mask] <= 1e-12 + 1e-8 * np.abs(i_d[mask]))

    for v_ds, v_gs in ((0.05, 18.0), (5.0, 12.0), (20.0, 18.0), (40.0, 14.0)):
        def h(i_d_, v_ds=v_ds, v_gs=v_gs):  # bind loop vars (ruff B023)
            return float(_channel_current(
                v_ds - r_tot * i_d_, v_gs - r_g * i_d_, p)) - i_d_

        root = brentq(h, 0.0, v_ds / r_tot, xtol=1e-15, rtol=8.9e-16)
        assert float(vdmos_id(v_ds, v_gs, p)) == pytest.approx(root, rel=1e-10)


def test_r_sk_lifts_gate_degeneration_monotonically():
    """At fixed r_d and r_s, raising r_sk (moving the kelvin tap toward the
    die) can only INCREASE the terminal current: the gate loop sees a smaller
    degeneration r_s - r_sk. Strictly increasing while r_sk < r_s."""
    base = GAINED
    v_ds, v_gs = 20.0, 18.0
    currents = [float(vdmos_id(v_ds, v_gs, base.replace(r_sk=r_sk)))
                for r_sk in (0.0, 0.25e-3, 0.5e-3, 0.75e-3, 1.0e-3)]
    assert np.all(np.diff(currents) > 0.0)
    # endpoint: r_sk = r_s equals the clean-gate reference
    clean = float(vdmos_id(v_ds, v_gs,
                           base.replace(r_s=0.0, r_sk=0.0,
                                        r_d=base.r_d + base.r_s)))
    assert currents[-1] == pytest.approx(clean, rel=1e-12)


def test_terminal_gm_degeneration_knob():
    """Terminal g_fs follows g_m,eff = g_m/(1 + (r_s - r_sk)*g_m): with the
    gate loop clean (r_sk = r_s) the model is bit-identical to moving the
    whole source path into the drain leg (same r_tot, same output-conductance
    feedback); with r_sk = 0 the legacy loop is clearly degenerated."""
    p_clean = GAINED.replace(r_d=2.0e-3, r_s=1.0e-3, r_sk=1.0e-3)
    p_drain = GAINED.replace(r_d=3.0e-3, r_s=0.0, r_sk=0.0)  # same r_tot
    p_leg = GAINED.replace(r_d=2.0e-3, r_s=1.0e-3, r_sk=0.0)

    def gfs(p):
        v_gs = 12.0
        h = 1e-3
        return (float(vdmos_id(20.0, v_gs + h, p))
                - float(vdmos_id(20.0, v_gs - h, p))) / (2 * h)

    gcl, gdr, gleg = gfs(p_clean), gfs(p_drain), gfs(p_leg)
    # clean gate loop: the r_s/r_d split is invisible -- bit-identical solve
    assert gcl == gdr
    # legacy loop (full r_s degeneration) must be clearly lower
    assert gleg < 0.97 * gcl


# ---------------------------------------------------------------- guards
def test_r_sk_out_of_range_raises():
    """0 <= r_sk <= r_s is a hard guard at the vdmos_id entry (raise, not
    clamp -- k_beta_trap precedent)."""
    with pytest.raises(ValueError, match="r_sk must satisfy"):
        vdmos_id(5.0, 10.0, DEFAULT.replace(r_s=1e-3, r_sk=2e-3))   # r_sk > r_s
    with pytest.raises(ValueError, match="r_sk must satisfy"):
        vdmos_id(5.0, 10.0, DEFAULT.replace(r_sk=-1e-6))            # negative
    # r_sk == r_s exactly is legal (ideal kelvin tap)
    assert np.isfinite(vdmos_id(5.0, 10.0, DEFAULT.replace(r_s=1e-3, r_sk=1e-3)))


def test_solver_rejects_negative_r_g():
    """Direct solver callers: r_g < 0 breaks the bracket proof -> ValueError."""
    from vdmos_fit.model import _solve_package_r

    z = np.zeros(1)
    with pytest.raises(ValueError, match="r_g"):
        _solve_package_r(z + 1.0, z + 5.0, 1e-3, 1e-3,
                         lambda d, g: d, r_g=-1e-4)


# ---------------------------------------------------------------- bracket + gain
def test_bracket_survives_gain_with_partial_degeneration():
    """The HANDOVER re-check: trap-filling gain active (gate gain up to ~12)
    with r_g = r_s - r_sk > 0. The proof needs only monotonicity:
    dg/di = -(r_d+r_s)*fd - r_g*fg - 1 <= -1. Dense grid stays finite,
    monotone in v_gs, and inside the current-loop bracket i_d <= v_ds/(r_d+r_s)."""
    p = GAINED
    r_tot = p.r_d + p.r_s
    vd2, vg2 = np.meshgrid(np.linspace(0.0, 60.0, 241),
                           np.linspace(0.0, 20.0, 101))
    i_d = vdmos_id(vd2, vg2, p)
    assert np.all(np.isfinite(i_d))
    assert np.all(i_d <= vd2 / r_tot + 1e-9)
    # column 0 is v_ds = 0 where f(0,.) = 0 exactly (all-zero slice); the
    # strict-increase check applies to the forward-polarity columns.
    for col in range(1, vd2.shape[1]):
        assert np.all(np.diff(i_d[:, col]) > 0.0)


# ---------------------------------------------------------------- fitting space
def test_r_sk_fitting_space_decisions():
    """Linear space (0 is the legitimate "no kelvin contact" value), bounds
    (0.0, 0.5) like r_d/r_s, and NOT in DEFAULT_FREE (freeing it silently
    changes existing fits)."""
    assert "r_sk" not in LOG_PARAMS
    assert "r_sk" not in DEFAULT_FREE
    assert BOUNDS["r_sk"] == (0.0, 0.5)
