"""Tests for the package drain/source resistance embedding (r_d, r_s).

The package resistances live OUTSIDE the intrinsic channel:

    v_ds' = v_ds - (r_d + r_s)*i_d     v_gs' = v_gs - r_s*i_d
    i_d   = f_intrinsic(v_ds', v_gs')   (implicit, solved per data point)

These tests validate the solver against independent root-finders and exact
limits, pin the fast path (r_d = r_s = 0) as bit-identical to the old model,
and check that the fitter can recover r_d/r_s from synthetic data.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.optimize import brentq

from vdmos_fit import (
    DEFAULT,
    Params,
    fit,
    identifiability,
    make_synthetic,
    vdmos_id,
    vdmos_id_level1,
)
from vdmos_fit.params import BOUNDS, DEFAULT_FREE

RTRUE = Params(
    k_p=30.0, v_th0=3.1, k_theta=0.02, v_vsat=1.2, k_eta=0.03, k_lambda=0.015,
    g_k=0.3, v_thk=1.5, r_d=2e-3, r_s=1e-3,
)


# ---------------------------------------------------------------- fast path
def test_zero_package_r_is_bit_identical_fast_path():
    """r_d = r_s = 0 must evaluate the OLD explicit model, bit for bit."""
    from vdmos_fit.model import _channel_current

    v_ds = np.linspace(0.0, 40.0, 97)
    v_gs = np.linspace(0.0, 14.0, 61)
    vd2, vg2 = np.meshgrid(v_ds, v_gs)
    assert np.array_equal(vdmos_id(vd2, vg2, DEFAULT), _channel_current(vd2, vg2, DEFAULT))


def test_level1_fast_path_unchanged():
    """The vdmos_id_level1 refactor must not move any reference number.

    Hand check: 0.5*Kp*v_ov^2*(1+k_lambda*v_ds) = 0.5*28*7^2*1.4 = 960.4 A.
    """
    assert float(vdmos_id_level1(20.0, 10.0)) == pytest.approx(960.4, rel=1e-9)


# ---------------------------------------------------------------- solver correctness
def test_matches_independent_brentq_root():
    """Cross-validate the internal solver against scipy brentq (independent code)."""
    from vdmos_fit.model import _channel_current

    p = DEFAULT.replace(k_theta=0.0, k_eta=0.0, k_lambda=0.0, g_k=0.0, v_vsat=1e9,
                        r_d=2e-3, r_s=1e-3)
    rtot = p.r_d + p.r_s
    v_ds, v_gs = 0.01, 10.0

    def h(i_d):
        return float(_channel_current(v_ds - rtot * i_d, v_gs - p.r_s * i_d, p)) - i_d

    root = brentq(h, 0.0, v_ds / rtot, xtol=1e-15, rtol=8.9e-16)
    got = float(vdmos_id(v_ds, v_gs, p))
    assert got == pytest.approx(root, rel=1e-12)


def test_solver_self_consistency_on_grid():
    """Returned i_d must satisfy the fixed-point equation to <1e-8 relative,
    OR below the solver's absolute tolerance (atol = 1e-12 A).

    Session 7 note: with the gds leak deleted, the deep off state sits at the
    softplus floor (~1e-14 A on this grid), below atol. There the solver's
    convergence criterion |g| <= atol + rtol*|i_d| bounds the residual only
    absolutely — demanding 1e-8 *relative* consistency at 1e-14 A would ask
    for sub-attoampere accuracy the solver never promised. The mixed rule
    below is exactly the solver contract plus the original relative bound.
    """
    from vdmos_fit.model import _channel_current

    p = DEFAULT.replace(r_d=2e-3, r_s=1e-3)
    rtot = p.r_d + p.r_s
    vd2, vg2 = np.meshgrid(np.linspace(0.0, 40.0, 101), np.linspace(0.0, 14.0, 71))
    i_d = vdmos_id(vd2, vg2, p)
    mask = i_d != 0.0
    resid = np.abs(_channel_current(vd2 - rtot * i_d, vg2 - p.r_s * i_d, p) - i_d)
    assert np.all(resid[mask] <= 1e-12 + 1e-8 * np.abs(i_d[mask]))


def test_off_state_solves_to_channel_floor():
    """v_gs = 0 with package R: i_d is the tiny gate-dependent channel floor.

    The old exactness check (i_d = gds*v_ds/(1 + gds*(r_d+r_s))) died with the gds
    term (session 7). The solver must still converge self-consistently at the
    softplus floor: ~1e-13 A at (v_ds, v_gs) = (20, 0), so the package-R voltage
    drops are ~1e-15 V and the terminal and internal currents agree.
    """
    from vdmos_fit.model import _channel_current

    p = DEFAULT.replace(g_k=0.0, r_d=5e-3, r_s=5e-3)
    i_d = float(vdmos_id(20.0, 0.0, p))
    assert 0.0 < i_d < 1e-9
    assert i_d == pytest.approx(float(_channel_current(20.0, 0.0, p)), rel=1e-6)


def test_package_dominated_limit():
    """k_p -> infinity: the channel becomes ideal, i_d -> v_ds/(r_d+r_s)."""
    p = DEFAULT.replace(k_p=1e9, r_d=2e-3, r_s=1e-3)
    assert float(vdmos_id(20.0, 10.0, p)) == pytest.approx(20.0 / 3e-3, rel=1e-4)


def test_terminal_ron_shifts_by_package_resistance():
    """Terminal Ron ~= intrinsic Ron + (r_d+r_s).

    The residual is slightly BELOW the naive sum (~70 uOhm here) because the
    intrinsic *chord* resistance is evaluated at the smaller internal v_ds,
    where the model's knee curvature is smaller. The bound pins that physics.
    """
    p = DEFAULT.replace(r_d=2e-3, r_s=2e-3)
    ron_int = 0.05 / float(vdmos_id(0.05, 10.0, DEFAULT))
    ron_pkg = 0.05 / float(vdmos_id(0.05, 10.0, p))
    assert ron_pkg > ron_int + 0.5 * (p.r_d + p.r_s)  # clearly shifted
    assert abs(ron_pkg - (ron_int + p.r_d + p.r_s)) < 1.5e-4


def test_level1_package_embedding_fair_comparison():
    """Level-1 with r_d/r_s gets the identical embedding: Ron shifts by r_d+r_s."""
    ron_int = 0.05 / float(vdmos_id_level1(0.05, 5.0))
    ron_pkg = 0.05 / float(vdmos_id_level1(0.05, 5.0, r_d=1e-3, r_s=1e-3))
    assert ron_pkg - ron_int == pytest.approx(2e-3, abs=2e-5)


# ---------------------------------------------------------------- guards
def test_negative_vds_with_package_r_raises():
    """Channel-only model has no defined reverse characteristic -- fail loudly."""
    p = DEFAULT.replace(r_d=2e-3, r_s=1e-3)
    with pytest.raises(ValueError, match="v_ds >= 0"):
        vdmos_id(np.linspace(-5.0, 5.0, 11), 10.0, p)
    with pytest.raises(ValueError, match="v_ds >= 0"):
        vdmos_id_level1(-1.0, 5.0, r_d=1e-3)


def test_negative_rd_or_rs_raises():
    from vdmos_fit.model import _solve_package_r

    z = np.zeros(1)
    with pytest.raises(ValueError, match="r_d and r_s must be"):
        _solve_package_r(z + 1.0, z + 5.0, -1e-3, 0.0, lambda d, g: d)


# ---------------------------------------------------------------- stall guard
def test_solver_survives_newton_two_cycle():
    """High-v_ds extrapolation of the fitted super-junction set used to stall.

    Around v_ds ~ 178 V with the fitted cap set (r_s = 4.64 mOhm, v_ovmax = 1.91 V,
    k_ovsharp = 9.87), g(i_d) has two nearly flat shoulders (slope ~ -1.2) separated by
    a steep gate-degeneration cliff: r_s*i_d pulls the internal gate below
    threshold once i_d passes ~1080 A. Newton maps each shoulder into the other,
    both endpoints stay in-bracket, and the bracket then shrinks by only ~1e-5
    per iteration -- 60 iterations were not enough and the solve raised
    RuntimeError with a ~800 A residual. The stall detector forces bisection
    after three consecutive non-halving iterations.

    The cliff is gate-degeneration physics, NOT the deleted gds leak: the same
    parameter set without gds still lands in the stall band (i_d(178) ~ 789 A)
    and only passes with the guard. The gds kwarg was stripped from the set in
    session 7; the pinned values below are the session-4 fit minus gds.

    The grid step matters: the stall band is ~177.8-178.7 V, so a coarse
    41-point sweep misses it entirely. 401 points (step 0.5 V) lands on it.
    """
    from vdmos_fit.model import _channel_current

    p = Params(
        k_p=91.9481, v_th0=4.7327, v_sub=0.06, v_ovmax=1.90909, k_ovsharp=9.86552,
        k_theta=0.0, v_vsat=1e4, k_eta=0.0, k_lambda=0.0383371,
        g_k=0.0, v_thk=1.6, v_kovsoft=0.1, v_dbsat=0.653951, r_d=26.03e-3, r_s=4.644e-3,
    )
    rtot = p.r_d + p.r_s
    v_ds = np.linspace(0.0, 200.0, 401)
    i_d = vdmos_id(v_ds, 10.0, p)  # the pre-guard solver raised here
    assert np.all(np.isfinite(i_d))
    resid = np.abs(_channel_current(v_ds - rtot * i_d, 10.0 - p.r_s * i_d, p) - i_d)
    assert (resid / np.maximum(np.abs(i_d), 1e-12)).max() < 1e-8


# ---------------------------------------------------------------- sanity sweeps
def test_monotone_and_finite_with_package_r():
    """Same convention as test_improved_model_has_no_negative_gds: the known
    ~6 % soft-knee rounding of the intrinsic model is tolerated (gradient > -2),
    and the package embedding must not make it worse."""
    p = DEFAULT.replace(r_d=2e-3, r_s=1e-3)
    v_ds = np.linspace(0.0, 40.0, 4000)
    for vg in (6.0, 10.0):
        g = np.gradient(vdmos_id(v_ds, vg, p), v_ds).min()
        assert g > -2.0
    i_d = vdmos_id(np.linspace(0.0, 100.0, 500), np.linspace(-5.0, 25.0, 500)[:, None], p)
    assert np.all(np.isfinite(i_d))


# ---------------------------------------------------------------- fitting
def test_rd_rs_are_fitted_by_default():
    assert "r_d" in DEFAULT_FREE and "r_s" in DEFAULT_FREE
    assert BOUNDS["r_d"] == (0.0, 0.5) and BOUNDS["r_s"] == (0.0, 0.5)


def test_fit_recovers_package_resistance():
    """Synthetic round trip: datasheet-style 2-terminal data pins r_d+r_s tightly;
    the r_d/r_s *split* leans on r_s's gm-degeneration signature and is looser.

    Measured (seed 3): r_d + r_s recovers to 0.3 %, r_d to -6.5 %, r_s to +14 %,
    with the identifiability report surfacing the r_d <-> r_s anti-correlation.
    """
    v_gs = np.linspace(2.0, 14.0, 25)
    v_ds = np.array([0.05, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0])
    data = make_synthetic(RTRUE, v_gs, v_ds, noise=0.005, seed=3)

    res = fit(data, start=DEFAULT, weight="relative")
    assert res.success
    assert res.rmse_rel < 0.02

    got = res.params
    assert got.r_d + got.r_s == pytest.approx(3e-3, rel=0.05)
    assert got.r_d == pytest.approx(2e-3, rel=0.25)
    assert got.r_s == pytest.approx(1e-3, rel=0.25)
    # package R must not corrupt the channel parameters
    assert got.k_p == pytest.approx(RTRUE.k_p, rel=0.05)
    assert got.v_th0 == pytest.approx(RTRUE.v_th0, rel=0.02)

    ident = identifiability(res)
    assert ident.rank == ident.n_free
    pairs = {(a, b): c for a, b, c in ident.worst_pairs}
    assert ("r_d", "r_s") in pairs
    assert abs(pairs[("r_d", "r_s")]) > 0.9  # only the sum is strongly pinned
