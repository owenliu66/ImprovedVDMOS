"""Tests for the package drain/source resistance embedding (rd, rs).

The package resistances live OUTSIDE the intrinsic channel:

    Vds' = Vds - (rd + rs)*Id     Vgs' = Vgs - rs*Id
    Id   = f_intrinsic(Vds', Vgs')   (implicit, solved per data point)

These tests validate the solver against independent root-finders and exact
limits, pin the fast path (rd = rs = 0) as bit-identical to the old model,
and check that the fitter can recover rd/rs from synthetic data.
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
    kp=30.0, vth0=3.1, theta=0.02, vsat=1.2, eta=0.03, lam=0.015,
    gds=2e-3, gk=0.3, vthk=1.5, rd=2e-3, rs=1e-3,
)


# ---------------------------------------------------------------- fast path
def test_zero_package_r_is_bit_identical_fast_path():
    """rd = rs = 0 must evaluate the OLD explicit model, bit for bit."""
    from vdmos_fit.model import _channel_current

    vds = np.linspace(0.0, 40.0, 97)
    vgs = np.linspace(0.0, 14.0, 61)
    vd2, vg2 = np.meshgrid(vds, vgs)
    assert np.array_equal(vdmos_id(vd2, vg2, DEFAULT), _channel_current(vd2, vg2, DEFAULT))


def test_level1_fast_path_unchanged():
    """The vdmos_id_level1 refactor must not move any reference number.

    Hand check: 0.5*Kp*Vov^2*(1+lam*Vds) = 0.5*28*7^2*1.4 = 960.4 A.
    """
    assert float(vdmos_id_level1(20.0, 10.0)) == pytest.approx(960.4, rel=1e-9)


# ---------------------------------------------------------------- solver correctness
def test_matches_independent_brentq_root():
    """Cross-validate the internal solver against scipy brentq (independent code)."""
    from vdmos_fit.model import _channel_current

    p = DEFAULT.replace(theta=0.0, eta=0.0, lam=0.0, gds=0.0, gk=0.0, vsat=1e9,
                        rd=2e-3, rs=1e-3)
    rtot = p.rd + p.rs
    vds, vgs = 0.01, 10.0

    def h(id_):
        return float(_channel_current(vds - rtot * id_, vgs - p.rs * id_, p)) - id_

    root = brentq(h, 0.0, vds / rtot, xtol=1e-15, rtol=8.9e-16)
    got = float(vdmos_id(vds, vgs, p))
    assert got == pytest.approx(root, rel=1e-12)


def test_solver_self_consistency_on_grid():
    """Returned Id must satisfy the fixed-point equation to <1e-8 relative."""
    from vdmos_fit.model import _channel_current

    p = DEFAULT.replace(rd=2e-3, rs=1e-3)
    rtot = p.rd + p.rs
    vd2, vg2 = np.meshgrid(np.linspace(0.0, 40.0, 101), np.linspace(0.0, 14.0, 71))
    id_ = vdmos_id(vd2, vg2, p)
    mask = id_ != 0.0
    resid = np.abs(_channel_current(vd2 - rtot * id_, vg2 - p.rs * id_, p) - id_)
    assert (resid[mask] / np.abs(id_[mask])).max() < 1e-8


def test_off_state_exact_series_fix():
    """Vgs = 0 (channel off, gk = 0): Id = gds*Vds/(1 + gds*(rd+rs)) exactly."""
    p = DEFAULT.replace(gk=0.0, rd=5e-3, rs=5e-3)
    exact = p.gds * 20.0 / (1.0 + p.gds * (p.rd + p.rs))
    assert float(vdmos_id(20.0, 0.0, p)) == pytest.approx(exact, rel=1e-6)


def test_package_dominated_limit():
    """kp -> infinity: the channel becomes ideal, Id -> Vds/(rd+rs)."""
    p = DEFAULT.replace(kp=1e9, rd=2e-3, rs=1e-3)
    assert float(vdmos_id(20.0, 10.0, p)) == pytest.approx(20.0 / 3e-3, rel=1e-4)


def test_terminal_ron_shifts_by_package_resistance():
    """Terminal Ron ~= intrinsic Ron + (rd+rs).

    The residual is slightly BELOW the naive sum (~70 uOhm here) because the
    intrinsic *chord* resistance is evaluated at the smaller internal Vds,
    where the model's knee curvature is smaller. The bound pins that physics.
    """
    p = DEFAULT.replace(rd=2e-3, rs=2e-3)
    ron_int = 0.05 / float(vdmos_id(0.05, 10.0, DEFAULT))
    ron_pkg = 0.05 / float(vdmos_id(0.05, 10.0, p))
    assert ron_pkg > ron_int + 0.5 * (p.rd + p.rs)  # clearly shifted
    assert abs(ron_pkg - (ron_int + p.rd + p.rs)) < 1.5e-4


def test_level1_package_embedding_fair_comparison():
    """Level-1 with rd/rs gets the identical embedding: Ron shifts by rd+rs."""
    ron_int = 0.05 / float(vdmos_id_level1(0.05, 5.0))
    ron_pkg = 0.05 / float(vdmos_id_level1(0.05, 5.0, rd=1e-3, rs=1e-3))
    assert ron_pkg - ron_int == pytest.approx(2e-3, abs=2e-5)


# ---------------------------------------------------------------- guards
def test_negative_vds_with_package_r_raises():
    """Channel-only model has no defined reverse characteristic -- fail loudly."""
    p = DEFAULT.replace(rd=2e-3, rs=1e-3)
    with pytest.raises(ValueError, match="Vds >= 0"):
        vdmos_id(np.linspace(-5.0, 5.0, 11), 10.0, p)
    with pytest.raises(ValueError, match="Vds >= 0"):
        vdmos_id_level1(-1.0, 5.0, rd=1e-3)


def test_negative_rd_or_rs_raises():
    from vdmos_fit.model import _solve_package_r

    z = np.zeros(1)
    with pytest.raises(ValueError, match="rd and rs must be"):
        _solve_package_r(z + 1.0, z + 5.0, -1e-3, 0.0, lambda d, g: d)


# ---------------------------------------------------------------- stall guard
def test_solver_survives_newton_two_cycle():
    """High-Vds extrapolation of the fitted super-junction set used to stall.

    Around Vds ~ 178 V with the fitted cap set (rs = 4.64 mOhm, kcap = 1.91 V,
    m = 9.87), g(Id) has two nearly flat shoulders (slope ~ -1.2) separated by
    a steep gate-degeneration cliff: rs*Id pulls the internal gate below
    threshold once Id passes ~1080 A. Newton maps each shoulder into the other,
    both endpoints stay in-bracket, and the bracket then shrinks by only ~1e-5
    per iteration -- 60 iterations were not enough and the solve raised
    RuntimeError with a ~800 A residual. The stall detector forces bisection
    after three consecutive non-halving iterations.

    The grid step matters: the stall band is ~177.8-178.7 V, so a coarse
    41-point sweep misses it entirely. 401 points (step 0.5 V) lands on it.
    """
    from vdmos_fit.model import _channel_current

    p = Params(
        kp=91.9481, vth0=4.7327, ksub=0.06, kcap=1.90909, m=9.86552,
        theta=0.0, vsat=1e4, eta=0.0, lam=0.0383371, gds=6.75611e-5,
        gk=0.0, vthk=1.6, kk=0.1, vk=0.653951, rd=26.03e-3, rs=4.644e-3,
    )
    rtot = p.rd + p.rs
    vds = np.linspace(0.0, 200.0, 401)
    ids = vdmos_id(vds, 10.0, p)  # the pre-guard solver raised here
    assert np.all(np.isfinite(ids))
    resid = np.abs(_channel_current(vds - rtot * ids, 10.0 - p.rs * ids, p) - ids)
    assert (resid / np.maximum(np.abs(ids), 1e-12)).max() < 1e-8


# ---------------------------------------------------------------- sanity sweeps
def test_monotone_and_finite_with_package_r():
    """Same convention as test_improved_model_has_no_negative_gds: the known
    ~6 % soft-knee rounding of the intrinsic model is tolerated (gradient > -2),
    and the package embedding must not make it worse."""
    p = DEFAULT.replace(rd=2e-3, rs=1e-3)
    vds = np.linspace(0.0, 40.0, 4000)
    for vg in (6.0, 10.0):
        g = np.gradient(vdmos_id(vds, vg, p), vds).min()
        assert g > -2.0
    ids = vdmos_id(np.linspace(0.0, 100.0, 500), np.linspace(-5.0, 25.0, 500)[:, None], p)
    assert np.all(np.isfinite(ids))


# ---------------------------------------------------------------- fitting
def test_rd_rs_are_fitted_by_default():
    assert "rd" in DEFAULT_FREE and "rs" in DEFAULT_FREE
    assert BOUNDS["rd"] == (0.0, 0.5) and BOUNDS["rs"] == (0.0, 0.5)


def test_fit_recovers_package_resistance():
    """Synthetic round trip: datasheet-style 2-terminal data pins rd+rs tightly;
    the rd/rs *split* leans on rs's gm-degeneration signature and is looser.

    Measured (seed 3): rd + rs recovers to 0.3 %, rd to -6.5 %, rs to +14 %,
    with the identifiability report surfacing the rd <-> rs anti-correlation.
    """
    vgs = np.linspace(2.0, 14.0, 25)
    vds = np.array([0.05, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0])
    data = make_synthetic(RTRUE, vgs, vds, noise=0.005, seed=3)

    res = fit(data, start=DEFAULT, weight="relative")
    assert res.success
    assert res.rmse_rel < 0.02

    got = res.params
    assert got.rd + got.rs == pytest.approx(3e-3, rel=0.05)
    assert got.rd == pytest.approx(2e-3, rel=0.25)
    assert got.rs == pytest.approx(1e-3, rel=0.25)
    # package R must not corrupt the channel parameters
    assert got.kp == pytest.approx(RTRUE.kp, rel=0.05)
    assert got.vth0 == pytest.approx(RTRUE.vth0, rel=0.02)

    ident = identifiability(res)
    assert ident.rank == ident.n_free
    pairs = {(a, b): c for a, b, c in ident.worst_pairs}
    assert ("rd", "rs") in pairs
    assert abs(pairs[("rd", "rs")]) > 0.9  # only the sum is strongly pinned
