"""Tests for parameter extraction: can the fitter recover known parameters?"""

from __future__ import annotations

import math

import numpy as np
import pytest

from vdmos_fit import (
    BOUNDS,
    DEFAULT,
    CurveSet,
    FitResult,
    Params,
    example_dataset,
    fit,
    identifiability,
    make_synthetic,
    parameter_uncertainty,
    sensitivity_report,
    vdmos_id,
)

TRUTH = Params(
    k_p=30.0, v_th0=3.1, k_theta=0.02, v_vsat=1.2, k_eta=0.03, k_lambda=0.015, g_k=0.3, v_thk=1.5
)


@pytest.fixture(scope="module")
def synthetic():
    v_gs = np.linspace(2.0, 14.0, 25)
    v_ds = np.array([0.05, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0])
    return make_synthetic(TRUTH, v_gs, v_ds, noise=0.005, seed=3)


def test_fit_recovers_synthetic_parameters(synthetic):
    """Noiseless-ish synthetic data must round-trip back to the truth.

    Only the strongly-sensed parameters are checked tightly; k_theta and v_thk are
    weakly identified from this grid (see the identifiability tests) and are
    deliberately not asserted here.
    """
    res = fit(synthetic, start=DEFAULT, weight="relative")
    assert res.success
    assert res.rmse_rel < 0.02
    got, want = res.params.as_dict(), TRUTH.as_dict()
    for name in ("k_p", "v_th0", "v_vsat", "k_eta", "k_lambda", "g_k"):
        assert got[name] == pytest.approx(want[name], rel=0.2), name


def test_fit_improves_on_default_start(synthetic):
    res = fit(synthetic, start=DEFAULT)
    start_err = np.sqrt(
        np.mean(
            (
                (vdmos_id(synthetic.v_ds, synthetic.v_gs, DEFAULT) - synthetic.i_d)
                / np.maximum(np.abs(synthetic.i_d), 1e-3)
            )
            ** 2
        )
    )
    assert res.rmse_rel < 0.2 * start_err


def test_bounds_are_respected(synthetic):
    """A deliberately absurd start must still land inside physical bounds.
    k_theta is now SIGNED (linear space, may be negative -> wider knee), so the
    check is 'inside the box', not 'positive' like the strictly-positive params."""
    bad = DEFAULT.replace(k_p=1e-4, k_theta=1e-4, v_vsat=1e-3)
    res = fit(synthetic, start=bad)
    lo, hi = BOUNDS["k_theta"]
    assert res.params.k_p > 0 and res.params.v_vsat > 0
    assert lo <= res.params.k_theta <= hi


def test_fixed_parameters_stay_fixed(synthetic):
    """Fitting only k_p/v_th0 must not move k_theta."""
    res = fit(synthetic, start=DEFAULT, free=("k_p", "v_th0"))
    assert res.params.k_theta == DEFAULT.k_theta
    assert res.params.v_vsat == DEFAULT.v_vsat


def test_uncertainty_report_shape(synthetic):
    res = fit(synthetic)
    unc = parameter_uncertainty(res)
    assert set(unc) == set(
        ("k_p", "v_th0", "k_theta", "v_vsat", "k_eta", "k_lambda", "g_k", "v_thk", "r_d", "r_s")
    )
    assert all(v >= 0 for v in unc.values())


def test_kink_parameters_are_weakly_determined(synthetic):
    """The kink onset is barely sampled by this grid, so v_thk is poorly determined.

    v_thk = 1.5 V with v_kovsoft = 0.1 V means the kink sigmoid is essentially saturated
    (=1) for every v_gs >= 2 V in this grid, so only the very lowest curves carry
    information about it. The fit lands at ~1.7 V against a truth of 1.5 V -- a
    16 % error -- and the report must show a correspondingly huge sigma rather
    than hiding it. Fix v_thk from a measurement, or add dense v_gs points through
    the kink region, before trusting it.
    """
    ident = identifiability(fit(synthetic))
    assert ident.sigma_rel["v_thk"] > 1.0  # order-unity or worse
    assert ident.sigma_rel["k_p"] < 1.0  # core shape params are fine
    assert ident.rank == ident.n_free  # not degenerate, just weak
    # Threshold raised 1e4 -> 1e5 (2026, package-R change): with r_d/r_s free and
    # truth r_d = r_s = 0, the r_d <-> r_s columns are near-parallel (corr -0.999,
    # only the sum is sensed), which alone pushes cond to ~1.4e4.
    assert ident.cond < 1e5


def test_sensitivity_report_ranks_parameters(synthetic):
    """Pre-fit diagnostic must agree with the post-fit Jacobian diagnosis.

    Sensitivity is measured as mean |d ln i_d / d ln p|, so it is directly
    comparable across parameters with different units.
    """
    sens = sensitivity_report(DEFAULT, synthetic)
    assert sens["k_p"] > 0.5  # strongly sensed
    assert sens["v_thk"] < 1e-3  # ~2000x weaker than k_p
    assert sens["v_thk"] < 1e-3 * sens["k_p"]
    # weakly sensed parameters must be the ones with the biggest error bars
    ident = identifiability(fit(synthetic))
    weakest = min(sens, key=sens.get)
    assert ident.sigma_rel[weakest] > ident.sigma_rel["v_th0"]


def test_theta_vsat_correlation_is_reported(synthetic):
    """k_theta and v_vsat both flatten the knee; output curves alone cannot separate
    them, so the report must surface the pairing.

    Evaluated with the classic 9-parameter free set: once r_d/r_s join the fit,
    the knee-flattening trade-off redistributes across more parameters and the
    single k_theta<->v_vsat pair drops below the 0.8 reporting threshold (their
    individual sigmas stay large -- see
    tests/test_package_r.py::test_fit_recovers_package_resistance).
    """
    legacy_free = ("k_p", "v_th0", "k_theta", "v_vsat", "k_eta", "k_lambda", "g_k", "v_thk")
    ident = identifiability(fit(synthetic, free=legacy_free), legacy_free)
    pairs = {(a, b) for a, b, _ in ident.worst_pairs}
    assert ("k_theta", "v_vsat") in pairs
    # and both carry large individual error bars as a consequence. k_theta is
    # fitted in LINEAR space (it may be negative), so its sigma_rel is an
    # ABSOLUTE uncertainty on a ~0.02 1/V parameter -- poorly determined once it
    # is comparable to the value scale itself (was "> 1" in the old log space).
    assert ident.sigma_rel["k_theta"] > 0.02


def test_dead_parameter_reports_inf_not_a_small_sigma():
    """A parameter with (near-)zero sensitivity must report sigma_rel = inf.

    This is the trap `identifiability()` exists to avoid: a bare
    pinv(J.T @ J) returns *zero* variance for a near-zero Jacobian column,
    which reads as "perfectly determined" when the truth is the exact
    opposite. Reproduced synthetically on the ST8L65 situation: a dataset
    whose lowest v_gs is 8 V (v_th0 = 3.1) contains no subthreshold points, so
    v_sub carries essentially no information -- its column measures ~1e-9 of
    the largest column, and pinning v_sub at 0.30 V leaves the fit at the same
    optimum. The old absolute 1e-10 * max tolerance classified that column
    LIVE and reported a confident small sigma; the relative 1e-5 tolerance
    classifies it dead. k_eta (sigma ~ 12, genuinely sensed but weak) stays
    live, so the tolerance does not over-trigger.
    """
    truth = Params(k_p=30.0, v_th0=3.1, v_sub=0.05, k_theta=0.02, v_vsat=1.2,
                   k_eta=0.03, k_lambda=0.015, g_k=0.0)
    v_gs = np.linspace(8.0, 14.0, 25)
    v_ds = np.array([0.5, 2.0, 5.0, 10.0, 20.0])
    data = make_synthetic(truth, v_gs, v_ds, noise=0.0, seed=0)
    free = ("k_p", "v_th0", "v_sub", "k_theta", "v_vsat", "k_eta", "k_lambda")

    ident = identifiability(fit(data, start=truth, free=free), free)
    assert math.isinf(ident.sigma_rel["v_sub"])
    assert ident.unidentifiable == ("v_sub",)
    # rank drops by at least the dead column; k_eta <-> others are also near-
    # degenerate on this grid, so rank is not pinned exactly.
    assert ident.rank <= len(free) - 1
    # the near-dead column must not fake a correlation with a real parameter
    assert not any("v_sub" in (a, b) for a, b, _ in ident.worst_pairs)
    # genuinely sensed, merely weak -> must NOT be called dead
    assert not math.isinf(ident.sigma_rel["k_eta"])

    # Independent confirmation: pinning the dead parameter cannot change the fit.
    pinned = fit(data, start=truth.replace(v_sub=0.30),
                 free=tuple(n for n in free if n != "v_sub"))
    assert pinned.params.k_p == pytest.approx(30.0, rel=1e-3)


# ---------------------------------------------------- on_bounds pinned test
def _fake_fit(p: Params, free: tuple[str, ...], bounds: dict) -> FitResult:
    """FitResult without running a fit -- on_bounds only reads these fields."""
    return FitResult(
        params=p, cost=0.0, rmse_rel=0.0, n_points=1, n_free=len(free),
        success=True, message="", free=free, bounds_used=bounds,
    )


def test_on_bounds_not_flagged_well_inside_wide_box():
    """Regression (session 6 false positive): the pinned test must scale with
    the parameter's own magnitude, NOT with the width of the search box.

    The deleted gds bound (1e-9, 10) S had span 10 S, so the old tolerance
    rel_tol*span = 1e-2 S called any gds below 0.01 S "pinned LOWER" -- even
    the physically meaningful 6.8e-5 S. g_k inherited the same wide box (its
    bound spans 14 decades), so it reproduces the bug exactly: 6.8e-5 S must
    read as an interior value.
    """
    p = DEFAULT.replace(g_k=6.8e-5)
    res = _fake_fit(p, ("g_k",), {"g_k": (1e-12, 100.0)})
    assert res.on_bounds() == {}


def test_on_bounds_flags_genuinely_pinned_parameters():
    """The magnitude fix must not lose true positives.

    - r_s = 0 and r_s = 3.8e-9 Ohm (session-6 fit A value) are numerically at
      the zero lower bound -> pinned (ZERO_FLOOR makes 'at zero' decidable
      when the bound itself has no magnitude).
    - r_s = 0.66 mOhm (session-7 fit B value) is a real resistance -> free.
    - k_lambda = 0.004 at its device bound 0.004 1/V -> pinned UPPER (the real
      open problem from the ST8L65 fits).
    """
    dev = (0.0, 0.036)
    p0 = DEFAULT.replace(r_s=0.0)
    assert _fake_fit(p0, ("r_s",), {"r_s": dev}).on_bounds() == {"r_s": "LOWER 0"}
    p1 = DEFAULT.replace(r_s=3.8e-9)
    assert _fake_fit(p1, ("r_s",), {"r_s": dev}).on_bounds() == {"r_s": "LOWER 0"}
    p2 = DEFAULT.replace(r_s=0.66e-3)
    assert _fake_fit(p2, ("r_s",), {"r_s": dev}).on_bounds() == {}
    p3 = DEFAULT.replace(k_lambda=0.004)
    assert _fake_fit(p3, ("k_lambda",), {"k_lambda": (1e-6, 0.004)}).on_bounds() == {
        "k_lambda": "UPPER 0.004"
    }
    # just below the bound is NOT pinned: 0.0039 is an interior optimum
    p4 = DEFAULT.replace(k_lambda=0.0039)
    assert _fake_fit(p4, ("k_lambda",), {"k_lambda": (1e-6, 0.004)}).on_bounds() == {}


def test_curveset_validation():
    with pytest.raises(ValueError):
        CurveSet(np.zeros(3), np.zeros(2), np.zeros(3))


def test_example_dataset_is_usable():
    d = example_dataset()
    assert len(d) == 25 * 8
    assert len(d.v_gs_grid()) == 25
    assert len(d.v_ds_grid()) == 8
