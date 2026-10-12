"""Tests for parameter extraction: can the fitter recover known parameters?"""

from __future__ import annotations

import math

import numpy as np
import pytest

from vdmos_fit import (
    DEFAULT,
    CurveSet,
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
    kp=30.0, vth0=3.1, theta=0.02, vsat=1.2, eta=0.03, lam=0.015, gk=0.3, vthk=1.5
)


@pytest.fixture(scope="module")
def synthetic():
    vgs = np.linspace(2.0, 14.0, 25)
    vds = np.array([0.05, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0])
    return make_synthetic(TRUTH, vgs, vds, noise=0.005, seed=3)


def test_fit_recovers_synthetic_parameters(synthetic):
    """Noiseless-ish synthetic data must round-trip back to the truth.

    Only the strongly-sensed parameters are checked tightly; theta and vthk are
    weakly identified from this grid (see the identifiability tests) and are
    deliberately not asserted here.
    """
    res = fit(synthetic, start=DEFAULT, weight="relative")
    assert res.success
    assert res.rmse_rel < 0.02
    got, want = res.params.as_dict(), TRUTH.as_dict()
    for name in ("kp", "vth0", "vsat", "eta", "lam", "gk"):
        assert got[name] == pytest.approx(want[name], rel=0.2), name


def test_fit_improves_on_default_start(synthetic):
    res = fit(synthetic, start=DEFAULT)
    start_err = np.sqrt(
        np.mean(
            (
                (vdmos_id(synthetic.vds, synthetic.vgs, DEFAULT) - synthetic.ids)
                / np.maximum(np.abs(synthetic.ids), 1e-3)
            )
            ** 2
        )
    )
    assert res.rmse_rel < 0.2 * start_err


def test_bounds_are_respected(synthetic):
    """A deliberately absurd start must still land inside physical bounds."""
    bad = DEFAULT.replace(kp=1e-4, theta=1e-4, vsat=1e-3)
    res = fit(synthetic, start=bad)
    assert res.params.kp > 0 and res.params.vsat > 0 and res.params.theta > 0


def test_fixed_parameters_stay_fixed(synthetic):
    """Fitting only kp/vth0 must not move theta."""
    res = fit(synthetic, start=DEFAULT, free=("kp", "vth0"))
    assert res.params.theta == DEFAULT.theta
    assert res.params.vsat == DEFAULT.vsat


def test_uncertainty_report_shape(synthetic):
    res = fit(synthetic)
    unc = parameter_uncertainty(res)
    assert set(unc) == set(
        ("kp", "vth0", "theta", "vsat", "eta", "lam", "gk", "vthk", "rd", "rs")
    )
    assert all(v >= 0 for v in unc.values())


def test_kink_parameters_are_weakly_determined(synthetic):
    """The kink onset is barely sampled by this grid, so vthk is poorly determined.

    vthk = 1.5 V with kk = 0.1 V means the kink sigmoid is essentially saturated
    (=1) for every Vgs >= 2 V in this grid, so only the very lowest curves carry
    information about it. The fit lands at ~1.7 V against a truth of 1.5 V -- a
    16 % error -- and the report must show a correspondingly huge sigma rather
    than hiding it. Fix vthk from a measurement, or add dense Vgs points through
    the kink region, before trusting it.
    """
    ident = identifiability(fit(synthetic))
    assert ident.sigma_rel["vthk"] > 1.0  # order-unity or worse
    assert ident.sigma_rel["kp"] < 1.0  # core shape params are fine
    assert ident.rank == ident.n_free  # not degenerate, just weak
    # Threshold raised 1e4 -> 1e5 (2026, package-R change): with rd/rs free and
    # truth rd = rs = 0, the rd <-> rs columns are near-parallel (corr -0.999,
    # only the sum is sensed), which alone pushes cond to ~1.4e4.
    assert ident.cond < 1e5


def test_sensitivity_report_ranks_parameters(synthetic):
    """Pre-fit diagnostic must agree with the post-fit Jacobian diagnosis.

    Sensitivity is measured as mean |d ln Id / d ln p|, so it is directly
    comparable across parameters with different units.
    """
    sens = sensitivity_report(DEFAULT, synthetic)
    assert sens["kp"] > 0.5  # strongly sensed
    assert sens["vthk"] < 1e-3  # ~2000x weaker than kp
    assert sens["vthk"] < 1e-3 * sens["kp"]
    # weakly sensed parameters must be the ones with the biggest error bars
    ident = identifiability(fit(synthetic))
    weakest = min(sens, key=sens.get)
    assert ident.sigma_rel[weakest] > ident.sigma_rel["vth0"]


def test_theta_vsat_correlation_is_reported(synthetic):
    """theta and vsat both flatten the knee; output curves alone cannot separate
    them, so the report must surface the pairing.

    Evaluated with the classic 9-parameter free set: once rd/rs join the fit,
    the knee-flattening trade-off redistributes across more parameters and the
    single theta<->vsat pair drops below the 0.8 reporting threshold (their
    individual sigmas stay order-unity or worse -- see
    tests/test_package_r.py::test_fit_recovers_package_resistance).
    """
    legacy_free = ("kp", "vth0", "theta", "vsat", "eta", "lam", "gk", "vthk")
    ident = identifiability(fit(synthetic, free=legacy_free), legacy_free)
    pairs = {(a, b) for a, b, _ in ident.worst_pairs}
    assert ("theta", "vsat") in pairs
    # and both carry large individual error bars as a consequence
    assert ident.sigma_rel["theta"] > 1.0


def test_dead_parameter_reports_inf_not_a_small_sigma():
    """A parameter with (near-)zero sensitivity must report sigma_rel = inf.

    This is the trap `identifiability()` exists to avoid: a bare
    pinv(J.T @ J) returns *zero* variance for a near-zero Jacobian column,
    which reads as "perfectly determined" when the truth is the exact
    opposite. Reproduced synthetically on the ST8L65 situation: a dataset
    whose lowest Vgs is 8 V (vth0 = 3.1) contains no subthreshold points, so
    ksub carries essentially no information -- its column measures ~1e-9 of
    the largest column, and pinning ksub at 0.30 V leaves the fit at the same
    optimum. The old absolute 1e-10 * max tolerance classified that column
    LIVE and reported a confident small sigma; the relative 1e-5 tolerance
    classifies it dead. eta (sigma ~ 12, genuinely sensed but weak) stays
    live, so the tolerance does not over-trigger.
    """
    truth = Params(kp=30.0, vth0=3.1, ksub=0.05, theta=0.02, vsat=1.2,
                   eta=0.03, lam=0.015, gk=0.0)
    vgs = np.linspace(8.0, 14.0, 25)
    vds = np.array([0.5, 2.0, 5.0, 10.0, 20.0])
    data = make_synthetic(truth, vgs, vds, noise=0.0, seed=0)
    free = ("kp", "vth0", "ksub", "theta", "vsat", "eta", "lam")

    ident = identifiability(fit(data, start=truth, free=free), free)
    assert math.isinf(ident.sigma_rel["ksub"])
    assert ident.unidentifiable == ("ksub",)
    # rank drops by at least the dead column; eta <-> others are also near-
    # degenerate on this grid, so rank is not pinned exactly.
    assert ident.rank <= len(free) - 1
    # the near-dead column must not fake a correlation with a real parameter
    assert not any("ksub" in (a, b) for a, b, _ in ident.worst_pairs)
    # genuinely sensed, merely weak -> must NOT be called dead
    assert not math.isinf(ident.sigma_rel["eta"])

    # Independent confirmation: pinning the dead parameter cannot change the fit.
    pinned = fit(data, start=truth.replace(ksub=0.30),
                 free=tuple(n for n in free if n != "ksub"))
    assert pinned.params.kp == pytest.approx(30.0, rel=1e-3)


def test_curveset_validation():
    with pytest.raises(ValueError):
        CurveSet(np.zeros(3), np.zeros(2), np.zeros(3))


def test_example_dataset_is_usable():
    d = example_dataset()
    assert len(d) == 25 * 8
    assert len(d.vgs_grid()) == 25
    assert len(d.vds_grid()) == 8
