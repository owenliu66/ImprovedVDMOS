"""Tests for the interface-trap-filling gain (k_beta_trap, v_ov_eff2).

    v_ov_eff2 = v_ov_eff * (1 + k_beta_trap * v_ov_eff)    (k_beta_trap [1/V], 0 -> identity)

The gain models 4H-SiC channel-mobility recovery as the surface potential
deepens and interface traps fill -- the WBG counterpart (same knob family,
opposite sign) of the super-junction charge cap. v_ov_eff2 is THE effective
overdrive downstream: v_knee = max(v_ov_eff2/(1 + k_theta*v_ov_eff2), v_kneemin)
and the i_d prefactor both use it (composition decision: k_theta multiplies
v_ov_eff2, not v_ov_eff -- the knee is set by the FINAL channel charge).

These tests pin: identity at k_beta_trap = 0 (every pre-gain anchor stays
bit-identical), monotonicity in v_gs for k_beta_trap >= 0 (the property the
package-R bracket proof actually needs -- the gain breaks the
dv_ov_eff/dv_ov <= 1 bound but not monotonicity), the composite cap+gain
asymptote v_ovmax*(1 + k_beta_trap*v_ovmax), the k_beta_trap >= 0 guard
(ValueError + BOUNDS lower bound 0.0), a package-R solve with r_s > 0 AND the
gain (the untested silicon regime per HANDOVER 1.4), the gain-derivative
formula 1 + 2*k_beta_trap*v_ov_eff, and synthetic round trips that recover
k_beta_trap (cap off) plus the joint cap+gain case whose expected
k_beta_trap <-> v_ovmax anti-correlation is REPORTED, not asserted.
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
    fit,
    identifiability,
    make_synthetic,
    vdmos_id,
)
from vdmos_fit.model import _channel_current


# ---------------------------------------------------------------- identity
def test_k_beta_trap_zero_bit_identical():
    """Params(k_beta_trap=0.0) must equal Params() and reproduce the pre-gain
    model bit for bit -- fast path and through the package-R solver, cap
    active or not (the gain identity is exact for any cap setting)."""
    assert Params(k_beta_trap=0.0) == Params()
    assert DEFAULT.k_beta_trap == 0.0

    v_ds = np.linspace(0.0, 60.0, 121)
    v_gs = np.linspace(-2.0, 18.0, 101)
    vd2, vg2 = np.meshgrid(v_ds, v_gs)

    ref = vdmos_id(vd2, vg2, DEFAULT)
    assert np.array_equal(vdmos_id(vd2, vg2, Params(k_beta_trap=0.0)), ref)

    capped = DEFAULT.replace(v_ovmax=1.9, k_ovsharp=9.9, r_d=2e-3, r_s=1e-3)
    capped_off_gain = capped.replace(k_beta_trap=0.0)
    assert np.array_equal(vdmos_id(vd2, vg2, capped), vdmos_id(vd2, vg2, capped_off_gain))


# ---------------------------------------------------------------- fitting-space decisions
def test_k_beta_trap_fitting_space_decisions():
    """Linear space (0 is the legitimate identity), bounds (0.0, 2.0) with the
    lower bound EXACTLY 0.0 (user decision: k_beta_trap >= 0), and NOT in
    DEFAULT_FREE (freeing it silently changes existing fits)."""
    assert "k_beta_trap" not in LOG_PARAMS
    assert "k_beta_trap" not in DEFAULT_FREE
    assert BOUNDS["k_beta_trap"] == (0.0, 2.0)


# ---------------------------------------------------------------- guard
def test_negative_k_beta_trap_guard():
    """k_beta_trap < 0 is rejected outright (HANDOVER 1.4(4) decision):
    v(1 + k*v) is non-monotone for k < 0 -- derivative 1 + 2*k*v crosses zero
    at v = -1/(2k) -- making i_d multi-valued in v_gs and killing the
    package-R bracket proof. Physics: trap filling only RAISES mobility.
    Pinned at BOTH the model entry points (vdmos_id and _channel_current)."""
    p = DEFAULT.replace(k_beta_trap=-1e-6)
    with pytest.raises(ValueError, match="k_beta_trap must be >= 0"):
        vdmos_id(5.0, 10.0, p)
    with pytest.raises(ValueError, match="k_beta_trap must be >= 0"):
        vdmos_id(np.linspace(-2.0, 20.0, 11), 10.0, p)  # reverse-only grid too
    with pytest.raises(ValueError, match="k_beta_trap must be >= 0"):
        _channel_current(np.array([5.0]), np.array([10.0]), p)
    # k_beta_trap = 0 exactly is legal (the identity sentinel).
    assert np.isfinite(vdmos_id(5.0, 10.0, DEFAULT.replace(k_beta_trap=0.0)))


# ---------------------------------------------------------------- monotonicity
def test_gain_monotone_in_vgs():
    """For k_beta_trap >= 0, i_d is strictly increasing in v_gs at fixed v_ds
    over the whole sampled grid (forward polarity only). The gain lifts
    dv_ov_eff2/dv_ov_eff = 1 + 2*k_beta_trap*v_ov_eff above 1 but keeps it
    positive -- exactly the property the package-R bracket proof needs."""
    v_gs = np.linspace(-2.0, 18.0, 81)
    for k_beta_trap in (0.0, 0.5, 2.0):
        p = DEFAULT.replace(k_beta_trap=k_beta_trap)
        for v_ds in (0.5, 5.0, 20.0):
            i_d = vdmos_id(v_ds, v_gs, p)
            assert np.all(np.isfinite(i_d))
            assert np.all(np.diff(i_d) > 0.0), f"i_d not increasing, k_beta_trap={k_beta_trap}"
    # cap active + gain: composite monotone too
    p = DEFAULT.replace(v_ovmax=2.0, k_ovsharp=10.0, k_beta_trap=0.5)
    i_d = vdmos_id(10.0, v_gs, p)
    assert np.all(np.diff(i_d) > 0.0)


def test_gain_derivative_formula():
    """d i_d/d v_gs in the ohmic limit must follow the chain rule with the
    QUADRATIC gain derivative 1 + 2*k_beta_trap*v_ov_eff (pins the factor 2,
    i.e. the mechanism is quadratic in v_ov, not a linear re-weighting)."""
    p = Params(k_p=30.0, v_th0=3.0, v_sub=0.05, k_beta_trap=0.4,
               k_theta=0.0, v_vsat=1e4, k_eta=0.0, k_lambda=0.0,
               g_k=0.0, v_dbsat=10.0)
    v_gs, v_ds = 10.0, 0.01
    h = 1e-4
    num = (float(vdmos_id(v_ds, v_gs + h, p)) - float(vdmos_id(v_ds, v_gs - h, p))) / (2 * h)
    # z = (10-3)/0.05 = 140: softplus slope sigmoid(z) = 1.0, v_ov = 7.0 exactly
    v_ov = 7.0
    ana = p.k_p * v_ds * (1.0 + 2.0 * p.k_beta_trap * v_ov)
    assert num == pytest.approx(ana, rel=5e-4)


# ---------------------------------------------------------------- composition
def test_gain_raises_cap_ceiling():
    """With the cap active the COMPOSITE asymptote of v_ov_eff2 is
    v_ovmax*(1 + k_beta_trap*v_ovmax), NOT v_ovmax (HANDOVER 1.5 -- the
    silicon SJ composition gotcha). Read it through the ohmic channel
    conductance i_d/(k_p*v_ds) at v_gs far above the cap."""
    p = Params(k_p=30.0, v_th0=3.0, v_sub=0.05, v_ovmax=2.0, k_ovsharp=10.0,
               k_theta=0.0, v_vsat=1e4, k_eta=0.0, k_lambda=0.0,
               g_k=0.0, v_dbsat=10.0)
    v_gs, v_ds = 30.0, 1e-4  # v_ov ~ 27 >> v_ovmax: cap error ~ (2/10)*(2/27)^10 ~ 1e-12
    # residual triode correction (1 - v_dseff/(2*v_knee)) ~ 1 - 2.5e-5 (pre) /
    # 1 - 1.25e-5 (gained, v_knee = 4 V) -- inside the tolerances below.

    pre = float(vdmos_id(v_ds, v_gs, p)) / (p.k_p * v_ds)
    assert pre == pytest.approx(2.0, rel=1e-4)  # pre-gain asymptote: v_ovmax

    gained = p.replace(k_beta_trap=0.5)
    post = float(vdmos_id(v_ds, v_gs, gained)) / (gained.k_p * v_ds)
    assert post == pytest.approx(2.0 * (1.0 + 0.5 * 2.0), rel=1e-4)  # = 4.0 V


# ---------------------------------------------------------------- package R
def test_package_r_bracket_with_gain_cap_off():
    """r_s > 0 AND k_beta_trap > 0 (the untested silicon regime per HANDOVER
    1.4): the solver converges on a dense grid, matches independent brentq
    roots, and stays finite/non-negative. The bracket proof survives because
    the gain is monotone: dg/di = -(r_d+r_s) f_vds - r_s f_vgs - 1 <= -1."""
    p = Params(k_p=30.0, v_th0=3.1, v_sub=0.05, k_beta_trap=0.3,
               k_theta=0.02, v_vsat=1.2, k_eta=0.03, k_lambda=0.015,
               g_k=0.0, v_dbsat=10.0, r_d=2e-3, r_s=1e-3)
    rtot = p.r_d + p.r_s
    vd2, vg2 = np.meshgrid(np.linspace(0.0, 40.0, 401), np.linspace(0.0, 14.0, 71))
    i_d = vdmos_id(vd2, vg2, p)
    assert np.all(np.isfinite(i_d))
    assert np.all(i_d >= 0.0)
    resid = np.abs(_channel_current(vd2 - rtot * i_d, vg2 - p.r_s * i_d, p) - i_d)
    mask = i_d != 0.0
    assert np.all(resid[mask] <= 1e-12 + 1e-8 * np.abs(i_d[mask]))

    for v_ds, v_gs in ((0.05, 14.0), (5.0, 10.0), (40.0, 14.0)):
        def h(i_d_, v_ds=v_ds, v_gs=v_gs):  # bind loop vars (ruff B023)
            return float(_channel_current(v_ds - rtot * i_d_, v_gs - p.r_s * i_d_, p)) - i_d_

        root = brentq(h, 0.0, v_ds / rtot, xtol=1e-15, rtol=8.9e-16)
        assert float(vdmos_id(v_ds, v_gs, p)) == pytest.approx(root, rel=1e-10)


def test_package_r_stall_guard_with_cap_and_gain():
    """Gain analogue of test_package_r.py::test_solver_survives_newton_two_cycle:
    the fitted SJ cap set + k_beta_trap = 0.06 + r_s > 0. The gain inflates the
    bracket upper bound i_max = f(v_ds, v_gs) AND steepens the gate dependence
    (Newton overshoots more), so the stall guard must carry the solve through
    the gate-degeneration cliff at high v_ds."""
    p = Params(
        k_p=91.9481, v_th0=4.7327, v_sub=0.06, v_ovmax=1.90909, k_ovsharp=9.86552,
        k_beta_trap=0.06,
        k_theta=0.0, v_vsat=1e4, k_eta=0.0, k_lambda=0.0383371,
        g_k=0.0, v_thk=1.6, v_kovsoft=0.1, v_dbsat=0.653951,
        r_d=26.03e-3, r_s=4.644e-3,
    )
    rtot = p.r_d + p.r_s
    v_ds = np.linspace(0.0, 200.0, 401)
    i_d = vdmos_id(v_ds, 10.0, p)
    assert np.all(np.isfinite(i_d))
    resid = np.abs(_channel_current(v_ds - rtot * i_d, 10.0 - p.r_s * i_d, p) - i_d)
    assert (resid / np.maximum(np.abs(i_d), 1e-12)).max() < 1e-8


# ---------------------------------------------------------------- fitting
GAINTRUTH = Params(
    k_p=30.0, v_th0=3.1, v_sub=0.05, k_beta_trap=0.3,
    k_theta=0.02, v_vsat=1.2, k_eta=0.03, k_lambda=0.015,
    g_k=0.0, v_thk=1.6, v_kovsoft=0.1, v_dbsat=10.0,
)

CAPGAINTRUTH = GAINTRUTH.replace(v_ovmax=2.0, k_ovsharp=8.0, k_beta_trap=0.15)


def test_synthetic_roundtrip_gain():
    """make_synthetic with k_beta_trap in the truth (cap off) on a grid with
    dense low-v_ds transfer points -> fit recovers it with finite sigma_rel.
    The low-v_ds points carry the super-linear conductance signature
    (dv_ov_eff2/dv_gs rising past 1) that a pre-gain model cannot produce."""
    v_gs = np.linspace(2.0, 10.0, 33)
    v_ds = np.array([0.05, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0])
    data = make_synthetic(GAINTRUTH, v_gs, v_ds, noise=0.005, seed=21)

    free = ("k_p", "v_th0", "v_sub", "k_beta_trap",
            "k_theta", "v_vsat", "k_eta", "k_lambda")
    start = DEFAULT.replace(k_p=30.0, g_k=0.0, k_beta_trap=0.1)
    res = fit(data, start=start, free=free, weight="relative", max_nfev=8000)

    assert res.success
    assert res.rmse_rel < 0.02
    assert res.params.k_beta_trap == pytest.approx(0.3, rel=0.15)
    assert res.params.k_p == pytest.approx(30.0, rel=0.05)
    sigma = identifiability(res, free=res.free).sigma_rel
    assert np.isfinite(sigma["k_beta_trap"])


def test_synthetic_roundtrip_cap_and_gain():
    """Joint cap+gain case (HANDOVER 1.5): with BOTH mechanisms free the
    composite asymptote v_ovmax*(1 + k_beta_trap*v_ovmax) is what the data
    pins, so k_beta_trap <-> v_ovmax anti-correlate and v_ovmax drops to
    compensate. The correlation is REPORTED, not asserted; the fit must still
    reproduce the curve."""
    v_gs = np.linspace(2.0, 10.0, 33)
    v_ds = np.array([0.05, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0])
    data = make_synthetic(CAPGAINTRUTH, v_gs, v_ds, noise=0.005, seed=22)

    free = ("k_p", "v_th0", "v_sub", "v_ovmax", "k_ovsharp", "k_beta_trap",
            "k_theta", "v_vsat", "k_eta", "k_lambda")
    start = DEFAULT.replace(k_p=30.0, g_k=0.0, v_ovmax=1.0, k_ovsharp=5.0,
                            k_beta_trap=0.1)
    res = fit(data, start=start, free=free, weight="relative", max_nfev=8000)

    assert res.success
    assert res.rmse_rel < 0.03
    # composite ceiling is the identifiable combination, not the individual knobs
    ceiling = res.params.v_ovmax * (1.0 + res.params.k_beta_trap * res.params.v_ovmax)
    truth_ceiling = CAPGAINTRUTH.v_ovmax * (1.0 + CAPGAINTRUTH.k_beta_trap * CAPGAINTRUTH.v_ovmax)
    assert ceiling == pytest.approx(truth_ceiling, rel=0.05)

    ident = identifiability(res, free=res.free)
    pairs = {(a, b): c for a, b, c in ident.worst_pairs}
    corr = next((c for (a, b), c in pairs.items()
                 if {a, b} == {"k_beta_trap", "v_ovmax"}), None)
    print(f"cap+gain round trip: corr(k_beta_trap, v_ovmax) = "
          f"{'not in worst pairs (<0.8)' if corr is None else f'{corr:+.3f}'}")
