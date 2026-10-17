"""Improved VDMOS channel current: one continuous expression, no triode/sat switch.

    i_d = k_p*v_ov_eff*(v_dseff - v_dseff**2/(2*v_knee)) / (1 + v_dseff/v_vsat)
            * (1 + k_lambda*v_ds*k_db)
          + g_k*v_ds*k_kink*k_db

    k_db     = tanh(v_ds/v_dbsat)                 drain-bias factor: 0 ohmic, ~1 above knee
    v_th     = v_th0 - k_eta*v_ds*k_db            DIBL (bounded by the same factor)
    v_ov     = v_sub*ln(1 + exp((v_gs - v_th)/v_sub))  smooth subthreshold -> strong inversion
    v_ov_eff = v_ov/(1 + (v_ov/v_ovmax)**k_ovsharp)**(1/k_ovsharp)
                                                soft charge cap (v_ovmax = 0 -> off)
    v_knee   = v_ov_eff/(1 + k_theta*v_ov_eff)    knee voltage (field-dependent mobility)
    v_dseff  = v_knee*tanh(v_ds/v_knee)           soft drain variable (k_kneesharp = 0, legacy)
             = v_ds/(1 + (v_ds/v_knee)**k_kneesharp)**(1/k_kneesharp)   power-law soft-min

Naming: <type prefix>_<suffix> with k_ for dimensionless factors and Greek letters
(k_lambda / k_eta / k_theta keep their industry-standard symbols; k_p is A/V^2).
`db` = drain **bias** (not drain-to-body: this model has no bulk terminal).
Convention table and unit rules: README.md, "Naming convention".

With k_theta = v_vsat = k_eta = g_k = 0 this is a smooth approximation of the Level-1
equations with m_triode = 1: exact in the ohmic limit and in deep saturation,
up to ~6 % low around the knee (tests/test_model.py::test_level1_limit_agreement).

The body diode (diode.py, LTspice diode DC branch) is assembled at the TERMINALS:

    i_d = i_ch(v_ds > 0) - i_di(-v_ds)

The diode bypasses r_d/r_s (only channel current flows through the package
resistances), so the package-R solve and the diode r_diode solve stay
independent one-dimensional problems. The channel contribution is 0 for
v_ds <= 0 (the channel formulation is forward-polarity only; sign-symmetric
channel conduction for synchronous rectification is a Queue item, not this
model). With the diode off (i_s = i_sr = 0, the default) vdmos_id() is
bit-identical to the pre-diode model.

Design rationale -- Level-1 comparison, charge-cap physics, package-R solver
derivation, body-diode topology and LTspice fidelity -- lives in the project
README.md, not in these docstrings.
"""

from __future__ import annotations

import numpy as np

from .diode import _solve_diode_series_r, vdmos_diode
from .params import DEFAULT, Params

__all__ = ["vdmos_id", "vdmos_id_level1", "vdmos_diode"]


def vdmos_id(v_ds, v_gs, p: Params = DEFAULT, *, t_c: float = 25.0) -> np.ndarray:
    """Drain current [A] of the improved VDMOS: channel + body diode + package R.

    Terminal assembly (drain-referenced, positive into the drain):

        i_d = i_ch(v_ds > 0) - i_di(v_sd = -v_ds)

    The channel is evaluated at the internal node voltages through the
    package r_d/r_s solve (forward polarity only); the body diode connects
    directly to the terminals and bypasses r_d/r_s, carrying all reverse
    current. With i_s = i_sr = 0 (the default) the diode is exactly zero and
    this function is bit-identical to the pre-diode model.

    Parameters
    ----------
    v_ds, v_gs : float or array-like, broadcast together
        Drain-source and gate-source *terminal* voltages [V].
    p : Params
        Model parameter set (see :mod:`vdmos_fit.params`). r_d = r_s = 0 selects
        the intrinsic fast path (bit-identical to the pre-package-R model).
    t_c : float
        Device temperature [degC] for the diode branch (ignored when the
        diode is off; the future Tj node supersedes it).

    Returns
    -------
    ndarray | float
        i_d [A], same broadcast shape as (v_ds, v_gs).
    """
    v_ds = np.asarray(v_ds, dtype=float)
    v_gs = np.asarray(v_gs, dtype=float)
    v_ds, v_gs = np.broadcast_arrays(v_ds, v_gs)

    diode_off = p.i_s == 0.0 and p.i_sr == 0.0
    if diode_off:
        i_di = 0.0
    else:
        t_k = np.full(v_ds.shape, t_c + 273.15, dtype=float)
        i_di = _solve_diode_series_r(-v_ds, p.r_diode, p, t_k)

    fwd = v_ds > 0.0
    if not np.any(fwd):
        return np.zeros_like(v_ds) if diode_off else -i_di

    if p.r_d == 0.0 and p.r_s == 0.0:
        if np.all(fwd):
            i_ch = _channel_current(v_ds, v_gs, p)
        else:
            i_ch = np.zeros_like(v_ds)
            i_ch[fwd] = _channel_current(v_ds[fwd], v_gs[fwd], p)
    else:
        # Solve on the forward-polarity points only: the channel model is
        # singular in reverse, and the diode carries v_ds <= 0 by itself.
        # The wrapper is 0 at v_ds = 0, so the package-R bracket proof
        # (g strictly decreasing, f(0) = 0) is unchanged.
        vd, vg = v_ds[fwd], v_gs[fwd]
        i_ch = np.zeros_like(v_ds)
        i_ch[fwd] = _solve_package_r(
            vd, vg, p.r_d, p.r_s, lambda d, g: _channel_current(d, g, p)
        )

    return i_ch - i_di


def _channel_current(v_ds, v_gs, p: Params) -> np.ndarray:
    """Intrinsic channel current at the *internal* node voltages (no package R)."""
    # Common drain-bias factor: 0 in the ohmic region, ~1 well above the knee.
    k_db = np.tanh(v_ds / p.v_dbsat)

    # 1) DIBL-shifted threshold (bounded by the same factor).
    v_th = p.v_th0 - p.k_eta * v_ds * k_db

    # 2) Smooth overdrive, overflow-safe softplus.
    z = np.clip((v_gs - v_th) / p.v_sub, -700.0, 700.0)
    v_ov = p.v_sub * np.log1p(np.exp(z))

    # 2b) Soft charge cap (v_ovmax = 0 is the identity; see README).
    v_ov_eff = _soft_cap(v_ov, p.v_ovmax, p.k_ovsharp)

    # 3) Knee voltage from vertical-field mobility degradation, built from the
    # CAPPED overdrive (charge and knee saturate together; see README).
    v_knee = np.maximum(v_ov_eff / (1.0 + p.k_theta * v_ov_eff), p.v_kneemin)

    # 4) Soft drain variable: v_dseff -> v_ds for v_ds << v_knee, -> v_knee for
    #    v_ds >> v_knee. k_kneesharp = 0 (default) keeps the legacy tanh
    #    soft-switch; k_kneesharp > 0 uses the same power-law soft-min as the
    #    charge cap, which decouples the knee WIDTH from the knee voltage:
    #    dv_dseff/dv_ds(0) = 1 and v_dseff -> v_knee hold for every
    #    k_kneesharp, so Ron and the i_d,sat plateau are untouched while the
    #    transition sharpens (k_kneesharp -> inf is a hard min) or softens
    #    (k_kneesharp = 1, harmonic). Sign-preserving so negative v_ds stays
    #    finite (the model is unphysical there until the body diode, but must
    #    not emit NaN).
    if p.k_kneesharp == 0.0:
        v_dseff = v_knee * np.tanh(v_ds / v_knee)
    else:
        v_dseff = np.sign(v_ds) * _soft_min(np.abs(v_ds), v_knee, p.k_kneesharp)

    # 5) Unified channel charge term (triode + saturation in one C^1 expression).
    i_ch = p.k_p * v_ov_eff * (v_dseff - 0.5 * v_dseff * v_dseff / v_knee)

    # 6) Velocity saturation + gated output conductance.
    i_ch /= 1.0 + v_dseff / p.v_vsat
    i_ch *= 1.0 + p.k_lambda * v_ds * k_db

    # 7) Trench kink: a gate-gated parallel drain conductance, NOT a Vgs-dependent
    #    CLM (it is additive and independent of the channel term -- README,
    #    "Trench kink: one-sided sigmoid"). One-sided sigmoid on purpose: bare
    #    tanh() is odd and injects non-physical negative current below the kink
    #    onset; the 0.5*(1 + tanh(.)) form matches tanh above onset and decays to
    #    zero below.
    k_kink = 0.5 * (1.0 + np.tanh((v_gs - p.v_thk) / p.v_kovsoft))
    i_kink = p.g_k * v_ds * k_kink * k_db

    return i_ch + i_kink


def _soft_min(x, k, n):
    """Power-law soft-min (x**-n + k**-n)**(-1/n) = x/(1 + (x/k)**n)**(1/n).

    Log-space (logaddexp) so x**n cannot overflow (n ~ 10-40 is a fitted
    regime). Requires x >= 0, k > 0, n > 0. Properties (shared by the charge
    cap and the drain soft-switch): identity for x << k, asymptote k for
    x >> k, derivative (1 + (x/k)**n)**(-1/n - 1) so 0 <= df/dx <= 1 --
    which is what keeps the package-R solver bracket invariant
    dg/di_d <= -1 intact. """
    with np.errstate(divide="ignore", invalid="ignore"):
        t = n * (np.log(np.maximum(x, 1e-300)) - np.log(k))
    return x * np.exp(-np.logaddexp(0.0, t) / n)


def _soft_cap(v_ov, v_ovmax: float, k_ovsharp: float):
    """Soft power-law cap v_ov/(1 + (v_ov/v_ovmax)**k_ovsharp)**(1/k_ovsharp), in log space.

    Log-space form prevents overflow of (v_ov/v_ovmax)**k_ovsharp (k_ovsharp ~
    10-40 is a fitted regime). v_ovmax <= 0 returns v_ov unchanged -- the exact
    pre-cap model, bit for bit. Properties and why 0 <= dv_ov_eff/dv_ov <= 1
    matters: README, "Charge-cap numerics".
    """
    if v_ovmax <= 0.0:
        return v_ov
    return _soft_min(v_ov, v_ovmax, k_ovsharp)


def _solve_package_r(
    v_ds: np.ndarray,
    v_gs: np.ndarray,
    r_d: float,
    r_s: float,
    f,
    *,
    atol: float = 1e-12,
    rtol: float = 1e-10,
    max_iter: int = 60,
) -> np.ndarray:
    """Solve i_d = f(v_ds - (r_d + r_s)*i_d, v_gs - r_s*i_d) for the terminal currents.

    Bracketed safeguarded Newton on i_d in [0, v_ds/(r_d + r_s)] with a stall
    guard and bisection fallback; fixed-point iteration is NOT used (it diverges
    -- see README, "Package resistance (r_d, r_s)" and "Package-R solver").
    Forward polarity only (v_ds >= 0). Raises RuntimeError if any grid point
    fails to converge; never silently returns an unconverged current.
    """
    if r_d < 0.0 or r_s < 0.0:
        raise ValueError(f"r_d and r_s must be >= 0 (got r_d={r_d}, r_s={r_s})")
    if np.any(v_ds < 0.0):
        raise ValueError(
            "package-R embedding requires v_ds >= 0: the channel-only model is "
            "singular/garbage for reverse polarity (needs the body diode)."
        )
    r_tot = r_d + r_s
    if r_tot == 0.0:
        return f(v_ds, v_gs)

    lo = np.zeros_like(v_ds)
    hi = v_ds / r_tot
    i_d = 0.5 * (lo + hi)
    stall = np.zeros(v_ds.shape, dtype=int)
    prev_width = None

    g = None
    for _ in range(max_iter):
        v_ds_int = v_ds - r_tot * i_d
        v_gs_int = v_gs - r_s * i_d
        g = f(v_ds_int, v_gs_int) - i_d

        # g > 0 -> root above i_d; g < 0 -> root below (g is strictly decreasing).
        lo = np.where(g > 0.0, i_d, lo)
        hi = np.where(g < 0.0, i_d, hi)

        # Stall detector: a Newton two-cycle can shrink the bracket only
        # ~1e-5 per iteration; force bisection until it halves again.
        width = hi - lo
        if prev_width is not None:
            stall = np.where(width > 0.9 * prev_width, stall + 1, 0)
        prev_width = width

        done = np.abs(g) <= atol + rtol * np.abs(i_d)
        if done.all():
            return i_d

        # Central-difference partials of f (model is C1).
        hd = 1e-6 * np.maximum(1.0, np.abs(v_ds_int))
        hg = 1e-6 * np.maximum(1.0, np.abs(v_gs_int))
        fd = (f(v_ds_int + hd, v_gs_int) - f(v_ds_int - hd, v_gs_int)) / (2.0 * hd)
        fg = (f(v_ds_int, v_gs_int + hg) - f(v_ds_int, v_gs_int - hg)) / (2.0 * hg)
        dg = -r_tot * fd - r_s * fg - 1.0  # <= -1 always

        cand = i_d - g / dg
        # Accept Newton only inside the bracket, with real progress, and when
        # the bracket has not stalled; otherwise bisect.
        inside = (cand > lo) & (cand < hi)
        moved = np.abs(cand - i_d) > 1e-12 * (hi - lo)
        cand = np.where(inside & moved & (stall < 3), cand, 0.5 * (lo + hi))
        i_d = np.where(done, i_d, cand)

    raise RuntimeError(
        f"package-R solve failed after {max_iter} iterations: "
        f"max |residual| = {np.abs(g).max():.3g} A"
    )


def vdmos_id_level1(
    v_ds,
    v_gs,
    *,
    k_p=28.0,
    v_th0=3.0,
    v_sub=0.05,
    m_triode=1.2,
    k_lambda=0.02,
    r_d=0.0,
    r_s=0.0,
):
    """LTspice Level-1 VDMOS drain current, as actually implemented.

    The triode/saturation knee is at **v_ds = v_ov / m_triode** (not v_ov), which
    makes Level-1 C1 continuous for any mtriode; m_triode is a pure shape
    knob. Verified against LTspice 26.0.1 .meas -- see README, "Why", and
    tests/test_model.py::test_level1_is_c1_continuous. r_d/r_s use the same
    external series embedding as :func:`vdmos_id` (shared solver); pass
    r_d = r_s = 0 for the bare equation reference.

    ``m_triode`` deliberately keeps the LTspice model-card spelling: this
    function is a foreign-model adapter, and the name is the mapping.
    """
    v_ds = np.asarray(v_ds, dtype=float)
    v_gs = np.asarray(v_gs, dtype=float)
    v_ds, v_gs = np.broadcast_arrays(v_ds, v_gs)
    if r_d == 0.0 and r_s == 0.0:
        return _level1_current(
            v_ds, v_gs, k_p=k_p, v_th0=v_th0, v_sub=v_sub, m_triode=m_triode, k_lambda=k_lambda
        )
    return _solve_package_r(
        v_ds,
        v_gs,
        r_d,
        r_s,
        lambda d, g: _level1_current(
            d, g, k_p=k_p, v_th0=v_th0, v_sub=v_sub, m_triode=m_triode, k_lambda=k_lambda
        ),
    )


def _level1_current(v_ds, v_gs, *, k_p, v_th0, v_sub, m_triode, k_lambda):
    z = np.clip((v_gs - v_th0) / v_sub, -700.0, 700.0)
    v_ov = v_sub * np.log1p(np.exp(z))
    v_knee = v_ov / m_triode
    i_tri = m_triode * k_p * (v_ov * v_ds - 0.5 * m_triode * v_ds**2) * (1.0 + k_lambda * v_ds)
    i_sat = 0.5 * k_p * v_ov**2 * (1.0 + k_lambda * v_ds)
    return np.where(v_ds < v_knee, i_tri, i_sat)
