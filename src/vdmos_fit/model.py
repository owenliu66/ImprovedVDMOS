"""Improved VDMOS channel current: one continuous expression, no triode/sat switch.

    I_D = kp*Vov_eff*(v - v^2/(2*Vk)) / (1 + v/vsat) * (1 + lam*Vds*tanh(Vds/Vk))
          + gk*Vds*tanh((Vgs - Vth_k)/kk)*tanh(Vds/Vk)

    Vth     = Vth0 - eta*Vds*tanh(Vds/Vk)          DIBL (bounded)
    Vov     = ksub*ln(1 + exp((Vgs - Vth)/ksub))   smooth subthreshold -> strong inversion
    Vov_eff = Vov/(1 + (Vov/kcap)**m)**(1/m)       soft charge cap (kcap = 0 -> off)
    Vk      = Vov_eff/(1 + theta*Vov_eff)          knee voltage (field-dependent mobility)
    v       = Vk*tanh(Vds/Vk)                      soft drain variable (replaces the knee switch)

With theta = vsat = eta = gk = 0 this is a smooth approximation of the Level-1
equations with m_triode = 1: exact in the ohmic limit and in deep saturation,
up to ~6 % low around the knee (tests/test_model.py::test_level1_limit_agreement).

Design rationale -- Level-1 comparison, gate-cap physics, package-R solver
derivation -- lives in the project README.md, not in these docstrings.
"""

from __future__ import annotations

import numpy as np

from .params import DEFAULT, Params

__all__ = ["vdmos_id", "vdmos_id_level1"]


def vdmos_id(vds, vgs, p: Params = DEFAULT) -> np.ndarray:
    """Drain current [A] of the improved VDMOS, including package rd/rs.

    Parameters
    ----------
    vds, vgs : float or array-like, broadcast together
        Drain-source and gate-source *terminal* voltages [V].
    p : Params
        Model parameter set (see :mod:`vdmos_fit.params`). rd = rs = 0 selects
        the intrinsic fast path (bit-identical to the pre-package-R model).

    Returns
    -------
    ndarray | float
        I_D [A], same broadcast shape as (vds, vgs).
    """
    vds = np.asarray(vds, dtype=float)
    vgs = np.asarray(vgs, dtype=float)
    vds, vgs = np.broadcast_arrays(vds, vgs)
    if p.rd == 0.0 and p.rs == 0.0:
        return _channel_current(vds, vgs, p)
    return _solve_package_r(
        vds, vgs, p.rd, p.rs, lambda d, g: _channel_current(d, g, p)
    )


def _channel_current(vds, vgs, p: Params) -> np.ndarray:
    """Intrinsic channel current at the *internal* node voltages (no package R)."""
    # Common drain-bias gate: 0 in the ohmic region, ~1 well above the knee.
    gate = np.tanh(vds / p.vk)

    # 1) DIBL-shifted threshold (bounded by the same gate).
    vth = p.vth0 - p.eta * vds * gate

    # 2) Smooth overdrive, overflow-safe softplus.
    z = np.clip((vgs - vth) / p.ksub, -700.0, 700.0)
    vov = p.ksub * np.log1p(np.exp(z))

    # 2b) Soft charge cap (kcap = 0 is the identity; see README).
    vov = _soft_cap(vov, p.kcap, p.m)

    # 3) Knee voltage from vertical-field mobility degradation, built from the
    # CAPPED overdrive (charge and knee saturate together; see README).
    vknee = np.maximum(vov / (1.0 + p.theta * vov), p.vov_floor)

    # 4) Soft drain variable: v -> Vds for Vds << Vk, v -> Vk for Vds >> Vk.
    v = vknee * np.tanh(vds / vknee)

    # 5) Unified channel charge term (triode + saturation in one C^1 expression).
    ich = p.kp * vov * (v - 0.5 * v * v / vknee)

    # 6) Velocity saturation + gated output conductance.
    ich /= 1.0 + v / p.vsat
    ich *= 1.0 + p.lam * vds * gate

    # 7) Trench kink. One-sided sigmoid on purpose: bare tanh() is odd and
    #    injects non-physical negative current below the kink onset; the
    #    0.5*(1 + tanh(.)) form matches tanh above onset and decays to zero
    #    below (see README, "Trench kink: one-sided sigmoid").
    kink_gate = 0.5 * (1.0 + np.tanh((vgs - p.vthk) / p.kk))
    ikink = p.gk * vds * kink_gate * gate

    return ich + ikink


def _soft_cap(vov, kcap: float, m: float):
    """Soft power-law cap vov/(1 + (vov/kcap)**m)**(1/m), evaluated in log space.

    Log-space form prevents overflow of (vov/kcap)**m (m ~ 10-40 is a fitted
    regime). kcap <= 0 returns vov unchanged -- the exact pre-cap model, bit
    for bit. Properties and why 0 <= dVov_eff/dVov <= 1 matters: README,
    "Gate-cap numerics".
    """
    if kcap <= 0.0:
        return vov
    with np.errstate(divide="ignore", invalid="ignore"):
        t = m * (np.log(np.maximum(vov, 1e-300)) - np.log(kcap))
    return vov * np.exp(-np.logaddexp(0.0, t) / m)


def _solve_package_r(
    vds: np.ndarray,
    vgs: np.ndarray,
    rd: float,
    rs: float,
    f,
    *,
    atol: float = 1e-12,
    rtol: float = 1e-10,
    max_iter: int = 60,
) -> np.ndarray:
    """Solve Id = f(Vds - (rd + rs)*Id, Vgs - rs*Id) for the terminal currents.

    Bracketed safeguarded Newton on Id in [0, Vds/(rd+rs)] with a stall guard
    and bisection fallback; fixed-point iteration is NOT used (it diverges --
    see README, "Package resistance (rd, rs)" and "Package-R solver").
    Forward polarity only (Vds >= 0). Raises RuntimeError if any grid point
    fails to converge; never silently returns an unconverged current.
    """
    if rd < 0.0 or rs < 0.0:
        raise ValueError(f"rd and rs must be >= 0 (got rd={rd}, rs={rs})")
    if np.any(vds < 0.0):
        raise ValueError(
            "package-R embedding requires Vds >= 0: the channel-only model is "
            "singular/garbage for reverse polarity (needs the body diode)."
        )
    rtot = rd + rs
    if rtot == 0.0:
        return f(vds, vgs)

    lo = np.zeros_like(vds)
    hi = vds / rtot
    idc = 0.5 * (lo + hi)
    stall = np.zeros(vds.shape, dtype=int)
    prev_width = None

    g = None
    for _ in range(max_iter):
        vd_i = vds - rtot * idc
        vg_i = vgs - rs * idc
        g = f(vd_i, vg_i) - idc

        # g > 0 -> root above idc; g < 0 -> root below (g is strictly decreasing).
        lo = np.where(g > 0.0, idc, lo)
        hi = np.where(g < 0.0, idc, hi)

        # Stall detector: a Newton two-cycle can shrink the bracket only
        # ~1e-5 per iteration; force bisection until it halves again.
        width = hi - lo
        if prev_width is not None:
            stall = np.where(width > 0.9 * prev_width, stall + 1, 0)
        prev_width = width

        done = np.abs(g) <= atol + rtol * np.abs(idc)
        if done.all():
            return idc

        # Central-difference partials of f (model is C1).
        hd = 1e-6 * np.maximum(1.0, np.abs(vd_i))
        hg = 1e-6 * np.maximum(1.0, np.abs(vg_i))
        fd = (f(vd_i + hd, vg_i) - f(vd_i - hd, vg_i)) / (2.0 * hd)
        fg = (f(vd_i, vg_i + hg) - f(vd_i, vg_i - hg)) / (2.0 * hg)
        dg = -rtot * fd - rs * fg - 1.0  # <= -1 always

        cand = idc - g / dg
        # Accept Newton only inside the bracket, with real progress, and when
        # the bracket has not stalled; otherwise bisect.
        inside = (cand > lo) & (cand < hi)
        moved = np.abs(cand - idc) > 1e-12 * (hi - lo)
        cand = np.where(inside & moved & (stall < 3), cand, 0.5 * (lo + hi))
        idc = np.where(done, idc, cand)

    raise RuntimeError(
        f"package-R solve failed after {max_iter} iterations: "
        f"max |residual| = {np.abs(g).max():.3g} A"
    )


def vdmos_id_level1(
    vds,
    vgs,
    *,
    kp=28.0,
    vth0=3.0,
    ksub=0.05,
    m_triode=1.2,
    lam=0.02,
    rd=0.0,
    rs=0.0,
):
    """LTspice Level-1 VDMOS drain current, as actually implemented.

    The triode/saturation knee is at **Vds = Vov / m_triode** (not Vov), which
    makes Level-1 C1 continuous for any mtriode; m_triode is a pure shape
    knob. Verified against LTspice 26.0.1 .meas -- see README, "Why", and
    tests/test_model.py::test_level1_is_c1_continuous. rd/rs use the same
    external series embedding as :func:`vdmos_id` (shared solver); pass
    rd = rs = 0 for the bare equation reference.
    """
    vds = np.asarray(vds, dtype=float)
    vgs = np.asarray(vgs, dtype=float)
    vds, vgs = np.broadcast_arrays(vds, vgs)
    if rd == 0.0 and rs == 0.0:
        return _level1_current(vds, vgs, kp=kp, vth0=vth0, ksub=ksub, m_triode=m_triode, lam=lam)
    return _solve_package_r(
        vds,
        vgs,
        rd,
        rs,
        lambda d, g: _level1_current(
            d, g, kp=kp, vth0=vth0, ksub=ksub, m_triode=m_triode, lam=lam
        ),
    )


def _level1_current(vds, vgs, *, kp, vth0, ksub, m_triode, lam):
    z = np.clip((vgs - vth0) / ksub, -700.0, 700.0)
    vov = ksub * np.log1p(np.exp(z))
    knee = vov / m_triode
    tri = m_triode * kp * (vov * vds - 0.5 * m_triode * vds**2) * (1.0 + lam * vds)
    sat = 0.5 * kp * vov**2 * (1.0 + lam * vds)
    return np.where(vds < knee, tri, sat)
