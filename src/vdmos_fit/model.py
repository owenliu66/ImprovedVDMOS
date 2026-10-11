"""Improved VDMOS channel current: one continuous expression, no triode/sat switch.

    I_D = kp*Vov*(v - v^2/(2*Vk)) / (1 + v/vsat) * (1 + lam*Vds*tanh(Vds/Vk))
          + gds*Vds
          + gk*Vds*tanh((Vgs - Vth_k)/kk)*tanh(Vds/Vk)

    Vth  = Vth0 - eta*Vds*tanh(Vds/Vk)          DIBL (bounded)
    Vov  = ksub*ln(1 + exp((Vgs - Vth)/ksub))   smooth subthreshold -> strong inversion
    Vk   = Vov/(1 + theta*Vov)                  knee voltage (field-dependent mobility)
    v    = Vk*tanh(Vds/Vk)                      soft drain variable (replaces the knee switch)

Limits
------
Vds << Vk :  I_D -> kp*Vov*Vds*(1 + ...)                     ohmic, m_triode = 1
Vds >> Vk :  I_D -> 0.5*kp*Vov^2*(0.5+th*Vov)/(1+th*Vov)^2 / (1+Vk/vsat)   saturated

With theta = vsat = eta = gds = gk = 0 the model becomes a *smooth approximation*
of the Level-1 equations with m_triode = 1: exact in the ohmic limit (Vds/Vov -> 0)
and in deep saturation (Vds/Vov >~ 5), but up to ~6 % low around the knee
(Vds/Vov ~ 1), because v = Vk*tanh(Vds/Vk) rounds off the hard corner rather than
reproducing it. See tests/test_model.py::test_level1_limit_agreement.
"""

from __future__ import annotations

import numpy as np

from .params import DEFAULT, Params

__all__ = ["vdmos_id", "vdmos_id_level1"]


def vdmos_id(vds, vgs, p: Params = DEFAULT) -> np.ndarray:
    """Drain current [A] of the improved VDMOS channel.

    Parameters
    ----------
    vds, vgs : float or array-like, broadcast together
        Drain-source and gate-source terminal voltages [V].
    p : Params
        Model parameter set (see :mod:`vdmos_fit.params`).

    Returns
    -------
    ndarray | float
        I_D [A], same broadcast shape as (vds, vgs).
    """
    vds = np.asarray(vds, dtype=float)
    vgs = np.asarray(vgs, dtype=float)
    vds, vgs = np.broadcast_arrays(vds, vgs)

    # Common drain-bias gate: 0 in the ohmic region, ~1 well above the knee.
    gate = np.tanh(vds / p.vk)

    # 1) DIBL-shifted threshold (bounded by the same gate).
    vth = p.vth0 - p.eta * vds * gate

    # 2) Smooth overdrive, overflow-safe softplus.
    z = np.clip((vgs - vth) / p.ksub, -700.0, 700.0)
    vov = p.ksub * np.log1p(np.exp(z))

    # 3) Knee voltage from vertical-field mobility degradation.
    vknee = np.maximum(vov / (1.0 + p.theta * vov), p.vov_floor)

    # 4) Soft drain variable: v -> Vds for Vds << Vk, v -> Vk for Vds >> Vk.
    v = vknee * np.tanh(vds / vknee)

    # 5) Unified channel charge term (triode + saturation in one C^1 expression).
    ich = p.kp * vov * (v - 0.5 * v * v / vknee)

    # 6) Velocity saturation + gated output conductance.
    ich /= 1.0 + v / p.vsat
    ich *= 1.0 + p.lam * vds * gate

    # 7) Drift conductance and trench kink.
    #    The kink gate must be a one-sided sigmoid, not tanh() itself: tanh is odd,
    #    so a bare tanh((Vgs-Vthk)/kk) injects NEGATIVE current below the kink onset
    #    (about -40 mA at Vgs = 0, Vds = 1 V for the nominal set), which is
    #    non-physical and can seed convergence problems in the off state.
    #    0.5*(1 + tanh(.)) is identical to tanh() once Vgs is a few kk above Vthk
    #    (tanh already saturates at 1.0 in float), and decays smoothly to zero below.
    idrift = p.gds * vds
    kink_gate = 0.5 * (1.0 + np.tanh((vgs - p.vthk) / p.kk))
    ikink = p.gk * vds * kink_gate * gate

    return ich + idrift + ikink


def vdmos_id_level1(vds, vgs, *, kp=28.0, vth0=3.0, ksub=0.05, m_triode=1.2, lam=0.02):
    """LTspice Level-1 VDMOS drain current, as actually implemented.

    IMPORTANT: the triode/saturation knee is at **Vds = Vov / m_triode**, not at
    Vds = Vov. LTspice places the knee there deliberately so the model is C1
    continuous for *any* m_triode: the triode slope
    m*Kp*(Vov - m*Vds) vanishes exactly at Vds = Vov/m, which equals the
    saturation slope. Verified against LTspice 26.0.1 .meas (see
    tests/test_model.py::test_level1_is_c1_continuous):

        Vgs=5, Vto=3 (Vov=2), Kp=2, Lambda=0
        m=1 knee 2.0 V: Id = 4.000000 at 1.9999 / 2.0000 / 2.0001 V
        m=2 knee 1.0 V: Id = 4.000000 at 0.9999 / 1.0000 / 1.0001 V
        m=4 knee 0.5 V: Id = 4.000000 at 0.4999 / 0.5000 / 0.5001 V
        min Id over 0..6 V = 4.0 for every m  ->  no negative-gds region.

    m_triode is therefore a pure *shape* knob: it sets Ron = 1/(m*Kp*Vov) and the
    knee position, but cannot change Id,sat = 0.5*Kp*Vov^2. That inability to
    set Ron, knee shape and saturation level independently -- together with the
    fixed square-law exponent -- is the real Level-1 limitation, not any
    discontinuity.
    """
    vds = np.asarray(vds, dtype=float)
    vgs = np.asarray(vgs, dtype=float)
    vds, vgs = np.broadcast_arrays(vds, vgs)
    z = np.clip((vgs - vth0) / ksub, -700.0, 700.0)
    vov = ksub * np.log1p(np.exp(z))
    knee = vov / m_triode
    tri = m_triode * kp * (vov * vds - 0.5 * m_triode * vds**2) * (1.0 + lam * vds)
    sat = 0.5 * kp * vov**2 * (1.0 + lam * vds)
    return np.where(vds < knee, tri, sat)
