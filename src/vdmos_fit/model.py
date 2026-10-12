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

Package resistance (rd, rs)
---------------------------
The package drain/source resistances sit OUTSIDE the intrinsic device, so the
internal transistor nodes see

    Vds' = Vds - (rd + rs)*Id        Vgs' = Vgs - rs*Id
    Id   = f_intrinsic(Vds', Vgs')   <- implicit scalar equation per data point

This is the same embedding LTspice applies to its Level-1 Rd/Rs model-card
parameters, so `vdmos_id_level1()` gets the identical treatment for fair A/B.
rs additionally *degenerates* the gate (gm_eff = gm/(1 + rs*gm)); it is not a
pure series drop.

Drain resistance that grows with drain bias (kappa)
---------------------------------------------------
`rd` alone cannot describe a super-junction / long-drift device: measured there
(Vgs 8-10 V, Vds ~ 19 V) the gate exponent collapses to alpha ~ 0.1 and the
output curves go nearly flat and gate-independent, which no setting of
theta/vsat/lam can produce. Setting `kappa` > 0 makes the drain-side series
resistance grow with the voltage *it* supports:

    u     = voltage across the drain-side resistance
    R_d(u)= rd*(1 + kappa*u)   ->  u = Id*rd / (1 - Id*rd*kappa)   (closed form)
    Vds'  = Vds - u - rs*Id     Vgs' = Vgs - rs*Id

With the channel strongly enhanced the terminal characteristic becomes
Id -> Vds/(rd*(1 + kappa*Vds)), i.e. it flattens toward the finite drift-limited
plateau Id_max = 1/(rd*kappa) -- the RESURF/JFET-depletion behaviour. `kappa = 0`
(the default) reproduces the plain rd/rs embedding exactly.

The implicit equation is solved with bracketed safeguarded Newton, NOT fixed-point
iteration: in the ohmic region the fixed-point contraction factor is
(rdp-rs)*kp*Vov, which exceeds 1 for realistic values (0.005 * 34.7 * 9 = 1.6 for
DEFAULT) and would diverge. With rd = rs = 0 (the default) the explicit fast path
is taken and results are bit-identical to the pre-package-R model.

The package-R path requires Vds >= 0 (forward polarity); the channel-only model
has no defined reverse characteristic until the body diode is added.
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
        vds, vgs, p.rd, p.rs, lambda d, g: _channel_current(d, g, p), kappa=p.kappa
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


def _solve_package_r(
    vds: np.ndarray,
    vgs: np.ndarray,
    rd: float,
    rs: float,
    f,
    *,
    kappa: float = 0.0,
    atol: float = 1e-12,
    rtol: float = 1e-10,
    max_iter: int = 60,
) -> np.ndarray:
    """Solve Id = f(Vds - R_d*Id - rs*Id, Vgs - rs*Id) for the terminal currents.

    R_d = rd for kappa = 0 (plain package resistance), else rd*(1 + kappa*u) with
    u the voltage the drain-side resistance supports, giving the closed form
    u = Id*rd/(1 - Id*rd*kappa). kappa = 0 takes the original code path, so
    existing results are bit-identical.

    Bracket: Id in [0, Id_hi]. Two independent upper bounds, take the min:
      (a) Id <= Vds/(rd+rs), because R_d >= rd so the terminal resistance is
          never below rd + rs.
      (b) With kappa > 0 the channel must support a non-negative internal Vds,
          so u = Id*rd/(1 - Id*rd*kappa) <= Vds, which rearranges to
          Id <= Vds/(rd*(1 + kappa*Vds)). This is much tighter than (a) at high
          Vds and, crucially, stays far below the u(Id) pole at 1/(rd*kappa).
          Without (b), Newton can sit near the pole where du/dId diverges, the
          Newton step becomes vanishingly small, and the solve stalls.
    g(Id) = f(.) - Id changes sign across the bracket: f has the sign of the
    internal Vds and f(0, .) = 0 exactly (every model term is gated by Vds or
    tanh(Vds/Vk)), which also makes clamping the internal Vds at 0 harmless when
    a trial point overshoots. dg/dId <= -1 keeps Newton well-conditioned;
    bisection is the fallback whenever a Newton step leaves the bracket or makes
    no progress. Fixed-point iteration is NOT used -- see module docstring.

    Raises RuntimeError if any grid point fails to converge; never silently
    returns an unconverged current.

    Scope: Vds >= 0 only. The channel model is forward-polarity by construction
    (the velocity-saturation divisor 1 + v/vsat is singular near Vds ~ -vsat,
    and the numerator is not odd), so the bracket argument fails for Vds < 0.
    Reverse operation needs the body diode (STATUS problem 7), not this solver.
    """
    if rd < 0.0 or rs < 0.0:
        raise ValueError(f"rd and rs must be >= 0 (got rd={rd}, rs={rs})")
    if kappa < 0.0:
        raise ValueError(f"kappa must be >= 0 (got {kappa})")
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
    if kappa > 0.0 and rd > 0.0:
        # Bound (b): keeps the iterate well below the u(Id) pole.
        hi = np.minimum(hi, vds / (rd * (1.0 + kappa * vds)))
    idc = 0.5 * (lo + hi)

    g = None
    for _ in range(max_iter):
        if kappa == 0.0:
            vd_i = vds - rtot * idc
        else:
            u = idc * rd / (1.0 - idc * rd * kappa)
            # A trial point may overshoot; the channel is then fully shorted and
            # f(0, .) = 0 exactly, which drives the bracket back down.
            vd_i = np.maximum(vds - u - rs * idc, 0.0)
        vg_i = vgs - rs * idc
        g = f(vd_i, vg_i) - idc

        # g > 0 -> root above idc; g < 0 -> root below (g is strictly decreasing).
        lo = np.where(g > 0.0, idc, lo)
        hi = np.where(g < 0.0, idc, hi)

        done = np.abs(g) <= atol + rtol * np.abs(idc)
        if done.all():
            return idc

        # Central-difference partials of f (model is C1).
        hd = 1e-6 * np.maximum(1.0, np.abs(vd_i))
        hg = 1e-6 * np.maximum(1.0, np.abs(vg_i))
        fd = (f(vd_i + hd, vg_i) - f(vd_i - hd, vg_i)) / (2.0 * hd)
        fg = (f(vd_i, vg_i + hg) - f(vd_i, vg_i - hg)) / (2.0 * hg)
        if kappa == 0.0:
            # Exactly the original expression -- keeps kappa=0 bit-identical.
            dg = -rtot * fd - rs * fg - 1.0  # <= -1 always
        else:
            # d(Vds_internal)/dId = -(du/dId + rs); du/dId = rd/(1-Id*rd*kappa)^2
            dudI = rd / (1.0 - idc * rd * kappa) ** 2  # >= rd
            dg = -((dudI + rs) * fd + rs * fg) - 1.0  # <= -1 always

        cand = idc - g / dg
        # Accept Newton only if it stays inside the bracket AND actually makes
        # progress; otherwise bisect. The progress test matters when kappa > 0:
        # near the u(Id) pole dg diverges, so a Newton step can be in-bracket
        # yet arbitrarily small and the iteration would stall.
        inside = (cand > lo) & (cand < hi)
        moved = np.abs(cand - idc) > 1e-12 * (hi - lo)
        cand = np.where(inside & moved, cand, 0.5 * (lo + hi))
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

    rd/rs: LTspice's Level-1 model card also has Rd and Rs, applied as the same
    external series embedding as :func:`vdmos_id` (implicit solve, shared
    solver). Pass rd = rs = 0 for the bare equation reference.
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
