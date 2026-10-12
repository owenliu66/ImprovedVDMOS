"""Improved VDMOS channel current: one continuous expression, no triode/sat switch.

    I_D = kp*Vov_eff*(v - v^2/(2*Vk)) / (1 + v/vsat) * (1 + lam*Vds*tanh(Vds/Vk))
          + gds*Vds
          + gk*Vds*tanh((Vgs - Vth_k)/kk)*tanh(Vds/Vk)

    Vth     = Vth0 - eta*Vds*tanh(Vds/Vk)          DIBL (bounded)
    Vov     = ksub*ln(1 + exp((Vgs - Vth)/ksub))   smooth subthreshold -> strong inversion
    Vov_eff = Vov/(1 + (Vov/kcap)**m)**(1/m)       soft charge cap (kcap = 0 -> off)
    Vk      = Vov_eff/(1 + theta*Vov_eff)          knee voltage (field-dependent mobility)
    v       = Vk*tanh(Vds/Vk)                      soft drain variable (replaces the knee switch)

Limits
------
Vds << Vk :  I_D -> kp*Vov_eff*Vds*(1 + ...)                   ohmic, m_triode = 1
Vds >> Vk :  I_D -> 0.5*kp*Vov_eff^2*(0.5+th*Vov_eff)/(1+th*Vov_eff)^2 / (1+Vk/vsat)

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

Gate-overdrive cap (kcap, m) -- finite charge in the SJ pillars
--------------------------------------------------------------
On a super-junction part the measured gate exponent alpha = d ln Id / d ln Vov is
strongly bimodal: ~1.6-1.8 at low overdrive, collapsing to ~0.1 above Vgs ~ 8 V
(ST8L65N044M9, 25 C pulsed: 1.68 / 1.61 / 1.00 / 0.16 / 0.07 across Vgs
5->6->7->8->9->10 at the top of each curve). The lightly doped N pillars under
the trench hold only a finite amount of free charge, so once the surface channel
is strongly enhanced the mobile charge in the current path stops growing with the
gate -- as if two channels in series had one of them saturate. No setting of
theta/vsat/lam can make alpha drop by an order of magnitude.

The previous fix for this was `kappa`, a drain-bias-dependent series resistance.
It is removed: it is a *voltage divider*, not a charge limit, and the fit made
that obvious -- rd had to grow to ~25 mOhm and then carried 75 % of the terminal
voltage already at the ohmic Ron anchor (Vds = 0.483 V, Id = 14 A), i.e. kappa
was acting as an implicit pseudo-cap while pretending to be a package resistance.
The cap is the honest mechanism:

    Vov_eff = Vov / (1 + (Vov/kcap)**m)**(1/m)

Identity for Vov << kcap (subthreshold and the Ron limit untouched), asymptote
kcap for Vov >> kcap, and 0 <= dVov_eff/dVov <= 1 everywhere -- which is what
keeps the package-R solver's bracket invariant dg/dId <= -1 intact. m sets the
sharpness; the measured collapse is abrupt, so m ~ 10 rather than the m = 1 of a
tanh. The cap is applied BEFORE Vk (variant A) so charge and knee saturate
together: with a depleted pillar the pinch-off voltage is set by the depleted
column, not by the surface charge. `kcap = 0` (default) is the exact identity and
leaves every pre-cap anchor bit-identical.

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

    # 2b) Soft charge cap: finite mobile charge in the super-junction pillars.
    # kcap = 0 is the identity, so the pre-cap model is reproduced bit for bit.
    vov = _soft_cap(vov, p.kcap, p.m)

    # 3) Knee voltage from vertical-field mobility degradation. Built from the
    # CAPPED overdrive (variant A): with a depleted pillar the pinch-off voltage
    # is set by the depleted column, not by the surface charge.
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


def _soft_cap(vov, kcap: float, m: float):
    """Soft power-law cap vov/(1 + (vov/kcap)**m)**(1/m), evaluated in log space.

    Identity for vov << kcap, asymptote kcap for vov >> kcap, and
    0 <= dVov_eff/dVov <= 1 everywhere (the derivative is (1+u**m)**(-1/m-1)
    with u = vov/kcap). That bound is what keeps the package-R Newton bracket
    valid: the channel's dI/dVds can only shrink, never blow up. The log-space
    form avoids overflowing (vov/kcap)**m, since m ~ 10-40 is a fitted regime,
    not a hypothetical one. kcap <= 0 returns vov unchanged -- the exact pre-cap
    model, bit for bit.
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

    Bracket: Id in [0, Vds/(rd+rs)]. Bound: Id <= Vds/(rd+rs), because the
    terminal resistance is never below rd + rs.
    g(Id) = f(.) - Id changes sign across the bracket: f has the sign of the
    internal Vds and f(0, .) = 0 exactly (every model term is gated by Vds or
    tanh(Vds/Vk)), which also makes clamping the internal Vds at 0 harmless when
    a trial point overshoots. dg/dId <= -1 keeps Newton well-conditioned; the
    gate-overdrive cap preserves this because 0 <= dVov_eff/dVov <= 1 can only
    shrink the channel's dI/dVds. Bisection is the fallback whenever a Newton
    step leaves the bracket or makes no progress. Fixed-point iteration is NOT
    used -- see module docstring.

    Raises RuntimeError if any grid point fails to converge; never silently
    returns an unconverged current.

    Scope: Vds >= 0 only. The channel model is forward-polarity by construction
    (the velocity-saturation divisor 1 + v/vsat is singular near Vds ~ -vsat,
    and the numerator is not odd), so the bracket argument fails for Vds < 0.
    Reverse operation needs the body diode (STATUS problem 7), not this solver.
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

        # Stall detector: a Newton two-cycle (e.g. gate-degeneration cliffs from
        # rs*Id near the cap knee) can keep both endpoints inside the bracket and
        # shrink it only by ~1e-5 per iteration. Force bisection until the bracket
        # halves again.
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
        # Accept Newton only if it stays inside the bracket, actually makes
        # progress, and the bracket has not stalled (guards a zero-step Newton
        # and the two-cycle above).
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
