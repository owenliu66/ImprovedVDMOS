"""Parameter container, physical bounds and log-space packing for the fitter."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import asdict, dataclass, fields

import numpy as np

__all__ = [
    "Params",
    "DEFAULT",
    "PARAM_NAMES",
    "LOG_PARAMS",
    "BOUNDS",
    "DEFAULT_FREE",
    "pack",
    "unpack",
]


@dataclass(frozen=True)
class Params:
    """Improved VDMOS drain-current parameters.

    Attributes
    ----------
    kp : float
        Process transconductance parameter [A/V^2]. NOTE: must be re-extracted
        when theta/vsat are added -- it is not the Level-1 kp (28 -> ~34.7).
    vth0 : float   Zero-bias threshold voltage [V].
    ksub : float   Subthreshold swing factor [V] (n*k*T/q).
    kcap : float   Soft cap on the effective gate overdrive [V]; 0 = off.

                   Vov_eff = Vov / (1 + (Vov/kcap)**m)**(1/m)

                   Physical picture: in a super-junction VDMOS the lightly doped
                   N pillars that carry the drain current hold only a finite
                   amount of free charge. Once the surface channel is strongly
                   enhanced, driving the gate harder cannot add much more mobile
                   charge in the pillar, so the *effective* overdrive saturates
                   even though Vgs keeps rising. The measured signature on the
                   650 V SJ part is a bimodal gate exponent: alpha = d ln Id /
                   d ln(Vgs - Vth) is 1.6-1.8 at low overdrive and collapses to
                   ~0.1 above Vgs ~ 8 V. No setting of theta/vsat/lam can
                   produce that; the cap can, and it does so with a physically
                   bounded output conductance.

                   Applied to Vov BEFORE the knee voltage (variant A), so the
                   charge factor and the pinch-off/knee voltage saturate
                   together. That is the self-consistent choice for a depleted
                   pillar: once the pillar charge caps, its pinch-off voltage is
                   set by the depleted column, not by the surface charge.

                   kcap = 0 (default) selects the exact identity, so every
                   pre-cap anchor stays bit-identical. "Off" is reached by
                   leaving kcap out of the free set -- NOT by driving it to the
                   lower bound, because kcap -> 0+ means "hard cap at zero
                   volts", which is the opposite of off.
    m : float      Sharpness of the gate-overdrive cap [-]. m = 1 is a gentle
                   roll-off; m >> 1 approaches a hard clip at kcap. Only
                   meaningful when kcap > 0. Fitted m ~ 10 on the ST8L65N044M9
                   data: the alpha collapse is abrupt (1.0 -> 0.15 -> 0.07 over
                   two gate steps), which a tanh-shaped cap (m = 1) cannot
                   reproduce.
    theta : float  Mobility-degradation coefficient [1/V] -> Vk = Vov/(1+theta*Vov).
    vsat : float   Velocity-saturation knee [V]; sets Id,sat ~ Vov**alpha, alpha < 2.
    eta : float    DIBL coefficient [1/V].
    lam : float    Channel/drift-length modulation [1/V].
    gds : float    Drift/substrate output conductance [S].
    gk : float     Trench kink (parasitic accumulation channel) conductance [S].
                   One-sided: contributes only for Vgs above vthk.
    vthk : float   Kink gate-onset voltage [V].
    kk : float     Kink gate transition softness [V].
    vk : float     Drain-bias softness [V] for the tanh() gates.
    vov_floor : float  Numerical floor on Vk [V]; keep fixed, not fitted.
    rd : float     Drain-side series resistance [Ohm] (bond wires, leadframe,
                   pins, plus the undepleted part of the drift region).
                   Embedded outside the intrinsic channel: the internal drain
                   sees Vds - (rd+rs)*Id. Fitted by default because datasheet
                   Rds(on) is a 2-terminal spec that already contains rd+rs.
                   NOTE: rd + rs is always BELOW the measured Rds(on) chord --
                   the rest of Rds(on) is the intrinsic channel resistance,
                   which the model computes separately. rd+rs is the *lower*
                   bound on Rds(on), approached only as Vgs -> large.
                   With the gate cap in place there is no bias-dependent drain
                   element left to distinguish "package" from "undepleted
                   drift": the two are exactly degenerate (correlation +1.000
                   when split), so they are kept as one parameter.
    rs : float     Source series resistance [Ohm]. Also *degenerates* the
                   gate drive: the internal gate-source voltage is
                   Vgs - rs*Id, so gm_eff = gm/(1 + rs*gm).
    """

    kp: float = 34.7
    vth0: float = 3.0
    ksub: float = 0.05
    kcap: float = 0.0
    m: float = 1.0
    theta: float = 0.015
    vsat: float = 1.5
    eta: float = 0.02
    lam: float = 0.02
    gds: float = 1e-3
    gk: float = 0.4
    vthk: float = 1.6
    kk: float = 0.1
    vk: float = 10.0
    vov_floor: float = 1e-3
    rd: float = 0.0
    rs: float = 0.0

    def as_dict(self) -> dict:
        return asdict(self)

    def replace(self, **kw) -> Params:
        return Params(**{**asdict(self), **kw})


DEFAULT = Params()

PARAM_NAMES: tuple[str, ...] = tuple(f.name for f in fields(Params))

# Parameters fitted in log space (strictly positive, span decades).
# kcap is NOT here: 0 is a legitimate value ("cap off") and cannot be logged.
LOG_PARAMS: frozenset[str] = frozenset(
    {"kp", "ksub", "theta", "vsat", "eta", "lam", "gds", "gk", "kk", "vk", "m"}
)

# Physically sensible search bounds, same units as Params.
BOUNDS: dict[str, tuple[float, float]] = {
    "kp": (1e-3, 1e3),
    "vth0": (0.0, 10.0),
    "ksub": (5e-3, 1.0),
    # Linear space: kcap = 0 means "off" and cannot be logged. The lower bound
    # is deliberately well above 0 -- between 0 and 0.05 V the cap is a hard
    # clip at nearly zero volts, so a free kcap must never wander there.
    "kcap": (0.05, 100.0),
    "m": (0.5, 40.0),
    "theta": (1e-4, 5.0),
    "vsat": (0.05, 100.0),
    "eta": (1e-5, 1.0),
    "lam": (1e-6, 5.0),
    "gds": (1e-9, 10.0),
    "gk": (1e-12, 100.0),
    "vthk": (-5.0, 10.0),
    "kk": (1e-3, 5.0),
    "vk": (0.1, 200.0),
    "vov_floor": (1e-6, 1e-1),
    # Linear space, not log: the physical answer can be exactly zero, and
    # datasheet package resistance never approaches 0.5 Ohm on a power FET.
    # 0.5 Ohm is a loose *global* ceiling only; per device, cap rd and rs at
    # the datasheet Rds(on) (see scripts/fit_st8l65.py DEVICE_BOUNDS).
    "rd": (0.0, 0.5),
    "rs": (0.0, 0.5),
}

# Sensible starting free set: everything that shapes the static I-V curves.
# `kcap`/`m` are deliberately NOT here: the cap is 0 (off) for the low-VT
# anchor device and only needed for super-junction / long-drift parts, where it
# must be freed explicitly (see scripts/fit_st8l65.py). Adding them to the
# default set would silently change every existing fit and anchor.
DEFAULT_FREE: tuple[str, ...] = (
    "kp",
    "vth0",
    "theta",
    "vsat",
    "eta",
    "lam",
    "gds",
    "gk",
    "vthk",
    "rd",
    "rs",
)


def pack(params: Params, free: Sequence[str]) -> np.ndarray:
    """Free parameters -> flat optimisation vector (log space where appropriate)."""
    out = []
    for name in free:
        val = getattr(params, name)
        out.append(np.log(val) if name in LOG_PARAMS else val)
    return np.asarray(out, dtype=float)


def unpack(x: np.ndarray, free: Sequence[str], base: Params = DEFAULT) -> Params:
    """Flat optimisation vector -> full :class:`Params` (non-free values from `base`)."""
    if len(x) != len(free):
        raise ValueError(f"expected {len(free)} values, got {len(x)}")
    values = asdict(base)
    for name, xi in zip(free, x, strict=True):
        values[name] = float(np.exp(xi)) if name in LOG_PARAMS else float(xi)
    return Params(**values)


def vector_bounds(
    free: Sequence[str],
    overrides: dict[str, tuple[float, float]] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Bounds for the flat vector, mapped through the log transform.

    `overrides` replaces individual entries of :data:`BOUNDS` for one call
    without mutating the module-level table (useful for exploratory fits that
    need a wider or tighter search box than the physical defaults, e.g. a
    per-device cap of rd/rs at the datasheet Rds(on)).
    """
    table = {**BOUNDS, **(overrides or {})}
    lo, hi = [], []
    for name in free:
        a, b = table[name]
        if name in LOG_PARAMS:
            lo.append(np.log(a))
            hi.append(np.log(b))
        else:
            lo.append(a)
            hi.append(b)
    return np.asarray(lo), np.asarray(hi)
