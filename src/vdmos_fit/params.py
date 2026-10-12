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
    rd : float     Package drain resistance [Ohm] (bond wires, leadframe, pins).
                   Embedded outside the intrinsic channel: the internal drain
                   sees Vds - (rd+rs)*Id. Fitted by default because datasheet
                   Rds(on) is a 2-terminal spec that already contains rd+rs.
                   NOTE: rd + rs is always BELOW the measured Rds(on) chord --
                   the rest of Rds(on) is the intrinsic channel resistance,
                   which the model computes separately. rd+rs is the *lower*
                   bound on Rds(on), approached only as Vgs -> large.
    rs : float     Package source resistance [Ohm]. Also *degenerates* the
                   gate drive: the internal gate-source voltage is
                   Vgs - rs*Id, so gm_eff = gm/(1 + rs*gm).
    kappa : float  Drain-bias coefficient of the drain-side series resistance
                   [1/V], 0 by default (plain rd/rs embedding).

                   With kappa > 0 the drain-side resistance grows with the
                   voltage it supports: R_d(u) = rd*(1 + kappa*u), so the
                   channel sees Vds - u - rs*Id with u = Id*rd/(1 - Id*rd*kappa).
                   Physically this is super-junction JFET/RESURF depletion: the
                   charge-balanced columns deplete as drain bias rises, so the
                   drift specific on-resistance climbs with Vds. With the
                   channel strongly enhanced the terminal characteristic becomes
                   Id -> Vds/(rd*(1 + kappa*Vds)), flattening toward the finite
                   drift-limited plateau Id_max = 1/(rd*kappa). That is the only
                   mechanism in this model that can produce the gate exponent
                   alpha << 1 (measured ~0.1 at Vgs 8-10 V, Vds ~ 19 V on the
                   650 V SJ part) that theta/vsat/lam cannot reach.

                   KNOWN MODELING TENSION: kappa multiplies `rd`, which is
                   defined as *package* resistance (wirebonds, pins) and should
                   not vary with drain bias. The drain-region resistance that
                   actually depletes is a separate physical element. A cleaner
                   form adds an explicit drift resistance (e.g. rdrift with its
                   own kappa) and keeps rd/rs purely parasitic -- see STATUS.
    """

    kp: float = 34.7
    vth0: float = 3.0
    ksub: float = 0.05
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
    kappa: float = 0.0

    def as_dict(self) -> dict:
        return asdict(self)

    def replace(self, **kw) -> Params:
        return Params(**{**asdict(self), **kw})


DEFAULT = Params()

PARAM_NAMES: tuple[str, ...] = tuple(f.name for f in fields(Params))

# Parameters fitted in log space (strictly positive, span decades).
LOG_PARAMS: frozenset[str] = frozenset(
    {"kp", "ksub", "theta", "vsat", "eta", "lam", "gds", "gk", "kk", "vk"}
)

# Physically sensible search bounds, same units as Params.
BOUNDS: dict[str, tuple[float, float]] = {
    "kp": (1e-3, 1e3),
    "vth0": (0.0, 10.0),
    "ksub": (5e-3, 1.0),
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
    "rd": (0.0, 0.5),
    "rs": (0.0, 0.5),
    # Linear space: kappa = 0 means "no drain-bias effect" and is a legitimate
    # answer (plain package resistance). Upper bound keeps rd*kappa*Vds within
    # a sane range over the fitted window.
    "kappa": (0.0, 2.0),
}

# Sensible starting free set: everything that shapes the static I-V curves.
# `kappa` is deliberately NOT here: it is 0 for the low-VT anchor device and
# only needed for super-junction / long-drift parts, where it must be freed
# explicitly (see scripts/fit_st8l65.py). Adding it to the default set would
# silently change every existing fit and anchor.
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
    need a wider search box than the physical defaults).
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
