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

    Physical rationale for the gate cap (kcap, m) and the package resistances
    (rd, rs) is in the project README.md ("Gate-overdrive charge cap",
    "Package resistance"); the notes here are units and gotchas only.

    Attributes
    ----------
    kp : float
        Process transconductance parameter [A/V^2]. NOTE: not the Level-1 kp
        (28 -> ~34.7); must be re-extracted when theta/vsat are added.
    vth0 : float   Zero-bias threshold voltage [V].
    ksub : float   Subthreshold swing factor [V] (n*k*T/q).
    kcap : float   Soft cap on the effective gate overdrive [V]; 0 = off.
                   "Off" is reached by leaving kcap out of the free set --
                   NOT by driving it to the lower bound (kcap -> 0+ means
                   "hard cap at zero volts", the opposite of off).
    m : float      Sharpness of the gate-overdrive cap [-]; m >> 1 approaches
                   a hard clip at kcap. Only meaningful when kcap > 0.
    theta : float  Mobility-degradation coefficient [1/V] -> Vk = Vov/(1+theta*Vov).
    vsat : float   Velocity-saturation knee [V]; sets Id,sat ~ Vov**alpha, alpha < 2.
    nd : float     Drain-knee softness exponent [-]; 0 = legacy tanh soft-switch.
                   nd > 0 replaces v = Vk*tanh(Vds/Vk) with the SAME power-law
                   soft-min used for the gate cap: v = Vds/(1 + (Vds/Vk)**nd)**(1/nd).
                   Shape-only knob: dv/dVds(0) = 1 and v -> Vk hold for every nd,
                   so Ron and the Id,sat plateau are untouched while the knee
                   WIDTH scales with nd (nd -> inf is a hard min, nd = 1 is the
                   harmonic soft-min; nd ~ 2.8 reproduces the legacy tanh).
                   "Off" is the default 0 -- freeing it changes the knee shape,
                   so it stays out of DEFAULT_FREE like kcap/m.
    eta : float    DIBL coefficient [1/V].
    lam : float    Channel/drift-length modulation [1/V].
    gk : float     Trench kink (parasitic accumulation channel) conductance [S].
                   One-sided: contributes only for Vgs above vthk.
    vthk : float   Kink gate-onset voltage [V].
    kk : float     Kink gate transition softness [V].
    vk : float     Drain-bias softness [V] for the tanh() gates.
    vov_floor : float  Numerical floor on Vk [V]; keep fixed, not fitted.
    rd : float     Drain-side series resistance [Ohm] (bond wires, leadframe,
                   pins, plus the undepleted drift). Embedded outside the
                   intrinsic channel; datasheet Rds(on) already contains
                   rd + rs, and rd + rs is a LOWER bound on the Rds(on)
                   chord (see README, "Package resistance (rd, rs)").
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
    nd: float = 0.0
    eta: float = 0.02
    lam: float = 0.02
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
# nd likewise: 0 is the legitimate "legacy tanh" value.
LOG_PARAMS: frozenset[str] = frozenset(
    {"kp", "ksub", "theta", "vsat", "eta", "lam", "gk", "kk", "vk", "m"}
)

# Physically sensible search bounds, same units as Params.
# Rationale for the non-obvious entries: README, "Bounds and fitting-space
# decisions".
BOUNDS: dict[str, tuple[float, float]] = {
    "kp": (1e-2, 1e4),
    "vth0": (0.0, 10.0),
    "ksub": (5e-3, 1.0),
    # Linear space; lower bound deliberately well above 0 (a free kcap must
    # never enter the "hard clip at ~0 V" region).
    "kcap": (0.05, 100.0),
    "m": (0.5, 40.0),
    "theta": (1e-4, 5.0),
    "vsat": (0.05, 10000.0),
    # Linear space; lower bound above 0 for the same reason as kcap -- a free
    # nd must never enter the region near the sentinel 0 ("legacy tanh")
    # where a fitted value and the sentinel are numerically indistinguishable.
    # Relaxed 0.5 -> 0.1 (session 10): the fit-B relaxed-bound probe landed
    # at nd = 0.47, just under the old guard, so the data wants a knee softer
    # than the old floor allowed.
    "nd": (0.1, 40.0),
    "eta": (1e-5, 1.0),
    "lam": (1e-12, 1.0),
    "gk": (1e-12, 100.0),
    "vthk": (-5.0, 10.0),
    "kk": (1e-3, 5.0),
    "vk": (0.1, 1300.0),
    "vov_floor": (1e-6, 1e-1),
    # Linear space (the physical answer can be exactly zero); 0.5 Ohm is a
    # loose global ceiling only -- per device cap at datasheet Rds(on).
    "rd": (0.0, 0.5),
    "rs": (0.0, 0.5),
}

# Sensible starting free set: everything that shapes the static I-V curves.
# kcap/m are deliberately absent -- freeing them would silently change every
# existing fit and anchor (see README, "Bounds and fitting-space decisions").
DEFAULT_FREE: tuple[str, ...] = (
    "kp",
    "vth0",
    "theta",
    "vsat",
    "eta",
    "lam",
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
    without mutating the module-level table (e.g. a per-device cap of rd/rs
    at the datasheet Rds(on)).
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
