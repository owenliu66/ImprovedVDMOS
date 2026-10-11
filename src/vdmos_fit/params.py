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
}

# Sensible starting free set: everything that shapes the static I-V curves.
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


def vector_bounds(free: Sequence[str]) -> tuple[np.ndarray, np.ndarray]:
    """Bounds for the flat vector, mapped through the log transform."""
    lo, hi = [], []
    for name in free:
        a, b = BOUNDS[name]
        if name in LOG_PARAMS:
            lo.append(np.log(a))
            hi.append(np.log(b))
        else:
            lo.append(a)
            hi.append(b)
    return np.asarray(lo), np.asarray(hi)
