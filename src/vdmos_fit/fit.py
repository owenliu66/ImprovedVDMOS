"""Least-squares parameter extraction for the improved VDMOS model."""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

from .data import CurveSet
from .model import vdmos_id
from .params import BOUNDS, DEFAULT, DEFAULT_FREE, Params, pack, unpack, vector_bounds

__all__ = [
    "FitResult",
    "Identifiability",
    "fit",
    "identifiability",
    "parameter_uncertainty",
    "residuals",
    "sensitivity_report",
]


@dataclass(frozen=True)
class FitResult:
    params: Params
    cost: float
    rmse_rel: float
    n_points: int
    n_free: int
    success: bool
    message: str
    jac: np.ndarray | None = None
    free: tuple[str, ...] = ()
    bounds_used: dict[str, tuple[float, float]] | None = None

    def on_bounds(self, rel_tol: float = 1e-3) -> dict[str, str]:
        """Parameters that ended on a search bound, and which one.

        A parameter pinned at a bound is not a physical result -- the data
        wanted to push it further than the box allowed (or it has no
        sensitivity there). Always check this before quoting fitted values.
        """
        if not self.bounds_used:
            return {}
        out: dict[str, str] = {}
        for name in self.free:
            lo, hi = self.bounds_used[name]
            val = getattr(self.params, name)
            span = hi - lo
            if abs(val - lo) <= rel_tol * max(span, abs(lo), 1e-12):
                out[name] = f"LOWER {lo:g}"
            elif abs(val - hi) <= rel_tol * max(span, abs(hi), 1e-12):
                out[name] = f"UPPER {hi:g}"
        return out

    def summary(self) -> str:
        lines = [
            f"fit: {'converged' if self.success else 'NOT converged'} "
            f"({self.n_points} pts, {self.n_free} free params)",
            f"cost = {self.cost:.4e}   relative RMSE = {self.rmse_rel * 100:.2f} %",
            "parameters:",
        ]
        pinned = self.on_bounds()
        for name, val in self.params.as_dict().items():
            tag = f"   <-- pinned at {pinned[name]}" if name in pinned else ""
            lines.append(f"  {name:<10s} = {val:.6g}{tag}")
        return "\n".join(lines)


def residuals(
    x: np.ndarray,
    free: Sequence[str],
    data: CurveSet,
    base: Params = DEFAULT,
    *,
    id_floor: float = 1e-3,
    weight: str = "relative",
) -> np.ndarray:
    """Residual vector for `least_squares`.

    weight = "relative" : (Im - Is)/max(|Is|, id_floor)  -> log-ish, treats
        subthreshold and high-current points as equally important.
    weight = "abs"      : Amperes; dominated by the high-current corner, which
        is what a Ron/Rds(on) spec cares about.
    """
    p = unpack(x, free, base)
    model = vdmos_id(data.vds, data.vgs, p)
    err = model - data.ids
    if weight == "relative":
        err = err / np.maximum(np.abs(data.ids), id_floor)
    elif weight != "abs":
        raise ValueError(f"unknown weight {weight!r}")
    return err


def fit(
    data: CurveSet,
    *,
    start: Params = DEFAULT,
    free: Sequence[str] = DEFAULT_FREE,
    weight: str = "relative",
    max_nfev: int = 4000,
    verbose: bool = False,
    bounds_overrides: dict[str, tuple[float, float]] | None = None,
) -> FitResult:
    """Extract model parameters from a static I-V dataset.

    Positive parameters are optimised in log space with physical bounds, which
    keeps the trust-region solver from wandering through non-physical regions
    (negative kp, zero ksub) where the model becomes singular. Pass
    `bounds_overrides` to widen individual bounds for an exploratory fit; the
    result records them so `FitResult.on_bounds()` can flag pinned parameters.
    """
    free = tuple(free)
    x0 = pack(start, free)
    lo, hi = vector_bounds(free, bounds_overrides)
    x0 = np.clip(x0, lo + 1e-9, hi - 1e-9)

    sol = least_squares(
        residuals,
        x0,
        bounds=(lo, hi),
        args=(free, data, start),
        kwargs={"weight": weight},
        method="trf",
        x_scale="jac",
        max_nfev=max_nfev,
        verbose=2 if verbose else 0,
    )

    p = unpack(sol.x, free, start)
    rel = residuals(sol.x, free, data, start, weight="relative")
    table = {**BOUNDS, **(bounds_overrides or {})}
    return FitResult(
        params=p,
        cost=float(sol.cost),
        rmse_rel=float(np.sqrt(np.mean(rel**2))),
        n_points=len(data),
        n_free=len(free),
        success=bool(sol.success),
        message=str(sol.message),
        jac=sol.jac,
        free=free,
        bounds_used={n: table[n] for n in free},
    )


def parameter_uncertainty(
    result: FitResult,
    free: Sequence[str] = DEFAULT_FREE,
) -> dict[str, float]:
    """1-sigma relative uncertainty per fitted parameter (see :func:`identifiability`).

    Returns `math.inf` for parameters the data has no sensitivity to.
    """
    return identifiability(result, free).sigma_rel


@dataclass(frozen=True)
class Identifiability:
    """How well the fitted data actually constrains each parameter.

    sigma_rel
        1-sigma relative uncertainty, in the same relative units the parameter
        was optimised in (log space -> relative). `inf` means the parameter has
        **zero sensitivity** over the sampled range: the fit is free to move it
        anywhere and the residual will not notice.
    rank, n_free
        Rank of the Jacobian. rank < n_free means the parameter set is
        under-determined by this dataset.
    cond
        Condition number over the non-trivial singular directions. Above ~1e6
        the linearised error bars are unreliable.
    worst_pairs
        Parameter pairs with |correlation| > 0.8 -- these trade off against each
        other, so their individual values are meaningless even when the fit is
        excellent. Fix one from independent data and refit.
    """

    sigma_rel: dict[str, float]
    rank: int
    n_free: int
    cond: float
    unidentifiable: tuple[str, ...]
    worst_pairs: tuple[tuple[str, str, float], ...]

    def summary(self) -> str:
        lines = [
            f"Jacobian rank {self.rank}/{self.n_free}   condition number {self.cond:.3g}"
            + ("   <-- RANK DEFICIENT" if self.rank < self.n_free else ""),
            "parameter            sigma_rel   note",
        ]
        for name, s in self.sigma_rel.items():
            note = "NOT SENSED by this data" if math.isinf(s) else ""
            lines.append(f"  {name:<18s} {s:9.3g}   {note}")
        if self.worst_pairs:
            lines.append("strongly correlated pairs (trade off against each other):")
            for a, b, c in self.worst_pairs:
                lines.append(f"  {a} <-> {b}: {c:+.3f}")
        return "\n".join(lines)


def identifiability(
    result: FitResult,
    free: Sequence[str] = DEFAULT_FREE,
) -> Identifiability:
    """Diagnose which parameters the dataset can actually determine.

    Built from the Jacobian at the optimum. Assumes unit-weight residuals (the
    relative-RMSE scaling), so these are order-of-magnitude identifiability
    numbers, not statistical confidence intervals.

    A zero-sensitivity column is reported as `inf`, never as 0. That distinction
    matters: a bare `pinv(J.T @ J)` returns *zero variance* for a column of
    zeros, which reads as "perfectly determined" when the truth is the exact
    opposite. For example vthk with kk = 0.1 V is fully saturated across a data
    set whose lowest Vgs is 2 V, so the kink transition is never sampled and
    vthk is unconstrained.

    The dead test is RELATIVE to the strongest Jacobian column (1e-5 of it),
    not an absolute threshold: a real datasheet fit has columns spanning many
    orders of magnitude, and an absolute cutoff either misses genuinely dead
    parameters (ST8L65 `ksub` measured 7.8e-7 of max yet was called live and
    given a confident small sigma) or mislabels merely-weak ones. Correlations
    and the condition number are computed over live columns only, so a dead
    column can no longer fake a +/-1.000 pair or inflate cond.
    """
    jac = result.jac
    free = tuple(free)
    if jac is None:
        return Identifiability({n: math.inf for n in free}, 0, len(free), math.inf, free, ())

    col_norm = np.linalg.norm(jac, axis=0)
    col_max = float(col_norm.max(initial=0.0))
    # Dead-column test is RELATIVE to the strongest column. The old absolute
    # 1e-10 * max was ~5 orders of magnitude too permissive: on the ST8L65
    # super-junction fit ksub's column measured 7.8e-7 of max -- provably dead
    # (pinning it at 0.007/0.02/0.04/0.06/0.09 gives identical RMSE and
    # identical other parameters) -- yet it was classified live and reported
    # sigma_rel = 0.038, i.e. "determined to 4 %" for a parameter carrying
    # zero information. That is precisely the pinv-returns-zero-variance trap
    # this function exists to avoid, and it also manufactured a spurious
    # ksub <-> vk = +1.000 correlation and inflated cond by ~3 orders.
    # 1e-5 sits in the empty gap between real and dead: the weakest genuinely
    # live columns measured in this project are ~3e-4 of max (theta/gds/vthk
    # on the synthetic grid), while dead ones land at ~1e-6.
    DEAD_REL_TOL = 1e-5
    tol = max(col_max * DEAD_REL_TOL, max(col_max, 1.0) * 1e-10)
    dead = col_norm < tol

    svals = np.linalg.svd(jac, compute_uv=False)
    live = svals[svals > tol]
    rank = int(live.size)
    cond = float(svals.max() / live.min()) if live.size else math.inf

    cov = np.linalg.pinv(jac.T @ jac)
    sigma: dict[str, float] = {}
    for i, name in enumerate(free):
        sigma[name] = math.inf if dead[i] else float(np.sqrt(max(cov[i, i], 0.0)))

    unident = tuple(n for n, d in zip(free, dead, strict=True) if d)

    # Correlations only among parameters that are actually sensed.
    ok = ~dead
    pairs: list[tuple[str, str, float]] = []
    if ok.any():
        sd = np.sqrt(np.clip(np.diag(cov), 0.0, None))
        sd[~ok] = 1.0
        corr = cov / np.outer(sd, sd)
        for i, a in enumerate(free):
            for j, b in enumerate(free):
                if j > i and ok[i] and ok[j] and abs(corr[i, j]) > 0.8:
                    pairs.append((a, b, float(corr[i, j])))
    return Identifiability(sigma, rank, len(free), cond, unident, tuple(pairs))


def sensitivity_report(
    params: Params,
    data: CurveSet,
    free: Sequence[str] = DEFAULT_FREE,
    *,
    rel_step: float = 1e-3,
) -> dict[str, float]:
    """Mean |d ln Id / d ln p| over the dataset, evaluated before fitting.

    Cheap diagnostic for deciding what to fit at all: a parameter whose current
    value changes Id by <1e-3 relative across the whole data range cannot be
    extracted from it. For instance vthk only matters where the data actually
    brackets the kink onset, and gds only matters in the off state.

    Parameters sitting at their (non-negative) lower bound -- rd, rs at 0 --
    get a one-sided probe of size `rel_step` in absolute units, and the
    reported number is the mean *relative* Id change that the probe produces.
    (The log-space formula d ln Id / d ln p is meaningless at p = 0, and a
    two-sided probe would request negative resistance.)
    """
    base = vdmos_id(data.vds, data.vgs, params)
    scale = np.maximum(np.abs(base), 1e-12)
    out: dict[str, float] = {}
    for name in free:
        p0 = getattr(params, name)
        lo_b = BOUNDS.get(name, (-np.inf, np.inf))[0]
        if abs(p0) > 1e-12 and p0 - rel_step * abs(p0) >= lo_b:
            step = rel_step * abs(p0)
            up = vdmos_id(data.vds, data.vgs, params.replace(**{name: p0 + step}))
            dn = vdmos_id(data.vds, data.vgs, params.replace(**{name: p0 - step}))
            sens = np.abs(up - dn) / (2.0 * step) / scale * abs(p0)
        else:
            step = rel_step
            up = vdmos_id(data.vds, data.vgs, params.replace(**{name: p0 + step}))
            sens = np.abs(up - base) / step / scale
        out[name] = float(np.mean(sens))
    return out
