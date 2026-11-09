"""Least-squares parameter extraction for the improved VDMOS model."""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

from .data import CurveSet, DiodeCurveSet
from .diode import vdmos_diode
from .model import vdmos_id
from .params import BOUNDS, DEFAULT, DEFAULT_FREE, Params, pack, unpack, vector_bounds

__all__ = [
    "FitResult",
    "Identifiability",
    "fit",
    "fit_diode",
    "diode_residuals",
    "diode_sensitivity_report",
    "identifiability",
    "parameter_uncertainty",
    "residuals",
    "sensitivity_report",
]


# "At zero" floor for the on_bounds test (parameters whose bound is 0, like
# r_d/r_s, count as pinned below rel_tol * ZERO_FLOOR); rationale in README,
# "Fitter diagnostics".
ZERO_FLOOR = 1e-5


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

        A parameter pinned at a bound is not a physical result -- always
        check this before quoting fitted values. The test is relative to the
        *bound's own magnitude* (via max(|bound|, |value|, ZERO_FLOOR)),
        never to the width of the search box; rationale in README,
        "Fitter diagnostics".
        """
        if not self.bounds_used:
            return {}
        out: dict[str, str] = {}
        for name in self.free:
            lo, hi = self.bounds_used[name]
            val = getattr(self.params, name)
            if abs(val - lo) <= rel_tol * max(abs(lo), abs(val), ZERO_FLOOR):
                out[name] = f"LOWER {lo:g}"
            elif abs(val - hi) <= rel_tol * max(abs(hi), abs(val), ZERO_FLOOR):
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
            lines.append(f"  {name:<12s} = {val:.6g}{tag}")
        return "\n".join(lines)


def residuals_p(
    p: Params,
    data: CurveSet,
    *,
    i_d_floor: float = 1e-3,
    weight: str = "relative",
) -> np.ndarray:
    """Residual vector for a ready-made Params (see :func:`residuals`)."""
    model = vdmos_id(data.v_ds, data.v_gs, p)
    err = model - data.i_d
    if weight == "relative":
        err = err / np.maximum(np.abs(data.i_d), i_d_floor)
    elif weight != "abs":
        raise ValueError(f"unknown weight {weight!r}")
    return err


def residuals(
    x: np.ndarray,
    free: Sequence[str],
    data: CurveSet,
    base: Params = DEFAULT,
    *,
    i_d_floor: float = 1e-3,
    weight: str = "relative",
) -> np.ndarray:
    """Residual vector for `least_squares`.

    weight = "relative" : (Im - Is)/max(|Is|, i_d_floor)  -> log-ish, treats
        subthreshold and high-current points as equally important.
    weight = "abs"      : Amperes; dominated by the high-current corner, which
        is what a Ron/Rds(on) spec cares about.
    """
    return residuals_p(unpack(x, free, base), data, i_d_floor=i_d_floor,
                       weight=weight)


def fit(
    data: CurveSet,
    *,
    start: Params = DEFAULT,
    free: Sequence[str] = DEFAULT_FREE,
    weight: str = "relative",
    max_nfev: int = 4000,
    verbose: bool = False,
    bounds_overrides: dict[str, tuple[float, float]] | None = None,
    progress: Callable[[str], None] | None = None,
) -> FitResult:
    """Extract model parameters from a static I-V dataset.

    Positive parameters are optimised in log space with physical bounds, which
    keeps the trust-region solver from wandering through non-physical regions
    (negative k_p, zero v_sub) where the model becomes singular. Pass
    `bounds_overrides` to widen individual bounds for an exploratory fit; the
    result records them so `FitResult.on_bounds()` can flag pinned parameters.

    `progress(msg)` -- optional callback invoked every 25 function
    evaluations with a convergence line ("N evals, cost X"), so a GUI can
    show live progress during long fits.

    With `r_sk` free the optimisation column is the RATIO f_sk = r_sk/r_s in
    [0, 1] (session-22 kelvin doctrine): r_sk = f_sk * r_s is the downstream
    segment of the source path, so the model guard 0 <= r_sk <= r_s can
    never fire mid-fit. FitResult still reports r_sk itself; f_sk pinned at
    its upper bound means the ideal kelvin tap (r_sk = r_s).
    """
    free = tuple(free)
    rk = free.index("r_sk") if "r_sk" in free else None
    if rk is not None:
        bounds_overrides = {**(bounds_overrides or {}), "r_sk": (0.0, 1.0)}
    x0 = pack(start, free)
    if rk is not None:
        x0[rk] = (start.r_sk / start.r_s) if start.r_s > 0 else 0.0
    lo, hi = vector_bounds(free, bounds_overrides)
    x0 = np.clip(x0, lo + 1e-9, hi - 1e-9)

    def expand(xv: np.ndarray) -> Params:
        p = unpack(xv, free, start)
        if rk is None:
            return p
        return p.replace(r_sk=float(np.clip(xv[rk], 0.0, 1.0)) * p.r_s)

    nfev = [0]

    def res_fun(xv: np.ndarray) -> np.ndarray:
        r = residuals_p(expand(xv), data, weight=weight)
        if progress is not None:
            nfev[0] += 1
            if nfev[0] % 25 == 0:
                progress(f"{nfev[0]} evals, cost "
                         f"{0.5 * float(np.dot(r, r)):.4g}")
        return r

    sol = least_squares(
        res_fun,
        x0,
        bounds=(lo, hi),
        method="trf",
        x_scale="jac",
        max_nfev=max_nfev,
        verbose=2 if verbose else 0,
    )

    p = expand(sol.x)
    rel = residuals_p(p, data, weight="relative")
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


def diode_residuals(
    x: np.ndarray,
    free: Sequence[str],
    data: DiodeCurveSet,
    base: Params = DEFAULT,
    *,
    i_d_floor: float = 1e-3,
    weight: str = "relative",
) -> np.ndarray:
    """Residual vector for the diode fit: model i_sd(v_sd, T) minus measured.

    Same weighting semantics as :func:`residuals` (relative by default --
    datasheet diode curves are read off log-ish axes; "abs" weights amperes).
    Each point is evaluated at its own measured temperature.
    """
    p = unpack(x, free, base)
    model = vdmos_diode(data.v_sd, p, t_c=data.temp_c)
    err = model - data.i_sd
    if weight == "relative":
        err = err / np.maximum(np.abs(data.i_sd), i_d_floor)
    elif weight != "abs":
        raise ValueError(f"unknown weight {weight!r}")
    return err


def fit_diode(
    data: DiodeCurveSet,
    *,
    start: Params,
    free: Sequence[str],
    weight: str = "relative",
    max_nfev: int = 4000,
    verbose: bool = False,
    bounds_overrides: dict[str, tuple[float, float]] | None = None,
    progress: Callable[[str], None] | None = None,
) -> FitResult:
    """Extract body-diode parameters from a multi-temperature diode dataset.

    `free` has no default on purpose: the diode free set IS the protocol
    (which mechanisms the data can support) -- see scripts/fit_st8l65_diode.py
    and README, "Body diode". Positive diode parameters are optimised in log
    space (params.LOG_PARAMS); r_diode is linear like r_d/r_s. `start` is
    required because the all-off DEFAULT is a singular starting point for a
    diode fit (zero sensitivity everywhere).
    """
    free = tuple(free)
    x0 = pack(start, free)
    lo, hi = vector_bounds(free, bounds_overrides)
    x0 = np.clip(x0, lo + 1e-9, hi - 1e-9)

    nfev = [0]

    def res_fun(xv: np.ndarray) -> np.ndarray:
        r = diode_residuals(xv, free, data, start, weight=weight)
        if progress is not None:
            nfev[0] += 1
            if nfev[0] % 25 == 0:
                progress(f"{nfev[0]} evals, cost "
                         f"{0.5 * float(np.dot(r, r)):.4g}")
        return r

    sol = least_squares(
        res_fun,
        x0,
        bounds=(lo, hi),
        method="trf",
        x_scale="jac",
        max_nfev=max_nfev,
        verbose=2 if verbose else 0,
    )

    p = unpack(sol.x, free, start)
    rel = diode_residuals(sol.x, free, data, start, weight="relative")
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


def diode_sensitivity_report(
    params: Params,
    data: DiodeCurveSet,
    free: Sequence[str],
    *,
    rel_step: float = 1e-3,
) -> dict[str, float]:
    """Mean |d ln i_sd / d ln p| over the diode dataset, evaluated before fitting.

    Same semantics as :func:`sensitivity_report` (one-sided absolute probe at
    zero-bounded parameters), but evaluates vdmos_diode at each point's own
    temperature.
    """
    base = vdmos_diode(data.v_sd, params, t_c=data.temp_c)
    scale = np.maximum(np.abs(base), 1e-12)
    out: dict[str, float] = {}
    for name in free:
        p0 = getattr(params, name)
        lo_b = BOUNDS.get(name, (-np.inf, np.inf))[0]
        if abs(p0) > 1e-12 and p0 - rel_step * abs(p0) >= lo_b:
            step = rel_step * abs(p0)
            up = vdmos_diode(data.v_sd, params.replace(**{name: p0 + step}), t_c=data.temp_c)
            dn = vdmos_diode(data.v_sd, params.replace(**{name: p0 - step}), t_c=data.temp_c)
            sens = np.abs(up - dn) / (2.0 * step) / scale * abs(p0)
        else:
            step = rel_step
            up = vdmos_diode(data.v_sd, params.replace(**{name: p0 + step}), t_c=data.temp_c)
            sens = np.abs(up - base) / step / scale
        out[name] = float(np.mean(sens))
    return out


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

    Built from the Jacobian at the optimum with unit-weight residuals: these
    are order-of-magnitude identifiability numbers, not statistical
    confidence intervals. Zero-sensitivity columns report `inf` (never 0) and
    are excluded from rank/cond/correlations; the dead-column test is
    relative to the strongest Jacobian column. Rationale in README,
    "Fitter diagnostics".
    """
    jac = result.jac
    free = tuple(free)
    if jac is None:
        return Identifiability({n: math.inf for n in free}, 0, len(free), math.inf, free, ())

    col_norm = np.linalg.norm(jac, axis=0)
    col_max = float(col_norm.max(initial=0.0))
    # Dead test RELATIVE to the strongest column; 1e-5 sits in the empty gap
    # between live (~3e-4 of max) and dead (~1e-6) columns. See README.
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
    """Mean |d ln i_d / d ln p| over the dataset, evaluated before fitting.

    Cheap diagnostic for deciding what to fit at all: a parameter whose
    current value changes i_d by <1e-3 relative across the whole data range
    cannot be extracted from it. Parameters sitting at a zero lower bound
    (r_d, r_s at 0) get a one-sided absolute probe of size `rel_step` -- the
    log-space formula is meaningless at p = 0 and a two-sided probe would
    request negative resistance.
    """
    base = vdmos_id(data.v_ds, data.v_gs, params)
    scale = np.maximum(np.abs(base), 1e-12)
    out: dict[str, float] = {}
    for name in free:
        p0 = getattr(params, name)
        lo_b = BOUNDS.get(name, (-np.inf, np.inf))[0]
        if abs(p0) > 1e-12 and p0 - rel_step * abs(p0) >= lo_b:
            step = rel_step * abs(p0)
            up = vdmos_id(data.v_ds, data.v_gs, params.replace(**{name: p0 + step}))
            dn = vdmos_id(data.v_ds, data.v_gs, params.replace(**{name: p0 - step}))
            sens = np.abs(up - dn) / (2.0 * step) / scale * abs(p0)
        else:
            step = rel_step
            up = vdmos_id(data.v_ds, data.v_gs, params.replace(**{name: p0 + step}))
            sens = np.abs(up - base) / step / scale
        out[name] = float(np.mean(sens))
    return out
