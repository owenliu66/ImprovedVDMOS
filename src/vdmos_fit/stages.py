"""Fitting stages for the release pipeline (session 30, Phase B).

The stage protocol agreed for releases (user, session 29/30):

    channel @ reference temp (all I_ds > 0 rows at T_ref)
    -> diode @ reference temp (I_ds < 0 with gate off)
    -> TC (everything at T != T_ref: Rds(on) chords + Vth specs)
    -> [recovery, caps, export: later phases]

Third-quadrant rows (I_ds < 0 with the gate ON) are NOT diode data --
they are mirror-channel data (session 29 doctrine). They are returned
separately for reporting/validation; feeding them into the channel fit
is a later enhancement (the model already captures them).

The TC stage generalizes the ISC040 recipe: Rds(on) chord residuals as
CURRENT at each row's own test current (resolved per row from the I
column, the condition text, or a user-supplied default), and Vth
residuals on the Vth VOLTAGE (solve the model gate voltage at
V_ds = V_gs giving I_test, vs the datasheet Vth) -- the trick that
avoids the exponential-current sensitivity the single-softplus
subthreshold cannot reproduce.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from scipy.optimize import brentq, least_squares

from .data import CurveSet, DiodeCurveSet
from .fit import FitResult, fit, fit_diode
from .ingest import IngestError, IVTable, RdsSpec, VthSpec, parse_condition
from .model import vdmos_id
from .params import BOUNDS, Params
from .presets import TC_HI, TC_LO, TC_START, Preset

T0K = 298.15


# ---------------------------------------------------------------------------
# table -> curve extraction
# ---------------------------------------------------------------------------

def channel_rows(iv: IVTable, ref_temp: float) -> CurveSet:
    """I_ds > 0 rows at the reference temperature."""
    m = (iv.ids > 0) & np.isclose(iv.tj, ref_temp)
    return CurveSet(iv.vgs[m], iv.vds[m], iv.ids[m],
                    f"channel@{ref_temp:g}C")


def diode_rows(iv: IVTable) -> DiodeCurveSet | None:
    """I_ds < 0 rows with the gate off (V_gs <= 0): body-diode data.

    Returned in the diode convention (I_sd, V_sd positive forward).
    """
    m = (iv.ids < 0) & (iv.vgs <= 0)
    if not np.any(m):
        return None
    return DiodeCurveSet(-iv.ids[m], -iv.vds[m], iv.tj[m],
                         label=f"diode[{iv.source}]")


def third_quadrant_rows(iv: IVTable) -> CurveSet | None:
    """I_ds < 0 rows with the gate ON: mirror-channel data (report/
    validate only for now)."""
    m = (iv.ids < 0) & (iv.vgs > 0)
    if not np.any(m):
        return None
    return CurveSet(iv.vgs[m], iv.vds[m], iv.ids[m],
                    f"third-quadrant[{iv.source}]")


def resolve_test_current(spec: RdsSpec, default: float | None) -> np.ndarray:
    """Per-row Rds(on) chord test current: column > condition > default."""
    it = spec.i_test.copy()
    for i, cond in enumerate(spec.condition):
        if it[i] == 0.0 and cond:
            c = parse_condition(str(cond))
            if "ids" in c:
                it[i] = abs(c["ids"])
    if default is not None:
        it = np.where(it == 0.0, default, it)
    if np.any(it <= 0):
        raise IngestError(spec.source, "R_dson",
                          "chord test current unknown: add an I_d column, "
                          "a condition text, or pass rdson_i_test=")
    return it


def rdson_chord_curves(spec: RdsSpec, itest: np.ndarray,
                       ref_temp: float) -> CurveSet:
    """Rds(on) chord rows at ref temp -> ohmic IV points (V = R*I)."""
    m = np.isclose(spec.tj, ref_temp)
    return CurveSet(spec.vgs[m], spec.rds[m] * itest[m], itest[m],
                    f"rdson-chords@{ref_temp:g}C")


# ---------------------------------------------------------------------------
# stage 1: channel at the reference temperature
# ---------------------------------------------------------------------------

def fit_channel_ref(iv: IVTable, rds_spec: RdsSpec | None, preset: Preset,
                    ref_temp: float = 25.0,
                    rdson_i_test: float | None = None,
                    r_d_upper: float | None = None,
                    multistart_limit: int | None = None,
                    max_nfev: int = 8000,
                    free_override: tuple[str, ...] | None = None,
                    bounds_override: dict[str, tuple[float, float]]
                    | None = None,
                    start_override: dict[str, float] | None = None,
                    progress: Callable[[str], None] | None = None
                    ) -> FitResult:
    free = (tuple(free_override) if free_override
            else preset.channel_free25)
    D = channel_rows(iv, ref_temp)
    bounds: dict[str, tuple[float, float]] = dict(bounds_override or {})
    if rds_spec is not None:
        itest = resolve_test_current(rds_spec, rdson_i_test)
        D = _combine(D, rdson_chord_curves(rds_spec, itest, ref_temp))
        if r_d_upper is None:
            # physical ceiling: the smallest measured chord (r_d is only
            # part of it) -- same discipline as the ISC040 datasheet bound.
            # A digitized Rds(on)-vs-T curve may have NO row exactly at
            # ref temp (IMCQ: closest 16.0/32.1 C) -- then no bound from
            # chords; the fit's default r_d bounds still apply.
            m = np.isclose(rds_spec.tj, ref_temp)
            if np.any(m):
                r_d_upper = float(np.min(rds_spec.rds[m] * itest[m]))
    if r_d_upper is not None and "r_d" in free:
        # physical ceiling: the smallest measured chord wins over any
        # user-set max (the chord is data, the box is a guess)
        lo, hi = bounds.get("r_d", BOUNDS["r_d"])
        bounds["r_d"] = (float(lo), min(float(hi), r_d_upper))
    starts = list(preset.channel_multistart)
    if multistart_limit is not None:
        starts = starts[:multistart_limit]
    if start_override is not None:
        # Interactive GUI: honour the user's live values as the FIRST start
        # verbatim (this is what carries a knob the preset never knew about --
        # e.g. k_kneeform -- through the fit instead of the preset default),
        # then keep the preset multistart deltas on top of that base for
        # robustness. Without start_override the official protocol is unchanged.
        starts = [dict()] + starts
    base = (start_override if start_override is not None
            else preset.channel_start)
    best: FitResult | None = None
    n_start = len(starts)
    for si, s in enumerate(starts):
        vals = {**preset.channel_fixed, **base, **s}
        if "v_ovmax" not in free:
            # user rule (session 30): the charge-cap term is OFF unless
            # v_ovmax is actually fitted -- never export a locked
            # start-value cap.
            vals["v_ovmax"] = 0.0
        start = Params().replace(**vals)
        prog = (None if progress is None else
                (lambda msg, _si=si, _n=n_start: progress(
                    f"start {_si + 1}/{_n}: {msg}")))
        r = fit(D, start=start, free=free,
                weight="relative", max_nfev=max_nfev, bounds_overrides=bounds,
                progress=prog)
        if best is None or r.rmse_rel < best.rmse_rel:
            best = r
    if best is None:
        raise ValueError(f"preset {preset.name!r} has no channel multistarts")
    return best


def _combine(*sets: CurveSet) -> CurveSet:
    return CurveSet(np.concatenate([s.v_gs for s in sets]),
                    np.concatenate([s.v_ds for s in sets]),
                    np.concatenate([s.i_d for s in sets]),
                    "+".join(s.label for s in sets))


# ---------------------------------------------------------------------------
# stage 2: channel temperature coefficients (reference fit locked)
# ---------------------------------------------------------------------------

def params_at(p25: Params, tc_c: float, tc: dict) -> Params:
    """Full Params at temperature tc_c [degC] through the record TC laws."""
    dt = tc_c - 25.0
    tk = tc_c + 273.15
    fac = 1.0 + tc["a_rd"] * dt + tc["b_rd"] * dt * dt
    return Params(**{**p25.as_dict(),
                     "k_p": p25.k_p * (tk / T0K) ** (-tc["beta_p"]),
                     "v_th0": p25.v_th0 + tc["s_vth"] * dt,
                     "r_d": p25.r_d * max(fac, 1e-6)})


def model_ron(p25: Params, tc_c: float, vg: float, itest: float,
              tc: dict) -> float:
    p = params_at(p25, tc_c, tc)
    v = brentq(lambda x: float(vdmos_id(x, vg, p)) - itest, 1e-5, 20.0)
    return v / itest


def model_vth(p25: Params, tc_c: float, itest: float, tc: dict) -> float:
    p = params_at(p25, tc_c, tc)
    return brentq(lambda vg: float(vdmos_id(vg, vg, p)) - itest, 0.5, 8.0)


def fit_channel_tc(p25: Params, rds_spec: RdsSpec | None,
                   vth_spec: VthSpec | None, preset: Preset,
                   rdson_i_test: float | None = None,
                   progress: Callable[[str], None] | None = None) -> dict:
    """TC-stage fit on Rds(on) chords (all T) + Vth specs. Names follow
    preset.tc_free (default s_vth, beta_p, a_rd, b_rd)."""
    if (rds_spec is None or not np.any(~np.isclose(rds_spec.tj, 25.0))) \
            and vth_spec is None:
        raise ValueError("TC stage needs Rds(on) specs at T != ref or "
                         "Vth specs")
    rv: list[tuple[float, float, float, float]] = []
    if rds_spec is not None:
        itest = resolve_test_current(rds_spec, rdson_i_test)
        rv = [(float(t), float(v), float(R), float(i))
              for t, v, R, i in zip(rds_spec.tj, rds_spec.vgs, rds_spec.rds,
                                    itest, strict=True)]
    vt = [(float(t), float(v), float(i))
          for t, v, i in zip(vth_spec.tj, vth_spec.vth, vth_spec.i_test,
                             strict=True)] \
        if vth_spec is not None else []
    names = preset.tc_free
    idx = {n: j for j, n in enumerate(names)}

    def resid(x: np.ndarray) -> np.ndarray:
        tcd = {n: float(x[j]) for n, j in idx.items()}
        for n in ("s_vth", "beta_p", "a_rd", "b_rd"):
            tcd.setdefault(n, 0.0)
        r1 = [(model_ron(p25, t, v, i, tcd) - R) / R for t, v, R, i in rv]
        r2 = [(model_vth(p25, t, i, tcd) - vth) / vth for t, vth, i in vt]
        return np.array(r1 + r2)

    pos = [names.index(n) for n in ("s_vth", "beta_p", "a_rd", "b_rd")
           if n in names]
    x0 = np.array([TC_START[names.index(n)] for n in
                   ("s_vth", "beta_p", "a_rd", "b_rd") if n in names])
    lo = np.array([TC_LO[names.index(n)] for n in
                   ("s_vth", "beta_p", "a_rd", "b_rd") if n in names])
    hi = np.array([TC_HI[names.index(n)] for n in
                   ("s_vth", "beta_p", "a_rd", "b_rd") if n in names])

    def resid4(x4: np.ndarray) -> np.ndarray:
        x = np.full(len(names), np.nan)
        # non-default TC names start at their TC_START value
        for j in range(len(names)):
            if j not in pos:
                x[j] = TC_START[j] if j < len(TC_START) else 0.0
        x[pos] = x4
        r = resid(x)
        if progress is not None:
            nfev[0] += 1
            if nfev[0] % 25 == 0:
                progress(f"{nfev[0]} evals, cost "
                         f"{0.5 * float(np.dot(r, r)):.4g}")
        return r

    nfev = [0]
    sol = least_squares(resid4, x0, bounds=(lo, hi), method="trf",
                        x_scale="jac", max_nfev=6000, ftol=1e-13, xtol=1e-13,
                        gtol=1e-13)
    x = np.full(len(names), 0.0)
    for j in range(len(names)):
        x[j] = TC_START[j] if j < len(TC_START) else 0.0
    x[pos] = sol.x
    rel = resid(x)
    n1 = len(rv)
    tcd = {n: float(x[j]) for n, j in idx.items()}
    return {**tcd,
            "rmse_rdson": float(np.sqrt(np.mean(rel[:n1] ** 2))) if n1 else 0,
            "rmse_vth": float(np.sqrt(np.mean(rel[n1:] ** 2))) if vt else 0,
            "rmse_overall": float(np.sqrt(np.mean(rel ** 2))),
            "success": bool(sol.success)}


# ---------------------------------------------------------------------------
# stage 3: body diode (joint DC + TC)
# ---------------------------------------------------------------------------

def fit_diode_stage(ivs: list[IVTable], preset: Preset,
                    multistart_limit: int | None = None,
                    max_nfev: int = 8000,
                    free_override: tuple[str, ...] | None = None,
                    bounds_override: dict[str, tuple[float, float]]
                    | None = None,
                    progress: Callable[[str], None] | None = None
                    ) -> FitResult | None:
    """Fit the diode from every IV table's gate-off reverse rows.

    Returns None when the preset has no diode free set (SiC: no body-diode
    data) or no diode rows exist -- the caller exports the placeholder
    card with a loud warning.
    """
    free = (tuple(free_override) if free_override
            else tuple(preset.diode_free))
    if not free:
        return None
    sets = [d for d in (diode_rows(iv) for iv in ivs) if d is not None]
    if not sets:
        return None
    data = sets[0] if len(sets) == 1 else _combine_diode(sets)
    starts = preset.diode_multistart
    if multistart_limit is not None:
        starts = starts[:multistart_limit]
    base = Params().replace(**{**preset.diode_fixed, **preset.diode_start})
    bounds = {**preset.diode_bounds, **(bounds_override or {})}
    best: FitResult | None = None
    n_start = len(starts)
    for si, s in enumerate(starts):
        # a start dict must never re-introduce a FIXED parameter
        prog = (None if progress is None else
                (lambda msg, _si=si, _n=n_start: progress(
                    f"start {_si + 1}/{_n}: {msg}")))
        r = fit_diode(data, start=base.replace(**{k: v for k, v in s.items()
                                                  if k in free}),
                      free=free, weight="relative",
                      max_nfev=max_nfev, bounds_overrides=bounds,
                      progress=prog)
        if best is None or r.rmse_rel < best.rmse_rel:
            best = r
    return best


def _combine_diode(sets: list[DiodeCurveSet]) -> DiodeCurveSet:
    return DiodeCurveSet(np.concatenate([s.i_sd for s in sets]),
                         np.concatenate([s.v_sd for s in sets]),
                         np.concatenate([s.temp_c for s in sets]),
                         label="+".join(s.label for s in sets))
